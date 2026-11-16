//! Seqair-backed base counting with Perbase's filtering and mate rules.
use super::{
    Position,
    mate_fix::{Base, MateRead, MateResolutionStrategy},
    pileup_position::PileupPosition,
};
use crate::read_filter::DefaultReadFilter;
use ::seqair::{
    bam::{
        pileup::{PileupAlignment, PileupEngine, PileupOp},
        record_idx::RecordIdx,
        record_store::{CustomizeRecordStore, FilterRawFields, RecordStore, SlimRecord},
    },
    reader::{AuxTags, IndexedReader, PartialIndexedReader, RecordFields, ResolveTid},
};
use anyhow::{Context, Result};
use seqair_types::Pos0;
use std::{cmp::Reverse, collections::BinaryHeap, path::Path, sync::Mutex};

/// Reads indexed regions with Seqair and counts their bases. Workers share
/// the parsed header and index, but own their file handles and read buffers.
pub struct SeqairCounter {
    template: Mutex<PartialIndexedReader>,
    filter: DefaultReadFilter,
    strategy: Option<MateResolutionStrategy>,
    max_depth: u32,
    min_base_quality: Option<u8>,
}

impl SeqairCounter {
    /// Open the alignment and its index. Supplying a reference also enables
    /// CRAM files whose bases must be reconstructed from that reference.
    pub fn open(
        path: &Path,
        reference: Option<&Path>,
        filter: DefaultReadFilter,
        strategy: Option<MateResolutionStrategy>,
        max_depth: u32,
        min_base_quality: Option<u8>,
    ) -> Result<Self> {
        let mut builder =
            IndexedReader::builder(path).fields(RecordFields::all().with_aux(AuxTags::None));
        if let Some(reference) = reference {
            builder = builder.reference(reference);
        }
        let reader = builder.open()?;
        Ok(Self {
            template: Mutex::new(reader),
            filter,
            strategy,
            max_depth,
            min_base_quality,
        })
    }

    /// Count the half-open interval `[start, stop)`. Positions returned here
    /// are zero-based; the caller adds coordinate conversion, uncovered rows
    /// and reference annotations just as it does for HTSlib results.
    pub fn count_region(&self, tid: u32, start: u32, stop: u32) -> Result<Vec<PileupPosition>> {
        let mut positions = Vec::new();
        self.worker()?
            .for_each_position(tid, start, stop, |position| {
                positions.push(position);
                Ok(())
            })?;
        Ok(positions)
    }

    /// Create a worker whose reader and allocations can be reused across regions.
    pub fn worker(&self) -> Result<SeqairWorker<'_>> {
        Ok(SeqairWorker {
            counter: self,
            reader: self.template.lock().expect("Seqair reader lock").fork()?,
            store: RecordStore::new(),
            ends: BinaryHeap::new(),
            info: Vec::new(),
            names: Vec::new(),
            heads: Vec::new(),
            links: Vec::new(),
            active: Vec::new(),
            members: Vec::new(),
            group_ends: Vec::new(),
            overlaps: Vec::new(),
        })
    }
}

/// A single counting worker. No output rows are retained between callbacks.
pub struct SeqairWorker<'a> {
    counter: &'a SeqairCounter,
    reader: PartialIndexedReader,
    store: RecordStore,
    ends: BinaryHeap<Reverse<i64>>,
    info: Vec<ReadInfo>,
    names: Vec<Option<RecordIdx>>,
    heads: Vec<usize>,
    links: Vec<usize>,
    active: Vec<usize>,
    members: Vec<usize>,
    group_ends: Vec<Option<u32>>,
    overlaps: Vec<bool>,
}

impl SeqairWorker<'_> {
    /// Count `[start, stop)` and immediately pass each zero-based row to `emit`.
    /// The callback may stop processing by returning an error, such as a failed
    /// output write. Read storage is reused on the next successful region.
    pub fn for_each_position(
        &mut self,
        tid: u32,
        start: u32,
        stop: u32,
        mut emit: impl FnMut(PileupPosition) -> Result<()>,
    ) -> Result<()> {
        self.for_each_position_ref(tid, start, stop, |position| emit(position.clone()))
    }

    /// Count `[start, stop)` without copying the reference name for every row.
    /// The callback borrows a row that is reused at the next position. Clone
    /// it only when retaining rows rather than immediately writing them.
    pub fn for_each_position_ref(
        &mut self,
        tid: u32,
        start: u32,
        stop: u32,
        mut emit: impl FnMut(&PileupPosition) -> Result<()>,
    ) -> Result<()> {
        anyhow::ensure!(start <= stop, "region start must not exceed its end");
        if start == stop {
            return Ok(());
        }
        let start = Pos0::new(start).context("region start exceeds Seqair's coordinate range")?;
        let last = Pos0::new(stop - 1).context("region end exceeds Seqair's coordinate range")?;
        let tid = tid.resolve_tid(self.reader.header())?;
        let name = self
            .reader
            .header()
            .target_name(tid)
            .context("unknown Seqair reference")?
            .to_owned();
        let mut reads = CountingReads;
        self.reader.fetch_into_customized(
            tid,
            (start..=last).into(),
            &mut self.store,
            &mut reads,
        )?;
        self.info.clear();
        self.info.extend(self.store.records().map(|read| {
            ReadInfo {
                group: 0,
                passes: self
                    .counter
                    .filter
                    .filter_fields(read.flags.raw(), read.mapq.as_byte()),
                admitted: true,
            }
        }));
        let restricted = self.admit()?;
        let input = std::mem::take(&mut self.store).prepare_for_pileup().input;
        let store = input.store();
        if self.counter.strategy.is_some() {
            // Store record indices in the name table, not copies of names.
            // Exact byte comparison resolves collisions and includes unpaired
            // and secondary reads, unlike primary-mate links.
            let capacity = (store.len() * 2).next_power_of_two();
            self.names.clear();
            self.names.resize(capacity, None);
            self.heads.clear();
            self.heads.resize(store.len(), usize::MAX);
            self.group_ends.clear();
            self.group_ends.resize(store.len(), None);
            self.overlaps.clear();
            self.overlaps.resize(store.len(), false);
            for index in store.indices() {
                let record = store.record(index).unwrap();
                let name = record.qname();
                let mut slot =
                    record.qname_hash().map_or(0, |hash| hash.get()) as usize & (capacity - 1);
                while let Some(other) = self.names[slot] {
                    if store.record(other).unwrap().qname() == name {
                        break;
                    }
                    slot = (slot + 1) & (capacity - 1);
                }
                let group = *self.names[slot].get_or_insert(index);
                let group = group.as_usize();
                self.info[index.as_usize()].group = group;
                let start = record.pos.as_u32();
                let end = record.end_pos.as_u32();
                // Sorted starts and inclusive ends let us prove that no two
                // records in this name group can contribute to one column.
                if let Some(previous_end) = self.group_ends[group] {
                    self.overlaps[group] |= start <= previous_end;
                    self.group_ends[group] = Some(previous_end.max(end));
                } else {
                    self.group_ends[group] = Some(end);
                }
            }
        }
        let mut engine = PileupEngine::new(input, (start..=last).into());
        // Seqair's per-position cap would discard different bases. Whole-read
        // admission above implements the HTSlib cap, so leave this uncapped.
        let cutoff = self.counter.max_depth - (self.counter.max_depth as f64 * 0.01) as u32;
        let mut position = PileupPosition::new(name.as_str().into(), start.as_u32());
        while let Some(column) = engine.pileups() {
            let depth = if restricted {
                column
                    .raw_alignments()
                    .filter(|read| self.info[read.record_idx().as_usize()].admitted)
                    .count()
            } else {
                column.depth()
            } as u32;
            if depth == 0 {
                continue;
            }
            position = PileupPosition {
                ref_seq: position.ref_seq,
                pos: column.pos().as_u32(),
                depth,
                near_max_depth: depth > cutoff,
                ..Default::default()
            };
            if let Some(strategy) = self.counter.strategy {
                self.active.clear();
                self.links.clear();
                self.members.clear();
                self.active.reserve(column.depth());
                self.links.reserve(column.depth());
                self.members.reserve(column.depth());
                for (index, alignment) in column.raw_alignments().enumerate() {
                    let info = self.info[alignment.record_idx().as_usize()];
                    if restricted && !info.admitted {
                        self.links.push(usize::MAX);
                        continue;
                    }
                    let group = info.group;
                    if !self.overlaps[group] {
                        self.links.push(usize::MAX);
                        update(
                            &mut position,
                            alignment,
                            info.passes,
                            self.counter.min_base_quality,
                            None,
                            false,
                        );
                        continue;
                    }
                    if self.heads[group] == usize::MAX {
                        self.active.push(group);
                    }
                    self.links.push(self.heads[group]);
                    self.heads[group] = index;
                }
                let alignment = |index| column.alignment_at(index).unwrap().alignment();
                let passes =
                    |read: &PileupAlignment| self.info[read.record_idx().as_usize()].passes;
                for &group in self.active.iter() {
                    let last = self.heads[group];
                    let previous = self.links[last];
                    if previous == usize::MAX {
                        let best = alignment(last);
                        update(
                            &mut position,
                            best,
                            passes(best),
                            self.counter.min_base_quality,
                            None,
                            false,
                        );
                        self.heads[group] = usize::MAX;
                        continue;
                    }
                    if self.links[previous] == usize::MAX {
                        // Preserve the two-element stable sort's comparison
                        // order, including ties and recommended-base selection.
                        let (a, b) = (alignment(previous), alignment(last));
                        let swap = strategy.compare(b, a, passes(b), passes(a)).ordering
                            == std::cmp::Ordering::Greater;
                        let (best, second) = if swap { (b, a) } else { (a, b) };
                        let base = strategy
                            .compare(best, second, passes(best), passes(second))
                            .recommended_base;
                        position.depth -= 1;
                        update(
                            &mut position,
                            best,
                            passes(best),
                            self.counter.min_base_quality,
                            base,
                            true,
                        );
                        self.heads[group] = usize::MAX;
                        continue;
                    }
                    self.members.clear();
                    let mut current = self.heads[group];
                    while current != usize::MAX {
                        self.members.push(current);
                        current = self.links[current];
                    }
                    self.members.reverse();
                    let reads = &mut self.members;
                    reads.sort_by(|&a, &b| {
                        let (a, b) = (alignment(a), alignment(b));
                        strategy
                            .compare(a, b, passes(a), passes(b))
                            .ordering
                            .reverse()
                    });
                    let best = alignment(reads[0]);
                    let base = reads.get(1).and_then(|&second| {
                        let second = alignment(second);
                        strategy
                            .compare(best, second, passes(best), passes(second))
                            .recommended_base
                    });
                    position.depth -= reads.len() as u32 - 1;
                    update(
                        &mut position,
                        best,
                        passes(best),
                        self.counter.min_base_quality,
                        base,
                        reads.len() > 1,
                    );
                    self.heads[group] = usize::MAX;
                }
            } else {
                for read in column.raw_alignments() {
                    if restricted && !self.info[read.record_idx().as_usize()].admitted {
                        continue;
                    }
                    update(
                        &mut position,
                        read,
                        self.info[read.record_idx().as_usize()].passes,
                        self.counter.min_base_quality,
                        None,
                        false,
                    );
                }
            }
            emit(&position)?;
        }
        self.store = engine.reclaim_allocation().unwrap_or_default();
        Ok(())
    }

    /// Mark whole reads excluded by HTSlib's depth-admission rule. Small
    /// regions cannot hit the limit and need no heap or CIGAR scan.
    /// Reads ending exactly at the next start still occupy a slot. The first
    /// read at a new start bypasses the cap, except at reference 0, position 0.
    /// Like Goose, this runs after loading the region: the depth limit does
    /// not limit how much read data is retained in memory.
    fn admit(&mut self) -> Result<bool> {
        let limit = self.counter.max_depth.min(i32::MAX as u32) as usize;
        if limit > 0 && self.store.len() <= limit {
            return Ok(false);
        }
        self.ends.clear();
        let mut previous = None;
        let mut restricted = false;
        for (read, info) in self.store.records().zip(&mut self.info) {
            let start = read.pos.as_i64();
            while self.ends.peek().is_some_and(|end| end.0 < start) {
                self.ends.pop();
            }
            let at_cursor = previous == Some(start)
                || (previous.is_none() && read.tid.unwrap().as_u32() == 0 && start == 0);
            if at_cursor && self.ends.len() >= limit {
                info.admitted = false;
                restricted = true;
                continue;
            }
            // Seqair normalizes both zero- and one-base spans to end == pos.
            // Only that ambiguous case needs the CIGAR; longer spans are cached.
            let span = if read.end_pos == read.pos {
                read.cigar(&self.store)?
                    .iter()
                    .filter(|op| op.consumes_ref())
                    .map(|op| i64::from(op.len()))
                    .sum()
            } else {
                read.end_pos.as_i64() - start + 1
            };
            if span > 0 || !at_cursor {
                self.ends.push(Reverse(start + span));
            }
            previous = Some(start);
        }
        Ok(restricted)
    }
}

#[derive(Clone, Copy)]
struct ReadInfo {
    group: usize,
    passes: bool,
    admitted: bool,
}

#[derive(Clone, Default)]
struct CountingReads;

impl CustomizeRecordStore for CountingReads {
    type Extra = ();

    fn filter_raw(&mut self, fields: &FilterRawFields<'_>) -> bool {
        !fields.flags.is_unmapped() && fields.tid.is_some()
    }

    fn compute(&mut self, _: &SlimRecord, _: &RecordStore) {}
}

impl MateRead for PileupAlignment {
    fn flags(&self) -> u16 {
        self.flags.raw()
    }
    fn mapq(&self) -> u8 {
        self.mapq.as_byte()
    }
    fn qpos(&self) -> Option<usize> {
        match self.op {
            PileupOp::Match { qpos, .. }
            | PileupOp::Insertion { qpos, .. }
            | PileupOp::SoftClip { qpos, .. } => Some(qpos.as_usize()),
            _ => None,
        }
    }
    fn base(&self) -> Base {
        let base = match self.op {
            PileupOp::Match { base, .. }
            | PileupOp::Insertion { base, .. }
            | PileupOp::SoftClip { base, .. } => base,
            _ => seqair_types::Base::Unknown,
        };
        match base {
            seqair_types::Base::A => Base::A,
            seqair_types::Base::C => Base::C,
            seqair_types::Base::G => Base::G,
            seqair_types::Base::T => Base::T,
            _ => Base::N,
        }
    }
    fn base_qual(&self) -> u8 {
        match self.op {
            PileupOp::Match { qual, .. }
            | PileupOp::Insertion { qual, .. }
            | PileupOp::SoftClip { qual, .. } => qual.as_byte(),
            _ => 255,
        }
    }
}

#[inline]
fn update(
    position: &mut PileupPosition,
    read: &PileupAlignment,
    passes: bool,
    min_quality: Option<u8>,
    recommended: Option<Base>,
    mates_resolved: bool,
) {
    if !passes {
        position.depth -= 1;
        position.fail += 1;
        return;
    }
    if matches!(
        read.op,
        PileupOp::RefSkip
            | PileupOp::ComplexIndel {
                is_refskip: true,
                ..
            }
    ) {
        position.ref_skip += 1;
        position.depth -= 1;
    } else if matches!(
        read.op,
        PileupOp::Deletion { .. } | PileupOp::ComplexIndel { .. }
    ) {
        position.del += 1;
    } else {
        let base = if min_quality.is_some_and(|q| read.seq_len == 0 || read.base_qual() < q) {
            Base::N
        } else {
            recommended.unwrap_or_else(|| MateRead::base(read))
        };
        match base {
            Base::A => position.a += 1,
            Base::C => position.c += 1,
            Base::G => position.g += 1,
            Base::T => position.t += 1,
            Base::R => position.r += 1,
            Base::Y => position.y += 1,
            Base::S => position.s += 1,
            Base::W => position.w += 1,
            Base::K => position.k += 1,
            Base::M => position.m += 1,
            _ => position.n += 1,
        }
        if matches!(read.op, PileupOp::Insertion { insert_len, .. } if insert_len > 0) {
            position.ins += 1;
        }
    }
    if mates_resolved {
        position.count_of_mate_resolutions += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::position::pileup_position::MateCounter;
    use proptest::prelude::*;
    use rust_htslib::bam::{self, Read};

    // name, contig, start, CIGAR shape, flags, MAPQ, sequence offset, quality
    type ReadSpec = (u8, u8, u8, u8, u8, u8, u8, u8);

    fn fixture(specs: &[ReadSpec]) -> (tempfile::TempDir, std::path::PathBuf) {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reads.bam");
        let mut header = bam::Header::new();
        for name in ["a", "b"] {
            header.push_record(
                bam::header::HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", name)
                    .push_tag(b"LN", 64),
            );
        }
        let view = bam::HeaderView::from_header(&header);
        let mut records = Vec::new();
        for &(name, tid, start, shape, flags, mapq, base, quality) in specs {
            let qname = if name == 0 {
                "*".to_owned()
            } else {
                format!("name-{name}")
            };
            let sequence: String = (0..5)
                .map(|i| b"ACGTN"[(usize::from(base) + i) % 5] as char)
                .collect();
            let qualities = if quality == 3 {
                "*".to_owned()
            } else {
                std::iter::repeat_n((33 + [0u8, 20, 40][quality as usize]) as char, 5).collect()
            };
            let cigar = ["5M", "2M2D3M", "2M3N3M", "1S2M1I1M", "5I"][shape as usize];
            let flags = [0, 4, 65, 129, 256, 2048][flags as usize];
            let mapq = [0, 20, 60, 255][mapq as usize];
            let sam = format!(
                "{qname}\t{flags}\t{}\t{}\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t{qualities}",
                ["a", "b"][tid as usize],
                start + 1
            );
            records.push(bam::Record::from_sam(&view, sam.as_bytes()).unwrap());
        }
        records.sort_by_key(|read| (read.tid(), read.pos()));
        let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
        for record in records {
            writer.write(&record).unwrap();
        }
        drop(writer);
        bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();
        (directory, path)
    }

    fn tsv(positions: &[PileupPosition]) -> Vec<u8> {
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(Vec::new());
        for position in positions {
            writer.serialize(position).unwrap();
        }
        writer.into_inner().unwrap()
    }

    fn compare(
        path: &Path,
        strategy: Option<MateResolutionStrategy>,
        depth: u32,
        filter: DefaultReadFilter,
        baseq: Option<u8>,
    ) {
        let counter = SeqairCounter::open(path, None, filter, strategy, depth, baseq).unwrap();
        let mut reader = bam::IndexedReader::from_path(path).unwrap();
        let header = reader.header().clone();
        let mut worker = counter.worker().unwrap();
        for (tid, start, stop) in [(0, 0, 64), (0, 5, 18), (1, 0, 64), (0, 5, 18), (1, 8, 8)] {
            reader.fetch((tid, start, stop)).unwrap();
            let mut pileups = reader.pileup();
            pileups.set_max_depth(depth);
            let mut mates = MateCounter::default();
            let mut expected = Vec::new();
            for pileup in pileups {
                let pileup = pileup.unwrap();
                if pileup.pos() < start || pileup.pos() >= stop {
                    continue;
                }
                let raw_depth = pileup.depth();
                let mut position = match strategy {
                    Some(strategy) => mates.count(pileup, &header, &filter, baseq, strategy),
                    None => PileupPosition::from_pileup(pileup, &header, &filter, baseq),
                };
                position.near_max_depth = raw_depth > depth - (depth as f64 * 0.01) as u32;
                expected.push(position);
            }
            let mut actual = Vec::new();
            worker
                .for_each_position(tid, start, stop, |position| {
                    actual.push(position);
                    Ok(())
                })
                .unwrap();
            assert_eq!(
                String::from_utf8(tsv(&actual)).unwrap(),
                String::from_utf8(tsv(&expected)).unwrap(),
                "strategy={strategy:?}, depth={depth}, region={tid}:{start}-{stop}"
            );
        }
    }

    const STRATEGIES: [MateResolutionStrategy; 9] = [
        MateResolutionStrategy::Original,
        MateResolutionStrategy::IUPAC,
        MateResolutionStrategy::N,
        MateResolutionStrategy::BaseQualMapQualFirstInPair,
        MateResolutionStrategy::BaseQualMapQualIUPAC,
        MateResolutionStrategy::BaseQualMapQualN,
        MateResolutionStrategy::MapQualBaseQualFirstInPair,
        MateResolutionStrategy::MapQualBaseQualIUPAC,
        MateResolutionStrategy::MapQualBaseQualN,
    ];

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(16))]
        #[test]
        fn matches_htslib_for_overlaps_filters_and_depth_caps(
            reads in prop::collection::vec((0u8..7, 0u8..2, 0u8..20, 0u8..5, 0u8..6, 0u8..4, 0u8..5, 0u8..4), 0..45)
        ) {
            let (_directory, path) = fixture(&reads);
            for strategy in [None].into_iter().chain(STRATEGIES.map(Some)) {
                for depth in [0, 1, 2, 7, 100000] {
                    compare(&path, strategy, depth, DefaultReadFilter::new(0, 0, 0), None);
                    compare(&path, strategy, depth, DefaultReadFilter::new(1, 256, 20), Some(20));
                }
            }
        }
    }

    #[test]
    fn disjoint_name_groups_and_end_boundaries_match_htslib() {
        // Abutting spans, a one-base overlap, and a later overlap after a
        // disjoint prefix. The same worker is reused by compare().
        for starts in [[0, 5, 10], [0, 4, 10], [0, 5, 9], [0, 0, 10]] {
            for shape in 0..5 {
                let reads: Vec<_> = starts
                    .into_iter()
                    .enumerate()
                    .map(|(i, start)| (1, 0, start, shape, i as u8 + 2, 2, i as u8, 2))
                    .collect();
                let (_directory, path) = fixture(&reads);
                for strategy in STRATEGIES {
                    for depth in [1, 7, 100000] {
                        compare(
                            &path,
                            Some(strategy),
                            depth,
                            DefaultReadFilter::new(0, 0, 0),
                            None,
                        );
                        compare(
                            &path,
                            Some(strategy),
                            depth,
                            DefaultReadFilter::new(0, 256, 20),
                            Some(20),
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn counts_long_cigars_instead_of_the_placeholder() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("long.bam");
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "a")
                .push_tag(b"LN", 100000),
        );
        let mut record = bam::Record::new();
        let cigar = bam::record::CigarString(
            [bam::record::Cigar::Match(1), bam::record::Cigar::Ins(1)].repeat(32768),
        );
        record.set(
            b"long",
            Some(&cigar),
            &b"AC".repeat(32768),
            &vec![40; 65536],
        );
        record.set_tid(0);
        record.set_pos(1);
        record.set_flags(0);
        let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
        writer.write(&record).unwrap();
        drop(writer);
        bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();
        let counter = SeqairCounter::open(
            &path,
            None,
            DefaultReadFilter::new(0, 0, 0),
            None,
            100000,
            None,
        )
        .unwrap();
        let actual = counter.count_region(0, 1, 8).unwrap();
        // Each 1M1I contributes an A match with a following C insertion.
        // The BAM placeholder would instead describe a reference skip.
        let expected: Vec<_> = (1..8)
            .map(|pos| PileupPosition {
                ref_seq: "a".into(),
                pos,
                depth: 1,
                a: 1,
                ins: 1,
                ..Default::default()
            })
            .collect();
        assert_eq!(tsv(&actual), tsv(&expected));
    }

    #[test]
    fn preserves_ties_in_large_same_name_groups() {
        for size in [2, 3, 16, 32, 80] {
            let reads: Vec<_> = (0..size).map(|i| (1, 0, 0, 0, 0, 2, i % 5, 2)).collect();
            let (_directory, path) = fixture(&reads);
            for strategy in STRATEGIES {
                compare(
                    &path,
                    Some(strategy),
                    100000,
                    DefaultReadFilter::new(0, 0, 0),
                    None,
                );
            }
        }
    }
}
