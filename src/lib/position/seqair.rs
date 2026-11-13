//! Seqair-backed base counting with Perbase's filtering and mate rules.
use super::{
    Position,
    mate_fix::{Base, MateRead, MateResolutionStrategy},
    pileup_position::PileupPosition,
};
use crate::read_filter::DefaultReadFilter;
use ::seqair::{
    bam::{
        pileup::{PileupAlignment, PileupEngine},
        record_store::{
            CustomizeRecordStore, FilterRawFields, RecordStore, SlimRecord, qname_hash,
        },
    },
    reader::{IndexedReader, ResolveTid},
};
use anyhow::{Context, Result};
use seqair_types::Pos0;
use std::{cmp::Reverse, collections::BinaryHeap, path::Path, sync::Mutex};

/// Reads indexed regions with Seqair and counts their bases. Workers share
/// the parsed header and index, but own their file handles and read buffers.
pub struct SeqairCounter {
    template: Mutex<IndexedReader>,
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
        let reader = match reference {
            Some(reference) => IndexedReader::open_with_reference(path, reference)?,
            None => IndexedReader::open(path)?,
        };
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
            admission: Admission {
                ends: BinaryHeap::new(),
                previous: None,
                limit: self.max_depth.min(i32::MAX as u32) as usize,
                filter: self.filter,
            },
            info: Vec::new(),
            names: Vec::new(),
            heads: Vec::new(),
            links: Vec::new(),
            active: Vec::new(),
            members: Vec::new(),
        })
    }
}

/// A single counting worker. No output rows are retained between callbacks.
pub struct SeqairWorker<'a> {
    counter: &'a SeqairCounter,
    reader: IndexedReader,
    store: RecordStore<ReadInfo>,
    admission: Admission,
    info: Vec<ReadInfo>,
    names: Vec<usize>,
    heads: Vec<usize>,
    links: Vec<usize>,
    active: Vec<usize>,
    members: Vec<usize>,
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
        self.admission.ends.clear();
        self.admission.previous = None;
        self.reader.fetch_into_customized(
            tid,
            (start..=last).into(),
            &mut self.store,
            &mut self.admission,
        )?;
        let prepared = std::mem::take(&mut self.store).prepare_for_pileup();
        let store = prepared.input.store();
        self.info.clear();
        self.info.extend(
            store
                .indices()
                .map(|index| *store.record(index).unwrap().extra()),
        );
        if self.counter.strategy.is_some() {
            // Store record indices in the name table, not copies of names.
            // Exact byte comparison resolves collisions and includes unpaired
            // and secondary reads, unlike primary-mate links.
            let capacity = (store.len() * 2).next_power_of_two();
            self.names.clear();
            self.names.resize(capacity, usize::MAX);
            self.heads.clear();
            self.heads.resize(store.len(), usize::MAX);
            let records: Vec<_> = store.indices().collect();
            for index in store.indices() {
                let name = store.record(index).unwrap().qname();
                let mut slot = qname_hash(name) as usize & (capacity - 1);
                while self.names[slot] != usize::MAX {
                    let other = records[self.names[slot]];
                    if store.record(other).unwrap().qname() == name {
                        break;
                    }
                    slot = (slot + 1) & (capacity - 1);
                }
                if self.names[slot] == usize::MAX {
                    self.names[slot] = index.as_usize();
                }
                self.info[index.as_usize()].group = self.names[slot];
            }
        }
        let mut engine = PileupEngine::new(prepared.input, (start..=last).into());
        // Seqair's per-position cap would discard different bases. Whole-read
        // admission above implements the HTSlib cap, so leave this uncapped.
        let cutoff = self.counter.max_depth - (self.counter.max_depth as f64 * 0.01) as u32;
        while let Some(column) = engine.pileups() {
            let mut position = PileupPosition::new(name.as_str().into(), column.pos().as_u32());
            position.depth = column.depth() as u32;
            position.near_max_depth = position.depth > cutoff;
            if let Some(strategy) = self.counter.strategy {
                self.active.clear();
                self.links.clear();
                for (index, alignment) in column.raw_alignments().enumerate() {
                    let group = self.info[alignment.record_idx().as_usize()].group;
                    if self.heads[group] == usize::MAX {
                        self.active.push(group);
                    }
                    self.links.push(self.heads[group]);
                    self.heads[group] = index;
                }
                let alignment = |index| column.alignment_at(index).unwrap().alignment();
                let passes =
                    |read: &PileupAlignment| self.info[read.record_idx().as_usize()].passes;
                for &group in &self.active {
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
            emit(position)?;
        }
        self.store = engine.reclaim_allocation().unwrap_or_default();
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct ReadInfo {
    group: usize,
    passes: bool,
}

// HTSlib limits admission of whole reads, not the number of bases in each
// column. Reads ending exactly at the next start still occupy a slot, and
// the first read at a new start bypasses the cap (except at tid=0, pos=0).
#[derive(Clone)]
struct Admission {
    ends: BinaryHeap<Reverse<i64>>,
    previous: Option<i64>,
    limit: usize,
    filter: DefaultReadFilter,
}

impl CustomizeRecordStore for Admission {
    type Extra = ReadInfo;

    fn filter_raw(&mut self, fields: &FilterRawFields<'_>) -> bool {
        if fields.flags.is_unmapped() || fields.tid.is_none() {
            return false;
        }
        let start = fields.pos.as_i64();
        while self.ends.peek().is_some_and(|end| end.0 < start) {
            self.ends.pop();
        }
        let at_cursor = self.previous == Some(start)
            || (self.previous.is_none() && fields.tid.unwrap().as_u32() == 0 && start == 0);
        if at_cursor && self.ends.len() >= self.limit {
            return false;
        }
        // Use the actual reference span, including zero-reference CIGARs,
        // rather than Seqair's normalized one-position end for those reads.
        let span: i64 = fields
            .cigar
            .chunks_exact(4)
            .map(|bytes| {
                let op = u32::from_le_bytes(bytes.try_into().unwrap());
                if matches!(op & 15, 0 | 2 | 3 | 7 | 8) {
                    i64::from(op >> 4)
                } else {
                    0
                }
            })
            .sum();
        if span > 0 || !at_cursor {
            self.ends.push(Reverse(start + span));
        }
        self.previous = Some(start);
        true
    }

    fn compute(&mut self, read: &SlimRecord, _store: &RecordStore<ReadInfo>) -> ReadInfo {
        ReadInfo {
            group: 0,
            passes: self
                .filter
                .filter_fields(read.flags.raw(), read.mapq.as_byte()),
        }
    }
}

impl MateRead for PileupAlignment {
    fn flags(&self) -> u16 {
        self.flags.raw()
    }
    fn mapq(&self) -> u8 {
        self.mapq.as_byte()
    }
    fn qpos(&self) -> Option<usize> {
        self.qpos().map(|pos| pos.as_usize())
    }
    fn base(&self) -> Base {
        match self.base() {
            Some(seqair_types::Base::A) => Base::A,
            Some(seqair_types::Base::C) => Base::C,
            Some(seqair_types::Base::G) => Base::G,
            Some(seqair_types::Base::T) => Base::T,
            _ => Base::N,
        }
    }
    fn base_qual(&self) -> u8 {
        self.qual().map_or(255, |qual| qual.as_byte())
    }
}

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
    if read.is_refskip() {
        position.ref_skip += 1;
        position.depth -= 1;
    } else if read.is_del() {
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
        if read.insert_len() > 0 {
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
