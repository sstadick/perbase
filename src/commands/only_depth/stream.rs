//! Reuse readers and coverage arrays; emit depth rows without retaining them.
use super::{OnlyDepth, OnlyDepthProcessor, ReadFilter};
use crate::commands::stream::{self, Config, Output};
#[cfg(feature = "seqair")]
use anyhow::Context;
use anyhow::Result;
use perbase_lib::{
    par_granges::Region,
    position::{Position, range_positions::RangePositions},
    utils,
};
use rust_htslib::bam::{self, Read, ext::BamRecordExtensions};
use std::{collections::HashMap, io::Write};

#[cfg(feature = "seqair")]
use seqair::reader::{
    IndexedReader, PartialIndexedReader, PartialRecord, RecordFields, ResolveTid,
};
#[cfg(feature = "seqair")]
use seqair_types::Pos0;

enum Reader {
    Htslib(bam::IndexedReader, bam::Record),
    #[cfg(feature = "seqair")]
    Seqair(PartialIndexedReader, PartialRecord),
}

#[derive(Clone, Copy)]
struct Block {
    start: usize,
    stop: usize,
    no_stop: bool,
}

#[derive(Default)]
struct Counter {
    delta: Vec<i64>,
    groups: HashMap<Vec<u8>, Vec<Block>>,
}

impl Counter {
    fn reset(&mut self, length: usize) {
        self.delta.clear();
        self.delta.resize(length, 0);
        self.groups.clear();
    }

    fn count(delta: &mut [i64], block: Block) {
        delta[block.start] += 1;
        if !block.no_stop {
            delta[block.stop] -= 1;
        }
    }

    fn add(&mut self, name: Option<&[u8]>, begin: i64, end: i64, range: (u32, u32)) {
        let (start, stop) = (i64::from(range.0), i64::from(range.1));
        if begin >= stop || end <= start {
            return;
        }
        let block = Block {
            start: (begin.max(start) - start) as usize,
            stop: (end - start).min(stop - start - 1) as usize,
            no_stop: end >= stop,
        };
        if let Some(name) = name {
            // Borrow for lookup so the second mate does not allocate a name.
            if let Some(intervals) = self.groups.get_mut(name) {
                intervals.push(block);
            } else {
                self.groups.insert(name.to_vec(), vec![block]);
            }
        } else {
            Self::count(&mut self.delta, block);
        }
    }

    fn merge_mates(&mut self) {
        for intervals in self.groups.values_mut() {
            // Stable order preserves the first interval's clipped-end flag,
            // matching rust-lapper even when later intervals extend further.
            if intervals.len() > 1 {
                intervals.sort_by_key(|block| (block.start, block.stop));
            }
            let mut current = intervals[0];
            for block in &intervals[1..] {
                if block.start <= current.stop {
                    current.stop = current.stop.max(block.stop);
                } else {
                    Self::count(&mut self.delta, current);
                    current = *block;
                }
            }
            Self::count(&mut self.delta, current);
        }
    }

    fn emit(
        &self,
        processor: &OnlyDepthProcessor,
        start: u32,
        stop: u32,
        mut row: impl FnMut(u32, u32, u32) -> Result<()>,
    ) -> Result<()> {
        let mut sum = 0;
        let mut begin = start;
        let mut depth = self.delta[0];
        for (offset, change) in self.delta.iter().enumerate() {
            sum += change;
            let pos = start + offset as u32;
            if processor.no_merge {
                if sum != 0 || processor.keep_zeros {
                    row(pos, pos + 1, u32::try_from(sum)?)?;
                }
            } else if depth != sum {
                row(begin, pos, u32::try_from(depth)?)?;
                begin = pos;
                depth = sum;
            }
        }
        if !processor.no_merge {
            row(begin, stop, u32::try_from(depth)?)?;
        }
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn run(
    processor: &OnlyDepthProcessor,
    regions: &[Region],
    workers: usize,
    chunksize: u32,
    modifier: f64,
    bed: bool,
    use_seqair: bool,
    output: &mut dyn Write,
) -> Result<()> {
    let names =
        utils::with_indexed_reader(&processor.reads, processor.ref_fasta.as_deref(), |r| {
            r.header()
                .target_names()
                .iter()
                .map(|name| String::from_utf8(name.to_vec()))
                .collect::<std::result::Result<Vec<_>, _>>()
        })?;
    let longest = names.iter().map(String::len).max().unwrap_or(0);
    #[cfg(feature = "seqair")]
    let template = if use_seqair {
        let fields = if processor.mate_fix {
            RecordFields::core().with_names()
        } else {
            RecordFields::core()
        };
        let mut builder = IndexedReader::builder(&processor.reads).fields(fields);
        if let Some(reference) = &processor.ref_fasta {
            builder = builder.reference(reference);
        }
        Some(builder.open()?)
    } else {
        None
    };
    #[cfg(not(feature = "seqair"))]
    let _ = use_seqair;
    stream::run(
        regions,
        Config {
            workers,
            chunksize,
            modifier,
            estimated_row_bytes: longest + 32,
        },
        if bed { b"" } else { b"REF\tPOS\tEND\tDEPTH\n" },
        output,
        || {
            #[cfg(feature = "seqair")]
            if let Some(template) = &template {
                return Ok((
                    Reader::Seqair(template.fork()?, PartialRecord::default()),
                    Counter::default(),
                ));
            }
            let mut reader = bam::IndexedReader::from_path(&processor.reads)?;
            if let Some(reference) = &processor.ref_fasta {
                reader.set_reference(reference)?;
            }
            Ok((
                Reader::Htslib(reader, bam::Record::new()),
                Counter::default(),
            ))
        },
        |(reader, counter), (tid, start, stop), output| {
            if start == stop {
                return Ok(());
            }
            counter.reset((stop - start) as usize);
            load(reader, counter, processor, (tid, start, stop))?;
            counter.merge_mates();
            let name = stream::quote_name(&names[tid as usize]);
            counter.emit(processor, start, stop, |begin, end, depth| {
                append_row(
                    output,
                    &name,
                    begin + processor.coord_base,
                    end + processor.coord_base,
                    depth,
                    bed,
                )
            })
        },
    )
}

fn load(
    reader: &mut Reader,
    counter: &mut Counter,
    processor: &OnlyDepthProcessor,
    (tid, start, stop): Region,
) -> Result<()> {
    match reader {
        Reader::Htslib(reader, record) => {
            reader.fetch((tid, start, stop))?;
            while let Some(result) = reader.read(record) {
                result?;
                if !processor.read_filter.filter_read(record, None) {
                    continue;
                }
                let mut pos = record.reference_start();
                let name = (processor.mate_fix && OnlyDepth::maybe_overlaps_mate(record))
                    .then(|| record.qname());
                if processor.fast_mode {
                    counter.add(name, pos, record.reference_end(), (start, stop));
                    continue;
                }
                for op in record.cigar().iter() {
                    let end = pos + i64::from(op.len());
                    match op.char() {
                        'M' | '=' | 'X' | 'D' => {
                            counter.add(name, pos, end, (start, stop));
                            pos = end;
                        }
                        'N' => pos = end,
                        _ => (),
                    }
                }
            }
        }
        #[cfg(feature = "seqair")]
        Reader::Seqair(reader, record) => {
            let tid = tid.resolve_tid(reader.header())?;
            let first =
                Pos0::new(start).context("region start exceeds Seqair's coordinate range")?;
            let last =
                Pos0::new(stop - 1).context("region end exceeds Seqair's coordinate range")?;
            let mut query = reader.query_raw(tid, (first..=last).into())?;
            while query.read_record(record)? {
                if !processor
                    .read_filter
                    .filter_fields(record.flags().raw(), record.mapq().as_byte())
                {
                    continue;
                }
                let mut pos = record
                    .pos()
                    .context("indexed record has no position")?
                    .as_i64();
                let end = record
                    .end_pos()
                    .context("indexed record has no end")?
                    .as_i64()
                    + 1;
                if processor.fast_mode {
                    counter.add(None, pos, end, (start, stop));
                    continue;
                }
                let mate = record.mate_fields();
                let overlaps = processor.mate_fix
                    && record.flags().raw() & 1 != 0
                    && record.tid() == mate.tid
                    && (i64::from(mate.template_len.map_or(0, |len| len.get())).abs() / 2
                        < end - pos
                        || (mate.pos.map_or(-1, |pos| pos.as_i64()) - pos).abs() < end - pos);
                let name = if overlaps {
                    Some(record.qname()?)
                } else {
                    None
                };
                for op in record.cigar() {
                    let end = pos + i64::from(op.len());
                    match op.op_code() {
                        0 | 2 | 7 | 8 => {
                            counter.add(name, pos, end, (start, stop));
                            pos = end;
                        }
                        3 => pos = end,
                        _ => (),
                    }
                }
            }
        }
    }
    Ok(())
}

pub(super) fn collect(
    processor: &OnlyDepthProcessor,
    region: Region,
) -> Result<Vec<RangePositions>> {
    let (tid, start, stop) = region;
    if start == stop {
        return Ok(Vec::new());
    }
    let mut reader = bam::IndexedReader::from_path(&processor.reads)?;
    if let Some(reference) = &processor.ref_fasta {
        reader.set_reference(reference)?;
    }
    let name = std::str::from_utf8(reader.header().tid2name(tid))?.to_owned();
    let mut counter = Counter::default();
    counter.reset((stop - start) as usize);
    load(
        &mut Reader::Htslib(reader, bam::Record::new()),
        &mut counter,
        processor,
        region,
    )?;
    counter.merge_mates();
    let mut rows = Vec::new();
    counter.emit(processor, start, stop, |begin, end, depth| {
        let mut row = RangePositions::new(name.as_str().into(), begin + processor.coord_base);
        row.end = end + processor.coord_base;
        row.depth = depth;
        rows.push(row);
        Ok(())
    })?;
    Ok(rows)
}

fn append_row(
    output: &mut Output<'_>,
    name: &[u8],
    start: u32,
    stop: u32,
    depth: u32,
    bed: bool,
) -> Result<()> {
    output.text.extend_from_slice(name);
    for value in [start, stop] {
        output.text.push(b'\t');
        stream::append_unsigned(output.text, value);
    }
    if bed {
        output.text.push(b'\t');
    }
    output.text.push(b'\t');
    stream::append_unsigned(output.text, depth);
    output.text.push(b'\n');
    output.flush()
}
