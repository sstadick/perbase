//! Count and format in workers; write bounded byte blocks in region order.
use super::{BaseProcessor, DefaultReadFilter, MateCounter, PileupPosition};
use crate::commands::stream::{self, Config, Output, append_unsigned};
use anyhow::{Context, Result};
#[cfg(feature = "seqair")]
use perbase_lib::position::seqair::SeqairWorker;
use perbase_lib::{par_granges::Region, utils};
use rust_htslib::bam::Read;
use std::io::Write;

#[cfg(not(feature = "seqair"))]
type Worker = ();
#[cfg(feature = "seqair")]
type Worker<'a> = Option<SeqairWorker<'a>>;

const COLUMNS: &[u8] = b"DEPTH\tA\tC\tG\tT\tN\tR\tY\tS\tW\tK\tM\tINS\tDEL\tREF_SKIP\tFAIL\tCOUNT_OF_MATE_RESOLUTIONS\tNEAR_MAX_DEPTH\n";

pub(super) fn run(
    processor: &BaseProcessor<DefaultReadFilter>,
    regions: &[Region],
    workers: usize,
    chunksize: u32,
    modifier: f64,
    output: &mut dyn Write,
) -> Result<()> {
    let names =
        utils::with_indexed_reader(&processor.reads, processor.ref_fasta.as_deref(), |reader| {
            reader
                .header()
                .target_names()
                .iter()
                .map(|name| String::from_utf8(name.to_vec()))
                .collect::<std::result::Result<Vec<_>, _>>()
        })?;
    let longest = names.iter().map(String::len).max().unwrap_or(0);
    let mut header = b"REF\tPOS\t".to_vec();
    if processor.ref_fasta.is_some() {
        header.extend_from_slice(b"REF_BASE\t");
    }
    header.extend_from_slice(COLUMNS);
    stream::run(
        regions,
        Config {
            workers,
            chunksize,
            modifier,
            estimated_row_bytes: longest + 64,
        },
        &header,
        output,
        || {
            #[cfg(feature = "seqair")]
            {
                processor
                    .seqair
                    .as_ref()
                    .map(|counter| counter.worker())
                    .transpose()
            }
            #[cfg(not(feature = "seqair"))]
            {
                Ok(())
            }
        },
        |worker, region, output| write_region(processor, worker, region, &names, output),
    )
}

fn write_region(
    processor: &BaseProcessor<DefaultReadFilter>,
    worker: &mut Worker,
    (tid, start, stop): Region,
    names: &[String],
    output: &mut Output<'_>,
) -> Result<()> {
    let name = &names[tid as usize];
    let quoted = stream::quote_name(name);
    // Fetch lazily: an empty region without --keep-zeros needs no reference.
    let mut reference = None;
    let mut emit = |pos: u32, row: &PileupPosition| -> Result<()> {
        let base = if let Some(buffer) = &processor.ref_buffer {
            if reference.is_none() {
                reference = Some(buffer.seq(name)?);
            }
            Some(char::from(
                *reference
                    .as_ref()
                    .unwrap()
                    .get(pos as usize)
                    .context("Input SAM does not match reference")?,
            ))
        } else {
            None
        };
        append_row(output.text, &quoted, pos + processor.coord_base, base, row);
        output.flush()
    };
    let mut next = start;
    let zero = PileupPosition::default();
    let mut count_row = |row: &PileupPosition| {
        if processor.keep_zeros {
            while next < row.pos {
                emit(next, &zero)?;
                next += 1;
            }
        }
        next = row.pos + 1;
        emit(row.pos, row)
    };
    #[cfg(feature = "seqair")]
    let used_seqair = if let Some(worker) = worker {
        worker.for_each_position_ref(tid, start, stop, &mut count_row)?;
        true
    } else {
        false
    };
    #[cfg(not(feature = "seqair"))]
    let used_seqair = {
        let _ = worker;
        false
    };
    if !used_seqair && start < stop {
        // Use the same scheduling, formatting, and output buffers for both
        // backends. Only reading, pileup generation, and counting differ.
        utils::with_indexed_reader(
            &processor.reads,
            processor.ref_fasta.as_deref(),
            |reader| -> Result<()> {
                let header = reader.header().to_owned();
                reader.fetch((tid, start, stop))?;
                let mut pileups = reader.pileup();
                pileups.set_max_depth(processor.max_depth.min(i32::MAX as u32));
                // Name groups must not accumulate across an entire genome.
                let mut mates = MateCounter::default();
                for pileup in pileups {
                    let pileup = pileup?;
                    if pileup.pos() < start || pileup.pos() >= stop {
                        continue;
                    }
                    let depth = pileup.depth();
                    let mut row = if processor.mate_fix {
                        mates.count(
                            pileup,
                            &header,
                            &processor.read_filter,
                            processor.min_base_quality_score,
                            processor.mate_fix_strategy,
                        )
                    } else {
                        PileupPosition::from_pileup(
                            pileup,
                            &header,
                            &processor.read_filter,
                            processor.min_base_quality_score,
                        )
                    };
                    row.near_max_depth = depth > processor.max_depth_warnings_cutoff;
                    count_row(&row)?;
                }
                Ok(())
            },
        )?;
    }
    if processor.keep_zeros {
        while next < stop {
            emit(next, &zero)?;
            next += 1;
        }
    }
    Ok(())
}

fn append_row(text: &mut Vec<u8>, name: &[u8], pos: u32, base: Option<char>, row: &PileupPosition) {
    text.extend_from_slice(name);
    text.push(b'\t');
    append_unsigned(text, pos);
    if let Some(base) = base {
        text.push(b'\t');
        match base {
            '"' => text.extend_from_slice(b"\"\"\"\""),
            '\t' | '\r' | '\n' => text.extend_from_slice(&[b'"', base as u8, b'"']),
            _ => text.extend_from_slice(base.encode_utf8(&mut [0; 4]).as_bytes()),
        }
    }
    for value in [
        row.depth,
        row.a,
        row.c,
        row.g,
        row.t,
        row.n,
        row.r,
        row.y,
        row.s,
        row.w,
        row.k,
        row.m,
        row.ins,
        row.del,
        row.ref_skip,
        row.fail,
        row.count_of_mate_resolutions,
    ] {
        text.push(b'\t');
        append_unsigned(text, value);
    }
    text.extend_from_slice(if row.near_max_depth {
        b"\ttrue\n"
    } else {
        b"\tfalse\n"
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rows_match_csv_at_decimal_boundaries() {
        let mut values: Vec<u32> = (0..=101).collect();
        for power in [1_000, 10_000, 100_000, 1_000_000, 1_000_000_000] {
            values.extend([power - 1, power, power + 1]);
        }
        values.push(u32::MAX);
        for value in values {
            for base in [None, Some('A'), Some('"'), Some('\t'), Some('é')] {
                let row = PileupPosition {
                    ref_seq: "chr20".into(),
                    pos: value,
                    ref_base: base,
                    depth: value,
                    a: 1,
                    c: 9,
                    g: 10,
                    t: 99,
                    n: 100,
                    r: 101,
                    y: 999,
                    s: 1_000,
                    w: 10_000,
                    k: 100_000,
                    m: 1_000_000,
                    ins: 2,
                    del: 3,
                    ref_skip: 4,
                    fail: 5,
                    count_of_mate_resolutions: 6,
                    near_max_depth: value % 2 == 0,
                };
                let mut actual = Vec::new();
                append_row(&mut actual, b"chr20", row.pos, base, &row);
                let mut csv = csv::WriterBuilder::new()
                    .delimiter(b'\t')
                    .has_headers(false)
                    .from_writer(Vec::new());
                csv.serialize(&row).unwrap();
                assert_eq!(
                    actual,
                    csv.into_inner().unwrap(),
                    "value={value}, base={base:?}"
                );
            }
        }
    }
}
