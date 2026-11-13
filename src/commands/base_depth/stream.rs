//! Count and format in workers; write bounded byte blocks in region order.
use super::{BaseProcessor, DefaultReadFilter, PileupPosition, Position};
use anyhow::{Context, Result, anyhow};
use crossbeam_channel::{Sender, bounded};
use perbase_lib::{par_granges::Region, position::seqair::SeqairWorker, utils};
use rust_htslib::bam::Read;
use std::{io::Write, panic::AssertUnwindSafe, thread};

const BLOCK: usize = 65536;
const COLUMNS: &[u8] = b"DEPTH\tA\tC\tG\tT\tN\tR\tY\tS\tW\tK\tM\tINS\tDEL\tREF_SKIP\tFAIL\tCOUNT_OF_MATE_RESOLUTIONS\tNEAR_MAX_DEPTH\n";

pub(super) fn run(
    processor: &BaseProcessor<DefaultReadFilter>,
    regions: &[Region],
    workers: usize,
    chunksize: u32,
    modifier: f64,
    output: &mut dyn Write,
) -> Result<()> {
    anyhow::ensure!(
        modifier.is_finite() && modifier >= 0.0,
        "invalid channel-size-modifier"
    );
    if regions.is_empty() {
        return Ok(());
    }
    let names =
        utils::with_indexed_reader(&processor.reads, processor.ref_fasta.as_deref(), |reader| {
            reader
                .header()
                .target_names()
                .iter()
                .map(|name| String::from_utf8(name.to_vec()))
                .collect::<std::result::Result<Vec<_>, _>>()
        })?;
    // Match Goose's scheduling window and output budget. Region boundaries
    // never change when the budget changes, even with only one block of space.
    let longest = names.iter().map(String::len).max().unwrap_or(0);
    let estimate = f64::from(chunksize) * (longest + 64) as f64;
    let budget = modifier * 1_000_000_000.0 * workers as f64;
    let window = ((budget / estimate).min(workers as f64) as usize)
        .max(1)
        .min(regions.len());
    let blocks = ((budget / window as f64).min(estimate) as usize / BLOCK).max(1);
    let counter = processor
        .seqair
        .as_ref()
        .expect("Seqair streaming processor");
    let states = (0..window)
        .map(|_| counter.worker())
        .collect::<Result<Vec<_>>>()?;
    let (jobs_tx, jobs_rx) = bounded::<usize>(window);
    let (senders, receivers): (Vec<_>, Vec<_>) = (0..window)
        .map(|_| bounded::<Result<Vec<u8>>>(blocks))
        .unzip();

    thread::scope(|scope| {
        for mut worker in states {
            let jobs = jobs_rx.clone();
            let senders = &senders;
            let names = &names;
            scope.spawn(move || {
                let mut text = Vec::with_capacity(BLOCK + 256);
                for job in jobs {
                    let sender = &senders[job % window];
                    // A failed producer must close its region, not leave the
                    // ordered writer waiting forever for its next block.
                    let result = std::panic::catch_unwind(AssertUnwindSafe(|| {
                        write_region(
                            processor,
                            &mut worker,
                            regions[job],
                            names,
                            sender,
                            &mut text,
                        )
                    }))
                    .unwrap_or_else(|_| Err(anyhow!("counting worker panicked")));
                    let failed = result.is_err();
                    if sender.send(result.map(|()| Vec::new())).is_err() || failed {
                        break;
                    }
                }
            });
        }
        let result = (|| -> Result<()> {
            for job in 0..window {
                jobs_tx.send(job)?;
            }
            let mut started = false;
            for job in 0..regions.len() {
                loop {
                    let bytes = receivers[job % window]
                        .recv()
                        .context("counting worker disconnected")??;
                    if bytes.is_empty() {
                        break;
                    }
                    if !started {
                        output.write_all(b"REF\tPOS\t")?;
                        if processor.ref_fasta.is_some() {
                            output.write_all(b"REF_BASE\t")?;
                        }
                        output.write_all(COLUMNS)?;
                        started = true;
                    }
                    output.write_all(&bytes)?;
                }
                if job + window < regions.len() {
                    jobs_tx.send(job + window)?;
                }
            }
            Ok(())
        })();
        // Also run on output failure: dropping receivers unblocks senders;
        // dropping the job queue wakes idle workers. Scope joins every worker.
        drop(receivers);
        drop(jobs_tx);
        result
    })
}

fn send_text(sender: &Sender<Result<Vec<u8>>>, text: &mut Vec<u8>) -> Result<()> {
    for block in text.chunks(BLOCK) {
        sender
            .send(Ok(block.to_vec()))
            .map_err(|_| anyhow!("output consumer stopped"))?;
    }
    text.clear();
    Ok(())
}

fn write_region(
    processor: &BaseProcessor<DefaultReadFilter>,
    worker: &mut SeqairWorker<'_>,
    (tid, start, stop): Region,
    names: &[String],
    sender: &Sender<Result<Vec<u8>>>,
    text: &mut Vec<u8>,
) -> Result<()> {
    text.clear();
    let name = &names[tid as usize];
    let mut quoted = Vec::new();
    if name
        .bytes()
        .any(|byte| matches!(byte, b'\t' | b'\r' | b'\n' | b'"'))
    {
        quoted.push(b'"');
        for byte in name.bytes() {
            if byte == b'"' {
                quoted.push(byte);
            }
            quoted.push(byte);
        }
        quoted.push(b'"');
    } else {
        quoted.extend_from_slice(name.as_bytes());
    }
    // Fetch lazily: an empty region without --keep-zeros needs no reference.
    let mut reference = None;
    let mut emit = |mut row: PileupPosition| -> Result<()> {
        if let Some(buffer) = &processor.ref_buffer {
            if reference.is_none() {
                reference = Some(buffer.seq(name)?);
            }
            row.ref_base = Some(char::from(
                *reference
                    .as_ref()
                    .unwrap()
                    .get(row.pos as usize)
                    .context("Input SAM does not match reference")?,
            ));
        }
        row.pos += processor.coord_base;
        append_row(text, &quoted, &row);
        if text.len() >= BLOCK {
            send_text(sender, text)?;
        }
        Ok(())
    };
    let mut next = start;
    worker.for_each_position(tid, start, stop, |row| {
        if processor.keep_zeros {
            while next < row.pos {
                emit(PileupPosition::new(name.as_str().into(), next))?;
                next += 1;
            }
        }
        next = row.pos + 1;
        emit(row)
    })?;
    if processor.keep_zeros {
        while next < stop {
            emit(PileupPosition::new(name.as_str().into(), next))?;
            next += 1;
        }
    }
    send_text(sender, text)
}

fn append_row(text: &mut Vec<u8>, name: &[u8], row: &PileupPosition) {
    let mut number = itoa::Buffer::new();
    text.extend_from_slice(name);
    text.push(b'\t');
    text.extend_from_slice(number.format(row.pos).as_bytes());
    if let Some(base) = row.ref_base {
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
        text.extend_from_slice(number.format(value).as_bytes());
    }
    text.extend_from_slice(if row.near_max_depth {
        b"\ttrue\n"
    } else {
        b"\tfalse\n"
    });
}
