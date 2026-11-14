//! Bounded, ordered output shared by the depth commands.
use anyhow::{Context, Result, anyhow};
use crossbeam_channel::{Receiver, Sender, bounded};
use perbase_lib::par_granges::Region;
use std::{io::Write, panic::AssertUnwindSafe, thread};

const BLOCK: usize = 65536;

pub(super) fn quote_name(name: &str) -> Vec<u8> {
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
    quoted
}

/// Most depth counters fit in two digits. Larger values use itoa's stack buffer.
#[inline]
pub(super) fn append_unsigned(text: &mut Vec<u8>, value: u32) {
    if value < 10 {
        text.push(b'0' + value as u8);
    } else if value < 100 {
        text.push(b'0' + (value / 10) as u8);
        text.push(b'0' + (value % 10) as u8);
    } else {
        text.extend_from_slice(itoa::Buffer::new().format(value).as_bytes());
    }
}

pub(super) struct Config {
    pub workers: usize,
    pub chunksize: u32,
    pub modifier: f64,
    pub estimated_row_bytes: usize,
}

/// A worker's reusable output buffer. Full blocks leave immediately; a partial
/// tail stays here until another row fills it or the region finishes.
pub(super) struct Output<'a> {
    pub text: &'a mut Vec<u8>,
    sender: &'a Sender<Result<Vec<u8>>>,
    recycled: &'a Receiver<Vec<u8>>,
}

impl Output<'_> {
    pub fn flush(&mut self) -> Result<()> {
        self.send(false)
    }

    fn send(&mut self, finish: bool) -> Result<()> {
        let length = if finish {
            self.text.len()
        } else {
            self.text.len() / BLOCK * BLOCK
        };
        if length == 0 {
            return Ok(());
        }
        for block in self.text[..length].chunks(BLOCK) {
            let mut bytes = if block.len() == BLOCK {
                self.recycled
                    .try_recv()
                    .unwrap_or_else(|_| Vec::with_capacity(BLOCK))
            } else {
                Vec::with_capacity(block.len())
            };
            bytes.extend_from_slice(block);
            self.sender
                .send(Ok(bytes))
                .map_err(|_| anyhow!("output consumer stopped"))?;
        }
        self.text.copy_within(length.., 0);
        self.text.truncate(self.text.len() - length);
        Ok(())
    }
}

/// Count and format in workers, then write their blocks in input-region order.
/// A failed worker or output write disconnects the queues before joining the
/// threads, so producers cannot remain blocked behind a stopped consumer.
pub(super) fn run<W: Send>(
    regions: &[Region],
    config: Config,
    header: &[u8],
    output: &mut dyn Write,
    make_worker: impl Fn() -> Result<W>,
    process: impl Fn(&mut W, Region, &mut Output<'_>) -> Result<()> + Sync,
) -> Result<()> {
    anyhow::ensure!(
        config.modifier.is_finite() && config.modifier >= 0.0,
        "invalid channel-size-modifier"
    );
    if regions.is_empty() {
        return Ok(());
    }
    let estimate = f64::from(config.chunksize) * config.estimated_row_bytes as f64;
    let budget = config.modifier * 1_000_000_000.0 * config.workers as f64;
    // Queues hold blocks, not whole regions. A small output budget must not
    // serialize counting: give every worker at least one queued block.
    let window = config.workers.max(1).min(regions.len());
    let blocks = ((budget / window as f64).min(estimate) as usize / BLOCK).max(1);
    let states = (0..window)
        .map(|_| make_worker())
        .collect::<Result<Vec<_>>>()?;
    let (jobs_tx, jobs_rx) = bounded::<usize>(window);
    let (senders, receivers): (Vec<_>, Vec<_>) = (0..window)
        .map(|_| bounded::<Result<Vec<u8>>>(blocks))
        .unzip();
    // Retain just one spare allocation per slot, not entire completed regions.
    let (returns, recycled): (Vec<_>, Vec<_>) = (0..window).map(|_| bounded::<Vec<u8>>(1)).unzip();

    thread::scope(|scope| {
        for mut worker in states {
            let jobs = jobs_rx.clone();
            let senders = &senders;
            let recycled = &recycled;
            let process = &process;
            scope.spawn(move || {
                let mut text = Vec::with_capacity(BLOCK + 256);
                for job in jobs {
                    text.clear();
                    let sender = &senders[job % window];
                    let mut output = Output {
                        text: &mut text,
                        sender,
                        recycled: &recycled[job % window],
                    };
                    let result = std::panic::catch_unwind(AssertUnwindSafe(|| {
                        process(&mut worker, regions[job], &mut output)?;
                        output.send(true)
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
                    let mut bytes = receivers[job % window]
                        .recv()
                        .context("counting worker disconnected")??;
                    if bytes.is_empty() {
                        break;
                    }
                    if !started {
                        output.write_all(header)?;
                        started = true;
                    }
                    output.write_all(&bytes)?;
                    if bytes.capacity() == BLOCK {
                        bytes.clear();
                        let _ = returns[job % window].try_send(bytes);
                    }
                }
                if job + window < regions.len() {
                    jobs_tx.send(job + window)?;
                }
            }
            Ok(())
        })();
        drop(receivers);
        drop(jobs_tx);
        result
    })
}
