//! # ParGranges
//!
//! Iterates over chunked genomic regions in parallel.
use anyhow::Result;
use bio::io::bed;
use crossbeam::channel::{unbounded, Receiver};
use log::*;
use num_cpus;
use rayon::prelude::*;
use rust_htslib::bam::{HeaderView, IndexedReader, Read};
use rust_lapper::{Interval, Lapper};
use serde::Serialize;
use std::{path::PathBuf, thread};

/// RegionProcessor defines the methods that must be implemented to process a region
pub trait RegionProcessor {
    /// A vector of P make up the output of [`process_region`] and
    /// are values associated with each position.
    ///
    /// [`process_region`]: #method.process_region
    type P: 'static + Send + Sync + Serialize;

    /// A function that takes the tid, start, and stop and returns something serializable.
    /// Note, a common use of this function will be a `fetch` -> `pileup`. The pileup must
    /// be bounds checked.
    fn process_region(&self, tid: u32, start: u64, stop: u64) -> Vec<Self::P>;
}

/// ParGranges holds all the information and configuration needed to launch the
/// [`ParGranges::process`].
///
/// [`ParGranges::process`]: #method.process
#[derive(Debug)]
pub struct ParGranges<R: 'static + RegionProcessor + Send + Sync> {
    /// Path to an indexed BAM / CRAM file
    reads: PathBuf,
    /// Optional reference file for CRAM
    ref_fasta: Option<PathBuf>,
    /// Optional path to a BED file to restrict the regions iterated over
    regions_bed: Option<PathBuf>,
    /// Number of threads this is allowed to use, uses all if None
    threads: usize,
    /// The ideal number of basepairs each worker will receive. Total bp in memory at one time = `threads` * `chunksize`
    chunksize: usize,
    /// The rayon threadpool to operate in
    pool: rayon::ThreadPool,
    /// The implementation of [RegionProcessor] that will be used to process regions
    processor: R,
}

impl<R: RegionProcessor + Send + Sync> ParGranges<R> {
    /// Create a ParIO object
    ///
    /// # Arguments
    ///
    /// * `reads`- path to an indexed BAM/CRAM
    /// * `ref_fasta`- path to an indexed reference file for CRAM
    /// * `regions_bed`- Optional BED file path restricting the regions to be examined
    /// * `threads`- Optional threads to restrict the number of threads this process will use, defaults to all
    /// * `chunksize`- optional argument to change the default chunksize of 1_000_000. `chunksize` determines the number of bases
    ///                each worker will get to work on at one time.
    /// * `processor`- Something that implements [`RegionProcessor`](RegionProcessor)
    pub fn new(
        reads: PathBuf,
        ref_fasta: Option<PathBuf>,
        regions_bed: Option<PathBuf>,
        threads: Option<usize>,
        chunksize: Option<usize>,
        processor: R,
    ) -> Self {
        let threads = if let Some(threads) = threads {
            threads
        } else {
            num_cpus::get()
        };

        // Keep two around for main thread and thread running the pool
        let threads = std::cmp::max(threads.checked_sub(2).unwrap_or(0), 1);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();

        let chunksize = if let Some(chunksize) = chunksize {
            chunksize
        } else {
            1_000_000
        };
        info!("Using {} worker threads.", threads);
        info!("Using chunksize of {}.", chunksize);
        Self {
            reads,
            ref_fasta,
            regions_bed,
            threads,
            chunksize,
            pool,
            processor,
        }
    }

    /// Process each region.
    ///
    /// This method splits the sequences in the BAM/CRAM header into `chunksize` * `self.threads` regions (aka 'super chunks').
    /// It then queries that 'super chunk' against the intervals (either the BED file, or the whole genome broken up into `chunksize`
    /// regions). The results of that query are then processed by a pool of workers that apply `process_region` to reach interval to
    /// do perbase analysis on. The collected result for each region is then sent back over the returned `Receiver<R::P>` channel
    /// for the caller to use. The results will be returned in order according to the order of the intervals used to drive this method.
    ///
    /// While one 'super chunk' is being worked on by all workers, the last 'super chunks' results are being printed to either to
    /// a file or to STDOUT, in order.
    ///
    /// Note, a common use case of this will be to fetch a region and do a pileup. The bounds of bases being looked at should still be
    /// checked since a fetch will pull all reads that overlap the region in question.
    pub fn process(self) -> Result<Receiver<R::P>> {
        // let mut writer = self.get_writer()?;

        let (snd, rxv) = unbounded();
        thread::spawn(move || {
            self.pool.install(|| {
                info!("Reading from {:?}", self.reads);
                let mut reader = IndexedReader::from_path(&self.reads).expect("Indexed BAM/CRAM");
                // If passed add ref_fasta
                if let Some(ref_fasta) = &self.ref_fasta {
                    reader.set_reference(ref_fasta).expect("Set ref");
                }
                // Get a copy of the header
                let header = reader.header().to_owned();

                let intervals = if let Some(regions_bed) = &self.regions_bed {
                    Self::bed_to_intervals(&header, regions_bed).expect("Parsed BED to intervals")
                } else {
                    Self::header_to_intervals(&header, self.chunksize)
                        .expect("Parsed BAM/CRAM header to intervals")
                };

                // The number positions to try to process in one batch
                let serial_step_size = self.chunksize * self.threads; // aka superchunk
                for (tid, intervals) in intervals.into_iter().enumerate() {
                    let tid: u32 = tid as u32;
                    let tid_end = header.target_len(tid).unwrap() as u64;
                    // Result holds the processed positions to be sent to writer
                    let mut result = vec![];
                    for chunk_start in (0..tid_end).step_by(serial_step_size) {
                        let chunk_end =
                            std::cmp::min(chunk_start + serial_step_size as u64, tid_end);
                        info!("Batch Processing {}:{}-{}", tid, chunk_start, chunk_end);
                        let (r, _) = rayon::join(
                            || {
                                // Must be a vec so that par_iter works and results stay in order
                                let ivs: Vec<Interval<u64, ()>> = intervals
                                    .find(chunk_start, chunk_end)
                                    // Truncate intervals that extend forward or backward of chunk in question
                                    .map(|iv| Interval {
                                        start: std::cmp::max(iv.start, chunk_start),
                                        stop: std::cmp::min(iv.stop, chunk_end),
                                        val: (),
                                    })
                                    .collect();
                                ivs.into_par_iter()
                                    .flat_map(|iv| {
                                        info!("Processing {}:{}-{}", tid, iv.start, iv.stop);
                                        self.processor.process_region(tid, iv.start, iv.stop)
                                    })
                                    .collect()
                            },
                            || {
                                result.into_iter().for_each(|p| {
                                    snd.send(p).expect("Sent a serializable to writer")
                                })
                            },
                        );
                        result = r;
                    }
                    // Send final set of results
                    result
                        .into_iter()
                        .for_each(|p| snd.send(p).expect("Sent a serializable to writer"));
                }
            });
        });
        // rxv.into_iter()
        //     .for_each(|pos| writer.serialize(pos).unwrap());
        // writer.flush()?;
        Ok(rxv)
    }

    // Convert the header into intervals of equally sized chunks. The last interval may be short.
    fn header_to_intervals(header: &HeaderView, chunksize: usize) -> Result<Vec<Lapper<u64, ()>>> {
        let mut intervals = vec![vec![]; header.target_count() as usize];
        for tid in 0..(header.target_count()) {
            let tid_len = header.target_len(tid).unwrap();
            for start in (0..tid_len).step_by(chunksize) {
                let stop = std::cmp::min(start + chunksize as u64, tid_len);
                intervals[tid as usize].push(Interval {
                    start: start,
                    stop: stop,
                    val: (),
                });
            }
        }
        Ok(intervals.into_iter().map(|ivs| Lapper::new(ivs)).collect())
    }

    /// Read a bed file into a vector of lappers with the index representing the TID
    // TODO add a proper error message
    fn bed_to_intervals(header: &HeaderView, bed_file: &PathBuf) -> Result<Vec<Lapper<u64, ()>>> {
        let mut bed_reader = bed::Reader::from_file(bed_file)?;
        let mut intervals = vec![vec![]; header.target_count() as usize];
        for record in bed_reader.records() {
            let record = record?;
            let tid = header
                .tid(record.chrom().as_bytes())
                .expect("Chromosome not found in BAM/CRAM header");
            intervals[tid as usize].push(Interval {
                start: record.start(),
                stop: record.end(),
                val: (),
            });
        }

        Ok(intervals
            .into_iter()
            .map(|ivs| {
                let mut lapper = Lapper::new(ivs);
                lapper.merge_overlaps();
                lapper
            })
            .collect())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use proptest::prelude::*;
    use rust_htslib::bam;
    use bio::io::bed;
    use rust_lapper::{Interval};
    // The purpose of these tests is to demonstrate that positions are covered once under a variety of circumstances

    // Variables
    // - chunksize
    // - number of cpus
    // - number of intervals

    prop_compose! {
        fn arb_iv_start(max_iv: usize)(start in 0..max_iv/2) -> usize { start }
    }
    prop_compose! {
        fn arb_iv_stop(max_iv: usize)(stop in max_iv/2..max_iv) -> usize { stop }
    }
    prop_compose! {
        // Create an arbitrary interval where the min size == max_iv / 2
        // TODO: make minsize smaller
        fn arb_iv(max_iv: usize)(start in arb_iv_start(max_iv), stop in arb_iv_stop(max_iv)) -> Interval<usize, ()> {
            Interval {start, stop, val: ()}
        }
    }
    // Create an arbitrary number of intervals along with the expected number of positions they cover
    fn arb_ivs(max_iv: usize, max_ivs: usize) -> impl Strategy<Value = (Vec<Interval<usize, ()>>, usize)> {
        prop::collection::vec(arb_iv(max_iv), 0..max_ivs).prop_map(|vec| {
            let mut expected = 0;
            for iv in vec.iter() {
                expected += iv.stop - iv.start;
            }
            (vec, expected)
        })
    }
    // Create arbitrary number of contigs with arbitrary intervals each
    fn arb_chrs(max_chr: usize, max_iv: usize, max_ivs: usize) -> impl Strategy<Value = Vec<(Vec<Interval<usize, ()>>, usize)>> {
        prop::collection::vec(arb_ivs(max_iv, max_ivs), 0..max_chr)
    }
    // An empty BAM with correct header
    // A BED file with the randomly generated intervals (with expected number of positions)
    // proptest generate random chunksize, cpus
    proptest! {
        #[test]
        // add random chunksize and random cpus
        fn interval_set(chromosmes in arb_chrs(50, 5_000_000, 10_000), chunksize in any::<usize>(), cpus in any::<usize>()) {
            // Write a bam file with info
            // Write a bed file with info
            // Create the processor with a dumb impl of procesing that just returns posiotns with no counting
            // Validate that for each chr we get the expected number of bases


            println!("{:?} - {:?}", ivs, expected );
            assert!(false);
        }
    }
}
