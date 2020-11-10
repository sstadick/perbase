//! # ParGranges
//!
//! Iterates over chunked genomic regions in parallel.
use anyhow::Result;
use bio::io::bed;
use crossbeam::channel::{unbounded, Receiver};
use log::*;
use num_cpus;
use rayon::prelude::*;
use rust_htslib::{
    bam::{HeaderView, IndexedReader, Read},
    bcf::{Reader, Read as bcfRead}
};
use rust_lapper::{Interval, Lapper};
use serde::Serialize;
use std::{path::PathBuf, thread, convert::TryInto};

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
    /// Optional path to a BCF/VCF file to restrict the regions iterated over
    regions_bcf: Option<PathBuf>,
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
    /// * `regions_bcf`- Optional BCF/VCF file path restricting the regions to be examined
    /// * `threads`- Optional threads to restrict the number of threads this process will use, defaults to all
    /// * `chunksize`- optional argument to change the default chunksize of 1_000_000. `chunksize` determines the number of bases
    ///                each worker will get to work on at one time.
    /// * `processor`- Something that implements [`RegionProcessor`](RegionProcessor)
    pub fn new(
        reads: PathBuf,
        ref_fasta: Option<PathBuf>,
        regions_bed: Option<PathBuf>,
        regions_bcf: Option<PathBuf>,
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
            regions_bcf,
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

                // Work out if we are restricted to a subset of sites
                let bed_intervals = if let Some(regions_bed) = &self.regions_bed {
                    Some(Self::bed_to_intervals(&header, regions_bed).expect("Parsed BED to intervals"))
                } else { None };
                let bcf_intervals = if let Some(regions_bcf) = &self.regions_bcf {
                    Some(Self::bcf_to_intervals(&header, regions_bcf).expect("Parsed BCF/VCF to intervals"))
                } else { None };
                let restricted_ivs = match (bed_intervals, bcf_intervals) {
                    (Some(bed_ivs), Some(bcf_ivs)) => Some(Self::merge_intervals(bed_ivs, bcf_ivs)),
                    (Some(bed_ivs), None) => Some(bed_ivs),
                    (None, Some(bcf_ivs)) => Some(bcf_ivs),
                    (None, None) => None
                };

                let intervals = if let Some(restricted) = restricted_ivs {
                    restricted
                } else {
                    Self::header_to_intervals(&header, self.chunksize)
                        .expect("Parsed BAM/CRAM header to intervals")
                };

                // The number positions to try to process in one batch
                let serial_step_size = self
                    .chunksize
                    .checked_mul(self.threads)
                    .unwrap_or(usize::MAX); // aka superchunk
                for (tid, intervals) in intervals.into_iter().enumerate() {
                    let tid: u32 = tid as u32;
                    let tid_end = header.target_len(tid).unwrap() as u64;
                    // Result holds the processed positions to be sent to writer
                    let mut result = vec![];
                    for chunk_start in (0..tid_end).step_by(serial_step_size) {
                        let tid_name = std::str::from_utf8(header.tid2name(tid)).unwrap();
                        let chunk_end =
                            std::cmp::min(chunk_start + serial_step_size as u64, tid_end);
                        info!("Batch Processing {}:{}-{}", tid_name, chunk_start, chunk_end);
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
                                        info!("Processing {}:{}-{}", tid_name, iv.start, iv.stop);
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

    /// Read a BCF/VCF file into a vector of lappers with index representing the TID
    fn bcf_to_intervals(header: &HeaderView, bcf_file: &PathBuf) -> Result<Vec<Lapper<u64, ()>>> {
        let mut bcf_reader = Reader::from_path(bcf_file).expect("Error opening BCF/VCF file.");
        let bcf_header_reader = Reader::from_path(bcf_file).expect("Error opening BCF/VCF file.");
        let bcf_header = bcf_header_reader.header();
        let mut intervals = vec![vec![]; header.target_count() as usize];
        // TODO: validate the headers against eachother
        for record in bcf_reader.records() {
            let record = record?;
            let record_rid = bcf_header.rid2name(record.rid().unwrap()).unwrap();
            let tid = header.tid(record_rid).expect("Chromosome not found in BAM/CRAM header");
            let pos: u64 = record.pos().try_into().expect("Got a negative value for pos");
            intervals[tid as usize].push(Interval {
                start: pos, stop: pos + 1, val: ()
            });
        }

        Ok(intervals.into_iter().map(|ivs| {
            let mut lapper = Lapper::new(ivs);
            lapper.merge_overlaps();
            lapper
        }).collect())
    }

    /// Merge two sets of restriction intervals together
    fn merge_intervals(a_ivs: Vec<Lapper<u64, ()>>, b_ivs: Vec<Lapper<u64, ()>>) -> Vec<Lapper<u64, ()>> {
        let mut intervals = vec![vec![]; a_ivs.len()];
        for (i, (a_lapper, b_lapper)) in a_ivs.into_iter().zip(b_ivs.into_iter()).enumerate() {
            intervals[i] = a_lapper.into_iter().chain(b_lapper.into_iter()).collect();
        }
        intervals.into_iter().map(|ivs| {
            let mut lapper = Lapper::new(ivs);
            lapper.merge_overlaps();
            lapper
        }).collect()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use bio::io::bed;
    use num_cpus;
    use proptest::prelude::*;
    use rust_htslib::{bam, bcf};
    use rust_lapper::{Interval, Lapper};
    use std::collections::{HashMap, HashSet};
    use tempfile::tempdir;
    // The purpose of these tests is to demonstrate that positions are covered once under a variety of circumstances

    prop_compose! {
        fn arb_iv_start(max_iv: u64)(start in 0..max_iv/2) -> u64 { start }
    }
    prop_compose! {
        fn arb_iv_size(max_iv: u64)(size in 0..max_iv/2) -> u64 { size }
    }
    prop_compose! {
        // Create an arbitrary interval where the min size == max_iv / 2
        fn arb_iv(max_iv: u64)(start in arb_iv_start(max_iv), size in arb_iv_size(max_iv)) -> Interval<u64, ()> {
            Interval {start, stop: start + size, val: ()}
        }
    }
    // Create an arbitrary number of intervals along with the expected number of positions they cover
    fn arb_ivs(
        max_iv: u64,    // max iv size
        max_ivs: usize, // max number of intervals
    ) -> impl Strategy<Value = (Vec<Interval<u64, ()>>, u64, u64)> {
        prop::collection::vec(arb_iv(max_iv), 0..max_ivs).prop_map(|vec| {
            let mut furthest_right = 0;
            let lapper = Lapper::new(vec.clone());
            let expected = lapper.cov();
            for iv in vec.iter() {
                if iv.stop > furthest_right {
                    furthest_right = iv.stop;
                }
            }
            (vec, expected, furthest_right)
        })
    }
    // Create arbitrary number of contigs with arbitrary intervals each
    fn arb_chrs(
        max_chr: usize, // number of chromosomes to use
        max_iv: u64,    // max interval size
        max_ivs: usize, // max number of intervals
    ) -> impl Strategy<Value = Vec<(Vec<Interval<u64, ()>>, u64, u64)>> {
        prop::collection::vec(arb_ivs(max_iv, max_ivs), 0..max_chr)
    }
    // An empty BAM with correct header
    // A BED file with the randomly generated intervals (with expected number of positions)
    // proptest generate random chunksize, cpus
    proptest! {
        #[test]
        // add random chunksize and random cpus
        // NB: using any larger numbers for this tends to blow up the test runtime
        fn interval_set(chromosomes in arb_chrs(4, 10_000, 1_000), chunksize in any::<usize>(), cpus in 0..num_cpus::get(), use_bed in any::<bool>(), use_vcf in any::<bool>()) {
            let tempdir = tempdir().unwrap();
            let bam_path = tempdir.path().join("test.bam");
            let bed_path = tempdir.path().join("test.bed");
            let vcf_path = tempdir.path().join("test.vcf");

            // Build a BAM
            let mut header = bam::header::Header::new();
            for (i,chr) in chromosomes.iter().enumerate() {
                let mut chr_rec = bam::header::HeaderRecord::new(b"SQ");
                chr_rec.push_tag(b"SN", &i.to_string());
                chr_rec.push_tag(b"LN", &chr.2.to_string()); // set len as max observed
                header.push_record(&chr_rec);
            }
            let writer = bam::Writer::from_path(&bam_path, &header, bam::Format::BAM).expect("Opened test.bam for writing");
            drop(writer); // force flush the writer so the header info is written
            bam::index::build(&bam_path, None, bam::index::Type::BAI, 1).unwrap();

            // Build a bed
            let mut writer = bed::Writer::to_file(&bed_path).expect("Opened test.bed for writing");
            for (i, chr) in chromosomes.iter().enumerate() {
                for iv in chr.0.iter() {
                    let mut record = bed::Record::new();
                    record.set_start(iv.start);
                    record.set_end(iv.stop);
                    record.set_chrom(&i.to_string());
                    record.set_score(&0.to_string());
                    writer.write(&record).expect("Wrote to test.bed");
                }
            }
            drop(writer); // force flush

            // Build a VCF file
            let mut vcf_truth = HashMap::new();
            let mut header = bcf::header::Header::new();
            for (i,chr) in chromosomes.iter().enumerate() {
                header.push_record(format!("##contig=<ID={},length={}>", &i.to_string(), &chr.2.to_string()).as_bytes());
            }
            let mut writer = bcf::Writer::from_path(&vcf_path, &header, true, bcf::Format::VCF).expect("Failed to open test.vcf for writing");
            let mut record = writer.empty_record();
            for (i, chr) in chromosomes.iter().enumerate() {
                record.set_rid(Some(i as u32));
                let counter = vcf_truth.entry(i).or_insert(0);
                let mut seen = HashSet::new();
                for iv in chr.0.iter() {
                    if !seen.contains(&iv.start) {
                        *counter += 1;
                        seen.insert(iv.start);
                    }
                    record.set_pos(iv.start as i64);
                    writer.write(&record).expect("Failed to write to test.vcf")
                }
            }

            drop(writer); // force flush
            // Create the processor with a dumb impl of processing that just returns positions with no counting
            let test_processor = TestProcessor {};
            let par_granges_runner = ParGranges::new(
                bam_path,
                None,
                if use_bed { Some(bed_path) } else { None }, // do one with regions
                if use_vcf { Some(vcf_path) } else { None }, // do one with vcf regions
                Some(cpus),
                Some(chunksize),
                test_processor
            );
            let receiver = par_granges_runner.process().expect("Launch ParGranges Process");
            let mut chrom_counts = HashMap::new();
            receiver.into_iter().for_each(|p: PileupPosition| {
                let positions = chrom_counts.entry(p.ref_seq.parse::<usize>().expect("parsed chr")).or_insert(0u64);
                *positions += 1
            });

            // Validate that for each chr we get the expected number of bases
            for (chrom, positions) in chrom_counts.iter() {
                if use_bed  && !use_vcf {
                    // if this was with bed, should be equal to .1
                    prop_assert_eq!(chromosomes[*chrom].1, *positions, "chr: {}, expected: {}, found: {}", chrom, chromosomes[*chrom].1, positions);
                } else if use_bed && use_vcf {
                    // if this was with bed, should be equal to .1, bed restrictions and vcf restrctions should overlap
                    prop_assert_eq!(chromosomes[*chrom].1, *positions, "chr: {}, expected: {}, found: {}", chrom, chromosomes[*chrom].1, positions);
                } else if use_vcf && !use_bed {
                    // total positions should be equal to the number of records for that chr in the vcf
                    prop_assert_eq!(vcf_truth.get(chrom).unwrap(), positions, "chr: {}, expected: {}, found: {}", chrom, chromosomes[*chrom].1, positions);
                } else {
                    // if this was bam only, should be equal to rightmost postion
                    prop_assert_eq!(chromosomes[*chrom].2, *positions, "chr: {}, expected: {}, found: {}", chrom, chromosomes[*chrom].2, positions);
                }
            }

        }
    }

    use crate::position::{pileup_position::PileupPosition, Position};
    use smartstring::SmartString;
    struct TestProcessor {}
    impl RegionProcessor for TestProcessor {
        type P = PileupPosition;

        fn process_region(&self, tid: u32, start: u64, stop: u64) -> Vec<Self::P> {
            let mut results = vec![];
            for i in start..stop {
                let chr = SmartString::from(&tid.to_string());
                let pos = PileupPosition::new(chr, i as usize);
                results.push(pos);
            }
            results
        }
    }
}
