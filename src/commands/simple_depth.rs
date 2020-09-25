//! # Simple Depth
//!
//! Simple single pass over bam file to calculate depth at each position
//! as well as depth per nucleotide. Additionally counts the number of
//! insertions / deletions at each position.
use anyhow::Result;
use argh::FromArgs;
use log::*;
use perbase_lib::{par_io, utils};
use rust_htslib::{bam, bam::Read};
use serde::Serialize;
use smartstring::alias::String;
use std::{default, path::PathBuf};

/// Hold all information about a position.
#[derive(Debug, Serialize, Default)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
struct Position {
    /// Reference sequence name.
    #[serde(rename = "REF")]
    ref_seq: String,
    /// 1-based position in the sequence.
    pos: usize,
    /// Total depth at this position.
    depth: usize,
    /// Number of A bases at this position.
    a: usize,
    /// Number of C bases at this position.
    c: usize,
    /// Number of G bases at this position.
    g: usize,
    /// Number of T bases at this position.
    t: usize,
    /// Number of N bases at this position. Any unrecognized base will be counted as an N.
    n: usize,
    /// Number of insertions that start to the right of this position.
    /// Does not count toward depth.
    ins: usize,
    /// Number of deletions at this position.
    del: usize,
    /// Number of reads failing filters at this position.
    fail: usize,
    /// Number of refskips at this position. Does count toward depth.
    ref_skip: usize,
}
impl Position {
    /// Create a new position for the given ref_seq name.
    fn new(ref_seq: String, pos: usize) -> Self {
        Position {
            ref_seq,
            pos,
            ..default::Default::default()
        }
    }

    /// Convert a pileup into a `Position`.
    ///
    /// This will walk over each of the alignments and count the number each nucleotide it finds.
    /// It will also count the number of Ins/Dels/Skips that are at each position. The output of this 1-based.
    ///
    /// # Arguments
    ///
    /// * `pileup` - a [bam::pileup::Pileup] at a genomic position
    /// * `header` - a [bam::HeaderView] for the bam file being read, to get the sequence name
    /// * `filter_fn` - a function to filter out reads, returning false will cause a read to be filtered
    // TODOs:
    // Write my own pilelup engine
    // Reference base
    // Average base quality
    // Average map quality
    // Average dist to ends
    // A calculated error rate based on mismatches?
    // Optionally display read names at the position?
    fn from_pileup<F>(pileup: bam::pileup::Pileup, header: &bam::HeaderView, filter_fn: F) -> Self
    where
        F: Fn(&bam::record::Record) -> bool,
    {
        let name = std::str::from_utf8(header.tid2name(pileup.tid())).unwrap();
        // make output 1-based
        let mut pos = Self::new(String::from(name), (pileup.pos() + 1) as usize);
        pos.depth = pileup.depth() as usize;

        for alignment in pileup.alignments() {
            let record = alignment.record();
            if !filter_fn(&record) {
                pos.depth -= 1;
                pos.fail += 1;
                continue;
            }
            // NB: Order matters here, a refskip is true for both is_del and is_refskip
            // while a true del is only true for is_del
            if alignment.is_refskip() {
                pos.ref_skip += 1
            } else if alignment.is_del() {
                pos.del += 1
            } else {
                // We have an actual base!
                match (record.seq()[alignment.qpos().unwrap()] as char).to_ascii_uppercase() {
                    'A' => pos.a += 1,
                    'C' => pos.c += 1,
                    'T' => pos.t += 1,
                    'G' => pos.g += 1,
                    _ => pos.n += 1,
                }
                // Check for insertions
                match alignment.indel() {
                    bam::pileup::Indel::Ins(_len) => {
                        pos.ins += 1;
                    }
                    _ => (),
                }
            }
        }
        pos
    }
}

/// Calculate the depth at each base, per-nucleotide. Takes an indexed BAM/CRAM as <reads>.
#[derive(FromArgs)]
#[argh(subcommand, name = "simple-depth")]
pub struct SimpleDepth {
    /// input BAM/CRAM to analyze, must be indexed
    #[argh(positional)]
    reads: PathBuf,

    /// indexed reference fasta, set if using CRAM
    #[argh(option, short = 'r')]
    ref_fasta: Option<PathBuf>,

    /// a BED file containing regions of interest. If specified, only bases from the given regions will be reported on
    #[argh(option, short = 'b')]
    bed_file: Option<PathBuf>,

    /// output path. DEFAULT: stdout
    #[argh(option, short = 'o')]
    output: Option<PathBuf>,

    /// the number of threads to use. DEFAULT: max_available
    #[argh(option, short = 't', default = "num_cpus::get()")]
    threads: usize,

    /// the ideal number of basepairs each worker receives. Total bp in memory at one time is (threads - 2) * chunksize
    #[argh(option, short = 'c')]
    chunksize: Option<usize>, // default set by par_io at 1_000_000

    /// SAM flags to include. DEFAULT: 0
    #[argh(option, short = 'f', default = "0")]
    include_flags: u16,

    /// SAM flags to exclude, recommended 3848. DEFAULT: 0
    #[argh(option, short = 'F', default = "0")]
    exclude_flags: u16,

    /// minimum mapq for a read to count toward depth. DEFAULT: 0
    #[argh(option, short = 'q', default = "0")]
    min_mapq: u8,
}

impl SimpleDepth {
    // TODO: Add mate detection like sambamba
    // TODO: move rayon pool setting to par_io
    pub fn run(self) -> Result<()> {
        info!("Running simple-depth on: {:?}", self.reads);

        let cpus = utils::determine_allowed_cpus(self.threads)?;
        // Keep two around for main thread and thread running the pool
        let usable_cpus = std::cmp::max(cpus.checked_sub(2).unwrap_or(0), 1);
        utils::set_rayon_global_pools_size(usable_cpus)?;

        let par_io_runner = par_io::ParIO::new(
            self.reads.clone(),
            self.ref_fasta.clone(),
            self.bed_file.clone(),
            self.output.clone(),
            Some(usable_cpus),
            self.chunksize.clone(),
        );

        // par_io_runner.process(Self::process_region)?;
        par_io_runner.process(move |tid, start, stop| {
            info!("Processing region {}:{}-{}", tid, start, stop);
            // Create a reader
            let mut reader =
                bam::IndexedReader::from_path(&self.reads).expect("Indexed Reader for region");

            // If passed add ref_fasta
            if let Some(ref_fasta) = &self.ref_fasta {
                reader.set_reference(ref_fasta).expect("Set ref");
            }

            let header = reader.header().to_owned();
            // fetch the region of interest
            reader
                .fetch(tid, start as u64, stop as u64)
                .expect("Fetched a region");
            // Walk over pileups
            let result: Vec<Position> = reader
                .pileup()
                .flat_map(|p| {
                    let pileup = p.expect("Extracted a pileup");
                    // Verify that we are within the bounds of the chunk we are iterating on
                    if (pileup.pos() as usize) >= start && (pileup.pos() as usize) < stop {
                        Some(Position::from_pileup(pileup, &header, |record| {
                            let flags = record.flags();
                            (!flags) & &self.include_flags == 0
                                && flags & &self.exclude_flags == 0
                                && &record.mapq() >= &self.min_mapq
                        }))
                    } else {
                        None
                    }
                })
                .collect();
            result
        })?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;
    use rust_htslib::{bam, bam::record::Record};
    use std::path::PathBuf;

    #[fixture]
    fn positions() -> Vec<Vec<Position>> {
        // This keep the test bam up to date
        let path = PathBuf::from("test/test.bam");
        // Build a header
        let mut header = bam::header::Header::new();
        let mut chr1 = bam::header::HeaderRecord::new(b"SQ");
        chr1.push_tag(b"SN", &"chr1".to_owned());
        chr1.push_tag(b"LN", &"100".to_owned());
        let mut chr2 = bam::header::HeaderRecord::new(b"SQ");
        chr2.push_tag(b"SN", &"chr2".to_owned());
        chr2.push_tag(b"LN", &"100".to_owned());
        header.push_record(&chr1);
        header.push_record(&chr2);
        let view = bam::HeaderView::from_header(&header);

        // Add records
        let records = vec![
            // Chr1 - Nice
            Record::from_sam(&view, b"ONE\t67\tchr1\t1\t40\t25M\tchr1\t50\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            Record::from_sam(&view, b"TWO\t67\tchr1\t5\t40\t25M\tchr1\t55\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            Record::from_sam(&view, b"THREE\t67\tchr1\t10\t40\t25M\tchr1\t60\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            Record::from_sam(&view, b"FOUR\t67\tchr1\t15\t40\t25M\tchr1\t65\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            Record::from_sam(&view, b"FIVE\t67\tchr1\t20\t40\t25M\tchr1\t70\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            Record::from_sam(&view, b"ONE_2\t147\tchr1\t50\t40\t25M\tchr1\t1\t75\tTTTTTTTTTTTTTTTTTTTTTTTTT\t#########################").unwrap(),
            Record::from_sam(&view, b"TWO_2\t147\tchr1\t55\t40\t25M\tchr1\t5\t75\tGGGGGGGGGGGGGGGGGGGGGGGGG\t#########################").unwrap(),
            Record::from_sam(&view, b"THREE_2\t147\tchr1\t60\t40\t25M\tchr1\t10\t75\tCCCCCCCCCCCCCCCCCCCCCCCCC\t#########################").unwrap(),
            Record::from_sam(&view, b"FOUR_2\t147\tchr1\t65\t40\t25M\tchr1\t15\t75\tNNNNNNNNNNNNNNNNNNNNNNNNN\t#########################").unwrap(),
            Record::from_sam(&view, b"FIVE_2\t147\tchr1\t70\t40\t25M\tchr1\t20\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),

            // Chr2 - Complex
            // Ins
            Record::from_sam(&view, b"ONE\t67\tchr2\t1\t40\t2M2I21M\tchr2\t50\t75\tAAGGGAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            // Del
            Record::from_sam(&view, b"TWO\t67\tchr2\t5\t40\t2M5D23M\tchr2\t55\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            // Skip
            Record::from_sam(&view, b"THREE\t67\tchr2\t10\t40\t3M5N22M\tchr2\t60\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            // Mismatch
            Record::from_sam(&view, b"FOUR\t67\tchr2\t15\t40\t25M\tchr2\t65\t75\tATAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            // Overlapping mates
            Record::from_sam(&view, b"FIVE\t67\tchr2\t20\t40\t25M\tchr2\t35\t40\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            Record::from_sam(&view, b"FIVE_2\t147\tchr2\t35\t40\t25M\tchr2\t20\t40\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),

            // Other base
            Record::from_sam(&view, b"ONE_2\t147\tchr2\t50\t40\t25M\tchr2\t1\t75\tAAAAAAAAAAAAAAAAAAAAAYAAA\t#########################").unwrap(),

            Record::from_sam(&view, b"TWO_2\t147\tchr2\t55\t40\t25M\tchr2\t5\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            Record::from_sam(&view, b"THREE_2\t147\tchr2\t60\t40\t25M\tchr2\t10\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            // A failure of QC
            Record::from_sam(&view, b"FOUR_2\t659\tchr2\t65\t40\t25M\tchr2\t15\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
        ];

        // Update the test/test.bam file
        let mut writer =
            bam::Writer::from_path(&path, &header, bam::Format::BAM).expect("Created writer");
        for record in records.iter() {
            writer.write(record).expect("Wrote record");
        }
        drop(writer); // Drop writer so filehandle closes

        // Extract bam into Positions
        let mut reader = bam::Reader::from_path(&path).expect("Opened bam for reading");
        let header = reader.header().to_owned();
        let mut positions = vec![vec![], vec![]];
        for p in reader.pileup() {
            let pileup = p.unwrap();
            let tid = pileup.tid();
            let pos = Position::from_pileup(pileup, &header, |record| {
                let flags = record.flags();
                (!flags) & (0 as u16) == 0
                    && flags & (3848 as u16) == 0
                    && &record.mapq() >= &(0 as u8)
            });
            positions[tid as usize].push(pos);
        }
        positions
    }

    #[rstest]
    fn check_insertions(positions: Vec<Vec<Position>>) {
        assert_eq!(positions[1][0].ins, 0);
        assert_eq!(positions[1][1].ins, 1);
        assert_eq!(positions[1][2].ins, 0);
    }

    #[rstest]
    fn check_deletions(positions: Vec<Vec<Position>>) {
        assert_eq!(positions[1][5].del, 0);
        assert_eq!(positions[1][6].del, 1);
        assert_eq!(positions[1][7].del, 1);
        assert_eq!(positions[1][8].del, 1);
        assert_eq!(positions[1][9].del, 1);
        assert_eq!(positions[1][10].del, 1);
        assert_eq!(positions[1][11].del, 0);
    }

    #[rstest]
    fn check_refskips(positions: Vec<Vec<Position>>) {
        assert_eq!(positions[1][11].ref_skip, 0);
        assert_eq!(positions[1][12].ref_skip, 1);
        assert_eq!(positions[1][13].ref_skip, 1);
        assert_eq!(positions[1][14].ref_skip, 1);
        assert_eq!(positions[1][15].ref_skip, 1);
        assert_eq!(positions[1][16].ref_skip, 1);
        assert_eq!(positions[1][17].ref_skip, 0);
    }

    #[rstest]
    // TODO: Make mate aware
    fn check_depths(positions: Vec<Vec<Position>>) {
        dbg!(&positions[0]);
        assert_eq!(positions[0][0].depth, 1);
        assert_eq!(positions[0][4].depth, 2);
        assert_eq!(positions[0][9].depth, 3);
        assert_eq!(positions[0][14].depth, 4);
        assert_eq!(positions[0][19].depth, 5);
        assert_eq!(positions[0][25].depth, 4);
        assert_eq!(positions[0][29].depth, 3);
        assert_eq!(positions[0][34].depth, 2);
        assert_eq!(positions[0][39].depth, 1);
        // NB: -6 bc there are 6 positions with no coverage from 44-50
        assert_eq!(positions[0][78 - 6].depth, 4);
    }

    #[rstest]
    fn check_filters(positions: Vec<Vec<Position>>) {
        // Verify that a read that has flags saying it failed QC got filtered out
        assert_eq!(positions[1][81].depth, 1);
        assert_eq!(positions[1][84].depth, 0);
    }

    #[rstest]
    fn check_depths_insertions(positions: Vec<Vec<Position>>) {
        assert_eq!(positions[1][0].depth, 1);
        assert_eq!(positions[1][1].depth, 1); // Insertion is here
        assert_eq!(positions[1][2].depth, 1);
    }

    #[rstest]
    fn check_depths_deletions(positions: Vec<Vec<Position>>) {
        assert_eq!(positions[1][5].depth, 2);
        assert_eq!(positions[1][6].depth, 2); // Del
        assert_eq!(positions[1][7].depth, 2); // Del
        assert_eq!(positions[1][8].depth, 2); // Del
        assert_eq!(positions[1][9].depth, 3); // Del
        assert_eq!(positions[1][10].depth, 3); // Del
        assert_eq!(positions[1][11].depth, 3);
    }

    #[rstest]
    fn check_depths_refskips(positions: Vec<Vec<Position>>) {
        assert_eq!(positions[1][11].depth, 3);
        assert_eq!(positions[1][12].depth, 3); // Skip
        assert_eq!(positions[1][13].depth, 3); // Skip
        assert_eq!(positions[1][14].depth, 4); // Skip
        assert_eq!(positions[1][15].depth, 4); // Skip
        assert_eq!(positions[1][16].depth, 4);
        assert_eq!(positions[1][17].depth, 4);
    }
}
