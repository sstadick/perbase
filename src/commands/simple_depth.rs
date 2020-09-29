//! # Simple Depth
//!
//! Simple single pass over bam file to calculate depth at each position
//! as well as depth per nucleotide. Additionally counts the number of
//! insertions / deletions at each position.
use anyhow::Result;
use log::*;
use perbase_lib::{
    par_granges::{self, RegionProcessor},
    position::{Position, ReadFilter},
    utils,
};
use rust_htslib::{bam, bam::record::Record, bam::Read};
use std::path::PathBuf;
use structopt::StructOpt;

/// Calculate the depth at each base, per-nucleotide.
#[derive(StructOpt)]
#[structopt(author)]
pub struct SimpleDepth {
    /// Input indexed BAM/CRAM to analyze.
    reads: PathBuf,

    /// Indexed reference fasta, set if using CRAM.
    #[structopt(long, short = "r")]
    ref_fasta: Option<PathBuf>,

    /// A BED file containing regions of interest. If specified, only bases from the given regions will be reported on.
    #[structopt(long, short = "b")]
    bed_file: Option<PathBuf>,

    /// Output path, defaults to stdout.
    #[structopt(long, short = "o")]
    output: Option<PathBuf>,

    /// The number of threads to use.
    #[structopt(long, short = "t", default_value = utils::NUM_CPU.as_str())]
    threads: usize,

    /// The ideal number of basepairs each worker receives. Total bp in memory at one time is (threads - 2) * chunksize.
    #[structopt(long, short = "c")]
    chunksize: Option<usize>, // default set by par_granges at 1_000_000

    /// SAM flags to include.
    #[structopt(long, short = "f", default_value = "0")]
    include_flags: u16,

    /// SAM flags to exclude, recommended 3848.
    #[structopt(long, short = "F", default_value = "0")]
    exclude_flags: u16,

    /// Fix overlapping mates counts, see docs for full details.
    #[structopt(long, short = "m")]
    mate_fix: bool,

    /// Minimum MAPQ for a read to count toward depth.
    #[structopt(long, short = "q", default_value = "0")]
    min_mapq: u8,
}

impl SimpleDepth {
    pub fn run(self) -> Result<()> {
        info!("Running simple-depth on: {:?}", self.reads);
        let cpus = utils::determine_allowed_cpus(self.threads)?;

        let read_filter =
            SimpleReadFilter::new(self.include_flags, self.exclude_flags, self.min_mapq);
        let simple_processor = SimpleProcessor::new(
            self.reads.clone(),
            self.ref_fasta.clone(),
            self.mate_fix,
            read_filter,
        );

        let par_granges_runner = par_granges::ParGranges::new(
            self.reads.clone(),
            self.ref_fasta.clone(),
            self.bed_file.clone(),
            self.output.clone(),
            Some(cpus),
            self.chunksize.clone(),
            simple_processor,
        );

        par_granges_runner.process()?;
        Ok(())
    }
}

/// A Simple impl of [ReadFilter]
pub struct SimpleReadFilter {
    include_flags: u16,
    exclude_flags: u16,
    min_mapq: u8,
}

impl SimpleReadFilter {
    /// Create a SimpleReadFilter
    fn new(include_flags: u16, exclude_flags: u16, min_mapq: u8) -> Self {
        Self {
            include_flags,
            exclude_flags,
            min_mapq,
        }
    }
}

impl ReadFilter for SimpleReadFilter {
    /// Filter reads based SAM flags and mapping quality
    #[inline]
    fn filter_read(&self, read: &Record) -> bool {
        let flags = read.flags();
        (!flags) & &self.include_flags == 0
            && flags & &self.exclude_flags == 0
            && &read.mapq() >= &self.min_mapq
    }
}

/// Holds the info needed for [par_io::RegionProcessor] implementation
struct SimpleProcessor<F: ReadFilter> {
    /// path to indexed BAM/CRAM
    reads: PathBuf,
    /// path to indexed ref file
    ref_fasta: Option<PathBuf>,
    /// Indicate whether or not to account for overlapping mates.
    mate_fix: bool,
    /// implementation of [position::ReadFilter] that will be used
    read_filter: F,
}

impl<F: ReadFilter> SimpleProcessor<F> {
    /// Create a new SimpleProcessor
    fn new(reads: PathBuf, ref_fasta: Option<PathBuf>, mate_fix: bool, read_filter: F) -> Self {
        Self {
            reads,
            ref_fasta,
            mate_fix,
            read_filter,
        }
    }
}

/// Implement [par_io::RegionProcessor] for [SimpleProcessor]
impl<F: ReadFilter> RegionProcessor for SimpleProcessor<F> {
    /// Objects of [position::Position] will be returned by each call to [SimpleProcessor::process_region]
    type P = Position;

    /// Process a region by fetching it from a BAM/CRAM, getting a pileup, and then
    /// walking the pileup (checking bounds) to create Position objects according to
    /// the defined filters
    fn process_region(&self, tid: u32, start: u64, stop: u64) -> Vec<Position> {
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
            .fetch(tid, start, stop)
            .expect("Fetched a region");
        // Walk over pileups
        let result: Vec<Position> = reader
            .pileup()
            .flat_map(|p| {
                let pileup = p.expect("Extracted a pileup");
                // Verify that we are within the bounds of the chunk we are iterating on
                if (pileup.pos() as u64) >= start && (pileup.pos() as u64) < stop {
                    if self.mate_fix {
                        Some(Position::from_pileup_mate_aware(
                            pileup,
                            &header,
                            &self.read_filter,
                        ))
                    } else {
                        Some(Position::from_pileup(pileup, &header, &self.read_filter))
                    }
                } else {
                    None
                }
            })
            .collect();
        result
    }
}

#[cfg(test)]
#[allow(unused)]
mod tests {
    use super::*;
    use perbase_lib::position::Position;
    use rstest::*;
    use rust_htslib::{bam, bam::record::Record};
    use std::path::PathBuf;
    use tempfile::{tempdir, TempDir};

    #[fixture]
    fn read_filter() -> SimpleReadFilter {
        SimpleReadFilter::new(0, 512, 0)
    }

    #[fixture]
    fn bamfile() -> (PathBuf, TempDir) {
        // No longer - This keep the test bam up to date
        // let test_path = PathBuf::from("test/test.bam");
        let tempdir = tempdir().unwrap();
        let path = tempdir.path().join("test.bam");

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
            Record::from_sam(&view, b"ONE\t147\tchr1\t50\t40\t25M\tchr1\t1\t75\tTTTTTTTTTTTTTTTTTTTTTTTTT\t#########################").unwrap(),
            Record::from_sam(&view, b"TWO\t147\tchr1\t55\t40\t25M\tchr1\t5\t75\tGGGGGGGGGGGGGGGGGGGGGGGGG\t#########################").unwrap(),
            Record::from_sam(&view, b"THREE\t147\tchr1\t60\t40\t25M\tchr1\t10\t75\tCCCCCCCCCCCCCCCCCCCCCCCCC\t#########################").unwrap(),
            Record::from_sam(&view, b"FOUR\t147\tchr1\t65\t40\t25M\tchr1\t15\t75\tNNNNNNNNNNNNNNNNNNNNNNNNN\t#########################").unwrap(),
            Record::from_sam(&view, b"FIVE\t147\tchr1\t70\t40\t25M\tchr1\t20\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),

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
            Record::from_sam(&view, b"FIVE\t147\tchr2\t35\t40\t25M\tchr2\t20\t40\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),

            // Other base
            Record::from_sam(&view, b"ONE\t147\tchr2\t50\t40\t25M\tchr2\t1\t75\tAAAAAAAAAAAAAAAAAAAAAYAAA\t#########################").unwrap(),

            Record::from_sam(&view, b"TWO\t147\tchr2\t55\t40\t25M\tchr2\t5\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            Record::from_sam(&view, b"THREE\t147\tchr2\t60\t40\t25M\tchr2\t10\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            // A failure of QC
            Record::from_sam(&view, b"FOUR\t659\tchr2\t65\t40\t25M\tchr2\t15\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
        ];

        // Update the test/test.bam file
        let mut writer =
            bam::Writer::from_path(&path, &header, bam::Format::BAM).expect("Created writer");
        for record in records.iter() {
            writer.write(record).expect("Wrote record");
        }
        (path, tempdir)
    }

    #[fixture]
    fn non_mate_aware_positions(
        bamfile: (PathBuf, TempDir),
        read_filter: SimpleReadFilter,
    ) -> Vec<Vec<Position>> {
        // Extract bam into Positions
        let mut reader = bam::Reader::from_path(&bamfile.0).expect("Opened bam for reading");
        let header = reader.header().to_owned();
        let mut positions = vec![vec![], vec![]];
        for p in reader.pileup() {
            let pileup = p.unwrap();
            let tid = pileup.tid();
            let pos = Position::from_pileup(pileup, &header, &read_filter);
            positions[tid as usize].push(pos);
        }
        positions
    }

    #[fixture]
    fn mate_aware_positions(
        bamfile: (PathBuf, TempDir),
        read_filter: SimpleReadFilter,
    ) -> Vec<Vec<Position>> {
        // Extract bam into Positions
        let mut reader = bam::Reader::from_path(&bamfile.0).expect("Opened bam for reading");
        let header = reader.header().to_owned();
        let mut positions = vec![vec![], vec![]];
        for p in reader.pileup() {
            let pileup = p.unwrap();
            let tid = pileup.tid();
            let pos = Position::from_pileup_mate_aware(pileup, &header, &read_filter);
            positions[tid as usize].push(pos);
        }
        positions
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0)
    )]
    fn check_insertions(positions: Vec<Vec<Position>>, awareness_modifier: usize) {
        assert_eq!(positions[1][0].ins, 0);
        assert_eq!(positions[1][1].ins, 1);
        assert_eq!(positions[1][2].ins, 0);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0)
    )]
    fn check_deletions(positions: Vec<Vec<Position>>, awareness_modifier: usize) {
        assert_eq!(positions[1][5].del, 0);
        assert_eq!(positions[1][6].del, 1);
        assert_eq!(positions[1][7].del, 1);
        assert_eq!(positions[1][8].del, 1);
        assert_eq!(positions[1][9].del, 1);
        assert_eq!(positions[1][10].del, 1);
        assert_eq!(positions[1][11].del, 0);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0)
    )]
    fn check_refskips(positions: Vec<Vec<Position>>, awareness_modifier: usize) {
        assert_eq!(positions[1][11].ref_skip, 0);
        assert_eq!(positions[1][12].ref_skip, 1);
        assert_eq!(positions[1][13].ref_skip, 1);
        assert_eq!(positions[1][14].ref_skip, 1);
        assert_eq!(positions[1][15].ref_skip, 1);
        assert_eq!(positions[1][16].ref_skip, 1);
        assert_eq!(positions[1][17].ref_skip, 0);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0)
    )]
    fn check_depths(positions: Vec<Vec<Position>>, awareness_modifier: usize) {
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

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0)
    )]
    fn check_filters(positions: Vec<Vec<Position>>, awareness_modifier: usize) {
        // Verify that a read that has flags saying it failed QC got filtered out
        assert_eq!(positions[1][81].depth, 1);
        assert_eq!(positions[1][84].depth, 0);
        assert_eq!(positions[1][81].fail, 1);
        assert_eq!(positions[1][84].fail, 1);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0)
    )]
    fn check_depths_insertions(positions: Vec<Vec<Position>>, awareness_modifier: usize) {
        assert_eq!(positions[1][0].depth, 1);
        assert_eq!(positions[1][1].depth, 1); // Insertion is here
        assert_eq!(positions[1][2].depth, 1);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0)
    )]
    fn check_depths_deletions(positions: Vec<Vec<Position>>, awareness_modifier: usize) {
        assert_eq!(positions[1][5].depth, 2);
        assert_eq!(positions[1][6].depth, 2); // Del
        assert_eq!(positions[1][7].depth, 2); // Del
        assert_eq!(positions[1][8].depth, 2); // Del
        assert_eq!(positions[1][9].depth, 3); // Del
        assert_eq!(positions[1][10].depth, 3); // Del
        assert_eq!(positions[1][11].depth, 3);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0)
    )]
    fn check_depths_refskips(positions: Vec<Vec<Position>>, awareness_modifier: usize) {
        assert_eq!(positions[1][11].depth, 3);
        assert_eq!(positions[1][12].depth, 2); // Skip
        assert_eq!(positions[1][13].depth, 2); // Skip
        assert_eq!(positions[1][14].depth, 3); // Skip
        assert_eq!(positions[1][15].depth, 3); // Skip
        assert_eq!(positions[1][16].depth, 3); // Skip
        assert_eq!(positions[1][17].depth, 4);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 1)
    )]
    fn check_mate_detection(positions: Vec<Vec<Position>>, awareness_modifier: usize) {
        assert_eq!(positions[1][33].depth, 4);
        assert_eq!(positions[1][34].depth, 4 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][35].depth, 4 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][36].depth, 4 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][37].depth, 4 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][38].depth, 4 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][39].depth, 2 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][40].depth, 2 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][41].depth, 2 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][42].depth, 2 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][43].depth, 2 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][44].depth, 1);

        assert_eq!(positions[1][33].a, 4);
        assert_eq!(positions[1][34].a, 4 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][35].a, 4 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][36].a, 4 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][37].a, 4 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][38].a, 4 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][39].a, 2 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][40].a, 2 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][41].a, 2 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][42].a, 2 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][43].a, 2 - awareness_modifier); // mate overlap
        assert_eq!(positions[1][44].a, 1);
    }
}
