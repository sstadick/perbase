//! # Only Depth
//!
//! Calculates the depth only at each position. This uses the same algorithm described in
//! the [`mosdepth`](https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btx699/4583630?guestAccessKey=35b55064-4566-4ab3-a769-32916fa1c6e6)
//! paper.
use anyhow::Result;
use csv;
use grep_cli::stdout;
use itertools::Itertools;
use log::*;
use perbase_lib::{
    par_granges::{self, RegionProcessor},
    position::ReadFilter,
    utils,
};
use rust_htslib::{
    bam,
    bam::pileup::{Alignment, Pileup},
    bam::record::Cigar,
    bam::record::Record,
    bam::Read,
};
use serde::Serialize;
use smartstring::alias::String;
use std::{
    cmp::Ordering,
    default,
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};
use structopt::StructOpt;
use termcolor::ColorChoice;

/// Calculate the only the depth at each base.
#[derive(StructOpt)]
#[structopt(author)]
pub struct OnlyDepth {
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

impl OnlyDepth {
    pub fn run(self) -> Result<()> {
        info!("Running only-depth on: {:?}", self.reads);
        let cpus = utils::determine_allowed_cpus(self.threads)?;

        let mut writer = self.get_writer()?;

        let read_filter =
            OnlyDepthReadFilter::new(self.include_flags, self.exclude_flags, self.min_mapq);
        let processor = OnlyDepthProcessor::new(
            self.reads.clone(),
            self.ref_fasta.clone(),
            self.mate_fix,
            read_filter,
        );

        let par_granges_runner = par_granges::ParGranges::new(
            self.reads.clone(),
            self.ref_fasta.clone(),
            self.bed_file.clone(),
            Some(cpus),
            self.chunksize.clone(),
            processor,
        );

        let receiver = par_granges_runner.process()?;

        receiver
            .into_iter()
            .for_each(|pos| writer.serialize(pos).unwrap());
        writer.flush()?;
        Ok(())
    }

    /// Open a CSV Writer to a file or stdout
    fn get_writer(&self) -> Result<csv::Writer<Box<dyn Write>>> {
        let raw_writer: Box<dyn Write> = match &self.output {
            Some(path) if path.to_str().unwrap() != "-" => {
                Box::new(BufWriter::new(File::open(path)?))
            }
            _ => Box::new(stdout(ColorChoice::Never)),
        };
        Ok(csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(raw_writer))
    }
}

/// A straightforward read filter.
pub struct OnlyDepthReadFilter {
    include_flags: u16,
    exclude_flags: u16,
    min_mapq: u8,
}

impl OnlyDepthReadFilter {
    /// Create an OnlyDepthReadFilter
    fn new(include_flags: u16, exclude_flags: u16, min_mapq: u8) -> Self {
        Self {
            include_flags,
            exclude_flags,
            min_mapq,
        }
    }
}

impl ReadFilter for OnlyDepthReadFilter {
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
struct OnlyDepthProcessor<F: ReadFilter> {
    /// path to indexed BAM/CRAM
    reads: PathBuf,
    /// path to indexed ref file
    ref_fasta: Option<PathBuf>,
    /// Indicate whether or not to account for overlapping mates.
    mate_fix: bool,
    /// implementation of [position::ReadFilter] that will be used
    read_filter: F,
}

impl<F: ReadFilter> OnlyDepthProcessor<F> {
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
impl<F: ReadFilter> RegionProcessor for OnlyDepthProcessor<F> {
    /// Objects of [position::Position] will be returned by each call to [SimpleProcessor::process_region]
    type P = OnlyDepthPosition;

    /// Process a region by fetching it from a BAM/CRAM, getting a pileup, and then
    /// walking the pileup (checking bounds) to create Position objects according to
    /// the defined filters
    // TODO: make `Position` a trait that can be implemented for each new depth tester
    fn process_region(&self, tid: u32, start: u64, stop: u64) -> Vec<OnlyDepthPosition> {
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
        reader.fetch(tid, start, stop).expect("Fetched a region");
        // Walk over pileups
        let result: Vec<OnlyDepthPosition> = reader
            .pileup()
            .flat_map(|p| {
                let pileup = p.expect("Extracted a pileup");
                // Verify that we are within the bounds of the chunk we are iterating on
                if (pileup.pos() as u64) >= start && (pileup.pos() as u64) < stop {
                    if self.mate_fix {
                        Some(OnlyDepthPosition::from_pileup_mate_aware(
                            pileup,
                            &header,
                            &self.read_filter,
                        ))
                    } else {
                        Some(OnlyDepthPosition::from_pileup(
                            pileup,
                            &header,
                            &self.read_filter,
                        ))
                    }
                } else {
                    None
                }
            })
            .collect();
        result
    }
}

pub struct IterStartStops {
    // offset into the read ref positions to start at
    offset: usize,
    // hardstop ref positon in reead to stop at
    hardstop: usize,
    // index into the cigar string
    cigar_index: usize,
    // vec of cigar ops
    cigar: Vec<Cigar>,
}

impl IterStartStops {
    fn new(cigar: Vec<Cigar>, offset: usize, hardstop: usize) -> Self {
        IterStartStops {
            cigar_index: 0,
            hardstop,
            offset,
            cigar,
        }
    }
}

/// Iterator over the start-stop pairs given a cigar string and a position offset
impl Iterator for IterStartStops {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        while self.cigar_index < self.cigar.len() {
            let mut start = self.offset;
            let mut ref_consumed = 0;
            let entry = self.cigar[self.cigar_index];
            match entry {
                // refconsuming and not breaking
                Cigar::Match(olen) | Cigar::Diff(olen) | Cigar::Equal(olen) | Cigar::Del(olen) => {
                    ref_consumed += olen as usize
                }
                // non-ref consuming and not breaking
                Cigar::Ins(_) | Cigar::SoftClip(_) => (),
                // ref-consuming and breaking
                Cigar::RefSkip(olen) => {
                    let stop = start + ref_consumed + olen as usize;
                    let result = (start, stop);
                    self.offset = stop;
                    self.cigar_index += 1;
                    return Some(result);
                }
                // no-ops
                Cigar::Pad(_) | Cigar::HardClip(_) => (),
                _ => unreachable!(),
            }
            self.cigar_index += 1;
        }
        None
    }
}

/// Hold all information about a position.
#[derive(Debug, Serialize, Default)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub struct OnlyDepthPosition {
    /// Reference sequence name.
    #[serde(rename = "REF")]
    pub ref_seq: String,
    /// 1-based position in the sequence.
    pub pos: usize,
    /// Total depth at this position.
    pub depth: usize,
}

impl OnlyDepthPosition {
    pub fn new(ref_seq: String, pos: usize) -> Self {
        Self {
            ref_seq,
            pos,
            ..default::Default::default()
        }
    }

    /// Given a record, update the counts at this position
    fn update<F: ReadFilter>(&mut self, alignment: &Alignment, record: Record, read_filter: &F) {
        if !read_filter.filter_read(&record) || alignment.is_refskip() {
            self.depth -= 1;
        }
    }

    /// Convert a pileup into a `Position`.
    ///
    /// This will walk over each of the alignments and count the number each nucleotide it finds.
    /// It will also count the number of Ins/Dels/Skips that are at each position. The output of this 1-based.
    ///
    /// # Arguments
    ///
    /// * `pileup` - a pileup at a genomic position
    /// * `header` - a headerview for the bam file being read, to get the sequence name
    /// * `read_filter` - a function to filter out reads, returning false will cause a read to be filtered
    pub fn from_pileup<F: ReadFilter>(
        pileup: Pileup,
        header: &bam::HeaderView,
        read_filter: &F,
    ) -> Self {
        let name = std::str::from_utf8(header.tid2name(pileup.tid())).unwrap();
        // make output 1-based
        let mut pos = Self::new(String::from(name), (pileup.pos() + 1) as usize);
        pos.depth = pileup.depth() as usize;

        for alignment in pileup.alignments() {
            let record = alignment.record();
            &pos.update(&alignment, record, read_filter);
        }
        pos
    }

    /// Convert a pileup into a `Position`.
    ///
    /// This will walk over each of the alignments and count the number each nucleotide it finds.
    /// It will also count the number of Ins/Dels/Skips that are at each position. The output of this 1-based.
    /// Additionally, this method is mate aware. Before processing a position it will scan the alignments for mates.
    /// If a mate is found, it will try to take use the mate that has the highest MAPQ, breaking ties by choosing the
    /// first in pair that passes filters. In the event of both failing filters or not being first in pair, the first
    /// read encountered is kept.
    ///
    /// # Arguments
    ///
    /// * `pileup` - a pileup at a genomic position
    /// * `header` - a headerview for the bam file being read, to get the sequence name
    /// * `read_filter` - a function to filter out reads, returning false will cause a read to be filtered
    pub fn from_pileup_mate_aware<F: ReadFilter>(
        pileup: Pileup,
        header: &bam::HeaderView,
        read_filter: &F,
    ) -> Self {
        let name = std::str::from_utf8(header.tid2name(pileup.tid())).unwrap();
        // make output 1-based
        let mut pos = Self::new(String::from(name), (pileup.pos() + 1) as usize);
        pos.depth = pileup.depth() as usize;

        // Group records by qname
        let grouped_by_qname = pileup
            .alignments()
            .map(|aln| {
                let record = aln.record();
                (aln, record)
            })
            .sorted_by(|a, b| Ord::cmp(a.1.qname(), b.1.qname()))
            // TODO: I'm not sure there is a good way to remove this allocation
            .group_by(|a| a.1.qname().to_owned());

        for (_qname, reads) in grouped_by_qname.into_iter() {
            // Choose the best of the reads based on mapq, if tied, check which is first and passes filters
            let mut total_reads = 0; // count how many reads there were
            let (alignment, record) = reads
                .into_iter()
                .map(|x| {
                    total_reads += 1;
                    x
                })
                .max_by(|a, b| match a.1.mapq().cmp(&b.1.mapq()) {
                    Ordering::Greater => Ordering::Greater,
                    Ordering::Less => Ordering::Less,
                    Ordering::Equal => {
                        // Check if a is first in pair
                        if a.1.flags() & 64 == 0 && read_filter.filter_read(&a.1) {
                            Ordering::Greater
                        } else if b.1.flags() & 64 == 0 && read_filter.filter_read(&b.1) {
                            Ordering::Less
                        } else {
                            // Default to `a` in the event that there is no first in pair for some reason
                            Ordering::Greater
                        }
                    }
                })
                .unwrap();
            // decrement depth for each read not used
            pos.depth -= total_reads - 1;
            pos.update(&alignment, record, read_filter);
        }
        pos
    }
}

#[cfg(test)]
#[allow(unused)]
mod tests {
    use super::*;
    use perbase_lib::position::Position;
    use rstest::*;
    use rust_htslib::{bam, bam::record::Record};
    use smartstring::alias::*;
    use std::{collections::HashMap, path::PathBuf};
    use tempfile::{tempdir, TempDir};

    #[fixture]
    fn read_filter() -> OnlyDepthReadFilter {
        OnlyDepthReadFilter::new(0, 512, 0)
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
        drop(writer); // force it to flush so indexing can happen
                      // build the index
        bam::index::build(&path, None, bam::index::Type::BAI, 1).unwrap();
        (path, tempdir)
    }

    #[fixture]
    fn non_mate_aware_positions(
        bamfile: (PathBuf, TempDir),
        read_filter: OnlyDepthReadFilter,
    ) -> HashMap<String, Vec<OnlyDepthPosition>> {
        let cpus = utils::determine_allowed_cpus(8).unwrap();

        let simple_processor = OnlyDepthProcessor::new(bamfile.0.clone(), None, false, read_filter);

        let par_granges_runner = par_granges::ParGranges::new(
            bamfile.0,
            None,
            None,       // TODO - make a test with befile
            Some(cpus), // TODO - parameterize over this
            None,       // TODO - parameterize over this
            simple_processor,
        );
        let mut positions = HashMap::new();
        par_granges_runner
            .process()
            .unwrap()
            .into_iter()
            .for_each(|p| {
                let pos = positions.entry(p.ref_seq.clone()).or_insert(vec![]);
                pos.push(p)
            });
        positions
    }

    #[fixture]
    fn mate_aware_positions(
        bamfile: (PathBuf, TempDir),
        read_filter: OnlyDepthReadFilter,
    ) -> HashMap<String, Vec<OnlyDepthPosition>> {
        let cpus = utils::determine_allowed_cpus(8).unwrap();

        let simple_processor = OnlyDepthProcessor::new(
            bamfile.0.clone(),
            None,
            true, // mate aware
            read_filter,
        );

        let par_granges_runner = par_granges::ParGranges::new(
            bamfile.0,
            None,
            None,       // TODO - make a test with befile
            Some(cpus), // TODO - parameterize over this
            None,       // TODO - parameterize over this
            simple_processor,
        );
        let mut positions = HashMap::new();
        par_granges_runner
            .process()
            .unwrap()
            .into_iter()
            .for_each(|p| {
                let pos = positions.entry(p.ref_seq.clone()).or_insert(vec![]);
                pos.push(p)
            });
        positions
    }

    // #[rstest(
    //     positions,
    //     awareness_modifier,
    //     case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
    //     case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0)
    // )]
    // fn check_insertions(positions: HashMap<String, Vec<Position>>, awareness_modifier: usize) {
    //     assert_eq!(positions.get("chr2").unwrap()[0].ins, 0);
    //     assert_eq!(positions.get("chr2").unwrap()[1].ins, 1);
    //     assert_eq!(positions.get("chr2").unwrap()[2].ins, 0);
    // }

    // #[rstest(
    //     positions,
    //     awareness_modifier,
    //     case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
    //     case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0)
    // )]
    // fn check_deletions(positions: HashMap<String, Vec<Position>>, awareness_modifier: usize) {
    //     assert_eq!(positions.get("chr2").unwrap()[5].del, 0);
    //     assert_eq!(positions.get("chr2").unwrap()[6].del, 1);
    //     assert_eq!(positions.get("chr2").unwrap()[7].del, 1);
    //     assert_eq!(positions.get("chr2").unwrap()[8].del, 1);
    //     assert_eq!(positions.get("chr2").unwrap()[9].del, 1);
    //     assert_eq!(positions.get("chr2").unwrap()[10].del, 1);
    //     assert_eq!(positions.get("chr2").unwrap()[11].del, 0);
    // }

    // #[rstest(
    //     positions,
    //     awareness_modifier,
    //     case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
    //     case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0)
    // )]
    // fn check_refskips(positions: HashMap<String, Vec<Position>>, awareness_modifier: usize) {
    //     assert_eq!(positions.get("chr2").unwrap()[11].ref_skip, 0);
    //     assert_eq!(positions.get("chr2").unwrap()[12].ref_skip, 1);
    //     assert_eq!(positions.get("chr2").unwrap()[13].ref_skip, 1);
    //     assert_eq!(positions.get("chr2").unwrap()[14].ref_skip, 1);
    //     assert_eq!(positions.get("chr2").unwrap()[15].ref_skip, 1);
    //     assert_eq!(positions.get("chr2").unwrap()[16].ref_skip, 1);
    //     assert_eq!(positions.get("chr2").unwrap()[17].ref_skip, 0);
    // }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0)
    )]
    fn check_depths(positions: HashMap<String, Vec<OnlyDepthPosition>>, awareness_modifier: usize) {
        assert_eq!(positions.get("chr1").unwrap()[0].depth, 1);
        assert_eq!(positions.get("chr1").unwrap()[4].depth, 2);
        assert_eq!(positions.get("chr1").unwrap()[9].depth, 3);
        assert_eq!(positions.get("chr1").unwrap()[14].depth, 4);
        assert_eq!(positions.get("chr1").unwrap()[19].depth, 5);
        assert_eq!(positions.get("chr1").unwrap()[25].depth, 4);
        assert_eq!(positions.get("chr1").unwrap()[29].depth, 3);
        assert_eq!(positions.get("chr1").unwrap()[34].depth, 2);
        assert_eq!(positions.get("chr1").unwrap()[39].depth, 1);
        // NB: -6 bc there are 6 positions with no coverage from 44-50
        assert_eq!(positions.get("chr1").unwrap()[78 - 6].depth, 4);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0)
    )]
    fn check_filters(
        positions: HashMap<String, Vec<OnlyDepthPosition>>,
        awareness_modifier: usize,
    ) {
        // Verify that a read that has flags saying it failed QC got filtered out
        assert_eq!(positions.get("chr2").unwrap()[81].depth, 1);
        assert_eq!(positions.get("chr2").unwrap()[84].depth, 0);
        // assert_eq!(positions.get("chr2").unwrap()[81].fail, 1);
        // assert_eq!(positions.get("chr2").unwrap()[84].fail, 1);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0)
    )]
    fn check_depths_insertions(
        positions: HashMap<String, Vec<OnlyDepthPosition>>,
        awareness_modifier: usize,
    ) {
        assert_eq!(positions.get("chr2").unwrap()[0].depth, 1);
        assert_eq!(positions.get("chr2").unwrap()[1].depth, 1); // Insertion is here
        assert_eq!(positions.get("chr2").unwrap()[2].depth, 1);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0)
    )]
    fn check_depths_deletions(
        positions: HashMap<String, Vec<OnlyDepthPosition>>,
        awareness_modifier: usize,
    ) {
        assert_eq!(positions.get("chr2").unwrap()[5].depth, 2);
        assert_eq!(positions.get("chr2").unwrap()[6].depth, 2); // Del
        assert_eq!(positions.get("chr2").unwrap()[7].depth, 2); // Del
        assert_eq!(positions.get("chr2").unwrap()[8].depth, 2); // Del
        assert_eq!(positions.get("chr2").unwrap()[9].depth, 3); // Del
        assert_eq!(positions.get("chr2").unwrap()[10].depth, 3); // Del
        assert_eq!(positions.get("chr2").unwrap()[11].depth, 3);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0)
    )]
    fn check_depths_refskips(
        positions: HashMap<String, Vec<OnlyDepthPosition>>,
        awareness_modifier: usize,
    ) {
        assert_eq!(positions.get("chr2").unwrap()[11].depth, 3);
        assert_eq!(positions.get("chr2").unwrap()[12].depth, 2); // Skip
        assert_eq!(positions.get("chr2").unwrap()[13].depth, 2); // Skip
        assert_eq!(positions.get("chr2").unwrap()[14].depth, 3); // Skip
        assert_eq!(positions.get("chr2").unwrap()[15].depth, 3); // Skip
        assert_eq!(positions.get("chr2").unwrap()[16].depth, 3); // Skip
        assert_eq!(positions.get("chr2").unwrap()[17].depth, 4);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 1)
    )]
    fn check_mate_detection(
        positions: HashMap<String, Vec<OnlyDepthPosition>>,
        awareness_modifier: usize,
    ) {
        assert_eq!(positions.get("chr2").unwrap()[33].depth, 4);
        assert_eq!(
            positions.get("chr2").unwrap()[34].depth,
            4 - awareness_modifier
        ); // mate overlap
        assert_eq!(
            positions.get("chr2").unwrap()[35].depth,
            4 - awareness_modifier
        ); // mate overlap
        assert_eq!(
            positions.get("chr2").unwrap()[36].depth,
            4 - awareness_modifier
        ); // mate overlap
        assert_eq!(
            positions.get("chr2").unwrap()[37].depth,
            4 - awareness_modifier
        ); // mate overlap
        assert_eq!(
            positions.get("chr2").unwrap()[38].depth,
            4 - awareness_modifier
        ); // mate overlap
        assert_eq!(
            positions.get("chr2").unwrap()[39].depth,
            2 - awareness_modifier
        ); // mate overlap
        assert_eq!(
            positions.get("chr2").unwrap()[40].depth,
            2 - awareness_modifier
        ); // mate overlap
        assert_eq!(
            positions.get("chr2").unwrap()[41].depth,
            2 - awareness_modifier
        ); // mate overlap
        assert_eq!(
            positions.get("chr2").unwrap()[42].depth,
            2 - awareness_modifier
        ); // mate overlap
        assert_eq!(
            positions.get("chr2").unwrap()[43].depth,
            2 - awareness_modifier
        ); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[44].depth, 1);

        // assert_eq!(positions.get("chr2").unwrap()[33].a, 4);
        // assert_eq!(positions.get("chr2").unwrap()[34].a, 4 - awareness_modifier); // mate overlap
        // assert_eq!(positions.get("chr2").unwrap()[35].a, 4 - awareness_modifier); // mate overlap
        // assert_eq!(positions.get("chr2").unwrap()[36].a, 4 - awareness_modifier); // mate overlap
        // assert_eq!(positions.get("chr2").unwrap()[37].a, 4 - awareness_modifier); // mate overlap
        // assert_eq!(positions.get("chr2").unwrap()[38].a, 4 - awareness_modifier); // mate overlap
        // assert_eq!(positions.get("chr2").unwrap()[39].a, 2 - awareness_modifier); // mate overlap
        // assert_eq!(positions.get("chr2").unwrap()[40].a, 2 - awareness_modifier); // mate overlap
        // assert_eq!(positions.get("chr2").unwrap()[41].a, 2 - awareness_modifier); // mate overlap
        // assert_eq!(positions.get("chr2").unwrap()[42].a, 2 - awareness_modifier); // mate overlap
        // assert_eq!(positions.get("chr2").unwrap()[43].a, 2 - awareness_modifier); // mate overlap
        // assert_eq!(positions.get("chr2").unwrap()[44].a, 1);
    }
}
