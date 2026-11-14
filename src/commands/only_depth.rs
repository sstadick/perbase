//! # Only Depth
//!
//! Calculates the depth only at each position. This uses the same algorithm described in
//! the [`mosdepth`](https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btx699/4583630?guestAccessKey=35b55064-4566-4ab3-a769-32916fa1c6e6)
//! paper.
use anyhow::Result;
use log::*;
use perbase_lib::{
    par_granges::{self, RegionProcessor},
    position::range_positions::RangePositions,
    read_filter::{DefaultReadFilter, ReadFilter},
    utils,
};
use rust_htslib::{bam, bam::ext::BamRecordExtensions};
use std::path::PathBuf;
use structopt::{StructOpt, clap::AppSettings};

mod stream;

/// Calculate only the depth at each base.
#[derive(StructOpt)]
#[structopt(author, setting = AppSettings::ArgRequiredElseHelp)]
pub struct OnlyDepth {
    /// Input indexed BAM/CRAM to analyze.
    reads: PathBuf,

    /// Indexed reference fasta, set if using CRAM.
    #[structopt(long, short = "r")]
    ref_fasta: Option<PathBuf>,

    /// A BED file containing regions of interest. If specified, only bases from the given regions will be reported on.
    #[structopt(long, short = "b")]
    bed_file: Option<PathBuf>,

    /// A BCF/VCF file containing positions of interest. If specified, only bases from the given positions will be reported on.
    /// Note that it may be more efficient to calculate depth over regions if your positions are clustered tightly together.
    #[structopt(long, short = "B")]
    bcf_file: Option<PathBuf>,

    /// Output path, defaults to stdout.
    #[structopt(long, short = "o")]
    output: Option<PathBuf>,

    /// Output BED-like output format with the depth in the 5th column. Note, `-z` can be used with this to change coordinates to
    /// 0-based to be more BED-like.
    #[structopt(long)]
    bed_format: bool,

    /// Optionally bgzip the output.
    #[structopt(long, short = "Z")]
    bgzip: bool,

    /// The number of threads to use.
    #[structopt(long, short = "t", default_value = utils::NUM_CPU.as_str())]
    threads: usize,

    /// The number of threads to use for compressing output (specified by --bgzip)
    #[structopt(long, short = "T", default_value = "4")]
    compression_threads: usize,

    /// The level to use for compressing output (specified by --bgzip)
    #[structopt(long, short = "L", default_value = "2")]
    compression_level: u32,

    /// The ideal number of basepairs each worker receives. Total bp in memory at one time is (threads - 2) * chunksize. With defaults, expect ~20 MB per worker thread.
    #[structopt(long, short = "c", default_value=par_granges::CHUNKSIZE_STR.as_str())]
    chunksize: u32,

    /// The fraction of a gigabyte to allocate per thread for message passing, can be greater than 1.0.
    #[structopt(long, short = "C", default_value = "0.001")]
    channel_size_modifier: f64,

    /// SAM flags to include.
    #[structopt(long, short = "f", default_value = "0")]
    include_flags: u16,

    /// SAM flags to exclude, recommended 3848.
    #[structopt(long, short = "F", default_value = "0")]
    exclude_flags: u16,

    /// Fix overlapping mates counts, see docs for full details.
    #[structopt(long, short = "m", conflicts_with = "fast_mode")]
    mate_fix: bool,

    /// Calculate depth based only on read starts/stops, see docs for full details.
    #[structopt(long, short = "x")]
    fast_mode: bool,

    /// Skip merging adjacent bases that have the same depth.
    #[structopt(long, short = "n")]
    no_merge: bool,

    /// Skip mergeing togther regions specified in the optional BED or BCF/VCF files.
    ///
    /// **NOTE** If this is set it could result in duplicate output entries for regions that overlap.
    /// **NOTE** This may cause issues with downstream tooling.
    #[structopt(long, short = "M")]
    skip_merging_intervals: bool,

    /// Keep positions even if they have 0 depth.
    #[structopt(long, short = "k")]
    keep_zeros: bool,

    /// Minimum MAPQ for a read to count toward depth.
    #[structopt(long, short = "q", default_value = "0")]
    min_mapq: u8,

    /// Output positions as 0-based instead of 1-based.
    #[structopt(long, short = "z")]
    zero_base: bool,

    /// Read alignments with Seqair instead of HTSlib.
    #[cfg(feature = "seqair")]
    #[structopt(long)]
    seqair: bool,
}

impl OnlyDepth {
    pub fn run(self) -> Result<()> {
        info!("Running only-depth on: {:?}", self.reads);
        let cpus = utils::determine_allowed_cpus(self.threads)?;

        let writer = utils::get_writer(
            &self.output,
            self.bgzip,
            !self.bed_format,
            self.compression_threads,
            self.compression_level,
        )?;

        let read_filter =
            DefaultReadFilter::new(self.include_flags, self.exclude_flags, self.min_mapq);
        let processor = OnlyDepthProcessor::new(
            self.reads.clone(),
            self.ref_fasta.clone(),
            self.mate_fix,
            self.fast_mode,
            self.no_merge,
            self.keep_zeros,
            if self.zero_base { 0 } else { 1 },
            read_filter,
        );

        let par_granges_runner = par_granges::ParGranges::new(
            self.reads.clone(),
            self.ref_fasta.clone(),
            self.bed_file.clone(),
            self.bcf_file.clone(),
            !self.skip_merging_intervals,
            Some(cpus),
            Some(self.chunksize),
            Some(self.channel_size_modifier),
            processor,
        );

        let (processor, regions) = par_granges_runner.into_regions()?;
        let mut output = writer.into_inner().map_err(|error| error.into_error())?;
        #[cfg(feature = "seqair")]
        let use_seqair = self.seqair;
        #[cfg(not(feature = "seqair"))]
        let use_seqair = false;
        stream::run(
            &processor,
            &regions,
            cpus.saturating_sub(2).max(1),
            self.chunksize,
            self.channel_size_modifier,
            self.bed_format,
            use_seqair,
            &mut output,
        )?;
        output.flush()?;
        Ok(())
    }

    /// Detect an overlap of mates
    /// TODO: Add more tests for this
    #[inline]
    fn maybe_overlaps_mate(record: &bam::Record) -> bool {
        if !record.is_paired() || record.tid() != record.mtid() {
            return false;
        }
        let tlen = record.insert_size().abs();
        let read_len = record.reference_end() - record.reference_start();
        // if the half the tlen is less than the readlen, or the mate start is less than one readlen away, concider for mate detection
        // (is there any scenario where this is not true?)
        // This is a heuristic, so it just have to never have a false negative and not have false positives > N% of the time
        (tlen / 2) < read_len || (record.mpos() - record.pos()).abs() < read_len
    }
}

/// Holds the info needed for [par_io::RegionProcessor] implementation
struct OnlyDepthProcessor {
    /// path to indexed BAM/CRAM
    reads: PathBuf,
    /// path to indexed ref file
    ref_fasta: Option<PathBuf>,
    /// Indicate whether or not to account for overlapping mates. Not checked in fastmode
    mate_fix: bool,
    /// Indicate whether or not to run in fastmode, only using read starts and stops
    fast_mode: bool,
    /// Indicate whether or not to merge adjacent positions that have the same depth
    no_merge: bool,
    /// Indicate whether or not to keep positions with 0 depth
    keep_zeros: bool,
    /// implementation of [position::ReadFilter] that will be used
    read_filter: DefaultReadFilter,
    /// 0-based or 1-based coordinate output
    coord_base: u32,
}

impl OnlyDepthProcessor {
    /// Create a new OnlyDepthProcessor
    #[allow(clippy::too_many_arguments)]
    fn new(
        reads: PathBuf,
        ref_fasta: Option<PathBuf>,
        mate_fix: bool,
        fast_mode: bool,
        no_merge: bool,
        keep_zeros: bool,
        coord_base: u32,
        read_filter: DefaultReadFilter,
    ) -> Self {
        Self {
            reads,
            ref_fasta,
            fast_mode,
            mate_fix,
            no_merge,
            keep_zeros,
            coord_base,
            read_filter,
        }
    }
}

/// Collect rows for callers of the generic region-processing interface.
/// The CLI uses the same counter but writes each row immediately.
impl RegionProcessor for OnlyDepthProcessor {
    type P = RangePositions;

    fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<RangePositions> {
        stream::collect(self, (tid, start, stop)).expect("count region")
    }
}

#[cfg(test)]
#[allow(unused)]
mod tests {
    use super::*;
    use perbase_lib::{
        position::{Position, range_positions::RangePositions},
        read_filter::DefaultReadFilter,
    };
    use rstest::*;
    use rust_htslib::{bam, bam::record::Record};
    use smartstring::alias::*;
    use std::{collections::HashMap, path::PathBuf};
    use tempfile::{TempDir, tempdir};

    #[test]
    fn streaming_depth_matches_cigar_coverage() {
        let directory = tempdir().unwrap();
        let reads = directory.path().join("reads.bam");
        let output = directory.path().join("depth.tsv");
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "a")
                .push_tag(b"LN", 12),
        );
        let view = bam::HeaderView::from_header(&header);
        let mut writer = bam::Writer::from_path(&reads, &header, bam::Format::Bam).unwrap();
        for sam in [
            "pair\t65\ta\t1\t60\t4M3N4M\t=\t3\t11\tAAAAAAAA\tIIIIIIII",
            "pair\t129\ta\t3\t60\t6M\t=\t1\t-11\tAAAAAA\tIIIIII",
            "filtered\t512\ta\t4\t60\t2M\t*\t0\t0\tAA\tII",
            "single\t0\ta\t6\t60\t2M2D2M\t*\t0\t0\tAAAA\tIIII",
        ] {
            writer
                .write(&Record::from_sam(&view, sam.as_bytes()).unwrap())
                .unwrap();
        }
        drop(writer);
        bam::index::build(&reads, None, bam::index::Type::Bai, 1).unwrap();
        for (mode, depths) in [
            (None, [1, 1, 2, 2, 1, 2, 2, 3, 2, 2, 2, 0]),
            (Some("-m"), [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0]),
            (Some("-x"), [1, 1, 2, 2, 2, 3, 3, 3, 2, 2, 2, 0]),
        ] {
            for bed in [false, true] {
                for seqair in [false, cfg!(feature = "seqair")] {
                    let mut args = vec![
                        "only-depth",
                        reads.to_str().unwrap(),
                        "-n",
                        "-k",
                        "-z",
                        "-F",
                        "512",
                        "-c",
                        "12",
                        "-t",
                        "4",
                        "-o",
                        output.to_str().unwrap(),
                    ];
                    args.extend(mode);
                    if bed {
                        args.push("--bed-format");
                    }
                    if seqair {
                        args.push("--seqair");
                    }
                    OnlyDepth::from_iter_safe(args).unwrap().run().unwrap();
                    let mut expected = if bed {
                        std::string::String::new()
                    } else {
                        "REF\tPOS\tEND\tDEPTH\n".to_owned()
                    };
                    for (pos, depth) in depths.iter().enumerate() {
                        expected.push_str(&format!(
                            "a\t{pos}\t{}\t{}{depth}\n",
                            pos + 1,
                            if bed { "\t" } else { "" }
                        ));
                    }
                    assert_eq!(
                        std::fs::read_to_string(&output).unwrap(),
                        expected,
                        "mode={mode:?}, bed={bed}, seqair={seqair}"
                    );
                }
            }
        }
    }

    #[fixture]
    fn read_filter() -> DefaultReadFilter {
        DefaultReadFilter::new(0, 512, 0)
    }

    #[fixture]
    fn bamfile() -> (PathBuf, TempDir) {
        // No longer - This keep the test bam up to date
        // let path = PathBuf::from("test/test.bam");
        let tempdir = tempdir().unwrap();
        let path = tempdir.path().join("test.bam");

        // Build a header
        let mut header = bam::header::Header::new();
        let mut chr1 = bam::header::HeaderRecord::new(b"SQ");
        chr1.push_tag(b"SN", &"chr1".to_owned());
        chr1.push_tag(b"LN", &"3000000".to_owned());
        let mut chr2 = bam::header::HeaderRecord::new(b"SQ");
        chr2.push_tag(b"SN", &"chr2".to_owned());
        chr2.push_tag(b"LN", &"3000000".to_owned());
        let mut chr3 = bam::header::HeaderRecord::new(b"SQ");
        chr3.push_tag(b"SN", &"chr3".to_owned());
        chr3.push_tag(b"LN", &"3000000".to_owned());
        header.push_record(&chr1);
        header.push_record(&chr2);
        header.push_record(&chr3);
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
            // Base on last position of chunk should verify that we aren't counting the last position by not throwing an index error
            Record::from_sam(&view, b"FIVE\t147\tchr1\t999976\t40\t25M\tchr1\t1000026\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            Record::from_sam(&view, b"FIVE\t147\tchr1\t1000026\t40\t25M\tchr1\t999976\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),

            // Chr2 - Complex
            // Ins
            Record::from_sam(&view, b"ONE\t67\tchr2\t1\t40\t2M2I21M\tchr2\t50\t75\tAAGGGAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            // Del
            Record::from_sam(&view, b"TWO\t67\tchr2\t5\t40\t2M5D23M\tchr2\t55\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            // Skip
            Record::from_sam(&view, b"THREE\t67\tchr2\t10\t40\t3M5N22M\tchr2\t60\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            // Mismatch
            Record::from_sam(&view, b"FOUR\t67\tchr2\t15\t40\t25M\tchr2\t65\t75\tATAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            // Overlapping mates a then b
            Record::from_sam(&view, b"FIVE\t67\tchr2\t20\t40\t25M\tchr2\t35\t40\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            Record::from_sam(&view, b"FIVE\t147\tchr2\t35\t40\t25M\tchr2\t20\t-40\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),


            // Other base
            Record::from_sam(&view, b"ONE\t147\tchr2\t50\t40\t25M\tchr2\t1\t75\tAAAAAAAAAAAAAAAAAAAAAYAAA\t#########################").unwrap(),

            Record::from_sam(&view, b"TWO\t147\tchr2\t55\t40\t25M\tchr2\t5\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            Record::from_sam(&view, b"THREE\t147\tchr2\t60\t40\t25M\tchr2\t10\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            // A failure of QC
            Record::from_sam(&view, b"FOUR\t659\tchr2\t65\t40\t25M\tchr2\t15\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),

            // Overlapping mates exact overlap - calling proper pair, dubious
            Record::from_sam(&view, b"SIX\t67\tchr2\t120\t40\t25M\tchr2\t120\t40\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            Record::from_sam(&view, b"SIX\t147\tchr2\t120\t40\t25M\tchr2\t120\t-40\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),

            // Overlapping mates a contains b - calling proper pair, dubious
            Record::from_sam(&view, b"SEVEN\t67\tchr2\t160\t40\t50M\tchr2\t150\t50\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t##################################################").unwrap(),
            Record::from_sam(&view, b"SEVEN\t147\tchr2\t175\t40\t25M\tchr2\t175\t-50\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),

            // Chr3, weird inter-seq breakpoints
            // Huge skip at start
            Record::from_sam(&view, b"ONE\t67\tchr3\t1\t40\t2M2000000N23M\tchr3\t50\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            Record::from_sam(&view, b"TWO\t67\tchr3\t5\t40\t25M\tchr3\t55\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            // Huge skip at end
            Record::from_sam(&view, b"THREE\t67\tchr3\t10\t40\t23M2000000N2M\tchr3\t60\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            Record::from_sam(&view, b"FOUR\t67\tchr3\t15\t40\t25M\tchr3\t65\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            Record::from_sam(&view, b"FIVE\t67\tchr3\t20\t40\t25M\tchr3\t70\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
            Record::from_sam(&view, b"ONE\t147\tchr3\t50\t40\t25M\tchr3\t1\t75\tTTTTTTTTTTTTTTTTTTTTTTTTT\t#########################").unwrap(),
            Record::from_sam(&view, b"TWO\t147\tchr3\t55\t40\t25M\tchr3\t5\t75\tGGGGGGGGGGGGGGGGGGGGGGGGG\t#########################").unwrap(),
            Record::from_sam(&view, b"THREE\t147\tchr3\t60\t40\t25M\tchr3\t10\t75\tCCCCCCCCCCCCCCCCCCCCCCCCC\t#########################").unwrap(),
            Record::from_sam(&view, b"FOUR\t147\tchr3\t65\t40\t25M\tchr3\t15\t75\tNNNNNNNNNNNNNNNNNNNNNNNNN\t#########################").unwrap(),
            Record::from_sam(&view, b"FIVE\t147\tchr3\t70\t40\t25M\tchr3\t20\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
        ];

        // Update the test/test.bam file
        let mut writer =
            bam::Writer::from_path(&path, &header, bam::Format::Bam).expect("Created writer");
        for record in records.iter() {
            writer.write(record).expect("Wrote record");
        }
        drop(writer); // force it to flush so indexing can happen
        // build the index
        bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();
        (path, tempdir)
    }

    // Test that all regions of the test bam can be read and don't panic.
    #[rstest(
        fast_mode => [true, false],
        mate_fix => [true, false]
    )]
    fn test_can_parse(
        fast_mode: bool,
        mate_fix: bool,
        bamfile: (PathBuf, TempDir),
        read_filter: DefaultReadFilter,
    ) {
        let cpus = utils::determine_allowed_cpus(8).unwrap();

        let onlydepth_processor = OnlyDepthProcessor::new(
            bamfile.0.clone(),
            None,
            mate_fix,
            fast_mode,
            false,
            false,
            0,
            read_filter,
        );

        let par_granges_runner = par_granges::ParGranges::new(
            bamfile.0,
            None,
            None,
            None,
            true,
            Some(cpus),
            Some(1_000_000),
            Some(0.001),
            onlydepth_processor,
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
        assert_eq!(positions.keys().len(), 3);
    }

    #[fixture]
    fn vanilla_positions(
        bamfile: (PathBuf, TempDir),
        read_filter: DefaultReadFilter,
    ) -> HashMap<String, Vec<RangePositions>> {
        let cpus = utils::determine_allowed_cpus(8).unwrap();

        let onlydepth_processor = OnlyDepthProcessor::new(
            bamfile.0.clone(),
            None,
            false,
            false,
            false,
            false,
            0,
            read_filter,
        );

        let par_granges_runner = par_granges::ParGranges::new(
            bamfile.0,
            None,
            None,
            None,
            true,
            Some(cpus),
            Some(1_000_000),
            Some(0.001),
            onlydepth_processor,
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
    fn vanilla_positions_mate_fix(
        bamfile: (PathBuf, TempDir),
        read_filter: DefaultReadFilter,
    ) -> HashMap<String, Vec<RangePositions>> {
        let cpus = utils::determine_allowed_cpus(8).unwrap();

        let onlydepth_processor = OnlyDepthProcessor::new(
            bamfile.0.clone(),
            None,
            true,
            false,
            false,
            false,
            0,
            read_filter,
        );

        let par_granges_runner = par_granges::ParGranges::new(
            bamfile.0,
            None,
            None,
            None,
            true,
            Some(cpus),
            Some(1_000_000),
            Some(0.001),
            onlydepth_processor,
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
    fn fast_mode_positions(
        bamfile: (PathBuf, TempDir),
        read_filter: DefaultReadFilter,
    ) -> HashMap<String, Vec<RangePositions>> {
        let cpus = utils::determine_allowed_cpus(8).unwrap();

        let onlydepth_processor = OnlyDepthProcessor::new(
            bamfile.0.clone(),
            None,
            false,
            true,
            false,
            false,
            0,
            read_filter,
        );

        let par_granges_runner = par_granges::ParGranges::new(
            bamfile.0,
            None,
            None,
            None,
            true,
            Some(cpus),
            None,
            Some(0.001),
            onlydepth_processor,
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
    fn fast_mode_positions_mate_fix(
        bamfile: (PathBuf, TempDir),
        read_filter: DefaultReadFilter,
    ) -> HashMap<String, Vec<RangePositions>> {
        let cpus = utils::determine_allowed_cpus(8).unwrap();

        let onlydepth_processor = OnlyDepthProcessor::new(
            bamfile.0.clone(),
            None,
            true,
            true,
            false,
            false,
            0,
            read_filter,
        );

        let par_granges_runner = par_granges::ParGranges::new(
            bamfile.0,
            None,
            None,
            None,
            true,
            Some(cpus),
            None,
            Some(0.001),
            onlydepth_processor,
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

    #[rstest(
        positions,
        awareness_modifier,
        case::vanilla(vanilla_positions(bamfile(), read_filter()), 0),
        case::fast_mode(fast_mode_positions(bamfile(), read_filter()), 0),
        case::vanilla_mate_fix(vanilla_positions_mate_fix(bamfile(), read_filter()), 0),
        case::fast_mode_mate_fix(fast_mode_positions_mate_fix(bamfile(), read_filter()), 0)
    )]
    fn check_depths(positions: HashMap<String, Vec<RangePositions>>, awareness_modifier: usize) {
        assert_eq!(positions.get("chr1").unwrap()[0].depth, 1);
        assert_eq!(positions.get("chr1").unwrap()[1].depth, 2);
        assert_eq!(positions.get("chr1").unwrap()[2].depth, 3);
        assert_eq!(positions.get("chr1").unwrap()[3].depth, 4);
        assert_eq!(positions.get("chr1").unwrap()[4].depth, 5);
        assert_eq!(positions.get("chr1").unwrap()[5].depth, 4);
        assert_eq!(positions.get("chr1").unwrap()[6].depth, 3);
        assert_eq!(positions.get("chr1").unwrap()[7].depth, 2);
        assert_eq!(positions.get("chr1").unwrap()[8].depth, 1);
        assert_eq!(positions.get("chr1").unwrap()[9].depth, 0);
        assert_eq!(positions.get("chr1").unwrap()[13].depth, 4);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::vanilla(vanilla_positions(bamfile(), read_filter()), 0),
        case::fast_mode(fast_mode_positions(bamfile(), read_filter()), 0),
        case::vanilla_mate_fix(vanilla_positions_mate_fix(bamfile(), read_filter()), 0),
        case::fast_mode_mate_fix(fast_mode_positions_mate_fix(bamfile(), read_filter()), 0)
    )]
    fn check_ranges(positions: HashMap<String, Vec<RangePositions>>, awareness_modifier: usize) {
        assert_eq!(positions.get("chr1").unwrap()[0].pos, 0);
        assert_eq!(positions.get("chr1").unwrap()[1].pos, 4);
        assert_eq!(positions.get("chr1").unwrap()[2].pos, 9);
        assert_eq!(positions.get("chr1").unwrap()[3].pos, 14);
        assert_eq!(positions.get("chr1").unwrap()[4].pos, 19);
        assert_eq!(positions.get("chr1").unwrap()[5].pos, 25);
        assert_eq!(positions.get("chr1").unwrap()[6].pos, 29);
        assert_eq!(positions.get("chr1").unwrap()[7].pos, 34);
        assert_eq!(positions.get("chr1").unwrap()[8].pos, 39);
        assert_eq!(positions.get("chr1").unwrap()[9].pos, 44);
        assert_eq!(positions.get("chr1").unwrap()[13].pos, 64);

        assert_eq!(positions.get("chr1").unwrap()[0].end, 4);
        assert_eq!(positions.get("chr1").unwrap()[1].end, 9);
        assert_eq!(positions.get("chr1").unwrap()[2].end, 14);
        assert_eq!(positions.get("chr1").unwrap()[3].end, 19);
        assert_eq!(positions.get("chr1").unwrap()[4].end, 25);
        assert_eq!(positions.get("chr1").unwrap()[5].end, 29);
        assert_eq!(positions.get("chr1").unwrap()[6].end, 34);
        assert_eq!(positions.get("chr1").unwrap()[7].end, 39);
        assert_eq!(positions.get("chr1").unwrap()[8].end, 44);
        assert_eq!(positions.get("chr1").unwrap()[9].end, 49);
        assert_eq!(positions.get("chr1").unwrap()[13].end, 69);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::vanilla(vanilla_positions(bamfile(), read_filter()), 2),
        case::fast_mode(fast_mode_positions(bamfile(), read_filter()), 0),
        case::vanilla_mate_fix(vanilla_positions_mate_fix(bamfile(), read_filter()), 2),
        case::fast_mode_mate_fix(fast_mode_positions_mate_fix(bamfile(), read_filter()), 0)
    )]
    fn check_filters(positions: HashMap<String, Vec<RangePositions>>, awareness_modifier: usize) {
        // Verify that a read that has flags saying it failed QC got filtered out
        assert_eq!(
            positions.get("chr2").unwrap()[11 + awareness_modifier].depth,
            1
        );
        assert_eq!(
            positions.get("chr2").unwrap()[12 + awareness_modifier].depth,
            0
        );
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::vanilla(vanilla_positions(bamfile(), read_filter()), 0),
        case::fast_mode(fast_mode_positions(bamfile(), read_filter()), 0),
        case::vanilla_mate_fix(vanilla_positions_mate_fix(bamfile(), read_filter()), 0),
        case::fast_mode_mate_fix(fast_mode_positions_mate_fix(bamfile(), read_filter()), 0)
    )]
    fn check_depths_insertions(
        positions: HashMap<String, Vec<RangePositions>>,
        awareness_modifier: usize,
    ) {
        // insertion in this range
        assert_eq!(positions.get("chr2").unwrap()[0].depth, 1);
        assert_eq!(positions.get("chr2").unwrap()[0].pos, 0);
        assert_eq!(positions.get("chr2").unwrap()[0].end, 4);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::vanilla(vanilla_positions(bamfile(), read_filter()), 0),
        case::fast_mode(fast_mode_positions(bamfile(), read_filter()), 0),
        case::vanilla_mate_fix(vanilla_positions_mate_fix(bamfile(), read_filter()), 0),
        case::fast_mode_mate_fix(fast_mode_positions_mate_fix(bamfile(), read_filter()), 0)
    )]
    fn check_depths_deletions(
        positions: HashMap<String, Vec<RangePositions>>,
        awareness_modifier: usize,
    ) {
        assert_eq!(positions.get("chr2").unwrap()[1].depth, 2); // dels in this range
        assert_eq!(positions.get("chr2").unwrap()[1].pos, 4);
        assert_eq!(positions.get("chr2").unwrap()[1].end, 9);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::vanilla(vanilla_positions(bamfile(), read_filter()), 1),
        case::fast_mode(fast_mode_positions(bamfile(), read_filter()), 0),
        case::vanilla_mate_fix(vanilla_positions_mate_fix(bamfile(), read_filter()), 1),
        case::fast_mode_mate_fix(fast_mode_positions_mate_fix(bamfile(), read_filter()), 0)
    )]
    fn check_depths_refskips(
        positions: HashMap<String, Vec<RangePositions>>,
        awareness_modifier: usize,
    ) {
        if awareness_modifier == 1 {
            // ref_skips are not being counted
            assert_eq!(positions.get("chr2").unwrap()[2].depth, 3); // noskips
            assert_eq!(positions.get("chr2").unwrap()[2].pos, 9);
            assert_eq!(positions.get("chr2").unwrap()[2].end, 12);
            assert_eq!(positions.get("chr2").unwrap()[3].depth, 2); // skips
            assert_eq!(positions.get("chr2").unwrap()[3].pos, 12);
            assert_eq!(positions.get("chr2").unwrap()[3].end, 14);
        } else {
            assert_eq!(positions.get("chr2").unwrap()[2].depth, 3); // Skips in this range
            assert_eq!(positions.get("chr2").unwrap()[2].pos, 9);
            assert_eq!(positions.get("chr2").unwrap()[2].end, 14);
        }
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::fast_mode(fast_mode_positions(bamfile(), read_filter()), 0),
        case::fast_mode_mate_fix(fast_mode_positions_mate_fix(bamfile(), read_filter()), 1)
    )]
    fn check_mate_detection_fast(
        positions: HashMap<String, Vec<RangePositions>>,
        awareness_modifier: usize,
    ) {
        assert_eq!(positions.get("chr2").unwrap()[5].depth, 4);
        assert_eq!(positions.get("chr2").unwrap()[5].pos, 23);
        if awareness_modifier == 1 {
            assert_eq!(positions.get("chr2").unwrap()[5].end, 34);

            assert_eq!(positions.get("chr2").unwrap()[6].pos, 34);
            assert_eq!(positions.get("chr2").unwrap()[6].end, 39);
            assert_eq!(positions.get("chr2").unwrap()[6].depth, 3);

            assert_eq!(positions.get("chr2").unwrap()[7].pos, 39);
            assert_eq!(positions.get("chr2").unwrap()[7].end, 49);
            assert_eq!(positions.get("chr2").unwrap()[7].depth, 1);

            assert_eq!(positions.get("chr2").unwrap()[13].pos, 119);
            assert_eq!(positions.get("chr2").unwrap()[13].end, 144);
            assert_eq!(positions.get("chr2").unwrap()[13].depth, 1);

            assert_eq!(positions.get("chr2").unwrap()[15].pos, 159);
            assert_eq!(positions.get("chr2").unwrap()[15].end, 209);
            assert_eq!(positions.get("chr2").unwrap()[15].depth, 1);
        } else {
            assert_eq!(positions.get("chr2").unwrap()[5].end, 39);

            assert_eq!(positions.get("chr2").unwrap()[6].pos, 39);
            assert_eq!(positions.get("chr2").unwrap()[6].end, 44);
            assert_eq!(positions.get("chr2").unwrap()[6].depth, 2);

            assert_eq!(positions.get("chr2").unwrap()[7].pos, 44);
            assert_eq!(positions.get("chr2").unwrap()[7].end, 49);
            assert_eq!(positions.get("chr2").unwrap()[7].depth, 1);

            assert_eq!(positions.get("chr2").unwrap()[13].pos, 119);
            assert_eq!(positions.get("chr2").unwrap()[13].end, 144);
            assert_eq!(positions.get("chr2").unwrap()[13].depth, 2);

            assert_eq!(positions.get("chr2").unwrap()[15].pos, 159);
            assert_eq!(positions.get("chr2").unwrap()[15].end, 174);
            assert_eq!(positions.get("chr2").unwrap()[15].depth, 1);

            assert_eq!(positions.get("chr2").unwrap()[16].pos, 174);
            assert_eq!(positions.get("chr2").unwrap()[16].end, 199);
            assert_eq!(positions.get("chr2").unwrap()[16].depth, 2);

            assert_eq!(positions.get("chr2").unwrap()[17].pos, 199);
            assert_eq!(positions.get("chr2").unwrap()[17].end, 209);
            assert_eq!(positions.get("chr2").unwrap()[17].depth, 1);
        }
        assert_eq!(positions.get("chr2").unwrap()[8].pos, 49);
        assert_eq!(positions.get("chr2").unwrap()[8].end, 54);
        assert_eq!(positions.get("chr2").unwrap()[8].depth, 2);
    }
    #[rstest(
        positions,
        awareness_modifier,
        case::vanilla(vanilla_positions(bamfile(), read_filter()), 0),
        case::vanilla_mate_fix(vanilla_positions_mate_fix(bamfile(), read_filter()), 1)
    )]
    fn check_mate_detection_vanilla(
        positions: HashMap<String, Vec<RangePositions>>,
        awareness_modifier: usize,
    ) {
        assert_eq!(positions.get("chr2").unwrap()[7].depth, 4);
        assert_eq!(positions.get("chr2").unwrap()[7].pos, 23);
        if awareness_modifier == 1 {
            assert_eq!(positions.get("chr2").unwrap()[7].end, 34);

            assert_eq!(positions.get("chr2").unwrap()[8].pos, 34);
            assert_eq!(positions.get("chr2").unwrap()[8].end, 39);
            assert_eq!(positions.get("chr2").unwrap()[8].depth, 3);

            assert_eq!(positions.get("chr2").unwrap()[9].pos, 39);
            assert_eq!(positions.get("chr2").unwrap()[9].end, 49);
            assert_eq!(positions.get("chr2").unwrap()[9].depth, 1);

            assert_eq!(positions.get("chr2").unwrap()[15].pos, 119);
            assert_eq!(positions.get("chr2").unwrap()[15].end, 144);
            assert_eq!(positions.get("chr2").unwrap()[15].depth, 1);

            assert_eq!(positions.get("chr2").unwrap()[17].pos, 159);
            assert_eq!(positions.get("chr2").unwrap()[17].end, 209);
            assert_eq!(positions.get("chr2").unwrap()[17].depth, 1)
        } else {
            assert_eq!(positions.get("chr2").unwrap()[7].end, 39);

            assert_eq!(positions.get("chr2").unwrap()[8].pos, 39);
            assert_eq!(positions.get("chr2").unwrap()[8].end, 44);
            assert_eq!(positions.get("chr2").unwrap()[8].depth, 2);

            assert_eq!(positions.get("chr2").unwrap()[9].pos, 44);
            assert_eq!(positions.get("chr2").unwrap()[9].end, 49);
            assert_eq!(positions.get("chr2").unwrap()[9].depth, 1);

            assert_eq!(positions.get("chr2").unwrap()[15].pos, 119);
            assert_eq!(positions.get("chr2").unwrap()[15].end, 144);
            assert_eq!(positions.get("chr2").unwrap()[15].depth, 2);

            assert_eq!(positions.get("chr2").unwrap()[17].pos, 159);
            assert_eq!(positions.get("chr2").unwrap()[17].end, 174);
            assert_eq!(positions.get("chr2").unwrap()[17].depth, 1);

            assert_eq!(positions.get("chr2").unwrap()[18].pos, 174);
            assert_eq!(positions.get("chr2").unwrap()[18].end, 199);
            assert_eq!(positions.get("chr2").unwrap()[18].depth, 2);

            assert_eq!(positions.get("chr2").unwrap()[19].pos, 199);
            assert_eq!(positions.get("chr2").unwrap()[19].end, 209);
            assert_eq!(positions.get("chr2").unwrap()[19].depth, 1)
        }
        assert_eq!(positions.get("chr2").unwrap()[10].pos, 49);
        assert_eq!(positions.get("chr2").unwrap()[10].end, 54);
        assert_eq!(positions.get("chr2").unwrap()[10].depth, 2);
    }

    // Test that reads that overlap chunks are counted correctly.
    #[rstest(
        positions,
        awareness_modifier,
        case::vanilla(vanilla_positions(bamfile(), read_filter()), 0),
        case::fast_mode(fast_mode_positions(bamfile(), read_filter()), 2),
        case::vanilla_mate_fix(vanilla_positions_mate_fix(bamfile(), read_filter()), 0),
        case::fast_mode_mate_fix(fast_mode_positions_mate_fix(bamfile(), read_filter()), 2)
    )]
    fn check_chunk_ends(
        positions: HashMap<String, Vec<RangePositions>>,
        awareness_modifier: usize,
    ) {
        if awareness_modifier == 0 {
            // Aware of refskips
            assert_eq!(positions.get("chr3").unwrap()[19].pos, 94);
            assert_eq!(positions.get("chr3").unwrap()[19].end, 1_000_000);
            assert_eq!(positions.get("chr3").unwrap()[19].depth, 0);
            assert_eq!(positions.get("chr3").unwrap()[20].pos, 1_000_000);
            assert_eq!(positions.get("chr3").unwrap()[20].end, 2_000_000);
            assert_eq!(positions.get("chr3").unwrap()[20].depth, 0);
            assert_eq!(positions.get("chr3").unwrap()[21].pos, 2_000_000);
            assert_eq!(positions.get("chr3").unwrap()[21].end, 2_000_002);
            assert_eq!(positions.get("chr3").unwrap()[21].depth, 0);
            assert_eq!(positions.get("chr3").unwrap()[22].pos, 2_000_002);
            assert_eq!(positions.get("chr3").unwrap()[22].end, 2_000_025);
            assert_eq!(positions.get("chr3").unwrap()[22].depth, 1);
            assert_eq!(positions.get("chr3").unwrap()[23].pos, 2_000_025);
            assert_eq!(positions.get("chr3").unwrap()[23].end, 2_000_032);
            assert_eq!(positions.get("chr3").unwrap()[23].depth, 0);
            assert_eq!(positions.get("chr3").unwrap()[24].pos, 2_000_032);
            assert_eq!(positions.get("chr3").unwrap()[24].end, 2_000_034);
            assert_eq!(positions.get("chr3").unwrap()[24].depth, 1);
            assert_eq!(positions.get("chr3").unwrap()[25].pos, 2_000_034);
            assert_eq!(positions.get("chr3").unwrap()[25].end, 3_000_000);
            assert_eq!(positions.get("chr3").unwrap()[25].depth, 0);
        } else {
            // only counting read start / ends
            assert_eq!(positions.get("chr3").unwrap()[17].pos, 94);
            assert_eq!(positions.get("chr3").unwrap()[17].end, 1_000_000);
            assert_eq!(positions.get("chr3").unwrap()[17].depth, 2);
            assert_eq!(positions.get("chr3").unwrap()[18].pos, 1_000_000);
            assert_eq!(positions.get("chr3").unwrap()[18].end, 2_000_000);
            assert_eq!(positions.get("chr3").unwrap()[18].depth, 2);
            assert_eq!(positions.get("chr3").unwrap()[19].pos, 2_000_000);
            assert_eq!(positions.get("chr3").unwrap()[19].end, 2_000_025);
            assert_eq!(positions.get("chr3").unwrap()[19].depth, 2);
            assert_eq!(positions.get("chr3").unwrap()[20].pos, 2_000_025);
            assert_eq!(positions.get("chr3").unwrap()[20].end, 2_000_034);
            assert_eq!(positions.get("chr3").unwrap()[20].depth, 1);
            assert_eq!(positions.get("chr3").unwrap()[21].pos, 2_000_034);
            assert_eq!(positions.get("chr3").unwrap()[21].end, 3_000_000);
            assert_eq!(positions.get("chr3").unwrap()[21].depth, 0);
        }
    }
}
