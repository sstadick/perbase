//! # Base Depth
//!
//! Base single pass over bam file to calculate depth at each position
//! as well as depth per nucleotide. Additionally counts the number of
//! insertions / deletions at each position.
use anyhow::Result;
use bio::io::fasta::IndexedReader;
use log::*;
use perbase_lib::{
    par_granges::{self, RegionProcessor},
    position::{Position, mate_fix::MateResolutionStrategy, pileup_position::PileupPosition},
    read_filter::{DefaultReadFilter, ReadFilter},
    reference, utils,
};
use rust_htslib::{bam, bam::Read};
use std::{convert::TryInto, path::PathBuf};
use structopt::StructOpt;

/// Calculate the depth at each base, per-nucleotide.
#[derive(StructOpt)]
#[structopt(author)]
pub struct BaseDepth {
    /// Input indexed BAM/CRAM to analyze.
    reads: PathBuf,

    /// Indexed reference fasta, set if using CRAM.
    #[structopt(long, short = "r")]
    ref_fasta: Option<PathBuf>,

    /// A BED file containing regions of interest. If specified, only bases from the given regions will be reported on.
    #[structopt(long, short = "b")]
    bed_file: Option<PathBuf>,

    /// A BCF/VCF file containing positions of interest. If specified, only bases from the given positions will be reported on.
    #[structopt(long, short = "B")]
    bcf_file: Option<PathBuf>,

    /// Output path, defaults to stdout.
    #[structopt(long, short = "o")]
    output: Option<PathBuf>,

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

    /// The ideal number of basepairs each worker receives. Total bp in memory at one time is (threads - 2) * chunksize.
    #[structopt(long, short = "c", default_value=par_granges::CHUNKSIZE_STR.as_str())]
    chunksize: u32,

    // NB: If there are large regions of low coverage, bumping this up may be helpful.
    /// The fraction of a gigabyte to allocate per thread for message passing, can be greater than 1.0.
    #[structopt(long, short = "C", default_value = "0.15")]
    channel_size_modifier: f64,

    /// SAM flags to include.
    #[structopt(long, short = "f", default_value = "0")]
    include_flags: u16,

    /// SAM flags to exclude, recommended 3848.
    #[structopt(long, short = "F", default_value = "0")]
    exclude_flags: u16,

    /// Fix overlapping mates counts, see docs for full details.
    #[structopt(long, short = "m")]
    mate_fix: bool,

    /// If `mate_fix` is true, select the method to use for mate fixing.
    #[structopt(long, short = "M", default_value = "original")]
    mate_resolution_strategy: MateResolutionStrategy,

    /// Keep positions even if they have 0 depth.
    #[structopt(long, short = "k")]
    keep_zeros: bool,

    /// Skip merging togther regions specified in the optional BED or BCF/VCF files.
    ///
    /// **NOTE** If this is set it could result in duplicate output entries for regions that overlap.
    /// **NOTE** This may cause issues with downstream tooling.
    #[structopt(long, short = "M")]
    skip_merging_intervals: bool,

    /// Minimum MAPQ for a read to count toward depth.
    #[structopt(long, short = "q", default_value = "0")]
    min_mapq: u8,

    /// Minium base quality for a base to be counted toward [A, C, T, G]. If the base is less than the specified
    /// quality score it will instead be counted as an `N`. If nothing is set for this no cutoff will be applied.
    #[structopt(long, short = "Q")]
    min_base_quality_score: Option<u8>,

    /// Output positions as 0-based instead of 1-based.
    #[structopt(long, short = "z")]
    zero_base: bool,

    /// Number of Reference Sequences to hold in memory at one time. Smaller will decrease mem usage.
    #[structopt(long, default_value = "10")]
    ref_cache_size: usize,

    /// Set the max depth for a pileup. If a positions depth is within 1% of max-depth the `NEAR_MAX_DEPTH`
    /// output field will be set to true and that position should be viewed as suspect.
    #[structopt(long, short = "D", default_value = "100000")]
    max_depth: u32,
}

impl BaseDepth {
    pub fn run(self) -> Result<()> {
        info!("Running base-depth on: {:?}", self.reads);
        let cpus = utils::determine_allowed_cpus(self.threads)?;

        let mut writer = utils::get_writer(
            &self.output,
            self.bgzip,
            true,
            self.compression_threads,
            self.compression_level,
        )?;

        let read_filter =
            DefaultReadFilter::new(self.include_flags, self.exclude_flags, self.min_mapq);
        let base_processor = BaseProcessor::new(
            self.reads.clone(),
            self.ref_fasta.clone(),
            self.mate_fix,
            self.mate_resolution_strategy,
            self.keep_zeros,
            if self.zero_base { 0 } else { 1 },
            read_filter,
            self.max_depth,
            self.ref_cache_size,
            self.min_base_quality_score,
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
            base_processor,
        );

        let receiver = par_granges_runner.process()?;

        for pos in receiver.into_iter() {
            writer.serialize(pos)?
        }

        writer.flush()?;
        Ok(())
    }
}

/// Holds the info needed for [par_io::RegionProcessor] implementation
struct BaseProcessor<F: ReadFilter> {
    /// path to indexed BAM/CRAM
    reads: PathBuf,
    /// path to indexed ref file
    ref_fasta: Option<PathBuf>,
    /// Cached access to the ref_fasta
    ref_buffer: Option<reference::Buffer>,
    /// Indicate whether or not to account for overlapping mates.
    mate_fix: bool,
    /// Strategy to use for mate fixing
    mate_fix_strategy: MateResolutionStrategy,
    /// Indicate whether or not to keep postitions that have zero depth
    keep_zeros: bool,
    /// 0-based or 1-based coordiante output
    coord_base: u32,
    /// implementation of [position::ReadFilter] that will be used
    read_filter: F,
    /// max depth to pass to htslib pileup engine, max value is MAX(i32)
    max_depth: u32,
    /// the cutoff at which we start logging warnings about depth being close to max depth
    max_depth_warnings_cutoff: u32,
    /// an optional base quality score. If Some(number) if the base quality is not >= that number the base is treated as an `N`
    min_base_quality_score: Option<u8>,
}

impl<F: ReadFilter> BaseProcessor<F> {
    /// Create a new BaseProcessor
    #[allow(clippy::too_many_arguments)]
    fn new(
        reads: PathBuf,
        ref_fasta: Option<PathBuf>,
        mate_fix: bool,
        mate_fix_strategy: MateResolutionStrategy,
        keep_zeros: bool,
        coord_base: u32,
        read_filter: F,
        max_depth: u32,
        ref_buffer_capacity: usize,
        min_base_quality_score: Option<u8>,
    ) -> Self {
        let ref_buffer = ref_fasta.as_ref().map(|ref_fasta| {
            reference::Buffer::new(
                IndexedReader::from_file(ref_fasta).expect("Reading Indexed FASTA"),
                ref_buffer_capacity,
            )
        });
        Self {
            reads,
            ref_fasta,
            ref_buffer,
            mate_fix,
            mate_fix_strategy,
            keep_zeros,
            coord_base,
            read_filter,
            max_depth,
            // Set cutoff to 1% of whatever max_depth is.
            max_depth_warnings_cutoff: max_depth - (max_depth as f64 * 0.01) as u32,
            min_base_quality_score,
        }
    }
}

/// Implement [par_io::RegionProcessor] for [BaseProcessor]
impl<F: ReadFilter> RegionProcessor for BaseProcessor<F> {
    /// Objects of [pipeup_position::PileupPosition] will be returned by each call to [BaseProcessor::process_region]
    type P = PileupPosition;

    /// Process a region by fetching it from a BAM/CRAM, getting a pileup, and then
    /// walking the pileup (checking bounds) to create Position objects according to
    /// the defined filters
    fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<PileupPosition> {
        trace!("Processing region {}(tid):{}-{}", tid, start, stop);
        // Create a reader
        let mut reader =
            bam::IndexedReader::from_path(&self.reads).expect("Indexed Reader for region");

        // If passed add ref_fasta
        if let Some(ref_fasta) = &self.ref_fasta {
            reader.set_reference(ref_fasta).expect("Set ref");
        }

        let header = reader.header().to_owned();
        // fetch the region of interest
        reader.fetch((tid, start, stop)).expect("Fetched a region");
        // Walk over pileups
        let mut pileup = reader.pileup();
        pileup.set_max_depth(std::cmp::min(i32::MAX.try_into().unwrap(), self.max_depth));
        let result: Vec<PileupPosition> = pileup
            .flat_map(|p| {
                let pileup = p.expect("Extracted a pileup");
                // Verify that we are within the bounds of the chunk we are iterating on
                let pileup_depth = pileup.depth();
                if pileup.pos() >= start && pileup.pos() < stop {
                    let mut pos = if self.mate_fix {
                        PileupPosition::from_pileup_mate_aware(
                            pileup,
                            &header,
                            &self.read_filter,
                            self.min_base_quality_score,
                            self.mate_fix_strategy,
                        )
                    } else {
                        PileupPosition::from_pileup(
                            pileup,
                            &header,
                            &self.read_filter,
                            self.min_base_quality_score,
                        )
                    };
                    // Add the ref base if reference is available
                    if let Some(buffer) = &self.ref_buffer {
                        let seq = buffer
                            .seq(&pos.ref_seq)
                            .expect("Fetched reference sequence");
                        pos.ref_base = Some(char::from(
                            *seq.get(pos.pos as usize)
                                .expect("Input SAM does not match reference"),
                        ));
                    }
                    if pileup_depth > (self.max_depth_warnings_cutoff) {
                        pos.near_max_depth = true;
                    }
                    pos.pos += self.coord_base;
                    Some(pos)
                } else {
                    None
                }
            })
            .collect();

        if self.keep_zeros {
            let mut new_result = vec![];
            let name = PileupPosition::compact_refseq(&header, tid);
            let mut pos = start;

            for position in result.into_iter() {
                while pos < (position.pos - self.coord_base) {
                    new_result.push(PileupPosition::new(name.clone(), pos + self.coord_base));
                    pos += 1;
                }
                new_result.push(position);
                pos += 1;
            }
            while pos < stop {
                new_result.push(PileupPosition::new(name.clone(), pos + self.coord_base));
                pos += 1;
            }
            new_result
        } else {
            result
        }
    }
}

#[cfg(test)]
#[allow(unused)]
mod tests {
    use super::*;
    use perbase_lib::position::{Position, pileup_position::PileupPosition};
    use proptest::bits::u8;
    use rstest::*;
    use rust_htslib::{bam, bam::record::Record};
    use smartstring::alias::*;
    use std::{collections::HashMap, path::PathBuf};
    use tempfile::{TempDir, tempdir};

    #[fixture]
    fn read_filter() -> DefaultReadFilter {
        DefaultReadFilter::new(0, 512, 0)
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
            bam::Writer::from_path(&path, &header, bam::Format::Bam).expect("Created writer");
        for record in records.iter() {
            writer.write(record).expect("Wrote record");
        }
        drop(writer); // force it to flush so indexing can happen
        // build the index
        bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();
        (path, tempdir)
    }

    #[fixture]
    fn non_mate_aware_positions(
        bamfile: (PathBuf, TempDir),
        read_filter: DefaultReadFilter,
    ) -> HashMap<String, Vec<PileupPosition>> {
        let cpus = utils::determine_allowed_cpus(8).unwrap();

        // Use the number of cpus available as a proxy for how may ref seqs to hold in memory at one time.
        let base_processor = BaseProcessor::new(
            bamfile.0.clone(),
            None,
            false,
            MateResolutionStrategy::Original,
            false,
            1,
            read_filter,
            500_000,
            cpus,
            None,
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
            base_processor,
        );
        let mut positions = HashMap::new();
        par_granges_runner
            .process()
            .unwrap()
            .into_iter()
            .for_each(|p| {
                let pos = positions.entry(p.ref_seq.clone()).or_insert_with(Vec::new);
                pos.push(p)
            });
        positions
    }

    #[fixture]
    fn non_mate_aware_keep_zeros_positions(
        bamfile: (PathBuf, TempDir),
        read_filter: DefaultReadFilter,
    ) -> HashMap<String, Vec<PileupPosition>> {
        let cpus = utils::determine_allowed_cpus(8).unwrap();

        // Use the number of cpus available as a proxy for how may ref seqs to hold in memory at one time.
        let base_processor = BaseProcessor::new(
            bamfile.0.clone(),
            None,
            false,
            MateResolutionStrategy::Original,
            true, // keep-zeros
            1,
            read_filter,
            500_000,
            cpus,
            None,
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
            base_processor,
        );
        let mut positions = HashMap::new();
        par_granges_runner
            .process()
            .unwrap()
            .into_iter()
            .for_each(|p| {
                let pos = positions.entry(p.ref_seq.clone()).or_insert_with(Vec::new);
                pos.push(p)
            });
        positions
    }

    #[fixture]
    fn mate_aware_positions(
        bamfile: (PathBuf, TempDir),
        read_filter: DefaultReadFilter,
    ) -> HashMap<String, Vec<PileupPosition>> {
        let cpus = utils::determine_allowed_cpus(8).unwrap();

        // Use the number of cpus available as a proxy for how may ref seqs to hold in memory at one time.
        let base_processor = BaseProcessor::new(
            bamfile.0.clone(),
            None,
            true, // mate aware
            MateResolutionStrategy::Original,
            false,
            1,
            read_filter,
            500_000,
            cpus,
            None,
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
            base_processor,
        );
        let mut positions = HashMap::new();
        par_granges_runner
            .process()
            .unwrap()
            .into_iter()
            .for_each(|p| {
                let pos = positions.entry(p.ref_seq.clone()).or_insert_with(Vec::new);
                pos.push(p)
            });
        positions
    }

    #[fixture]
    fn mate_aware_keep_zeros_positions(
        bamfile: (PathBuf, TempDir),
        read_filter: DefaultReadFilter,
    ) -> HashMap<String, Vec<PileupPosition>> {
        let cpus = utils::determine_allowed_cpus(8).unwrap();

        // Use the number of cpus available as a proxy for how may ref seqs to hold in memory at one time.
        let base_processor = BaseProcessor::new(
            bamfile.0.clone(),
            None,
            false, // mate aware
            MateResolutionStrategy::Original,
            true, // keep-zeros
            1,
            read_filter,
            500_000,
            cpus,
            None,
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
            base_processor,
        );
        let mut positions = HashMap::new();
        par_granges_runner
            .process()
            .unwrap()
            .into_iter()
            .for_each(|p| {
                let pos = positions.entry(p.ref_seq.clone()).or_insert_with(Vec::new);
                pos.push(p)
            });
        positions
    }

    fn non_mate_aware_positions_base_qual(
        bamfile: (PathBuf, TempDir),
        read_filter: DefaultReadFilter,
        base_quality: u8,
    ) -> HashMap<String, Vec<PileupPosition>> {
        let cpus = utils::determine_allowed_cpus(8).unwrap();

        // Use the number of cpus available as a proxy for how may ref seqs to hold in memory at one time.
        let base_processor = BaseProcessor::new(
            bamfile.0.clone(),
            None,
            false,
            MateResolutionStrategy::Original,
            false,
            1,
            read_filter,
            500_000,
            cpus,
            Some(base_quality),
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
            base_processor,
        );
        let mut positions = HashMap::new();
        par_granges_runner
            .process()
            .unwrap()
            .into_iter()
            .for_each(|p| {
                let pos = positions.entry(p.ref_seq.clone()).or_insert_with(Vec::new);
                pos.push(p)
            });
        positions
    }

    fn mate_aware_positions_base_qual(
        bamfile: (PathBuf, TempDir),
        read_filter: DefaultReadFilter,
        base_quality: u8,
    ) -> HashMap<String, Vec<PileupPosition>> {
        let cpus = utils::determine_allowed_cpus(8).unwrap();

        // Use the number of cpus available as a proxy for how may ref seqs to hold in memory at one time.
        let base_processor = BaseProcessor::new(
            bamfile.0.clone(),
            None,
            true, // mate aware
            MateResolutionStrategy::Original,
            false,
            1,
            read_filter,
            500_000,
            cpus,
            Some(base_quality),
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
            base_processor,
        );
        let mut positions = HashMap::new();
        par_granges_runner
            .process()
            .unwrap()
            .into_iter()
            .for_each(|p| {
                let pos = positions.entry(p.ref_seq.clone()).or_insert_with(Vec::new);
                pos.push(p)
            });
        positions
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_unaware_bq(non_mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 0),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 0),
        case::mate_unaware_bq(non_mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 0),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 0)
    )]
    fn check_insertions(positions: HashMap<String, Vec<PileupPosition>>, awareness_modifier: u32) {
        assert_eq!(positions.get("chr2").unwrap()[0].ins, 0);
        assert_eq!(positions.get("chr2").unwrap()[1].ins, 1);
        assert_eq!(positions.get("chr2").unwrap()[2].ins, 0);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_unaware_bq(non_mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 0),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 0),
        case::mate_unaware_bq(non_mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 0),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 0)
    )]
    fn check_deletions(positions: HashMap<String, Vec<PileupPosition>>, awareness_modifier: u32) {
        assert_eq!(positions.get("chr2").unwrap()[5].del, 0);
        assert_eq!(positions.get("chr2").unwrap()[6].del, 1);
        assert_eq!(positions.get("chr2").unwrap()[7].del, 1);
        assert_eq!(positions.get("chr2").unwrap()[8].del, 1);
        assert_eq!(positions.get("chr2").unwrap()[9].del, 1);
        assert_eq!(positions.get("chr2").unwrap()[10].del, 1);
        assert_eq!(positions.get("chr2").unwrap()[11].del, 0);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_unaware_bq(non_mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 0),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 0),
        case::mate_unaware_bq(non_mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 0),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 0)
    )]
    fn check_refskips(positions: HashMap<String, Vec<PileupPosition>>, awareness_modifier: u32) {
        assert_eq!(positions.get("chr2").unwrap()[11].ref_skip, 0);
        assert_eq!(positions.get("chr2").unwrap()[12].ref_skip, 1);
        assert_eq!(positions.get("chr2").unwrap()[13].ref_skip, 1);
        assert_eq!(positions.get("chr2").unwrap()[14].ref_skip, 1);
        assert_eq!(positions.get("chr2").unwrap()[15].ref_skip, 1);
        assert_eq!(positions.get("chr2").unwrap()[16].ref_skip, 1);
        assert_eq!(positions.get("chr2").unwrap()[17].ref_skip, 0);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_unaware_bq(non_mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 0),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 0),
        case::mate_unaware_bq(non_mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 0),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 0)
    )]
    fn check_start(positions: HashMap<String, Vec<PileupPosition>>, awareness_modifier: u32) {
        assert_eq!(positions.get("chr1").unwrap()[0].pos, 1);
        assert_eq!(positions.get("chr1").unwrap()[4].pos, 5);
        assert_eq!(positions.get("chr1").unwrap()[9].pos, 10);
        assert_eq!(positions.get("chr1").unwrap()[14].pos, 15);
        assert_eq!(positions.get("chr1").unwrap()[19].pos, 20);
        assert_eq!(positions.get("chr1").unwrap()[25].pos, 26);
        assert_eq!(positions.get("chr1").unwrap()[29].pos, 30);
        assert_eq!(positions.get("chr1").unwrap()[34].pos, 35);
        assert_eq!(positions.get("chr1").unwrap()[39].pos, 40);
        // NB: -6 bc there are 6 positions with no coverage from 44-50
        assert_eq!(positions.get("chr1").unwrap()[72].pos, 78);
    }

    #[rstest(
        positions,
        keep_zeros,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), false),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), false),
        case::mate_unaware_bq(
            non_mate_aware_positions_base_qual(bamfile(), read_filter(), 1),
            false
        ),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 1), false),
        case::mate_unaware_bq(
            non_mate_aware_positions_base_qual(bamfile(), read_filter(), 3),
            false
        ),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 3), false),
        case::mate_unaware_keep_zeros(
            non_mate_aware_keep_zeros_positions(bamfile(), read_filter()),
            true
        ),
        case::mate_aware_keep_zeros(
            mate_aware_keep_zeros_positions(bamfile(), read_filter()),
            true
        )
    )]
    fn check_depths(positions: HashMap<String, Vec<PileupPosition>>, keep_zeros: bool) {
        assert_eq!(positions.get("chr1").unwrap()[0].depth, 1);
        assert_eq!(positions.get("chr1").unwrap()[1].depth, 1);
        assert_eq!(positions.get("chr1").unwrap()[4].depth, 2);
        assert_eq!(positions.get("chr1").unwrap()[9].depth, 3);
        assert_eq!(positions.get("chr1").unwrap()[14].depth, 4);
        assert_eq!(positions.get("chr1").unwrap()[19].depth, 5);
        assert_eq!(positions.get("chr1").unwrap()[25].depth, 4);
        assert_eq!(positions.get("chr1").unwrap()[29].depth, 3);
        assert_eq!(positions.get("chr1").unwrap()[34].depth, 2);
        assert_eq!(positions.get("chr1").unwrap()[39].depth, 1);

        if !keep_zeros {
            // NB: -6 bc there are 6 positions with no coverage from 44-50
            assert_eq!(positions.get("chr1").unwrap()[78 - 6].depth, 4);
        } else {
            // NB: with --keep-zeros then there should be no skipped loci
            assert_eq!(positions.get("chr1").unwrap().len(), 100);

            // NB: Since positions[0].pos == 1 for this test data: positions[i].pos == i+1
            for (i, p) in positions.get("chr1").unwrap().iter().enumerate() {
                assert_eq!(p.pos, (i + 1).try_into().unwrap());
            }

            assert_eq!(positions.get("chr1").unwrap()[43].depth, 1);
            assert_eq!(positions.get("chr1").unwrap()[44].depth, 0);
            assert_eq!(positions.get("chr1").unwrap()[45].depth, 0);
            assert_eq!(positions.get("chr1").unwrap()[46].depth, 0);
            assert_eq!(positions.get("chr1").unwrap()[47].depth, 0);
            assert_eq!(positions.get("chr1").unwrap()[48].depth, 0);
            assert_eq!(positions.get("chr1").unwrap()[49].depth, 1);
            assert_eq!(positions.get("chr1").unwrap()[53].depth, 1);
            assert_eq!(positions.get("chr1").unwrap()[54].depth, 2);
            assert_eq!(positions.get("chr1").unwrap()[55].depth, 2);
            assert_eq!(positions.get("chr1").unwrap()[56].depth, 2);
            assert_eq!(positions.get("chr1").unwrap()[57].depth, 2);
            assert_eq!(positions.get("chr1").unwrap()[58].depth, 2);
            assert_eq!(positions.get("chr1").unwrap()[59].depth, 3);
            assert_eq!(positions.get("chr1").unwrap()[60].depth, 3);
            assert_eq!(positions.get("chr1").unwrap()[61].depth, 3);
            assert_eq!(positions.get("chr1").unwrap()[62].depth, 3);
            assert_eq!(positions.get("chr1").unwrap()[63].depth, 3);
            assert_eq!(positions.get("chr1").unwrap()[64].depth, 4);
            assert_eq!(positions.get("chr1").unwrap()[68].depth, 4);
            assert_eq!(positions.get("chr1").unwrap()[93].depth, 1);
            assert_eq!(positions.get("chr1").unwrap()[94].depth, 0);
            assert_eq!(positions.get("chr1").unwrap()[95].depth, 0);
            assert_eq!(positions.get("chr1").unwrap()[96].depth, 0);
            assert_eq!(positions.get("chr1").unwrap()[99].depth, 0);
        }
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_unaware_bq(non_mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 0),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 0),
        case::mate_unaware_bq(non_mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 0),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 0)
    )]
    fn check_filters(positions: HashMap<String, Vec<PileupPosition>>, awareness_modifier: u32) {
        // Verify that a read that has flags saying it failed QC got filtered out
        assert_eq!(positions.get("chr2").unwrap()[81].depth, 1);
        assert_eq!(positions.get("chr2").unwrap()[84].depth, 0);
        assert_eq!(positions.get("chr2").unwrap()[81].fail, 1);
        assert_eq!(positions.get("chr2").unwrap()[84].fail, 1);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_unaware_bq(non_mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 0),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 0),
        case::mate_unaware_bq(non_mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 0),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 0)
    )]
    fn check_depths_insertions(
        positions: HashMap<String, Vec<PileupPosition>>,
        awareness_modifier: u32,
    ) {
        assert_eq!(positions.get("chr2").unwrap()[0].depth, 1);
        assert_eq!(positions.get("chr2").unwrap()[1].depth, 1); // Insertion is here
        assert_eq!(positions.get("chr2").unwrap()[2].depth, 1);
    }

    #[rstest(
        positions,
        awareness_modifier,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_unaware_bq(non_mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 0),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 0),
        case::mate_unaware_bq(non_mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 0),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 0)
    )]
    fn check_depths_deletions(
        positions: HashMap<String, Vec<PileupPosition>>,
        awareness_modifier: u32,
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
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 0),
        case::mate_unaware_bq(non_mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 0),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 0),
        case::mate_unaware_bq(non_mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 0),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 0)
    )]
    fn check_depths_refskips(
        positions: HashMap<String, Vec<PileupPosition>>,
        awareness_modifier: u32,
    ) {
        assert_eq!(positions.get("chr2").unwrap()[11].depth, 3);
        assert_eq!(positions.get("chr2").unwrap()[12].depth, 2); // Skip
        assert_eq!(positions.get("chr2").unwrap()[13].depth, 2); // Skip
        assert_eq!(positions.get("chr2").unwrap()[14].depth, 3); // Skip
        assert_eq!(positions.get("chr2").unwrap()[15].depth, 3); // Skip
        assert_eq!(positions.get("chr2").unwrap()[16].depth, 3); // Skip
        assert_eq!(positions.get("chr2").unwrap()[17].depth, 4);
    }

    #[rustfmt::skip]
    #[rstest(
        positions,
        awareness_modifier,
        base_qual_filtered_all,
        case::mate_unaware(non_mate_aware_positions(bamfile(), read_filter()), 0, false),
        case::mate_aware(mate_aware_positions(bamfile(), read_filter()), 1, false),
        case::mate_unaware_bq(
            non_mate_aware_positions_base_qual(bamfile(), read_filter(), 1),
            0,
            false
        ),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 1), 1, false),
        case::mate_unaware_bq(
            non_mate_aware_positions_base_qual(bamfile(), read_filter(), 3),
            0,
            true
        ),
        case::mate_aware_bq(mate_aware_positions_base_qual(bamfile(), read_filter(), 3), 1, true)
    )]
    fn check_mate_detection(
        positions: HashMap<String, Vec<PileupPosition>>,
        awareness_modifier: u32,
        base_qual_filtered_all: bool,
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

        assert_eq!(positions.get("chr2").unwrap()[33].a, if base_qual_filtered_all { 0 } else { 4 } );
        assert_eq!(positions.get("chr2").unwrap()[34].a, if base_qual_filtered_all { 0 } else { 4 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[35].a, if base_qual_filtered_all { 0 } else { 4 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[36].a, if base_qual_filtered_all { 0 } else { 4 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[37].a, if base_qual_filtered_all { 0 } else { 4 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[38].a, if base_qual_filtered_all { 0 } else { 4 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[39].a, if base_qual_filtered_all { 0 } else { 2 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[40].a, if base_qual_filtered_all { 0 } else { 2 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[41].a, if base_qual_filtered_all { 0 } else { 2 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[42].a, if base_qual_filtered_all { 0 } else { 2 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[43].a, if base_qual_filtered_all { 0 } else { 2 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[44].a, if base_qual_filtered_all { 0 } else { 1 });

        assert_eq!(positions.get("chr2").unwrap()[33].n, if !base_qual_filtered_all { 0 } else { 4 } );
        assert_eq!(positions.get("chr2").unwrap()[34].n, if !base_qual_filtered_all { 0 } else { 4 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[35].n, if !base_qual_filtered_all { 0 } else { 4 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[36].n, if !base_qual_filtered_all { 0 } else { 4 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[37].n, if !base_qual_filtered_all { 0 } else { 4 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[38].n, if !base_qual_filtered_all { 0 } else { 4 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[39].n, if !base_qual_filtered_all { 0 } else { 2 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[40].n, if !base_qual_filtered_all { 0 } else { 2 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[41].n, if !base_qual_filtered_all { 0 } else { 2 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[42].n, if !base_qual_filtered_all { 0 } else { 2 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[43].n, if !base_qual_filtered_all { 0 } else { 2 - awareness_modifier }); // mate overlap
        assert_eq!(positions.get("chr2").unwrap()[44].n, if !base_qual_filtered_all { 0 } else { 1 });
    }
}
