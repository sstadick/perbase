//! # Simple Depth
//!
//! Simple single pass over bam file to calculate depth at each position
//! as well as depth per nucleotide. Additionally counts the number of
//! insertions / deletions at each position.
use anyhow::Result;
use argh::FromArgs;
use grep_cli::stdout;
use log::*;
use perbase_lib::utils;
use rust_htslib::{bam, bam::Read};
use serde::Serialize;
use smartstring::alias::String;
use std::{
    default,
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};
use termcolor::ColorChoice;

/// Hold all information about a position.
#[derive(Debug, Serialize, Default)]
struct Position {
    /// Reference sequence name.
    #[serde(alias = "REF")]
    ref_seq: String,
    /// 1-based position in the sequence.
    #[serde(alias = "POS")]
    pos: usize,
    /// Total depth at this position.
    #[serde(alias = "DEPTH")]
    depth: usize,
    /// Number of A bases at this position.
    #[serde(alias = "A")]
    a: usize,
    /// Number of C bases at this position.
    #[serde(alias = "C")]
    c: usize,
    /// Number of G bases at this position.
    #[serde(alias = "G")]
    g: usize,
    /// Number of T bases at this position.
    #[serde(alias = "T")]
    t: usize,
    /// Number of N bases at this position. Any unrecognized base will be counted as an N.
    #[serde(alias = "N")]
    n: usize,
    /// Number of insertions that start to the right of this position.
    /// Does not count toward depth.
    #[serde(alias = "INS")]
    ins: usize,
    /// Number of deletions at this position.
    #[serde(alias = "DEL")]
    del: usize,
    /// Number of refskips at this position. Does count toward depth.
    #[serde(alias = "REF_SKIP")]
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
    /// Reference base
    /// Average base quality
    /// Average map quality
    /// Average dist to ends
    /// A calculated error rate based on mismatches?
    /// Optionally display read names at the position?
    fn from_pileup(pileup: bam::pileup::Pileup, header: &bam::HeaderView) -> Self {
        let name = std::str::from_utf8(header.tid2name(pileup.tid())).unwrap();
        // make output 1-based
        let mut pos = Self::new(String::from(name), (pileup.pos() + 1) as usize);
        pos.depth = pileup.depth() as usize;

        for alignment in pileup.alignments() {
            // NB: Order matters here, a refskip is true for both is_del and is_refskip
            // while a true del is only true for is_del
            if alignment.is_refskip() {
                pos.ref_skip += 1
            } else if alignment.is_del() {
                pos.del += 1
            } else {
                // We have an actual base!
                match (alignment.record().seq()[alignment.qpos().unwrap()] as char)
                    .to_ascii_uppercase()
                {
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

/// Calculate the depth at each base, per-base. Use samtools to prefilter
/// and pipe into this command.
#[derive(FromArgs)]
#[argh(subcommand, name = "simple-depth")]
pub struct SimpleDepth {
    /// input BAM/CRAM to analyze
    #[argh(positional)]
    reads: Option<PathBuf>,

    /// indexed reference fasta, set if using CRAM
    #[argh(option, short = 'r')]
    ref_fasta: Option<PathBuf>,

    /// output path, defaults to stdout
    #[argh(option, short = 'o')]
    output: Option<PathBuf>,

    /// the number of threads to use
    #[argh(option, short = 't', default = "num_cpus::get()")]
    threads: usize,
}

impl SimpleDepth {
    // TODO: Add special handling when the bam is indexed and read it multi-threaded mode, fetching on chr at a time
    // TODO: Allow specifying a region, multithreaded over regions
    // TODO: Add sam flag support as found here: https://github.com/rust-bio/rust-bio-tools/blob/master/src/bam/depth.rs
    // TODO: Add mate detection like sambamba
    // TODO: Update README and add a table explaining the output and options
    // TODO: Add a cached MD tag parser to rust_htslib
    pub fn run(self) -> Result<()> {
        info!("Running simple-depth on: {:?}", self.reads);

        let cpus = utils::determine_allowed_cpus(self.threads)?;
        let reader_num = if cpus > 1 {
            std::cmp::min(cpus - 1, 4)
        } else {
            0
        };
        info!("Using threads: {}", reader_num + 1);

        // Set up output writer
        let raw_writer: Box<dyn Write> = match self.output {
            Some(path) if path.to_str().unwrap() != "-" => {
                Box::new(BufWriter::new(File::open(path)?))
            }
            _ => Box::new(stdout(ColorChoice::Never)),
        };
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(raw_writer);

        // Set up input reader
        let mut reader = if let Some(path) = self.reads {
            info!("Reading from {:?}", path);
            bam::Reader::from_path(&path)?
        } else {
            info!("Reading from STDIN");
            bam::Reader::from_stdin()?
        };
        if reader_num >= 1 {
            reader.set_threads(reader_num)?;
        }
        // If passed add ref_fasta
        if let Some(ref_fasta) = self.ref_fasta {
            reader.set_reference(ref_fasta)?;
        }
        // Get a copy of the header
        let header = reader.header().to_owned();

        // Walk over pileups
        for p in reader.pileup() {
            let pileup = p?;
            let pos = Position::from_pileup(pileup, &header);
            writer.serialize(pos)?;
        }
        writer.flush()?;
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
            Record::from_sam(&view, b"FOUR_2\t147\tchr2\t65\t40\t25M\tchr2\t15\t75\tAAAAAAAAAAAAAAAAAAAAAAAAA\t#########################").unwrap(),
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
            let pos = Position::from_pileup(pileup, &header);
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
