/// Simple single pass over bam file to calculate depth at each position
/// as well as depth per nucleotide. Additionally counts the number of
/// insertions / deletions at each position.
use anyhow::Result;
use argh::FromArgs;
use grep_cli::stdout;
use log::*;
use perbase_lib::utils;
use rust_htslib::{bam, bam::Read};
use serde::Serialize;
use smartstring::alias::String;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};
use termcolor::ColorChoice;

#[derive(Debug, Serialize)]
struct Position {
    /// contig name
    chr: String,
    /// 1-based position
    pos: usize,
    /// number of A bases at this position
    #[serde(alias = "A")]
    a: usize,
    /// number of C bases at this position
    #[serde(alias = "C")]
    c: usize,
    /// number of T bases at this position
    #[serde(alias = "T")]
    t: usize,
    /// number of G bases at this position
    #[serde(alias = "G")]
    g: usize,
    /// number of N bases at this position
    #[serde(alias = "N")]
    n: usize,
    /// number of insertions at this position
    ins: usize,
    /// number of deletions at this position
    del: usize,
    /// total depth at this position
    depth: usize,
}
impl Position {
    fn new(chr: String, pos: usize) -> Self {
        Position {
            chr,
            pos,
            a: 0,
            c: 0,
            t: 0,
            g: 0,
            n: 0,
            ins: 0,
            del: 0,
            depth: 0,
        }
    }

    /// Convert a pileup into a `Position`.
    ///
    /// This will walk over each of the alignments and count the number each nucleotide it finds.
    /// It will also count the number of Ins/Dels that are at each position. The output of this 1-based.
    fn from_pileup(pileup: bam::pileup::Pileup, header: &bam::HeaderView) -> Self {
        let name = std::str::from_utf8(header.tid2name(pileup.tid())).unwrap();
        // make output 1-based
        let mut pos = Self::new(String::from(name), (pileup.pos() + 1) as usize);
        pos.depth = pileup.depth() as usize;

        for alignment in pileup.alignments() {
            if !alignment.is_del() && !alignment.is_refskip() {
                match (alignment.record().seq()[alignment.qpos().unwrap()] as char)
                    .to_ascii_uppercase()
                {
                    'A' => pos.a += 1,
                    'C' => pos.c += 1,
                    'T' => pos.t += 1,
                    'G' => pos.g += 1,
                    _ => pos.n += 1,
                }

                // Only report indel starts in the same way a VCF does
                match alignment.indel() {
                    bam::pileup::Indel::Ins(_len) => {
                        pos.ins += 1;
                    }
                    bam::pileup::Indel::Del(_len) => {
                        pos.del += 1;
                    }
                    bam::pileup::Indel::None => (),
                }
            } else {
                // htslib is counting these toward depth. Don't
                pos.depth -= 1;
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
    use csv::ReaderBuilder;
    use rust_htslib::{bam, bam::record::Record};
    use serde::Deserialize;
    use std::path::PathBuf;
    use tempfile::tempdir;

    fn generate_data(path: &PathBuf) {
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

        let mut writer =
            bam::Writer::from_path(path, &header, bam::Format::BAM).expect("Created writer");
        // Keep the `test/test.bam` up to date
        let mut test_writer = bam::Writer::from_path("test/test.bam", &header, bam::Format::BAM)
            .expect("Created writer");
        for record in records.iter() {
            writer.write(record).expect("Wrote record");
            test_writer.write(record).expect("Wrote to test dir");
        }
    }

    #[derive(Debug, Deserialize)]
    struct Expected {
        chr: String,
        pos: usize,
        #[serde(alias = "A")]
        a: usize,
        #[serde(alias = "C")]
        c: usize,
        #[serde(alias = "T")]
        t: usize,
        #[serde(alias = "G")]
        g: usize,
        #[serde(alias = "N")]
        n: usize,
        total: usize,
    }

    /// Expected results for the test data based on bam-readcounts
    fn expected_results() -> Vec<Expected> {
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path("./test/expected.tsv")
            .unwrap();
        let mut results = vec![];
        for record in reader.deserialize() {
            let record: Expected = record.unwrap();
            results.push(record);
        }
        results
    }

    #[test]
    fn sanity() {
        let tmp = tempdir().expect("Made tempdir");
        let bam = tmp.path().join("test.bam");
        generate_data(&bam);
        let expected = expected_results();
        let mut reader = bam::Reader::from_path(bam).expect("Opened bam for reading");
        let header = reader.header().to_owned();
        let mut positions = vec![vec![], vec![]];
        for (i, p) in reader.pileup().enumerate() {
            let pileup = p.unwrap();
            let tid = pileup.tid();
            let pos = Position::from_pileup(pileup, &header);
            assert_eq!(&expected[i].a, &pos.a);
            assert_eq!(&expected[i].c, &pos.c);
            assert_eq!(&expected[i].g, &pos.g);
            assert_eq!(&expected[i].t, &pos.t);
            assert_eq!(&expected[i].n, &pos.n);
            assert_eq!(&expected[i].total, &pos.depth);
            positions[tid as usize].push(pos);
        }

        // Confirm select Ins / Dels
        // NB: Using String::new because String is aliased by Smartstring
        assert_eq!(positions[1][0].ins, 0);
        assert_eq!(positions[1][1].ins, 1);
        assert_eq!(positions[1][2].ins, 0);
        assert_eq!(positions[1][4].del, 0);
        assert_eq!(positions[1][5].del, 1);
        assert_eq!(positions[1][6].del, 0);
    }
}
