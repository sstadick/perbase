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
    /// any base not covered above
    other_base: usize,
    /// number of insertions at this position
    ins: usize,
    /// number of deletions at this position
    del: usize,
    /// number of skips at this position
    skip: usize,
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
            other_base: 0,
            ins: 0,
            del: 0,
            skip: 0,
            depth: 0,
        }
    }
}

/// Calculate the depth at each base, per-base. Use samtools to prefilter
/// and pipe into this command.
#[derive(FromArgs)]
#[argh(subcommand, name = "simple-depth")]
pub struct SimpleDepth {
    /// input BAM/CRAM to analyze
    #[argh(positional)]
    reads: PathBuf,

    /// indexed reference simplea, set if using CRAM
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
        let mut reader = if self.reads.to_str().unwrap() == "-" {
            bam::Reader::from_stdin()?
        } else {
            bam::Reader::from_path(&self.reads)?
        };
        if reader_num >= 1 {
            reader.set_threads(reader_num)?;
        }
        // If passed add ref_simplea
        if let Some(ref_fasta) = self.ref_fasta {
            reader.set_reference(ref_fasta)?;
        }
        // Get a copy of the header
        let header = bam::HeaderView::from_header(&bam::Header::from_template(reader.header()));

        // Walk over pileups
        for p in reader.pileup() {
            let pileup = p?;
            let name = std::str::from_utf8(header.tid2name(pileup.tid())).unwrap();
            // make output 1-based
            let mut pos = Position::new(String::from(name), (pileup.pos() + 1) as usize);
            pos.depth = pileup.depth() as usize;

            for alignment in pileup.alignments() {
                if alignment.is_del() {
                    pos.del += 1;
                } else if alignment.is_refskip() {
                    // TODO: Should skips be counted as del's?
                    // pos.skip += 1;
                } else {
                    match (alignment.record().seq()[alignment.qpos().unwrap()] as char)
                        .to_ascii_uppercase()
                    {
                        'A' => pos.a += 1,
                        'C' => pos.c += 1,
                        'T' => pos.t += 1,
                        'G' => pos.g += 1,
                        'N' => pos.n += 1,
                        _ => pos.other_base += 1,
                    }

                    match alignment.indel() {
                        bam::pileup::Indel::Ins(_len) => pos.ins += 1,
                        bam::pileup::Indel::Del(_len) => pos.del += 1,
                        bam::pileup::Indel::None => (),
                    }
                }
            }
            writer.serialize(pos)?;
        }
        writer.flush()?;
        Ok(())
    }
}
