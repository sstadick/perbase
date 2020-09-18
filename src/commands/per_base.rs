use anyhow::Result;
use log::*;
use rust_htslib::{bam, bam::Read};
use serde::Serialize;
use smartstring::alias::String;
use std::io;
use std::path::PathBuf;
/// Simple single pass over bam file to calculate depth at each position
/// as well as depth per nucleotide. Additionally counts the number of
/// insertions / deletions at each position.

#[derive(Debug, Serialize)]
struct Position {
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
    /// Any base not covered above
    other_base: usize,
    ins: usize,
    del: usize,
    skip: usize,
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

// TODO: filter out none-primary alignments
// TODO: Add rayon? is that even a bottleneck
// TODO: Add faster printer / output file option
// TODO: Add reference base / mismatch score??
pub fn parse(reads: PathBuf, ref_fasta: PathBuf) -> Result<()> {
    info!("Parsing reads: {:?}", reads);

    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(io::stdout());

    let header_reader = bam::Reader::from_path(&reads)?;
    let mut reads = bam::Reader::from_path(&reads)?;
    let header = header_reader.header();
    for p in reads.pileup() {
        let pileup = p?;
        let name = std::str::from_utf8(header.tid2name(pileup.tid())).unwrap();
        let mut pos = Position::new(String::from(name), pileup.pos() as usize);
        pos.depth = pileup.depth() as usize;

        for alignment in pileup.alignments() {
            if alignment.is_del() {
                pos.del += 1;
            } else if alignment.is_refskip() {
                pos.skip += 1;
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
