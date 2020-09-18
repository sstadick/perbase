use anyhow::Result;
use log::*;
use rust_htslib::{bam, bam::Read};
use std::path::PathBuf;
/// Simple single pass over bam file to calculate depth at each position
/// as well as depth per nucleotide. Additionally counts the number of
/// insertions / deletions at each position.

pub fn parse(reads: PathBuf, ref_fasta: PathBuf) -> Result<()> {
    info!("Parsing reads: {:?}", reads);

    let mut reads = bam::Reader::from_path(reads)?;
    for p in reads.pileup() {
        let pileup = p?;
        println!("{}:{} depth {}", pileup.tid(), pileup.pos(), pileup.depth());
        for alignment in pileup.alignments() {
            if !alignment.is_del() && !alignment.is_refskip() {
                println!(
                    "Base {}",
                    alignment.record().seq()[alignment.qpos().unwrap()]
                );
            }
            // mark indel start
            match alignment.indel() {
                bam::pileup::Indel::Ins(len) => println!(
                    "Insertion of length {} between this and next position.",
                    len
                ),
                bam::pileup::Indel::Del(len) => {
                    println!("Deletion of length {} between this and next position.", len)
                }
                bam::pileup::Indel::None => (),
            }
        }
    }

    Ok(())
}
