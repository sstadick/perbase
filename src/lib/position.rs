//! # Position
//!
//! A module for holding all the information about a given genomic position.
use itertools::Itertools;
use rust_htslib::bam::{
    self,
    pileup::{Alignment, Pileup},
    record::Record,
};
use serde::Serialize;
use smartstring::alias::String;
use std::{cmp::Ordering, default};

// TODOs:
// Write my own pilelup engine
// Reference base
// Average base quality
// Average map quality
// Average dist to ends
// A calculated error rate based on mismatches?
// Optionally display read names at the position?

/// Anything that implements ReadFilter can apply a filter set to read
pub trait ReadFilter {
    /// filters a read, true is pass, false if fail
    fn filter_read(&self, read: &Record) -> bool;
}
/// Hold all information about a position.
#[derive(Debug, Serialize, Default)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub struct Position {
    /// Reference sequence name.
    #[serde(rename = "REF")]
    pub ref_seq: String,
    /// 1-based position in the sequence.
    pub pos: usize,
    /// Total depth at this position.
    pub depth: usize,
    /// Number of A bases at this position.
    pub a: usize,
    /// Number of C bases at this position.
    pub c: usize,
    /// Number of G bases at this position.
    pub g: usize,
    /// Number of T bases at this position.
    pub t: usize,
    /// Number of N bases at this position. Any unrecognized base will be counted as an N.
    pub n: usize,
    /// Number of insertions that start to the right of this position.
    /// Does not count toward depth.
    pub ins: usize,
    /// Number of deletions at this position.
    pub del: usize,
    /// Number of refskips at this position. Does not count toward depth.
    pub ref_skip: usize,
    /// Number of reads failing filters at this position.
    pub fail: usize,
}

impl Position {
    /// Create a new position for the given ref_seq name.
    pub fn new(ref_seq: String, pos: usize) -> Self {
        Position {
            ref_seq,
            pos,
            ..default::Default::default()
        }
    }

    /// Given a record, update the counts at this position
    fn update<F: ReadFilter>(&mut self, alignment: &Alignment, record: Record, read_filter: &F) {
        if !read_filter.filter_read(&record) {
            self.depth -= 1;
            self.fail += 1;
            return;
        }
        // NB: Order matters here, a refskip is true for both is_del and is_refskip
        // while a true del is only true for is_del
        if alignment.is_refskip() {
            self.ref_skip += 1;
            self.depth -= 1;
        } else if alignment.is_del() {
            self.del += 1;
        } else {
            // We have an actual base!
            match (record.seq()[alignment.qpos().unwrap()] as char).to_ascii_uppercase() {
                'A' => self.a += 1,
                'C' => self.c += 1,
                'T' => self.t += 1,
                'G' => self.g += 1,
                _ => self.n += 1,
            }
            // Check for insertions
            match alignment.indel() {
                bam::pileup::Indel::Ins(_len) => {
                    self.ins += 1;
                }
                _ => (),
            }
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
    /// * `pileup` - a [bam::pileup::Pileup] at a genomic position
    /// * `header` - a [bam::HeaderView] for the bam file being read, to get the sequence name
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
