//! An implementation of `Position` for dealing with pileups.
use crate::position::Position;
use crate::read_filter::ReadFilter;
use rust_htslib::bam::{
    self,
    pileup::{Alignment, Pileup},
    record::Record,
    HeaderView,
};
use serde::Serialize;
use smartstring::{alias::String, LazyCompact, SmartString};
use std::collections::HashMap;
use std::{cmp::Ordering, default};

/// Hold all information about a position.
// NB: The max depth that htslib will return is i32::MAX, and the type of pos for htlib is u32
// There is no reason to go bigger, for now at least
#[derive(Debug, Serialize, Default)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub struct PileupPosition {
    /// Reference sequence name.
    #[serde(rename = "REF")]
    pub ref_seq: String,
    /// 1-based position in the sequence.
    pub pos: u32,
    /// The reference base at this position.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ref_base: Option<char>,
    /// Total depth at this position.
    pub depth: u32,
    /// Number of A bases at this position.
    pub a: u32,
    /// Number of C bases at this position.
    pub c: u32,
    /// Number of G bases at this position.
    pub g: u32,
    /// Number of T bases at this position.
    pub t: u32,
    /// Number of N bases at this position. Any unrecognized base will be counted as an N.
    pub n: u32,
    /// Number of insertions that start to the right of this position.
    /// Does not count toward depth.
    pub ins: u32,
    /// Number of deletions at this position.
    pub del: u32,
    /// Number of refskips at this position. Does not count toward depth.
    pub ref_skip: u32,
    /// Number of reads failing filters at this position.
    pub fail: u32,
    /// Depth is within 1% of max_depth
    pub near_max_depth: bool,
}

impl Position for PileupPosition {
    /// Create a new position for the given ref_seq name.
    fn new(ref_seq: String, pos: u32) -> Self {
        PileupPosition {
            ref_seq,
            pos,
            ..default::Default::default()
        }
    }
}

impl PileupPosition {
    /// Given a record, update the counts at this position
    #[inline(always)]
    fn update<F: ReadFilter>(
        &mut self,
        alignment: &Alignment,
        record: Record,
        read_filter: &F,
        base_filter: Option<u8>,
    ) {
        if !read_filter.filter_read(&record, Some(alignment)) {
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

            // Check if we are checking the base quality score
            if let Some(base_qual_filter) = base_filter {
                // Check if the base quality score is greater or equal to than the cutoff
                // TODO: When `if let` + && / || stabilizes clean this up.
                if record.qual()[alignment.qpos().unwrap()] < base_qual_filter {
                    self.n += 1
                } else {
                    match (record.seq()[alignment.qpos().unwrap()] as char).to_ascii_uppercase() {
                        'A' => self.a += 1,
                        'C' => self.c += 1,
                        'T' => self.t += 1,
                        'G' => self.g += 1,
                        _ => self.n += 1,
                    }
                }
            } else {
                match (record.seq()[alignment.qpos().unwrap()] as char).to_ascii_uppercase() {
                    'A' => self.a += 1,
                    'C' => self.c += 1,
                    'T' => self.t += 1,
                    'G' => self.g += 1,
                    _ => self.n += 1,
                }
            }
            // Check for insertions
            if let bam::pileup::Indel::Ins(_len) = alignment.indel() {
                self.ins += 1;
            }
        }
    }

    /// Convert a pileup into a `Position`.
    ///
    /// This will walk over each of the alignments and count the number each nucleotide it finds.
    /// It will also count the number of Ins/Dels/Skips that are at each position.
    ///
    /// # Arguments
    ///
    /// * `pileup` - a pileup at a genomic position
    /// * `header` - a headerview for the bam file being read, to get the sequence name
    /// * `read_filter` - a function to filter out reads, returning false will cause a read to be filtered
    /// * `base_filter` - an optional base quality score. If Some(number) if the base quality is not >= that number the base is treated as an `N`
    #[inline]
    pub fn from_pileup<F: ReadFilter>(
        pileup: Pileup,
        header: &bam::HeaderView,
        read_filter: &F,
        base_filter: Option<u8>,
    ) -> Self {
        let name = Self::compact_refseq(header, pileup.tid());
        // make output 1-based
        let mut pos = Self::new(name, pileup.pos());
        pos.depth = pileup.depth();

        for alignment in pileup.alignments() {
            let record = alignment.record();
            Self::update(&mut pos, &alignment, record, read_filter, base_filter);
        }
        pos
    }

    /// Convert a pileup into a `Position`.
    ///
    /// This will walk over each of the alignments and count the number each nucleotide it finds.
    /// It will also count the number of Ins/Dels/Skips that are at each position.
    ///
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
    /// * `base_filter` - an optional base quality score. If Some(number) if the base quality is not >= that number the base is treated as an `N`
    #[inline]
    pub fn from_pileup_mate_aware<F: ReadFilter>(
        pileup: Pileup,
        header: &bam::HeaderView,
        read_filter: &F,
        base_filter: Option<u8>,
    ) -> Self {
        let name = Self::compact_refseq(header, pileup.tid());
        // make output 1-based
        let mut pos = Self::new(name, pileup.pos());
        pos.depth = pileup.depth();

        let mut by_name = HashMap::new();

        pileup.alignments().for_each(|aln| {
            let record = aln.record();
            if record.tid() != record.mtid() || record.cigar().end_pos() < record.mpos() {
                // we know this doesn't overlap and it's likely left-most
                Self::update(&mut pos, &aln, record, read_filter, base_filter);
            } else {
                let name = record.qname().to_owned();
                by_name
                    .entry(name)
                    .or_insert_with(Vec::new)
                    .push((aln, record));
            }
        });

        // Group records by qname
        for (_qname, alignments) in by_name.into_iter() {
            // Choose the best of the reads based on mapq, if tied, check which is first and passes filters
            let mut total_reads = 0; // count how many reads there were
            if alignments.len() == 1 {
                // unwrap is safe here as we know we have one element
                let (aln, record) = alignments.into_iter().next().unwrap();
                Self::update(&mut pos, &aln, record, read_filter, base_filter);
                continue;
            }
            let (alignment, record) = alignments
                .into_iter()
                .map(|(aln, record)| {
                    total_reads += 1;
                    (aln, record)
                })
                .max_by(|a, b| match a.1.mapq().cmp(&b.1.mapq()) {
                    Ordering::Greater => Ordering::Greater,
                    Ordering::Less => Ordering::Less,
                    Ordering::Equal => {
                        // Check if a is first in pair
                        if a.1.flags() & 64 == 0 && read_filter.filter_read(&a.1, Some(&a.0)) {
                            Ordering::Greater
                        } else if b.1.flags() & 64 == 0 && read_filter.filter_read(&b.1, Some(&b.0))
                        {
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
            Self::update(&mut pos, &alignment, record, read_filter, base_filter);
        }
        pos
    }

    /// Convert a tid to a [`SmartString<LazyCompact>`].
    #[inline]
    pub fn compact_refseq(header: &HeaderView, tid: u32) -> SmartString<LazyCompact> {
        let name = std::str::from_utf8(header.tid2name(tid)).unwrap();
        String::from(name)
    }
}
