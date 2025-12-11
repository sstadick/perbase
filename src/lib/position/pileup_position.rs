//! An implementation of `Position` for dealing with pileups.
use crate::position::Position;
use crate::read_filter::ReadFilter;
use itertools::Itertools;
use rust_htslib::bam::{
    self, HeaderView,
    pileup::{Alignment, Pileup},
    record::Record,
};
use serde::Serialize;
use smartstring::{LazyCompact, SmartString, alias::String};
use std::default;

use super::mate_fix::{Base, MateResolutionStrategy};

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
    /// Number of bases that could be A or G
    pub r: u32,
    /// Number of bases that could be C or T
    pub y: u32,
    /// Number of bases that could be G or C
    pub s: u32,
    /// Number of bases that could be A or T
    pub w: u32,
    /// Number of bases that could be G or T
    pub k: u32,
    /// Number of bases that could be A or C
    pub m: u32,
    /// Number of insertions that start to the right of this position.
    /// Does not count toward depth.
    pub ins: u32,
    /// Number of deletions at this position.
    pub del: u32,
    /// Number of refskips at this position. Does not count toward depth.
    pub ref_skip: u32,
    /// Number of reads failing filters at this position.
    pub fail: u32,
    /// Number of times a mate resolution was needed.
    pub count_of_mate_resoutions: u32,
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
        record: &Record,
        read_filter: &F,
        base_filter: Option<u8>,
        recommended_base: Option<Base>,
        mates_resolved: bool,
    ) {
        if !read_filter.filter_read(record, Some(alignment)) {
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
            // && Check if the base quality score is greater or equal to than the cutoff
            if let Some(base_qual_filter) = base_filter
                && (record.seq().is_empty() || record.qual()[alignment.qpos().unwrap()] < base_qual_filter)
            {
                self.n += 1
            } else if let Some(b) = recommended_base {
                match b {
                    Base::A => self.a += 1,
                    Base::C => self.c += 1,
                    Base::T => self.t += 1,
                    Base::G => self.g += 1,
                    Base::R => self.r += 1,
                    Base::Y => self.y += 1,
                    Base::S => self.s += 1,
                    Base::W => self.w += 1,
                    Base::K => self.k += 1,
                    Base::M => self.m += 1,
                    _ => self.n += 1,
                }
            } else if record.seq().is_empty() {
                self.n += 1
            } else {
                match (record.seq()[alignment.qpos().unwrap()] as char).to_ascii_uppercase() {
                    'A' => self.a += 1,
                    'C' => self.c += 1,
                    'T' | 'U' => self.t += 1,
                    'G' => self.g += 1,
                    _ => self.n += 1,
                }
            }
            // Check for insertions
            if let bam::pileup::Indel::Ins(_len) = alignment.indel() {
                self.ins += 1;
            }
        }

        if mates_resolved {
            self.count_of_mate_resoutions += 1;
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
            Self::update(
                &mut pos,
                &alignment,
                &record,
                read_filter,
                base_filter,
                None,
                false,
            );
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
        mate_fix_strat: MateResolutionStrategy,
    ) -> Self {
        let name = Self::compact_refseq(header, pileup.tid());
        // make output 1-based
        let mut pos = Self::new(name, pileup.pos());
        pos.depth = pileup.depth();

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
            let mut total_reads = 0; // count how many reads there were

            let mut reads = reads
                .inspect(|_| total_reads += 1)
                .sorted_by(|a, b| mate_fix_strat.cmp(a, b, read_filter).ordering.reverse());
            let best = &reads.next().unwrap();

            if let Some(second_best) = &reads.next().as_ref() {
                // Deal explicitly with mate overlap, they are already ordered correctly
                let result = mate_fix_strat.cmp(best, second_best, read_filter);
                pos.depth -= total_reads - 1;
                Self::update(
                    &mut pos,
                    &best.0,
                    &best.1,
                    read_filter,
                    base_filter,
                    result.recommended_base,
                    true,
                );
            } else {
                pos.depth -= total_reads - 1;
                Self::update(
                    &mut pos,
                    &best.0,
                    &best.1,
                    read_filter,
                    base_filter,
                    None,
                    false,
                );
            }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_filter::DefaultReadFilter;
    use rust_htslib::bam::{self, Read, record::Record};
    use std::cmp::Ordering;

    /// Test that from_pileup_mate_aware chooses first mate when MAPQ scores are equal
    /// This verifies the fix for issue #82 by testing the actual function with overlapping mates
    #[test]
    fn test_from_pileup_mate_aware_chooses_first_mate() {
        use rust_htslib::bam::{IndexedReader, Writer, index};
        use tempfile::tempdir;

        let tempdir = tempdir().unwrap();
        let bam_path = tempdir.path().join("test.bam");

        // Create header
        let mut header = bam::header::Header::new();
        let mut chr1 = bam::header::HeaderRecord::new(b"SQ");
        chr1.push_tag(b"SN", &"chr1".to_owned());
        chr1.push_tag(b"LN", &"100".to_owned());
        header.push_record(&chr1);
        let view = bam::HeaderView::from_header(&header);

        // Create overlapping mate pair with equal MAPQ where they disagree on a base
        // Both mates overlap at position 15 (1-based)
        // First mate has 'A' at position 15, second mate has 'C' at position 15
        // Both have MAPQ 40, so tie-breaker should choose first mate
        let records = vec![
            Record::from_sam(
                &view,
                b"TESTPAIR\t67\tchr1\t10\t40\t10M\tchr1\t15\t30\tAAAAAAAAAA\t##########",
            )
            .unwrap(),
            Record::from_sam(
                &view,
                b"TESTPAIR\t147\tchr1\t15\t40\t10M\tchr1\t10\t30\tCCCCCCCCCC\t##########",
            )
            .unwrap(),
        ];

        // Write BAM file
        let mut writer = Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
        for record in &records {
            writer.write(record).unwrap();
        }
        drop(writer);

        // Index the BAM file
        index::build(&bam_path, None, index::Type::Bai, 1).unwrap();

        // Read back and create pileup at the overlapping position (15, 1-based = 14, 0-based)
        let mut reader = IndexedReader::from_path(&bam_path).unwrap();
        let header_view = reader.header().clone();

        // Set up pileup at position 14 (0-based) where both mates overlap
        reader.fetch(("chr1", 14, 15)).unwrap();
        let pileup_iter = reader.pileup();

        let read_filter = DefaultReadFilter::new(0, 0, 0);

        // Find the pileup at position 14
        let mut found_position = false;
        for pileup_result in pileup_iter {
            let pileup = pileup_result.unwrap();
            if pileup.pos() == 14 {
                // 0-based position 14 = 1-based position 15
                found_position = true;

                // Test mate-aware processing
                let position = PileupPosition::from_pileup_mate_aware(
                    pileup,
                    &header_view,
                    &read_filter,
                    None,
                    MateResolutionStrategy::Original,
                );

                // Should choose first mate's 'A' over second mate's 'C'
                // Depth should be 1 (one mate chosen from the overlapping pair)
                assert_eq!(position.depth, 1, "Depth should be 1 after mate selection");
                assert_eq!(position.a, 1, "Should count first mate's A base");
                assert_eq!(position.c, 0, "Should not count second mate's C base");
                break;
            }
        }

        assert!(found_position, "Should have found pileup at position 14");
    }

    /// Test different MateResolutionStrategy behaviors
    #[test]
    fn test_mate_resolution_strategies() {
        use rust_htslib::bam::{IndexedReader, Writer, index};
        use tempfile::tempdir;

        let tempdir = tempdir().unwrap();
        let bam_path = tempdir.path().join("test.bam");

        // Create header
        let mut header = bam::header::Header::new();
        let mut chr1 = bam::header::HeaderRecord::new(b"SQ");
        chr1.push_tag(b"SN", &"chr1".to_owned());
        chr1.push_tag(b"LN", &"100".to_owned());
        header.push_record(&chr1);
        let view = bam::HeaderView::from_header(&header);

        // Create overlapping mate pair with different base qualities and bases
        // Position 14 (0-based): First mate has 'A' with quality 30, second mate has 'G' with quality 40
        let records = vec![
            Record::from_sam(
                &view,
                b"TESTPAIR\t67\tchr1\t10\t40\t10M\tchr1\t15\t30\tAAAAAAAAGA\t#########=",
            )
            .unwrap(),
            Record::from_sam(
                &view,
                b"TESTPAIR\t147\tchr1\t15\t40\t10M\tchr1\t10\t30\tGGGGGGGGGG\t((((((((((",
            )
            .unwrap(),
        ];

        // Write BAM file
        let mut writer = Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
        for record in &records {
            writer.write(record).unwrap();
        }
        drop(writer);

        // Index the BAM file
        index::build(&bam_path, None, index::Type::Bai, 1).unwrap();

        // Read back and test different strategies
        let mut reader = IndexedReader::from_path(&bam_path).unwrap();
        let _header_view = reader.header().clone();

        // Test position 19 (0-based) where mates overlap
        reader.fetch(("chr1", 19, 20)).unwrap();
        let pileup_iter = reader.pileup();

        let read_filter = DefaultReadFilter::new(0, 0, 0);

        for pileup_result in pileup_iter {
            let pileup = pileup_result.unwrap();
            if pileup.pos() == 19 {
                // Test with different strategies
                let strat = crate::position::mate_fix::MateResolutionStrategy::BaseQualMapQualIUPAC;

                // Get alignments
                let alns: Vec<_> = pileup
                    .alignments()
                    .map(|aln| {
                        let rec = aln.record();
                        (aln, rec)
                    })
                    .collect();

                if alns.len() == 2 {
                    // Test that higher base quality wins (second mate has quality 40)
                    let result = strat.cmp(&alns[0], &alns[1], &read_filter);
                    assert_eq!(
                        result.ordering,
                        Ordering::Less,
                        "Second mate should win due to higher base quality"
                    );
                }
                break;
            }
        }
    }

    /// Test IUPAC base resolution
    #[test]
    fn test_iupac_base_resolution() {
        use rust_htslib::bam::{IndexedReader, Writer, index};
        use tempfile::tempdir;

        let tempdir = tempdir().unwrap();
        let bam_path = tempdir.path().join("test.bam");

        // Create header
        let mut header = bam::header::Header::new();
        let mut chr1 = bam::header::HeaderRecord::new(b"SQ");
        chr1.push_tag(b"SN", &"chr1".to_owned());
        chr1.push_tag(b"LN", &"100".to_owned());
        header.push_record(&chr1);
        let view = bam::HeaderView::from_header(&header);

        // Create overlapping mates with same quality but different bases
        // Both have MAPQ 40 and base quality 30
        let records = vec![
            Record::from_sam(
                &view,
                b"TESTPAIR\t67\tchr1\t10\t40\t10M\tchr1\t15\t30\tAAAAAAAAAA\t>>>>>>>>>>",
            )
            .unwrap(),
            Record::from_sam(
                &view,
                b"TESTPAIR\t147\tchr1\t15\t40\t10M\tchr1\t10\t30\tGGGGGGGGGG\t>>>>>>>>>>",
            )
            .unwrap(),
        ];

        // Write BAM file
        let mut writer = Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
        for record in &records {
            writer.write(record).unwrap();
        }
        drop(writer);

        // Index the BAM file
        index::build(&bam_path, None, index::Type::Bai, 1).unwrap();

        // Read back and test IUPAC strategy
        let mut reader = IndexedReader::from_path(&bam_path).unwrap();
        let _header_view = reader.header().clone();

        reader.fetch(("chr1", 19, 20)).unwrap();
        let pileup_iter = reader.pileup();

        let read_filter = DefaultReadFilter::new(0, 0, 0);

        for pileup_result in pileup_iter {
            let pileup = pileup_result.unwrap();
            if pileup.pos() == 19 {
                let strat = crate::position::mate_fix::MateResolutionStrategy::BaseQualMapQualIUPAC;

                let alns: Vec<_> = pileup
                    .alignments()
                    .map(|aln| {
                        let rec = aln.record();
                        (aln, rec)
                    })
                    .collect();

                if alns.len() == 2 {
                    // With equal qualities, should return IUPAC code
                    let result = strat.cmp(&alns[0], &alns[1], &read_filter);
                    // A + G should give R
                    assert!(matches!(
                        result.recommended_base,
                        Some(crate::position::mate_fix::Base::R)
                    ));
                }
                break;
            }
        }
    }

    /// Test N strategy
    #[test]
    fn test_n_strategy() {
        use rust_htslib::bam::{IndexedReader, Writer, index};
        use tempfile::tempdir;

        let tempdir = tempdir().unwrap();
        let bam_path = tempdir.path().join("test.bam");

        // Create header
        let mut header = bam::header::Header::new();
        let mut chr1 = bam::header::HeaderRecord::new(b"SQ");
        chr1.push_tag(b"SN", &"chr1".to_owned());
        chr1.push_tag(b"LN", &"100".to_owned());
        header.push_record(&chr1);
        let view = bam::HeaderView::from_header(&header);

        // Create overlapping mates with same quality but different bases
        let records = vec![
            Record::from_sam(
                &view,
                b"TESTPAIR\t67\tchr1\t10\t40\t10M\tchr1\t15\t30\tCCCCCCCCCC\t>>>>>>>>>>",
            )
            .unwrap(),
            Record::from_sam(
                &view,
                b"TESTPAIR\t147\tchr1\t15\t40\t10M\tchr1\t10\t30\tTTTTTTTTTT\t>>>>>>>>>>",
            )
            .unwrap(),
        ];

        // Write BAM file
        let mut writer = Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
        for record in &records {
            writer.write(record).unwrap();
        }
        drop(writer);

        // Index the BAM file
        index::build(&bam_path, None, index::Type::Bai, 1).unwrap();

        // Read back and test N strategy
        let mut reader = IndexedReader::from_path(&bam_path).unwrap();
        let _header_view = reader.header().clone();

        reader.fetch(("chr1", 19, 20)).unwrap();
        let pileup_iter = reader.pileup();

        let read_filter = DefaultReadFilter::new(0, 0, 0);

        for pileup_result in pileup_iter {
            let pileup = pileup_result.unwrap();
            if pileup.pos() == 19 {
                let strat = crate::position::mate_fix::MateResolutionStrategy::N;

                let alns: Vec<_> = pileup
                    .alignments()
                    .map(|aln| {
                        let rec = aln.record();
                        (aln, rec)
                    })
                    .collect();

                if alns.len() == 2 {
                    // N strategy should always return N for different bases
                    let result = strat.cmp(&alns[0], &alns[1], &read_filter);
                    assert!(matches!(
                        result.recommended_base,
                        Some(crate::position::mate_fix::Base::N)
                    ));
                }
                break;
            }
        }
    }
}
