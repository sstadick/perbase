//! Check for mate fixing.
//!
//! The strategy here is as follows:
//! - Always check first for the user-based filter status of each mate, if one fails, go with the other mate, if both fail, default to `a`
//! - Try one of the [`MateResolutionStrategy`] variants
//!
//! ## MateResolutionStrategy::BaseQualMapQualFirstInPair
//! If either is deletion / ref skip -> MAPQ -> first in pair
//! Else BASEQ -> MAPQ -> first in pair
//!
//! ## MateResolutionStrategy::BaseQualMapQualIUPAC
//! If either is deletion / ref skip -> MAPQ -> first in pair
//! Else BASEQ -> MAPQ -> IAUPAC
//!
//! ## MateResolutionStrategy::BaseQualMapQualN
//! If either is deletion / ref skip -> MAPQ -> first in pair
//! Else BASEQ -> MAPQ -> N
//!
//! ## MateResolutionStrategy::MapQualBaseQualFirstInPair
//! If either is deletion / ref skip -> MAPQ -> first in pair
//! Else MAPQ -> BASEQ-> first in pair
//!
//! ## MateResolutionStrategy::MapQualBaseQualIUPAC
//! If either is deletion / ref skip -> MAPQ -> first in pair
//! Else MAPQ -> BASEQ -> IUPAC
//!
//! ## MateResolutionStrategy::MapQualBaseQualN
//! If either is deletion / ref skip -> MAPQ -> first in pair
//! Else MAPQ -> BASEQ -> N
//!
//! ## MateResolutionStrategy::IUPAC
//! If either is deletion / ref skip -> MAPQ -> first in pair
//! Else IUPAC
//!
//! ## MateResolutionStrategy::N
//! If either is deletion / ref skip -> MAPQ -> first in pair
//! Else N
//!
//! ## MateResolutionStrategy::Original
//! MAPQ -> first in pair
use crate::read_filter::ReadFilter;
use rust_htslib::bam::{pileup::Alignment, record::Record};
use strum::EnumString;

use std::cmp::Ordering;

/// IUPAC Bases
#[derive(Copy, Clone, Debug)]
pub enum Base {
    /// A
    A,
    /// C
    C,
    /// T or U
    T,
    /// G
    G,
    /// N (any base)
    N,
    /// A or G
    R,
    /// C or T
    Y,
    /// G or C
    S,
    /// A or T
    W,
    /// G or T
    K,
    /// A or C
    M,
    /// C or G or T
    B,
    /// A or G or T
    D,
    /// A or C or T
    H,
    /// A or C or G
    V,
    //Gap, // Technically a gap but not needed for us
}

impl Base {
    /// Maps uncertainty of bases
    pub(crate) fn either_or<const DEFAULT_TO_N: bool>(a: Self, b: Self) -> Self {
        match (a, b) {
            (Self::A, Self::A) => Self::A,
            (Self::C, Self::C) => Self::C,
            (Self::G, Self::G) => Self::G,
            (Self::T, Self::T) => Self::T,
            (Self::A, Self::G) | (Self::G, Self::A) if !DEFAULT_TO_N => Self::R,
            (Self::C, Self::T) | (Self::T, Self::C) if !DEFAULT_TO_N => Self::Y,
            (Self::G, Self::C) | (Self::C, Self::G) if !DEFAULT_TO_N => Self::S,
            (Self::A, Self::T) | (Self::T, Self::A) if !DEFAULT_TO_N => Self::W,
            (Self::G, Self::T) | (Self::T, Self::G) if !DEFAULT_TO_N => Self::K,
            (Self::A, Self::C) | (Self::C, Self::A) if !DEFAULT_TO_N => Self::M,
            (_, _) => Self::N,
        }
    }
}

impl From<char> for Base {
    /// Lossy conversion, we don't support reading IUPAC bases, only outputting them for now
    ///
    /// If we wanted to support them, then we get crazy combos of either_or
    #[inline]
    fn from(value: char) -> Self {
        match value.to_ascii_uppercase() {
            'A' => Self::A,
            'C' => Self::C,
            'T' | 'U' => Self::T,
            'G' => Self::G,
            _ => Self::N,
        }
    }
}

/// Result of mate resolution comparison containing ordering and optional base recommendation.
///
/// This struct represents the decision made when comparing two overlapping mate reads
/// at the same genomic position. It contains both the ordering decision (which mate
/// to prefer) and an optional base recommendation for ambiguous cases.
#[derive(Debug, Copy, Clone)]
pub struct MateResolution {
    /// The ordering of mate a relative to mate b (Greater = choose a, Less = choose b)
    pub(crate) ordering: Ordering,
    /// Optional base recommendation for cases where both mates contribute to the final call
    pub(crate) recommended_base: Option<Base>,
}

impl MateResolution {
    pub(crate) fn new(ordering: Ordering, recommended_base: Option<Base>) -> Self {
        Self {
            ordering,
            recommended_base,
        }
    }
}

/// Strategies for resolving which mate to use when overlapping mates disagree on a base call.
///
/// All strategies first check user-based read filters. If one mate fails filters, the other
/// is chosen. If both fail, the first mate is chosen by default.
///
/// For reads that are deletions / ref skips or lack a base call, all strategies fall back to the Original
/// strategy (MAPQ → first in pair).
#[derive(Debug, Copy, Clone, EnumString)]
#[strum(ascii_case_insensitive)]
pub enum MateResolutionStrategy {
    /// Base quality → MAPQ → first in pair
    ///
    /// Priority order:
    /// 1. Higher base quality score wins
    /// 2. If base qualities are equal, higher MAPQ wins  
    /// 3. If both are equal, first mate in pair wins
    BaseQualMapQualFirstInPair,

    /// Base quality → MAPQ → IUPAC ambiguity code
    ///
    /// Priority order:
    /// 1. Higher base quality score wins
    /// 2. If base qualities are equal, higher MAPQ wins
    /// 3. If both are equal, return IUPAC ambiguity code (e.g., A+G→R)
    BaseQualMapQualIUPAC,

    /// Base quality → MAPQ → N
    ///
    /// Priority order:
    /// 1. Higher base quality score wins
    /// 2. If base qualities are equal, higher MAPQ wins
    /// 3. If both are equal, return N (unknown base)
    BaseQualMapQualN,

    /// MAPQ → base quality → first in pair
    ///
    /// Priority order:
    /// 1. Higher MAPQ wins
    /// 2. If MAPQ is equal, higher base quality wins
    /// 3. If both are equal, first mate in pair wins
    MapQualBaseQualFirstInPair,

    /// MAPQ → base quality → IUPAC ambiguity code
    ///
    /// Priority order:
    /// 1. Higher MAPQ wins
    /// 2. If MAPQ is equal, higher base quality wins
    /// 3. If both are equal, return IUPAC ambiguity code (e.g., A+G→R)
    MapQualBaseQualIUPAC,

    /// MAPQ → base quality → N
    ///
    /// Priority order:
    /// 1. Higher MAPQ wins
    /// 2. If MAPQ is equal, higher base quality wins
    /// 3. If both are equal, return N (unknown base)
    MapQualBaseQualN,

    /// Always return IUPAC ambiguity code for different bases
    ///
    /// Returns the appropriate IUPAC code for the base combination:
    /// - Same bases: return the base (A+A→A)
    /// - Different bases: return IUPAC code (A+G→R, C+T→Y, etc.)
    IUPAC,

    /// Always return N for different bases
    ///
    /// Returns N for different bases, but identical bases return themselves:
    /// - Same bases: return the base (A+A→A)
    /// - Different bases: return N
    N,

    /// Original strategy: MAPQ → first in pair
    ///
    /// Priority order:
    /// 1. Higher MAPQ wins
    /// 2. If MAPQ is equal, first mate in pair wins
    /// 3. If neither is marked as first in pair, choose first mate by default
    Original,
}

impl MateResolutionStrategy {
    pub(crate) fn cmp<F: ReadFilter>(
        &self,
        a: &(Alignment<'_>, Record),
        b: &(Alignment<'_>, Record),
        read_filter: &F,
    ) -> MateResolution {
        // Handle user-set filters.
        let a_pass = read_filter.filter_read(&a.1, Some(&a.0));
        let b_pass = read_filter.filter_read(&b.1, Some(&b.0));
        if a_pass && !b_pass {
            return MateResolution::new(Ordering::Greater, None);
        } else if b_pass && !a_pass {
            return MateResolution::new(Ordering::Less, None);
        } else if !a_pass && !b_pass {
            return MateResolution::new(Ordering::Greater, None);
        }

        match self {
            MateResolutionStrategy::BaseQualMapQualFirstInPair => {
                Self::base_qual_map_qual_first_in_pair(a, b)
            }
            MateResolutionStrategy::BaseQualMapQualIUPAC => Self::base_qual_map_qual_iupac(a, b),
            MateResolutionStrategy::BaseQualMapQualN => Self::base_qual_map_qual_n(a, b),
            MateResolutionStrategy::MapQualBaseQualFirstInPair => {
                Self::map_qual_base_qual_first_in_pair(a, b)
            }
            MateResolutionStrategy::MapQualBaseQualIUPAC => Self::map_qual_base_qual_iupac(a, b),
            MateResolutionStrategy::MapQualBaseQualN => Self::map_qual_base_qual_n(a, b),
            MateResolutionStrategy::IUPAC => Self::resolve_base::<false>(a, b),
            MateResolutionStrategy::N => Self::resolve_base::<true>(a, b),
            MateResolutionStrategy::Original => Self::original(a, b),
        }
    }

    pub(crate) fn resolve_base<const DEFAULT_TO_N: bool>(
        a: &(Alignment<'_>, Record),
        b: &(Alignment<'_>, Record),
    ) -> MateResolution {
        // First check that we have a base
        let a_is_not_base = a.0.qpos().is_none();
        let b_is_not_base = b.0.qpos().is_none();
        if a_is_not_base || b_is_not_base {
            return Self::original(a, b);
        }

        let a_base = Base::from(a.1.seq()[a.0.qpos().unwrap()] as char);
        let b_base = Base::from(b.1.seq()[b.0.qpos().unwrap()] as char);
        MateResolution::new(
            Ordering::Greater,
            Some(Base::either_or::<DEFAULT_TO_N>(a_base, b_base)),
        )
    }

    pub(crate) fn base_qual_map_qual_n(
        a: &(Alignment<'_>, Record),
        b: &(Alignment<'_>, Record),
    ) -> MateResolution {
        // First check that we have a base
        let a_is_not_base = a.0.qpos().is_none();
        let b_is_not_base = b.0.qpos().is_none();
        if a_is_not_base || b_is_not_base {
            return Self::original(a, b);
        }

        // compare baseq -> mapq -> default to first in pair
        let a_qual = a.1.qual()[a.0.qpos().unwrap()];
        let b_qual = b.1.qual()[b.0.qpos().unwrap()];

        match a_qual.cmp(&b_qual) {
            Ordering::Greater => MateResolution::new(Ordering::Greater, None),
            Ordering::Less => MateResolution::new(Ordering::Less, None),
            Ordering::Equal => match a.1.mapq().cmp(&b.1.mapq()) {
                Ordering::Greater => MateResolution::new(Ordering::Greater, None),
                Ordering::Less => MateResolution::new(Ordering::Less, None),
                Ordering::Equal => MateResolution::new(Ordering::Greater, Some(Base::N)), // how to pass this back?
            },
        }
    }

    pub(crate) fn base_qual_map_qual_iupac(
        a: &(Alignment<'_>, Record),
        b: &(Alignment<'_>, Record),
    ) -> MateResolution {
        // First check that we have a base
        let a_is_not_base = a.0.qpos().is_none();
        let b_is_not_base = b.0.qpos().is_none();
        if a_is_not_base || b_is_not_base {
            return Self::original(a, b);
        }

        // compare baseq -> mapq -> default to first in pair
        let a_qual = a.1.qual()[a.0.qpos().unwrap()];
        let b_qual = b.1.qual()[b.0.qpos().unwrap()];

        match a_qual.cmp(&b_qual) {
            Ordering::Greater => MateResolution::new(Ordering::Greater, None),
            Ordering::Less => MateResolution::new(Ordering::Less, None),
            Ordering::Equal => match a.1.mapq().cmp(&b.1.mapq()) {
                Ordering::Greater => MateResolution::new(Ordering::Greater, None),
                Ordering::Less => MateResolution::new(Ordering::Less, None),
                Ordering::Equal => {
                    let a_base = Base::from(a.1.seq()[a.0.qpos().unwrap()] as char);
                    let b_base = Base::from(b.1.seq()[b.0.qpos().unwrap()] as char);
                    MateResolution::new(
                        Ordering::Greater,
                        Some(Base::either_or::<false>(a_base, b_base)),
                    )
                }
            },
        }
    }

    pub(crate) fn base_qual_map_qual_first_in_pair(
        a: &(Alignment<'_>, Record),
        b: &(Alignment<'_>, Record),
    ) -> MateResolution {
        // First check that we have a base
        let a_is_not_base = a.0.qpos().is_none();
        let b_is_not_base = b.0.qpos().is_none();
        if a_is_not_base || b_is_not_base {
            return Self::original(a, b);
        }

        // compare baseq -> mapq -> default to first in pair
        let a_qual = a.1.qual()[a.0.qpos().unwrap()];
        let b_qual = b.1.qual()[b.0.qpos().unwrap()];

        match a_qual.cmp(&b_qual) {
            Ordering::Greater => MateResolution::new(Ordering::Greater, None),
            Ordering::Less => MateResolution::new(Ordering::Less, None),
            Ordering::Equal => Self::original(a, b),
        }
    }

    pub(crate) fn map_qual_base_qual_n(
        a: &(Alignment<'_>, Record),
        b: &(Alignment<'_>, Record),
    ) -> MateResolution {
        // First check that we have a base
        let a_is_not_base = a.0.qpos().is_none();
        let b_is_not_base = b.0.qpos().is_none();
        if a_is_not_base || b_is_not_base {
            return Self::original(a, b);
        }

        match a.1.mapq().cmp(&b.1.mapq()) {
            Ordering::Greater => MateResolution::new(Ordering::Greater, None),
            Ordering::Less => MateResolution::new(Ordering::Less, None),
            Ordering::Equal => {
                let a_qual = a.1.qual()[a.0.qpos().unwrap()];
                let b_qual = b.1.qual()[b.0.qpos().unwrap()];

                match a_qual.cmp(&b_qual) {
                    Ordering::Greater => MateResolution::new(Ordering::Greater, None),
                    Ordering::Less => MateResolution::new(Ordering::Less, None),
                    Ordering::Equal => MateResolution::new(Ordering::Greater, Some(Base::N)), // how to pass this back?
                }
            }
        }
    }

    pub(crate) fn map_qual_base_qual_iupac(
        a: &(Alignment<'_>, Record),
        b: &(Alignment<'_>, Record),
    ) -> MateResolution {
        // First check that we have a base
        let a_is_not_base = a.0.qpos().is_none();
        let b_is_not_base = b.0.qpos().is_none();
        if a_is_not_base || b_is_not_base {
            return Self::original(a, b);
        }

        match a.1.mapq().cmp(&b.1.mapq()) {
            Ordering::Greater => MateResolution::new(Ordering::Greater, None),
            Ordering::Less => MateResolution::new(Ordering::Less, None),
            Ordering::Equal => {
                let a_qual = a.1.qual()[a.0.qpos().unwrap()];
                let b_qual = b.1.qual()[b.0.qpos().unwrap()];

                match a_qual.cmp(&b_qual) {
                    Ordering::Greater => MateResolution::new(Ordering::Greater, None),
                    Ordering::Less => MateResolution::new(Ordering::Less, None),
                    Ordering::Equal => {
                        let a_base = Base::from(a.1.seq()[a.0.qpos().unwrap()] as char);
                        let b_base = Base::from(b.1.seq()[b.0.qpos().unwrap()] as char);
                        MateResolution::new(
                            Ordering::Greater,
                            Some(Base::either_or::<false>(a_base, b_base)),
                        )
                    }
                }
            }
        }
    }

    pub(crate) fn map_qual_base_qual_first_in_pair(
        a: &(Alignment<'_>, Record),
        b: &(Alignment<'_>, Record),
    ) -> MateResolution {
        // First check that we have a base
        let a_is_not_base = a.0.qpos().is_none();
        let b_is_not_base = b.0.qpos().is_none();
        if a_is_not_base || b_is_not_base {
            return Self::original(a, b);
        }

        match a.1.mapq().cmp(&b.1.mapq()) {
            Ordering::Greater => MateResolution::new(Ordering::Greater, None),
            Ordering::Less => MateResolution::new(Ordering::Less, None),
            Ordering::Equal => {
                let a_qual = a.1.qual()[a.0.qpos().unwrap()];
                let b_qual = b.1.qual()[b.0.qpos().unwrap()];

                match a_qual.cmp(&b_qual) {
                    Ordering::Greater => MateResolution::new(Ordering::Greater, None),
                    Ordering::Less => MateResolution::new(Ordering::Less, None),
                    Ordering::Equal => {
                        if a.1.flags() & 64 != 0 {
                            MateResolution::new(Ordering::Greater, None)
                        } else if b.1.flags() & 64 != 0 {
                            MateResolution::new(Ordering::Less, None)
                        } else {
                            // Default to `a` in the event that there is no first in pair for some reason
                            MateResolution::new(Ordering::Greater, None)
                        }
                    }
                }
            }
        }
    }

    /// Whichever has higher MAPQ, or if equal, whichever is first in pair
    pub(crate) fn original(
        a: &(Alignment<'_>, Record),
        b: &(Alignment<'_>, Record),
    ) -> MateResolution {
        match a.1.mapq().cmp(&b.1.mapq()) {
            Ordering::Greater => MateResolution::new(Ordering::Greater, None),
            Ordering::Less => MateResolution::new(Ordering::Less, None),
            Ordering::Equal => {
                // Check if a is first in pair
                if a.1.flags() & 64 != 0 {
                    MateResolution::new(Ordering::Greater, None)
                } else if b.1.flags() & 64 != 0 {
                    MateResolution::new(Ordering::Less, None)
                } else {
                    // Default to `a` in the event that there is no first in pair for some reason
                    MateResolution::new(Ordering::Greater, None)
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{self, record::Record};

    // Since we can't directly create Alignment objects, we need to work around this
    // by testing the internal functions directly

    // Base enum tests
    #[test]
    fn test_base_either_or_identical() {
        assert!(matches!(
            Base::either_or::<false>(Base::A, Base::A),
            Base::A
        ));
        assert!(matches!(
            Base::either_or::<false>(Base::C, Base::C),
            Base::C
        ));
        assert!(matches!(
            Base::either_or::<false>(Base::G, Base::G),
            Base::G
        ));
        assert!(matches!(
            Base::either_or::<false>(Base::T, Base::T),
            Base::T
        ));
    }

    #[test]
    fn test_base_either_or_iupac() {
        // Test IUPAC combinations
        assert!(matches!(
            Base::either_or::<false>(Base::A, Base::G),
            Base::R
        ));
        assert!(matches!(
            Base::either_or::<false>(Base::G, Base::A),
            Base::R
        ));

        assert!(matches!(
            Base::either_or::<false>(Base::C, Base::T),
            Base::Y
        ));
        assert!(matches!(
            Base::either_or::<false>(Base::T, Base::C),
            Base::Y
        ));

        assert!(matches!(
            Base::either_or::<false>(Base::G, Base::C),
            Base::S
        ));
        assert!(matches!(
            Base::either_or::<false>(Base::C, Base::G),
            Base::S
        ));

        assert!(matches!(
            Base::either_or::<false>(Base::A, Base::T),
            Base::W
        ));
        assert!(matches!(
            Base::either_or::<false>(Base::T, Base::A),
            Base::W
        ));

        assert!(matches!(
            Base::either_or::<false>(Base::G, Base::T),
            Base::K
        ));
        assert!(matches!(
            Base::either_or::<false>(Base::T, Base::G),
            Base::K
        ));

        assert!(matches!(
            Base::either_or::<false>(Base::A, Base::C),
            Base::M
        ));
        assert!(matches!(
            Base::either_or::<false>(Base::C, Base::A),
            Base::M
        ));
    }

    #[test]
    fn test_base_either_or_default_to_n() {
        // When DEFAULT_TO_N is true, all non-identical should return N
        assert!(matches!(Base::either_or::<true>(Base::A, Base::G), Base::N));
        assert!(matches!(Base::either_or::<true>(Base::C, Base::T), Base::N));
        assert!(matches!(Base::either_or::<true>(Base::G, Base::C), Base::N));

        // Identical bases should still return themselves
        assert!(matches!(Base::either_or::<true>(Base::A, Base::A), Base::A));
        assert!(matches!(Base::either_or::<true>(Base::C, Base::C), Base::C));
    }

    #[test]
    fn test_base_from_char() {
        // Uppercase
        assert!(matches!(Base::from('A'), Base::A));
        assert!(matches!(Base::from('C'), Base::C));
        assert!(matches!(Base::from('T'), Base::T));
        assert!(matches!(Base::from('G'), Base::G));

        // Lowercase
        assert!(matches!(Base::from('a'), Base::A));
        assert!(matches!(Base::from('c'), Base::C));
        assert!(matches!(Base::from('t'), Base::T));
        assert!(matches!(Base::from('g'), Base::G));

        // U -> T
        assert!(matches!(Base::from('U'), Base::T));
        assert!(matches!(Base::from('u'), Base::T));

        // Unknown -> N
        assert!(matches!(Base::from('X'), Base::N));
        assert!(matches!(Base::from('?'), Base::N));
        assert!(matches!(Base::from('1'), Base::N));
    }

    // Since we can't easily test the strategy functions that take Alignment parameters,
    // we'll add integration tests with the actual BAM reading in pileup_position.rs

    // Test that we can at least construct the types
    #[test]
    fn test_mate_resolution_construction() {
        let res = MateResolution::new(Ordering::Greater, Some(Base::N));
        assert_eq!(res.ordering, Ordering::Greater);
        assert!(matches!(res.recommended_base, Some(Base::N)));
    }

    // Helper function to create test BAM data and test a strategy
    #[allow(clippy::too_many_arguments)]
    fn test_strategy_with_bam(
        strategy: MateResolutionStrategy,
        read1_seq: &[u8],
        read1_qual: &[u8],
        read1_mapq: u8,
        read1_flags: u16,
        read2_seq: &[u8],
        read2_qual: &[u8],
        read2_mapq: u8,
        read2_flags: u16,
    ) -> MateResolution {
        use rust_htslib::bam::{HeaderView, IndexedReader, Read, Writer, index};
        use tempfile::tempdir;

        let tempdir = tempdir().unwrap();
        let bam_path = tempdir.path().join("test.bam");

        // Create header
        let mut header = bam::header::Header::new();
        let mut chr1 = bam::header::HeaderRecord::new(b"SQ");
        chr1.push_tag(b"SN", &"chr1".to_owned());
        chr1.push_tag(b"LN", &"100".to_owned());
        header.push_record(&chr1);
        let view = HeaderView::from_header(&header);

        // Create overlapping mate pair
        // Read 1: positions 10-19 (1-based: 11-20)
        // Read 2: positions 15-24 (1-based: 16-25)
        // They overlap at positions 15-19 (1-based: 16-20)
        let records = vec![
            Record::from_sam(
                &view,
                format!(
                    "TESTPAIR\t{}\tchr1\t11\t{}\t10M\tchr1\t16\t30\t{}\t{}",
                    read1_flags,
                    read1_mapq,
                    std::str::from_utf8(read1_seq).unwrap(),
                    std::str::from_utf8(read1_qual).unwrap()
                )
                .as_bytes(),
            )
            .unwrap(),
            Record::from_sam(
                &view,
                format!(
                    "TESTPAIR\t{}\tchr1\t16\t{}\t10M\tchr1\t11\t30\t{}\t{}",
                    read2_flags,
                    read2_mapq,
                    std::str::from_utf8(read2_seq).unwrap(),
                    std::str::from_utf8(read2_qual).unwrap()
                )
                .as_bytes(),
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

        // Read back and test strategy at overlapping position
        let mut reader = IndexedReader::from_path(&bam_path).unwrap();
        reader.fetch(("chr1", 15, 20)).unwrap(); // Fetch the overlapping region
        let pileup_iter = reader.pileup();

        let read_filter = crate::read_filter::DefaultReadFilter::new(0, 0, 0);

        for pileup_result in pileup_iter {
            let pileup = pileup_result.unwrap();
            if pileup.pos() >= 15 && pileup.pos() < 20 {
                // Check overlapping positions
                let alns: Vec<_> = pileup
                    .alignments()
                    .map(|aln| {
                        let rec = aln.record();
                        (aln, rec)
                    })
                    .collect();

                if alns.len() == 2 {
                    return strategy.cmp(&alns[0], &alns[1], &read_filter);
                }
            }
        }

        panic!("Failed to find overlapping reads at test position");
    }

    #[test]
    fn test_base_qual_map_qual_first_in_pair() {
        // Test different base qualities - higher base quality should win
        let result = test_strategy_with_bam(
            MateResolutionStrategy::BaseQualMapQualFirstInPair,
            b"AAAAAAAAAA",
            b"##########",
            30,
            67, // First mate: A, qual 2, MAPQ 30, first in pair
            b"CCCCCCCCCC",
            b">>>>>>>>>>",
            30,
            147, // Second mate: C, qual 30, MAPQ 30, second in pair
        );
        assert_eq!(result.ordering, Ordering::Less); // Second mate wins due to higher base quality
        assert!(result.recommended_base.is_none());

        // Test equal base qualities, different MAPQ - higher MAPQ should win
        let result = test_strategy_with_bam(
            MateResolutionStrategy::BaseQualMapQualFirstInPair,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            40,
            67, // First mate: A, qual 30, MAPQ 40, first in pair
            b"CCCCCCCCCC",
            b">>>>>>>>>>",
            30,
            147, // Second mate: C, qual 30, MAPQ 30, second in pair
        );
        assert_eq!(result.ordering, Ordering::Greater); // First mate wins due to higher MAPQ
        assert!(result.recommended_base.is_none());

        // Test equal base qualities and MAPQ - first in pair should win
        let result = test_strategy_with_bam(
            MateResolutionStrategy::BaseQualMapQualFirstInPair,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            30,
            67, // First mate: A, qual 30, MAPQ 30, first in pair
            b"CCCCCCCCCC",
            b">>>>>>>>>>",
            30,
            147, // Second mate: C, qual 30, MAPQ 30, second in pair
        );
        assert_eq!(result.ordering, Ordering::Greater); // First mate wins due to being first in pair
        assert!(result.recommended_base.is_none());
    }

    #[test]
    fn test_base_qual_map_qual_iupac() {
        // Test different base qualities - higher base quality should win
        let result = test_strategy_with_bam(
            MateResolutionStrategy::BaseQualMapQualIUPAC,
            b"AAAAAAAAAA",
            b"##########",
            30,
            67, // First mate: A, qual 2, MAPQ 30
            b"CCCCCCCCCC",
            b">>>>>>>>>>",
            30,
            147, // Second mate: C, qual 30, MAPQ 30
        );
        assert_eq!(result.ordering, Ordering::Less); // Second mate wins due to higher base quality
        assert!(result.recommended_base.is_none());

        // Test equal base qualities and MAPQ - should return IUPAC code
        let result = test_strategy_with_bam(
            MateResolutionStrategy::BaseQualMapQualIUPAC,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            30,
            67, // First mate: A, qual 30, MAPQ 30
            b"GGGGGGGGGG",
            b">>>>>>>>>>",
            30,
            147, // Second mate: G, qual 30, MAPQ 30
        );
        assert_eq!(result.ordering, Ordering::Greater); // Chooses first by default
        assert!(matches!(result.recommended_base, Some(Base::R))); // A + G = R
    }

    #[test]
    fn test_base_qual_map_qual_n() {
        // Test different base qualities - higher base quality should win
        let result = test_strategy_with_bam(
            MateResolutionStrategy::BaseQualMapQualN,
            b"AAAAAAAAAA",
            b"##########",
            30,
            67, // First mate: A, qual 2, MAPQ 30
            b"CCCCCCCCCC",
            b">>>>>>>>>>",
            30,
            147, // Second mate: C, qual 30, MAPQ 30
        );
        assert_eq!(result.ordering, Ordering::Less); // Second mate wins due to higher base quality
        assert!(result.recommended_base.is_none());

        // Test equal base qualities and MAPQ - should return N
        let result = test_strategy_with_bam(
            MateResolutionStrategy::BaseQualMapQualN,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            30,
            67, // First mate: A, qual 30, MAPQ 30
            b"CCCCCCCCCC",
            b">>>>>>>>>>",
            30,
            147, // Second mate: C, qual 30, MAPQ 30
        );
        assert_eq!(result.ordering, Ordering::Greater); // Chooses first by default
        assert!(matches!(result.recommended_base, Some(Base::N))); // Returns N
    }

    #[test]
    fn test_map_qual_base_qual_first_in_pair() {
        // Test different MAPQ - higher MAPQ should win
        let result = test_strategy_with_bam(
            MateResolutionStrategy::MapQualBaseQualFirstInPair,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            40,
            67, // First mate: A, qual 30, MAPQ 40
            b"CCCCCCCCCC",
            b">>>>>>>>>>",
            30,
            147, // Second mate: C, qual 30, MAPQ 30
        );
        assert_eq!(result.ordering, Ordering::Greater); // First mate wins due to higher MAPQ
        assert!(result.recommended_base.is_none());

        // Test equal MAPQ, different base qualities - higher base quality should win
        let result = test_strategy_with_bam(
            MateResolutionStrategy::MapQualBaseQualFirstInPair,
            b"AAAAAAAAAA",
            b"##########",
            30,
            67, // First mate: A, qual 2, MAPQ 30
            b"CCCCCCCCCC",
            b">>>>>>>>>>",
            30,
            147, // Second mate: C, qual 30, MAPQ 30
        );
        assert_eq!(result.ordering, Ordering::Less); // Second mate wins due to higher base quality
        assert!(result.recommended_base.is_none());

        // Test equal MAPQ and base qualities - first in pair should win
        let result = test_strategy_with_bam(
            MateResolutionStrategy::MapQualBaseQualFirstInPair,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            30,
            67, // First mate: A, qual 30, MAPQ 30, first in pair
            b"CCCCCCCCCC",
            b">>>>>>>>>>",
            30,
            147, // Second mate: C, qual 30, MAPQ 30, second in pair
        );
        assert_eq!(result.ordering, Ordering::Greater); // First mate wins due to being first in pair
        assert!(result.recommended_base.is_none());
    }

    #[test]
    fn test_map_qual_base_qual_iupac() {
        // Test different MAPQ - higher MAPQ should win
        let result = test_strategy_with_bam(
            MateResolutionStrategy::MapQualBaseQualIUPAC,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            40,
            67, // First mate: A, qual 30, MAPQ 40
            b"CCCCCCCCCC",
            b">>>>>>>>>>",
            30,
            147, // Second mate: C, qual 30, MAPQ 30
        );
        assert_eq!(result.ordering, Ordering::Greater); // First mate wins due to higher MAPQ
        assert!(result.recommended_base.is_none());

        // Test equal MAPQ and base qualities - should return IUPAC code
        let result = test_strategy_with_bam(
            MateResolutionStrategy::MapQualBaseQualIUPAC,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            30,
            67, // First mate: A, qual 30, MAPQ 30
            b"GGGGGGGGGG",
            b">>>>>>>>>>",
            30,
            147, // Second mate: G, qual 30, MAPQ 30
        );
        assert_eq!(result.ordering, Ordering::Greater); // Chooses first by default
        assert!(matches!(result.recommended_base, Some(Base::R))); // A + G = R
    }

    #[test]
    fn test_map_qual_base_qual_n() {
        // Test different MAPQ - higher MAPQ should win
        let result = test_strategy_with_bam(
            MateResolutionStrategy::MapQualBaseQualN,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            40,
            67, // First mate: A, qual 30, MAPQ 40
            b"CCCCCCCCCC",
            b">>>>>>>>>>",
            30,
            147, // Second mate: C, qual 30, MAPQ 30
        );
        assert_eq!(result.ordering, Ordering::Greater); // First mate wins due to higher MAPQ
        assert!(result.recommended_base.is_none());

        // Test equal MAPQ and base qualities - should return N
        let result = test_strategy_with_bam(
            MateResolutionStrategy::MapQualBaseQualN,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            30,
            67, // First mate: A, qual 30, MAPQ 30
            b"CCCCCCCCCC",
            b">>>>>>>>>>",
            30,
            147, // Second mate: C, qual 30, MAPQ 30
        );
        assert_eq!(result.ordering, Ordering::Greater); // Chooses first by default
        assert!(matches!(result.recommended_base, Some(Base::N))); // Returns N
    }

    #[test]
    fn test_iupac_strategy() {
        // Test with different bases - should return IUPAC code
        let result = test_strategy_with_bam(
            MateResolutionStrategy::IUPAC,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            30,
            67, // First mate: A
            b"GGGGGGGGGG",
            b">>>>>>>>>>",
            40,
            147, // Second mate: G
        );
        assert_eq!(result.ordering, Ordering::Greater); // Always chooses first read for ordering
        assert!(matches!(result.recommended_base, Some(Base::R))); // A + G = R

        // Test with same bases - should return that base
        let result = test_strategy_with_bam(
            MateResolutionStrategy::IUPAC,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            30,
            67, // First mate: A
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            40,
            147, // Second mate: A
        );
        assert_eq!(result.ordering, Ordering::Greater); // Always chooses first read for ordering
        assert!(matches!(result.recommended_base, Some(Base::A))); // A + A = A
    }

    #[test]
    fn test_n_strategy() {
        // Test with different bases - should always return N
        let result = test_strategy_with_bam(
            MateResolutionStrategy::N,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            30,
            67, // First mate: A
            b"CCCCCCCCCC",
            b">>>>>>>>>>",
            40,
            147, // Second mate: C
        );
        assert_eq!(result.ordering, Ordering::Greater); // Always chooses first read for ordering
        assert!(matches!(result.recommended_base, Some(Base::N))); // Always N

        // Test with same bases - should return the actual base since they're identical
        let result = test_strategy_with_bam(
            MateResolutionStrategy::N,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            30,
            67, // First mate: A
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            40,
            147, // Second mate: A
        );
        assert_eq!(result.ordering, Ordering::Greater); // Always chooses first read for ordering
        assert!(matches!(result.recommended_base, Some(Base::A))); // Identical bases return themselves
    }

    #[test]
    fn test_original_strategy() {
        // Test different MAPQ - higher MAPQ should win
        let result = test_strategy_with_bam(
            MateResolutionStrategy::Original,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            40,
            67, // First mate: A, MAPQ 40, first in pair
            b"CCCCCCCCCC",
            b">>>>>>>>>>",
            30,
            147, // Second mate: C, MAPQ 30, second in pair
        );
        assert_eq!(result.ordering, Ordering::Greater); // First mate wins due to higher MAPQ
        assert!(result.recommended_base.is_none());

        // Test equal MAPQ - first in pair should win
        let result = test_strategy_with_bam(
            MateResolutionStrategy::Original,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            30,
            67, // First mate: A, MAPQ 30, first in pair
            b"CCCCCCCCCC",
            b">>>>>>>>>>",
            30,
            147, // Second mate: C, MAPQ 30, second in pair
        );
        assert_eq!(result.ordering, Ordering::Greater); // First mate wins due to being first in pair
        assert!(result.recommended_base.is_none());

        // Test equal MAPQ, second read is first in pair
        let result = test_strategy_with_bam(
            MateResolutionStrategy::Original,
            b"AAAAAAAAAA",
            b">>>>>>>>>>",
            30,
            147, // First mate: A, MAPQ 30, second in pair
            b"CCCCCCCCCC",
            b">>>>>>>>>>",
            30,
            67, // Second mate: C, MAPQ 30, first in pair
        );
        assert_eq!(result.ordering, Ordering::Less); // Second mate wins due to being first in pair
        assert!(result.recommended_base.is_none());
    }
}
