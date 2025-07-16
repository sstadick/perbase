//! Check for mate fixing.
//!
//! The strategy here is as follows:
//! - Always check first for the user-based filter status of each mate, if one fails, go with the other mate, if both fail, default to `a`
//! - Try one of the [`MateResolutionStrategies`]
//!
//! ## MateResolutionStrategy::BaseQualMapQualFirstInPair
//! If either is indel / not a base -> MAPQ -> first in pair
//! Else BASEQ -> MAPQ -> first in pair
//!
//! ## MateResolutionStrategy::BaseQualMapQualIUPAC
//! If either is indel / not a base -> MAPQ -> first in pair
//! Else BASEQ -> MAPQ -> IAUPAC
//!
//! ## MateResolutionStrategy::BaseQualMapQualN
//! If either is indel / not a base -> MAPQ -> first in pair
//! Else BASEQ -> MAPQ -> N
//!
//! ## MateResolutionStrategy::MapQualBaseQualFirstInPair
//! If either is indel / not a base -> MAPQ -> first in pair
//! Else MAPQ -> BASEQ-> first in pair
//!
//! ## MateResolutionStrategy::MapQualBaseQualIUPAC
//! If either is indel / not a base -> MAPQ -> first in pair
//! Else MAPQ -> BASEQ -> IAUPAC
//!
//! ## MateResolutionStrategy::BaseQualMapQualN
//! If either is indel / not a base -> MAPQ -> first in pair
//! Else MAPQ -> BASEQ -> N
//!
//! ## MateResolutionStrategy::IUPAC
//! If either is indel / not a base -> MAPQ -> first in pair
//! Else IUPAC
//!
//! ## MateResolutionStrategy::N
//! If either is indel / not a base -> MAPQ -> first in pair
//! Else N
//!
//! ## MateResolutionStrategy::Original
//! MAPQ -> first in pair
use crate::read_filter::ReadFilter;
use rust_htslib::bam::{self, pileup::Alignment, record::Record};

use std::cmp::Ordering;

#[inline]
pub(crate) fn is_indel(indel: bam::pileup::Indel) -> bool {
    matches!(
        indel,
        bam::pileup::Indel::Ins(_) | bam::pileup::Indel::Del(_)
    )
}

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

#[derive(Debug, Copy, Clone)]
pub struct MateResolution {
    /// The ordering of a and b
    pub(crate) ordering: Ordering,
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

#[derive(Debug, Copy, Clone)]
pub enum MateResolutionStrategy {
    BaseQualMapQualFirstInPair,
    BaseQualMapQualIUPAC,
    BaseQualMapQualN,
    MapQualBaseQualFirstInPair,
    MapQualBaseQualIUPAC,
    MapQualBaseQualN,
    IUPAC,
    N,
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
        let a_is_not_base = a.0.is_refskip() || is_indel(a.0.indel());
        let b_is_not_base = b.0.is_refskip() || is_indel(b.0.indel());
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
        let a_is_not_base = a.0.is_refskip() || is_indel(a.0.indel());
        let b_is_not_base = b.0.is_refskip() || is_indel(b.0.indel());
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
        let a_is_not_base = a.0.is_refskip() || is_indel(a.0.indel());
        let b_is_not_base = b.0.is_refskip() || is_indel(b.0.indel());
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
        let a_is_not_base = a.0.is_refskip() || is_indel(a.0.indel());
        let b_is_not_base = b.0.is_refskip() || is_indel(b.0.indel());
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
        let a_is_not_base = a.0.is_refskip() || is_indel(a.0.indel());
        let b_is_not_base = b.0.is_refskip() || is_indel(b.0.indel());
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
        let a_is_not_base = a.0.is_refskip() || is_indel(a.0.indel());
        let b_is_not_base = b.0.is_refskip() || is_indel(b.0.indel());
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
        let a_is_not_base = a.0.is_refskip() || is_indel(a.0.indel());
        let b_is_not_base = b.0.is_refskip() || is_indel(b.0.indel());
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
