//! An implementation of `Position` that covers a range from `pos` to `end`.
//!
//! Whether this is 0-based or 1-based is up to the caller.
use crate::position::Position;
use serde::Serialize;
use smartstring::alias::String;
use std::default;

/// Hold all information about a range of positions.
#[derive(Debug, Serialize, Default)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub struct RangePositions {
    /// Reference sequence name.
    #[serde(rename = "REF")]
    pub ref_seq: String,
    /// 1-based position in the sequence.
    pub pos: u32,
    /// The point at which this depth ends, non-inclusive
    pub end: u32,
    /// Total depth at this position.
    pub depth: u32,
}

impl Position for RangePositions {
    /// Create a new position for the given ref_seq name.
    fn new(ref_seq: String, pos: u32) -> Self {
        RangePositions {
            ref_seq,
            pos,
            ..default::Default::default()
        }
    }
}

/// Bed-like formatting of [`RangePositions`]
#[derive(Debug, Serialize, Default)]
#[serde(rename_all = "camelCase")]
pub struct BedFormatRangePositions {
    /// Reference sequence name.
    pub chrom: String,
    /// 1-based position in the sequence.
    pub chrom_start: u32,
    /// The point at which this depth ends, non-inclusive
    pub chrom_end: u32,
    /// The name column of the bed file
    pub name: String,
    /// Total depth at this position.
    pub score: u32,
}

impl From<RangePositions> for BedFormatRangePositions {
    fn from(orig: RangePositions) -> Self {
        BedFormatRangePositions {
            chrom: orig.ref_seq,
            chrom_start: orig.pos,
            chrom_end: orig.end,
            name: String::from(""),
            score: orig.depth,
        }
    }
}
