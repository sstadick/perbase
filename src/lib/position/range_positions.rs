use crate::position::Position;
use serde::Serialize;
use smartstring::alias::String;
use std::{default};

/// Hold all information about a range of positions.
#[derive(Debug, Serialize, Default)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub struct RangePositions {
    /// Reference sequence name.
    #[serde(rename = "REF")]
    pub ref_seq: String,
    /// 1-based position in the sequence.
    pub pos: usize,
    /// The point at which this depth ends, non-inclusive
    pub end: usize,
    /// Total depth at this position.
    pub depth: usize,
}

impl Position for RangePositions {
    /// Create a new position for the given ref_seq name.
    fn new(ref_seq: String, pos: usize) -> Self {
        RangePositions {
            ref_seq,
            pos,
            ..default::Default::default()
        }
    }
}