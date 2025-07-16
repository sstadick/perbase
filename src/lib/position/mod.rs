//! A set of implementations of `Position` for different use cases
mod mate_fix;
pub mod pileup_position;
pub mod range_positions;

use serde::Serialize;
use smartstring::alias::String;
/// A serializable object meant to hold all information about a position.
pub trait Position: Default + Serialize {
    /// Create a new position with all other values zeroed
    fn new(ref_seq: String, pos: u32) -> Self;
}
