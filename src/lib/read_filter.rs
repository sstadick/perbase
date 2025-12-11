//! A trait and default implementation of a read filter.
use rust_htslib::bam::pileup::Alignment;
use rust_htslib::bam::record::Record;

/// Anything that implements ReadFilter can apply a filter set to read.
pub trait ReadFilter {
    /// filters a read, true is pass, false if fail
    fn filter_read(&self, read: &Record, alignment: Option<&Alignment>) -> bool;
}

/// A straightforward read filter.
pub struct DefaultReadFilter {
    include_flags: u16,
    exclude_flags: u16,
    min_mapq: u8,
}

impl DefaultReadFilter {
    /// Create an OnlyDepthReadFilter
    pub fn new(include_flags: u16, exclude_flags: u16, min_mapq: u8) -> Self {
        Self {
            include_flags,
            exclude_flags,
            min_mapq,
        }
    }
}

impl ReadFilter for DefaultReadFilter {
    /// Filter reads based SAM flags and mapping quality
    #[inline(always)]
    fn filter_read(&self, read: &Record, _alignment: Option<&Alignment>) -> bool {
        let flags = read.flags();
        (!flags) & self.include_flags == 0
            && flags & self.exclude_flags == 0
            && read.mapq() >= self.min_mapq
    }
}
