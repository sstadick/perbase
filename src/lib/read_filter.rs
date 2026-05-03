//! A trait and default implementation of a read filter.
use rust_htslib::bam::record::Record;

/// Minimal read-level data needed by [`ReadFilter`] implementations.
pub trait ReadView {
    /// SAM flags for the read.
    fn flags(&self) -> u16;

    /// Mapping quality for the read.
    fn mapq(&self) -> u8;
}

impl ReadView for Record {
    #[inline(always)]
    fn flags(&self) -> u16 {
        self.flags()
    }

    #[inline(always)]
    fn mapq(&self) -> u8 {
        self.mapq()
    }
}

impl<T: ReadView + ?Sized> ReadView for std::rc::Rc<T> {
    #[inline(always)]
    fn flags(&self) -> u16 {
        self.as_ref().flags()
    }

    #[inline(always)]
    fn mapq(&self) -> u8 {
        self.as_ref().mapq()
    }
}

/// Anything that implements ReadFilter can apply a filter set to read.
pub trait ReadFilter {
    /// filters a read, true is pass, false if fail
    fn filter_read<R: ReadView + ?Sized>(&self, read: &R) -> bool;
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
    fn filter_read<R: ReadView + ?Sized>(&self, read: &R) -> bool {
        let flags = read.flags();
        (!flags) & self.include_flags == 0
            && flags & self.exclude_flags == 0
            && read.mapq() >= self.min_mapq
    }
}
