//! General utility methods.
use anyhow::{Error, Result};
use grep_cli::stdout;
use gzp::{BgzfSyncReader, Compression, ZBuilder, deflate::Bgzf};
use lazy_static::lazy_static;
use log::{error, warn};
use rust_htslib::bam;
use std::{
    cell::RefCell,
    ffi::OsStr,
    fs::File,
    io::{self, BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
};
use termcolor::ColorChoice;

thread_local! {
    static INDEXED_READER: RefCell<Option<(PathBuf, Option<PathBuf>, bam::IndexedReader)>> = const { RefCell::new(None) };
}

/// Use a lazily opened indexed BAM/CRAM reader local to the calling thread.
///
/// Changing either input path replaces the cached reader. Callers must fetch their
/// region and finish consuming records or pileups within the callback. The reader
/// is closed when replaced or when the thread exits; callbacks must not reenter
/// this function on the same thread.
///
/// # Panics
///
/// Panics if the input cannot be opened, the reference cannot be set, or the
/// callback reenters this function.
pub fn with_indexed_reader<T>(
    reads: &Path,
    ref_fasta: Option<&Path>,
    f: impl FnOnce(&mut bam::IndexedReader) -> T,
) -> T {
    INDEXED_READER.with(|cell| {
        let mut slot = cell.borrow_mut();
        if slot
            .as_ref()
            .is_none_or(|(path, reference, _)| path != reads || reference.as_deref() != ref_fasta)
        {
            let mut reader =
                bam::IndexedReader::from_path(reads).expect("Indexed Reader for region");
            if let Some(reference) = ref_fasta {
                reader.set_reference(reference).expect("Set ref");
            }
            *slot = Some((
                reads.to_path_buf(),
                ref_fasta.map(Path::to_path_buf),
                reader,
            ));
        }
        f(&mut slot.as_mut().unwrap().2)
    })
}

/// Set rayon global thread pool size.
///
/// # Errors
///
/// - [`Error`] if an issue is encountered determening the allowed cpus.
pub fn set_rayon_global_pools_size(size: usize) -> Result<()> {
    let cpus = determine_allowed_cpus(size)?;
    rayon::ThreadPoolBuilder::new()
        .num_threads(cpus)
        .build_global()?;
    Ok(())
}

/// Check if err is a broken pipe.
/// Check if err is a broken pipe.
#[inline]
pub fn is_broken_pipe(err: &Error) -> bool {
    if let Some(io_err) = err.root_cause().downcast_ref::<io::Error>()
        && io_err.kind() == io::ErrorKind::BrokenPipe
    {
        return true;
    }
    false
}

/// Check that specified `desired` is valid.
///
/// If more threads are specified than available, the max available are used.
///
/// # Errors
///
/// - [`Error`] if less than or equal to 0 threads are selected
pub fn determine_allowed_cpus(desired: usize) -> Result<usize> {
    if desired == 0 {
        error!("Must select > 0 threads");
        Err(Error::msg("Too few threads selected. Min 4"))
    } else if desired > num_cpus::get() {
        let cpus = num_cpus::get();
        warn!("Specified more threads than are available, using {}", cpus);
        Ok(cpus)
    } else {
        Ok(desired)
    }
}

/// Detect if a path path ends with the usual bgzip extension.
pub fn is_bgzipped<P: AsRef<Path>>(path: P) -> bool {
    let ext = path.as_ref().extension().unwrap_or_else(|| OsStr::new(""));
    ext == "gz" || ext == "gzip" || ext == "bgzf"
}

/// Get a CSV Reader
///
/// # Errors
///
/// - If unable to open file
///
/// # Panics
///
/// - If unable to parse input path
pub fn get_reader<P: AsRef<Path>>(
    path: &Option<P>,
    has_headers: bool,
    bgzipped: bool,
) -> Result<csv::Reader<Box<dyn Read>>> {
    let raw_reader: Box<dyn Read> = match &path {
        Some(path) if path.as_ref().to_str().unwrap() != "-" => {
            let reader = BufReader::new(File::open(path)?);
            if bgzipped {
                Box::new(BgzfSyncReader::new(reader))
            } else {
                Box::new(reader)
            }
        }
        _ => {
            let reader = std::io::stdin();
            if bgzipped {
                Box::new(BgzfSyncReader::new(reader))
            } else {
                Box::new(reader)
            }
        }
    };

    Ok(csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(has_headers)
        .from_reader(raw_reader))
}

/// Open a CSV Writer to a file or stdout.
///
/// # Errors
///
/// - If file can't be created
///
/// # Panics
///
/// - If unable to parse input path
pub fn get_writer<P: AsRef<Path>>(
    path: &Option<P>,
    bgzipped: bool,
    write_headers: bool,
    threads: usize,
    compression_level: u32,
) -> Result<csv::Writer<Box<dyn Write>>> {
    let raw_writer: Box<dyn Write> = match &path {
        Some(path) if path.as_ref().to_str().unwrap() != "-" => {
            let writer = BufWriter::new(File::create(path)?);
            if bgzipped {
                Box::new(
                    ZBuilder::<Bgzf, _>::new()
                        .num_threads(threads)
                        .compression_level(Compression::new(compression_level))
                        .from_writer(writer),
                )
            } else {
                Box::new(writer)
            }
        }
        _ => {
            let writer = stdout(ColorChoice::Never);
            if bgzipped {
                Box::new(
                    ZBuilder::<Bgzf, _>::new()
                        .num_threads(threads)
                        .compression_level(Compression::new(compression_level))
                        .from_writer(writer),
                )
            } else {
                Box::new(writer)
            }
        }
    };
    Ok(csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(write_headers)
        .from_writer(raw_writer))
}

lazy_static! {
    /// Return the number of cpus as an &str
    pub static ref NUM_CPU: String = num_cpus::get().to_string();
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::Read;

    #[test]
    fn reused_cram_reader_switches_references() {
        let dir = tempfile::tempdir().unwrap();
        let references = [dir.path().join("a.fa"), dir.path().join("c.fa")];
        for (path, base) in references.iter().zip(["A", "C"]) {
            std::fs::write(path, format!(">chr1\n{}\n", base.repeat(100))).unwrap();
            rust_htslib::faidx::build(path).unwrap();
        }
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 100)
                .push_tag(b"UR", references[0].to_str().unwrap()),
        );
        let path = dir.path().join("reads.cram");
        let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Cram).unwrap();
        writer.set_reference(&references[0]).unwrap();
        let record = bam::Record::from_sam(
            &bam::HeaderView::from_header(&header),
            b"read\t0\tchr1\t11\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
        )
        .unwrap();
        writer.write(&record).unwrap();
        drop(writer);
        bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();

        for (reference, expected) in [
            (Some(references[0].as_path()), b'A'),
            (Some(references[1].as_path()), b'C'),
            (Some(references[0].as_path()), b'A'),
        ] {
            let sequences = with_indexed_reader(&path, reference, |reader| {
                // Deliberately decode against different reference bases to make
                // reference switching observable, rather than just changing paths
                // to identical FASTAs. Only this test disables MD5 validation.
                reader
                    .set_cram_options(rust_htslib::htslib::hts_fmt_option_CRAM_OPT_IGNORE_MD5, 1)
                    .unwrap();
                reader.fetch((0, 10, 20)).unwrap();
                reader
                    .records()
                    .map(|record| record.unwrap().seq().as_bytes())
                    .collect::<Vec<_>>()
            });
            assert_eq!(sequences, vec![vec![expected; 10]]);
        }
    }

    #[test]
    fn reused_reader_fetches_regions_and_switches_inputs() {
        let dir = tempfile::tempdir().unwrap();
        let mut header = bam::Header::new();
        for name in ["chr1", "chr2"] {
            header.push_record(
                bam::header::HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", name)
                    .push_tag(b"LN", 100),
            );
        }
        let view = bam::HeaderView::from_header(&header);
        let paths = [dir.path().join("a.bam"), dir.path().join("b.bam")];
        for (path, base) in paths.iter().zip(["A", "T"]) {
            let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
            for (name, pos, seq) in [("chr1", 11, base), ("chr1", 31, "C"), ("chr2", 6, "G")] {
                let record = bam::Record::from_sam(
                    &view,
                    format!("read\t0\t{name}\t{pos}\t60\t1M\t*\t0\t0\t{seq}\tI").as_bytes(),
                )
                .unwrap();
                writer.write(&record).unwrap();
            }
            drop(writer);
            bam::index::build(path, None, bam::index::Type::Bai, 1).unwrap();
        }

        // Revisit regions backwards, after empty fetches and across contigs,
        // then switch files and return to the original input on the same thread.
        for (file, base) in [(0, b'A'), (1, b'T'), (0, b'A')] {
            for (tid, start, stop, expected) in [
                (1, 0, 10, vec![b'G']),
                (0, 30, 40, vec![b'C']),
                (0, 0, 5, vec![]),
                (0, 10, 20, vec![base]),
                (0, 10, 20, vec![base]),
            ] {
                let actual = with_indexed_reader(&paths[file], None, |reader| {
                    reader.fetch((tid, start, stop)).unwrap();
                    reader
                        .records()
                        .map(|record| record.unwrap().seq()[0])
                        .collect::<Vec<_>>()
                });
                assert_eq!(actual, expected);
                let positions = with_indexed_reader(&paths[file], None, |reader| {
                    reader.fetch((tid, start, stop)).unwrap();
                    reader
                        .pileup()
                        .map(|p| p.unwrap().pos())
                        .collect::<Vec<_>>()
                });
                let expected_positions = if expected.is_empty() {
                    vec![]
                } else {
                    vec![if tid == 1 { 5 } else { start }]
                };
                assert_eq!(positions, expected_positions);
            }
        }
    }
}
