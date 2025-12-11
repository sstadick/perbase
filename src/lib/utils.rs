//! General utility methods.
use anyhow::{Error, Result};
use grep_cli::stdout;
use gzp::{BgzfSyncReader, Compression, ZBuilder, deflate::Bgzf};
use lazy_static::lazy_static;
use log::{error, warn};
use std::{
    ffi::OsStr,
    fs::File,
    io::{self, BufReader, BufWriter, Read, Write},
    path::Path,
};
use termcolor::ColorChoice;

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
