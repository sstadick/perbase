//! General utility methods.
use anyhow::{Error, Result};
use lazy_static::lazy_static;
use log::*;
use num_cpus;
use rayon;

/// Set rayon global thread pool size
pub fn set_rayon_global_pools_size(size: usize) -> Result<()> {
    let cpus = determine_allowed_cpus(size)?;
    rayon::ThreadPoolBuilder::new()
        .num_threads(cpus)
        .build_global()?;
    Ok(())
}

/// Check that specified `desired` is valid
pub fn determine_allowed_cpus(desired: usize) -> Result<usize> {
    if desired <= 0 {
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

lazy_static! {
    /// Return the number of cpus as an &str
    pub static ref NUM_CPU: String = num_cpus::get().to_string();
}
