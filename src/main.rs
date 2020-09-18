pub mod commands;
use anyhow::Result;
use argh::FromArgs;
use commands::per_base;
use env_logger::Env;
use log::*;
use std::path::PathBuf;

/// Calculate the depth at each base, per-base
#[derive(FromArgs)]
struct PerBase {
    /// input BAM/CRAM to analyze
    #[argh(positional)]
    reads: PathBuf,
    /// indexed reference fasta
    #[argh(option)]
    ref_fasta: PathBuf,
}

impl PerBase {
    fn run(self) -> Result<()> {
        per_base::parse(self.reads, self.ref_fasta)?;
        Ok(())
    }
}

fn main() -> Result<()> {
    env_logger::from_env(Env::default().default_filter_or("info")).init();
    if let Err(err) = argh::from_env::<PerBase>().run() {
        error!("{}", err);
        std::process::exit(1);
    }
    Ok(())
}
