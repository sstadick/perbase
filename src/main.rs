extern crate perbase_lib;
pub mod commands;
use anyhow::Result;
use argh::FromArgs;
use commands::*;
use env_logger::Env;
use log::*;

#[derive(FromArgs)]
/// Commands for generating per-base analysis
struct Args {
    #[argh(subcommand)]
    subcommand: Subcommand,
}

#[derive(FromArgs)]
#[argh(subcommand)]
enum Subcommand {
    SimpleDepth(simple_depth::SimpleDepth),
}

impl Subcommand {
    fn run(self) -> Result<()> {
        match self {
            Subcommand::SimpleDepth(x) => x.run()?,
        }
        Ok(())
    }
}

fn main() -> Result<()> {
    env_logger::from_env(Env::default().default_filter_or("info")).init();
    if let Err(err) = argh::from_env::<Args>().subcommand.run() {
        error!("{}", err);
        std::process::exit(1);
    }
    Ok(())
}
