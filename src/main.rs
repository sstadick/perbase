extern crate perbase_lib;
pub mod commands;
use anyhow::Result;
use commands::*;
use env_logger::Env;
use log::*;
use structopt::StructOpt;

#[derive(StructOpt)]
#[structopt(rename_all = "kebab-case", author, about)]
/// Commands for generating per-base analysis
struct Args {
    #[structopt(subcommand)]
    subcommand: Subcommand,
}

#[derive(StructOpt)]
enum Subcommand {
    BaseDepth(base_depth::BaseDepth),
    OnlyDepth(only_depth::OnlyDepth),
    MergeAdjacent(merge_adjacent::MergeAdjacent),
}

impl Subcommand {
    fn run(self) -> Result<()> {
        match self {
            Subcommand::BaseDepth(x) => x.run()?,
            Subcommand::OnlyDepth(x) => x.run()?,
            Subcommand::MergeAdjacent(x) => x.run()?,
        }
        Ok(())
    }
}

fn main() -> Result<()> {
    env_logger::from_env(Env::default().default_filter_or("info")).init();
    if let Err(err) = Args::from_args().subcommand.run() {
        error!("{}", err);
        std::process::exit(1);
    }
    Ok(())
}
