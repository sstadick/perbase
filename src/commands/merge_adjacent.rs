use anyhow::Result;
use log::*;
use perbase_lib::utils;
use serde::{Deserialize, Serialize};
use smartstring::alias::*;
use std::path::PathBuf;
use structopt::StructOpt;

/// Merge adjacent intervals that have the same depth. Input must be sorted like:
/// `sort -k1,1 -k2,2n in.bed > in.sorted.bed`
///
/// Generally accepts any file with no header tha is <chrom>\t<start>\t<stop>\t<depth>.
/// The <stop> is optional. See documentation for explaination of headers that are accepted.
#[derive(StructOpt)]
#[structopt(author)]
pub struct MergeAdjacent {
    /// Input bed-like file, defaults to STDIN.
    in_file: Option<PathBuf>,

    /// The number of threads to use for compressing output (specified by --bgzip)
    #[structopt(long, short = "T", default_value = utils::NUM_CPU.as_str())]
    compression_threads: usize,

    /// The level to use for compressing output (specified by --bgzip)
    #[structopt(long, short = "L", default_value = "2")]
    compression_level: u32,

    /// Indicate if the input file does not have a header
    #[structopt(long, short = "n")]
    no_header: bool,

    /// Optionally bgzip the output.
    #[structopt(long, short = "Z")]
    bgzip: bool,

    /// The output location, defaults to STDOUT
    #[structopt(long, short = "o")]
    output: Option<PathBuf>,
}

impl MergeAdjacent {
    pub fn run(&self) -> Result<()> {
        let mut reader = utils::get_reader(
            &self.in_file,
            !self.no_header,
            self.in_file
                .as_ref()
                .map(utils::is_bgzipped)
                .unwrap_or(false),
        )?;
        let mut writer = utils::get_writer(
            &self.output,
            self.bgzip,
            true,
            self.compression_threads,
            self.compression_level,
        )?;
        let mut iter = reader.deserialize().map(|r| {
            let rec: BedLike = r.expect("Deserialzied record");
            rec
        });

        // iterate over input records
        if let Some(first) = iter.next() {
            let mut curr_seq = first.ref_seq.clone();
            let mut curr_start = first.pos;
            let mut curr_end = if let Some(end) = first.end {
                end
            } else {
                first.pos + 1
            };
            let mut curr_depth = first.depth;

            for bedlike in iter {
                if bedlike.depth == curr_depth
                    && curr_end >= bedlike.pos
                    && bedlike.ref_seq == curr_seq
                {
                    curr_end = if let Some(end) = bedlike.end {
                        end
                    } else {
                        bedlike.pos + 1
                    };
                } else {
                    writer.serialize(BedLike {
                        ref_seq: curr_seq.clone(),
                        pos: curr_start,
                        end: Some(curr_end),
                        depth: curr_depth,
                    })?;
                    curr_seq = bedlike.ref_seq;
                    curr_start = bedlike.pos;
                    curr_end = if let Some(end) = bedlike.end {
                        end
                    } else {
                        bedlike.pos + 1
                    };
                    curr_depth = bedlike.depth;
                }
            }
            // Write last record
            writer.serialize(BedLike {
                ref_seq: curr_seq.clone(),
                pos: curr_start,
                end: Some(curr_end),
                depth: curr_depth,
            })?;
        } else {
            info!("No input records");
        }

        Ok(())
    }
}

#[derive(Serialize, Deserialize, Debug)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
struct BedLike {
    /// The contig name
    #[serde(alias = "REF", alias = "chrom", rename = "REF")]
    ref_seq: String, // Smartstring
    /// The start position
    #[serde(alias = "chromStart")]
    pos: usize,
    /// The stop position, optional, will be converted to start+1 if missing
    #[serde(alias = "chromEnd")]
    end: Option<usize>,
    /// The depth of this interval
    #[serde(alias = "COV")]
    depth: usize,
}
