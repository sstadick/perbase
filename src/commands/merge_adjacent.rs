use anyhow::Result;
use csv;
use grep_cli::stdout;
use log::*;
use serde::{Deserialize, Serialize};
use smartstring::alias::*;
use std::{
    fs::File,
    io::{BufReader, BufWriter, Read, Write},
    path::PathBuf,
};
use structopt::StructOpt;
use termcolor::ColorChoice;

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

    /// Indicate if the input file does not have a header
    #[structopt(long, short = "n")]
    no_header: bool,

    /// The output location, defaults to STDOUT
    #[structopt(long, short = "o")]
    output: Option<PathBuf>,
}

impl MergeAdjacent {
    pub fn run(&self) -> Result<()> {
        let mut reader = self.get_reader()?;
        let mut writer = self.get_writer()?;
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

            while let Some(bedlike) = iter.next() {
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

    /// Open a CSV Reader from file or stdin
    fn get_reader(&self) -> Result<csv::Reader<Box<dyn Read>>> {
        let raw_reader: Box<dyn Read> = match &self.in_file {
            Some(path) if path.to_str().unwrap() != "-" => {
                Box::new(BufReader::new(File::open(path)?))
            }
            _ => Box::new(std::io::stdin()),
        };

        Ok(csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(!self.no_header)
            .from_reader(raw_reader))
    }

    /// Open a CSV Writer to a file or stdout
    fn get_writer(&self) -> Result<csv::Writer<Box<dyn Write>>> {
        let raw_writer: Box<dyn Write> = match &self.output {
            Some(path) if path.to_str().unwrap() != "-" => {
                Box::new(BufWriter::new(File::open(path)?))
            }
            _ => Box::new(stdout(ColorChoice::Never)),
        };
        Ok(csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(raw_writer))
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
