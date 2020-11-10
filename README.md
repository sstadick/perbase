<p align="center">
<img width="500" height="250" src="./perbase.png">
</p>

![Publish](https://github.com/sstadick/perbase/workflows/Publish/badge.svg)
![Rust](https://github.com/sstadick/perbase/workflows/Rust/badge.svg)
[![API docs](https://img.shields.io/badge/API-documentation-blue.svg)](https://docs.rs/perbase)
[![Crates.io](https://img.shields.io/crates/v/perbase.svg)](https://crates.io/crates/perbase)
[![Conda](https://anaconda.org/anaconda/anaconda/badges/installer/conda.svg)](https://anaconda.org/bioconda/perbase)

A highly parallelized utility for analyzing metrics at a per-base level.

If a metric is missing, or performance is lacking. Please file a bug/feature ticket in issues.

## Why?

Why `perbase` when so many other tools are out there? `perbase` leverages Rust's concurrency system to automagically parallelize over your input regions. This leads to orders of magnitude faster runtimes that scale with the compute that you have available. Additionally, `perbase` aims to be more accurate than other tools. E.g.: `perbase` counts DELs toward depth, `bam-readcount` does not, `perbase` does not count REF_SKIPs toward depth, `sambamba` does.

## Installation

```bash
conda install -c bioconda perbase
# or via the rust toolchain
cargo install perbase
```

\* Check version, conda lags behind what you will find in the release page

You can also download a binary from the [releases](https://github.com/sstadick/perbase/releases) page.

## Tools

### base-depth

The `base-depth` tool walks over every position in the BAM/CRAM file and calculates the depth, as well as the number of each nucleotide at the given position. Additionally, it counts the numbers of Ins/Dels at each position.

The output columns are as follows:

| Column         | Description                                                                                        |
| -------------- | -------------------------------------------------------------------------------------------------- |
| REF            | The reference sequence name                                                                        |
| POS            | The position on the reference sequence                                                             |
| REF_BASE       | The reference base at the position, column excluded if no reference was supplied                   |
| DEPTH          | The total depth at the position SUM(A, C, T, G, DEL)                                               |
| A              | Total A nucleotides seen at this position                                                          |
| C              | Total C nucleotides seen at this position                                                          |
| G              | Total G nucleotides seen at this position                                                          |
| T              | Total T nucleotides seen at this position                                                          |
| N              | Total N nucleotides seen at this position                                                          |
| INS            | Total insertions that start at the base to the right of this position                              |
| DEL            | Total deletions covering this position                                                             |
| REF_SKIP       | Total reference skip operations covering this position                                             |
| FAIL           | Total reads failing filters that covered this position (their bases were not counted toward depth) |
| NEAR_MAX_DEPTH | Flag to indicate if this position came within 1% of the max depth specified                        |

```bash
perbase base-depth ./test/test.bam
```

Example output

```text
REF     POS     REF_BASE        DEPTH   A       C       G       T       N       INS     DEL     REF_SKIP        FAIL    NEAR_MAX_DEPTH
chr1    709636  T       16      0       0       0       16      0       0       0       0       0   false
chr1    709637  T       16      0       4       0       12      0       0       0       0       0   false
chr1    709638  A       16      16      0       0       0       0       0       0       0       0   false
chr1    709639  G       16      0       0       16      0       0       0       0       0       0   false
chr1    709640  A       16      16      0       0       0       0       0       0       0       0   false
chr1    709641  A       16      16      0       0       0       0       0       0       0       0   false
chr1    709642  G       16      0       0       16      0       0       0       0       0       0   false
chr1    709643  G       16      0       0       16      0       0       0       0       0       0   false
chr1    709644  T       16      0       0       0       16      0       0       0       0       0   false
chr1    709645  G       16      0       0       16      0       0       0       0       0       0   false
```

If the `--mate-fix` flag is passed, each position will first check if there are any mate overlaps and choose the mate with the hightest MAPQ, breaking ties by choosing the first mate that passes filters. Mates that are discarded are not counted toward `FAIL` or `DEPTH`.

If the `--reference-fasta` is supplied, the `REF_BASE` field will be filled in. The reference must be indexed an match the BAM/CRAM header of the input.

The output can be compressed and indexed as follows:

```bash
perbase base-depth ./test/test.bam | bgzip > output.tsv.gz
tabix -S 1 -s 1 -b 2 -e 2 ./output.tsv.gz
# Query all positions overlapping region
tabix output.tsv.gz chr1:5-10
```

Usage:

```text
Calculate the depth at each base, per-nucleotide

USAGE:
    perbase base-depth [FLAGS] [OPTIONS] <reads>

FLAGS:
    -h, --help         Prints help information
    -m, --mate-fix     Fix overlapping mates counts, see docs for full details
    -V, --version      Prints version information
    -z, --zero-base    Output positions as 0-based instead of 1-based

OPTIONS:
    -B, --bcf-file <bcf-file>                A BCF/VCF file containing positions of interest. If specified, only bases
                                             from the given positions will be reported on
    -b, --bed-file <bed-file>                A BED file containing regions of interest. If specified, only bases from
                                             the given regions will be reported on
    -c, --chunksize <chunksize>              The ideal number of basepairs each worker receives. Total bp in memory at
                                             one time is (threads - 2) * chunksize
    -F, --exclude-flags <exclude-flags>      SAM flags to exclude, recommended 3848 [default: 0]
    -f, --include-flags <include-flags>      SAM flags to include [default: 0]
    -D, --max-depth <max-depth>              Set the max depth for a pileup. If a positions depth is within 1% of max-
                                             depth the `NEAR_MAX_DEPTH` output field will be set to true and that
                                             position should be viewed as suspect [default: 100000]
    -q, --min-mapq <min-mapq>                Minimum MAPQ for a read to count toward depth [default: 0]
    -o, --output <output>                    Output path, defaults to stdout
        --ref-cache-size <ref-cache-size>    Number of Reference Sequences to hold in memory at one time. Smaller will
                                             decrease mem usage [default: 10]
    -r, --ref-fasta <ref-fasta>              Indexed reference fasta, set if using CRAM
    -t, --threads <threads>                  The number of threads to use [default: 16]

ARGS:
    <reads>    Input indexed BAM/CRAM to analyze

```

### only-depth

The `only-depth` tool walks over the input BAM/CRAM file and calculates the depth over all positions specified by either a BED file or in the BAM/CRAM header. Adjacent positions that have the same depth will be merged together to form a non-inclusive range (see example output).

There are two distinct modes that `only-depth` can run in, gated by the `--fast-mode` flag. When running in fast-mode, only depth over the area a read covers is only determined by the reads start and end postions, and no cigar related info is taken into account. `--mate-fix` may still be used in this mode, and areas where mates overlap will not be counted twice.

Without the `--fast-mode` flag, the depth at each position is determined in a manner similar to `base-depth` where `DEL` will count toward depth, but `REF_SKIP` will not. Additionally, any reads that fail the `--exclude-flags` will not be counted toward depth. Lastly, `--mate-fix` can be applied to avoid counting regions twice where mates may overlap.

Regarding mate fixes, `perbase` will make "fixes" based only on the counted regions in a read. For example, if you have a read that goes from "chr1:0-1000" with a CIGAR of "25M974N1M", and the mate aligns nicely at "chr1:45-70" with CIGAR "25M", the mate will count toward the depth over "chr1:45-74". This is in contrast to other tools that will reject the mate even though it overlaps a region of R1 that is not counted toward depth.

For the fastest possible output, use `only-depth --fast-mode`.

**Note** that it is possible that two adjacent positions may not merge if they fall at a `--chunksize` boundary. If this is an issue you can set the `--chunksize` to the size of the largest contig in question. At a future date this may be fixed or a post processing tool may be provided to fix it. For most use cases this should not be a problem. Additionally, you can pipe into `merge-adjcent` which will fix it as well. EX: `perbase only-depth -m file.bam | perbase merge-adjacent > out.tsv`.

Example output of `perbase only-depth --mate-fix --zero-base  ./test/test.bam`:

```text
REF     POS     END     DEPTH
chr2    0       4       1
chr2    4       9       2
chr2    9       12      3
chr2    12      14      2
chr2    14      17      3
chr2    17      19      4
chr2    19      23      5
chr2    23      34      4
chr2    34      39      3
chr2    39      49      1
chr2    49      54      2
chr2    54      64      3
chr2    64      74      4
chr2    74      79      3
chr2    79      84      2
chr2    84      89      1
```

Usage:

```text
Calculate the only the depth at each base

USAGE:
    perbase only-depth [FLAGS] [OPTIONS] <reads>

FLAGS:
    -x, --fast-mode    Calculate depth based only on read starts/stops, see docs for full details
    -h, --help         Prints help information
    -m, --mate-fix     Fix overlapping mates counts, see docs for full details
    -n, --no-merge     Skip merging adjacent bases that have the same depth
    -V, --version      Prints version information
    -z, --zero-base    Output positions as 0-based instead of 1-based

OPTIONS:
    -B, --bcf-file <bcf-file>              A BCF/VCF file containing positions of interest. If specified, only bases
                                           from the given positions will be reported on. Note that it may be more
                                           efficient to calculate depth over regions if your positions are clustered
                                           tightly together
    -b, --bed-file <bed-file>              A BED file containing regions of interest. If specified, only bases from the
                                           given regions will be reported on
    -c, --chunksize <chunksize>            The ideal number of basepairs each worker receives. Total bp in memory at one
                                           time is (threads - 2) * chunksize
    -F, --exclude-flags <exclude-flags>    SAM flags to exclude, recommended 3848 [default: 0]
    -f, --include-flags <include-flags>    SAM flags to include [default: 0]
    -q, --min-mapq <min-mapq>              Minimum MAPQ for a read to count toward depth [default: 0]
    -o, --output <output>                  Output path, defaults to stdout
    -r, --ref-fasta <ref-fasta>            Indexed reference fasta, set if using CRAM
    -t, --threads <threads>                The number of threads to use [default: 16]

ARGS:
    <reads>    Input indexed BAM/CRAM to analyze
e
```

## merge-adjacent

`merge-adjacent` is a utility to merge overlapping regions in a BED-like file.

It will take a file with four columns and no header as long as the columns are like:

```text
<contig>\t<start>\t<stop>\t<depth>\n
```

Or it can take files with three columns with headers that are like

```text
<REF|chrom>\t<POS|chromStart>\t<END|chromEnd>\t<DEPTH|COV>
```

The `END|chromEnd` column is optional.

```text
Merge adjacent intervals that have the same depth. Input must be sorted like: `sort -k1,1 -k2,2n in.bed > in.sorted.bed`

Generally accepts any file with no header tha is <chrom>\t<start>\t<stop>\t<depth>. The <stop> is optional. See
documentation for explaination of headers that are accepted.

USAGE:
    perbase merge-adjacent [FLAGS] [OPTIONS] [in-file]

FLAGS:
    -h, --help
            Prints help information

    -n, --no-header
            Indicate if the input file does not have a header

    -V, --version
            Prints version information


OPTIONS:
    -o, --output <output>
            The output location, defaults to STDOUT


ARGS:
    <in-file>
            Input bed-like file, defaults to STDIN
```

EX:

```bash
perbase only-depth indexed.bam | perbase merge-adjacent > out.tsv
```

## Similar Projects

- [`sambamba depth`](https://github.com/biod/sambamba/wiki/%5Bsambamba-depth%5D-documentation)
- [`samtools depth`](http://www.htslib.org/doc/samtools-depth.html)
- [`mosdepth`](https://github.com/brentp/mosdepth)
- [`bam-readcount`](https://github.com/genome/bam-readcount)
