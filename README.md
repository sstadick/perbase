![Publish](https://github.com/sstadick/perbase/workflows/Publish/badge.svg)
![Rust](https://github.com/sstadick/perbase/workflows/Rust/badge.svg)
# perbase

A highly parallelized utility for analyzing metrics at a per-base level.

If a metric is missing, or performance is lacking. Please file a bug/feature issue. This tool aims to improve upon [`bam-readcount`](https://github.com/genome/bam-readcount).

## Installation

```bash
conda install -c bioconda -c conda-forge perbase
```

\* Check version, conda lags behind what you will find in the release page

## Tools

### simple-depth

The `simple-depth` tool walks over every position in the BAM/CRAM file and calculates the depth, as well as the number of each nucleotide at the given position. Additionally, it counts the numbers of Ins/Dels at each position.

The output columns are as follows:

| Column   | Description                                                                                        |
| -------- | -------------------------------------------------------------------------------------------------- |
| REF      | The reference sequence name                                                                        |
| POS      | The 1-based position on the reference sequence                                                     |
| DEPTH    | The total depth at the position SUM(A, C, T, G, DEL)                                               |
| A        | Total A nucleotides seen at this position                                                          |
| C        | Total C nucleotides seen at this position                                                          |
| G        | Total G nucleotides seen at this position                                                          |
| T        | Total T nucleotides seen at this position                                                          |
| N        | Total N nucleotides seen at this position                                                          |
| INS      | Total insertions that start at the base to the right of this position                              |
| DEL      | Total deletions covering this position                                                             |
| REF_SKIP | Total reference skip operations covering this position                                             |
| FAIL     | Total reads failing filters that covered this position (their bases were not counted toward depth) |

```bash
perbase simple-depth ./test/test.bam
```

Example output

```text
REF     POS     DEPTH   A       C       G       T       N       INS     DEL     FAIL    REF_SKIP
chr2    1       1       1       0       0       0       0       0       0       0       0
chr2    2       1       1       0       0       0       0       1       0       0       0
chr2    3       1       0       0       1       0       0       0       0       0       0
chr2    4       1       1       0       0       0       0       0       0       0       0
chr2    5       2       2       0       0       0       0       0       0       0       0
chr2    6       2       2       0       0       0       0       0       0       0       0
chr2    7       2       1       0       0       0       0       0       1       0       0
chr2    8       2       1       0       0       0       0       0       1       0       0
chr2    9       2       1       0       0       0       0       0       1       0       0
chr2    10      3       2       0       0       0       0       0       1       0       0
chr2    11      3       2       0       0       0       0       0       1       0       0
chr2    12      3       3       0       0       0       0       0       0       0       0
chr2    13      3       2       0       0       0       0       0       0       0       1
chr2    14      3       2       0       0       0       0       0       0       0       1
chr2    15      4       3       0       0       0       0       0       0       0       1
chr2    16      4       2       0       0       1       0       0       0       0       1
chr2    17      4       3       0       0       0       0       0       0       0       1
```

If the `--mate-fix` flag is passed, each position will first check if there are any mate overlaps and choose the mate with the hightest MAPQ, breaking ties by choosing the first mate that passes filters. Using `--mate-fix` is slower by several orders of magnitude.

The output can be compressed and indexed as follows:

```bash
perbase simple-depth ./test/test.bam | bgzip > output.tsv.gz
tabix -S 1 -s 1 -b 2 -e 2 ./output.tsv.gz
# Query all positions overlapping region
tabix output.tsv.gz chr1:5-10
```

Usage:

```text
Usage: perbase simple-depth <reads> [-r <ref-fasta>] [-b <bed-file>] [-o <output>] [-t <threads>] [-c <chunksize>] [-f <include-flags>] [-F <exclude-flags>] [-m] [-q <min-mapq>]

Calculate the depth at each base, per-nucleotide. Takes an indexed BAM/CRAM as <reads>.

Options:
  -r, --ref-fasta   indexed reference fasta, set if using CRAM
  -b, --bed-file    a BED file containing regions of interest. If specified,
                    only bases from the given regions will be reported on
  -o, --output      output path. DEFAULT: stdout
  -t, --threads     the number of threads to use. DEFAULT: max_available
  -c, --chunksize   the ideal number of basepairs each worker receives. Total bp
                    in memory at one time is (threads - 2) * chunksize
  -f, --include-flags
                    SAM flags to include. DEFAULT: 0
  -F, --exclude-flags
                    SAM flags to exclude, recommended 3848. DEFAULT: 0
  -m, --mate-fix    fix overlapping mates counts, see docs for full details.
                    DEAFAULT: off
  -q, --min-mapq    minimum mapq for a read to count toward depth. DEFAULT: 0
  --help            display usage information
```

## TODOs

- [ ] Add more metrics to match `bam-readcount` as an `indepth` tool
- [ ] Add a strictly depth calculation a la `mosdepth` as an `onlydepth`
- [ ] Add bgzip output / auto tabix indexing support
