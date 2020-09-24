![Publish](https://github.com/sstadick/perbase/workflows/Publish/badge.svg)
![Rust](https://github.com/sstadick/perbase/workflows/Rust/badge.svg)
# perbase

A highly parallelized utility for analyzing metrics at a per-base level.

If a metric is missing, or performance is lacking. Please file a bug/feature issue.

## Installation

```bash
conda install -c bioconda perbase
```

## Tools

### simple-depth

The `simple-depth` tool walks over every position in the BAM/CRAM file and calculates the depth, as well as the number of each nucleotide at the given position. Additionally, it counts the numbers of Ins/Dels at each position.

| Column   | Description                                                                                        |
| -------- | -------------------------------------------------------------------------------------------------- |
| REF      | The reference sequence name                                                                        |
| POS      | The 1-based position on the reference sequence                                                     |
| DEPTH    | The total depth at the position SUM(A, C, T, G, DEL, REF_SKIP)                                     |
| A        | Total A nucleotides seen at this position                                                          |
| C        | Total C nucleotides seen at this position                                                          |
| G        | Total G nucleotides seen at this position                                                          |
| T        | Total T nucleotides seen at this position                                                          |
| N        | Total N nucleotides seen at this position                                                          |
| INS      | Total insertions that start at the base to the right of this position                              |
| DEL      | Total deletions covering this position                                                             |
| REF_SKIP | Total reference skip operations covering this position                                             |
| FAIL     | Total reads failing filters that covered this position (their bases were not counted toward depth) |

The output is 1-based.

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

The output can be compressed and indexed as follows:

```bash
perbase simple-depth ./test/test.bam | bgzip > output.tsv.gz
tabix -S 1 -s 1 -b 2 -e 2 ./output.tsv.gz
# Query all positions overlapping region
tabix output.tsv.gz chr1:5-10
```

## TODOs

- [ ] Add more metrics to match `bam-readcount` as an `indepth` tool
- [ ] Add a strictly depth calculation a la `mosdepth` as an `onlydepth`
- [ ] Add bgzip output / auto tabix indexing support
- [ ] Support limiting inputs with an interval_list file
