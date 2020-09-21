# perbase

A utility for analyzing metrics at a per-base level.

This includes many small tools optimized for their particular tasks. The goal is to provide fast tools for all per-base level calculations. If a metric is missing, or performance is lacking. Please file a bug/feature issue.

## Tools

### simple-depth

The `simple-depth` tool walks over every position in the BAM/CRAM file and calculates the depth, as well as the number of each nucleotide at the given position. Additionally, it counts the numbers of Ins/Dels at each position. These counts have been validated against the output of `bam-readcount` and should be an easy replacement.

If you need to filter the reads based on SAM flags, the current method is to use samtools and pipe the output of a `samtools view` command into this tool. EX:

```bash
samtools view -hu -F 3844 test/test.bam | perbase simple-bam
```

The same method can be used to filter by region.

Indels are reported in the same manner as a VCF.

The output is 1-based.

```bash
perbase simple-depth ./test/test.bam
```

Example output

```text
chr     pos     a       c       t       g       n       ins     del     depth
chr2    1       1       0       0       0       0       0       0       1
chr2    2       1       0       0       0       0       1       0       1
chr2    3       1       0       0       0       0       0       0       1
chr2    4       1       0       0       0       0       0       0       1
chr2    5       2       0       0       0       0       0       0       2
chr2    6       2       0       0       0       0       0       1       2
chr2    7       1       0       0       0       0       0       0       1
chr2    8       1       0       0       0       0       0       0       1
chr2    9       1       0       0       0       0       0       0       1
chr2    10      2       0       0       0       0       0       0       2
chr2    11      2       0       0       0       0       0       0       2
chr2    12      3       0       0       0       0       0       0       3
chr2    13      2       0       0       0       0       0       1       2
```

The output can be compressed and indexed as follows:

```bash
perbase simple-depth ./test/test.bam | bgzip > output.tsv.gz
tabix -S 1 -s 1 -b 2 -e 2 ./output.tsv.gz
# Query all positions overlapping region
tabix output.tsv.gz chr1:5-10
```

## TODOs

- [ ] Paralleleize pileup walk
- [ ] Add more metrics to match `bam-readcount` as a `comprehensive` tool
- [ ] Add a strictly depth calculation a la `mosdepth`
- [ ] Explore other BAM/CRAM parsers that have a tidier API
- [ ] Add bgzip output / auto tabix indexing support
- [ ] Support SAM flags for pre-filtering
- [ ] Support limiting inputs with an interval_list file

## Tests

### Generate Expected

```bash
bam-readcount test/test.bam 2>/dev/null \
    | awk -F'\t' \
        'BEGIN{OFS=FS} \
        { \
            for (i=5; i <= NF; i++) { \
                split($i, a, ":"); \
                pos[a[1]]=a[2] \
            } \
            asorti(pos, pos_sorted); \
            if (NR==1) { \
                printf("%s\t%s", "chr", "pos"); \
                for (k in pos_sorted) { \
                    printf("\t%s", pos_sorted[k]) \
                } \
                printf("\ttotal\n");\
            } \
            printf "%s\t%s", $1, $2; \
            for (k in pos_sorted) { \
                printf "\t%s",  pos[pos_sorted[k]] \
            } \
            printf "\t%s\n", $4; \
    }' | cut -f1,2,4-9 > test/expected.tsv
```

If you update any of the SAM records in a test, you will need to rerun ^.
