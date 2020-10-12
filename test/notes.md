# Notes on testing

* Check output is sorted for single chromosme: `awk 'NR != 1' output.tsv | sort -n --check`
* Check against sambamba (mate unaware, no filtering):

  ```bash
  sambamba depth base -t (nproc) -F ""  indexed.bam | awk -F '\t' 'BEGIN{OFS=FS} NR == 1 {print $0} NR != 1 {$2 += 1; $3 -= $9; print $0}' | cut -f1-9 > /tmp/sambamba_out.tsv

  ./target/release/perbase simple-depth indexed.bam | cut -f1-7,10-11 > /tmp/perbase_out.tsv

  vimdiff /tmp/perbase_out.tsv /tmp/sambamba_out.tsv
  ```

Add the `-m` flag to each to check with mate algorithms.

## Compare only-depth

Note, this has -m flag set

```bash
sambamba depth base -m -t (nproc) -F ""  ./test/sorted.bam | awk -F '\t' 'BEGIN{OFS=FS} NR == 1 {print $0} NR != 1 {$2 += 1; $3 -= $9; print $0}' | cut -f1-3 | ./target/release/perbase merge-adjacent > /tmp/sambamba_out.tsv

./target/release/perbase only-depth -m  ./test/sorted.bam | ./target/release/perbase merge-adjacent > /tmp/perbase_out.tsv

vimdiff /tmp/perbase_out.tsv /tmp/sambamba_out.tsv
```

There will be a difference around "chr3:65" where sambamba is removing the mate from the depth count but we keep it since the first in pair is ref_skip over that region

```bash
./target/release/perbase only-depth -x -z -F 1796 ./test/test.bam > /tmp/sambamba_out.tsv
mosdepth -x /tmp/tiny ./test/test.bam
vimdiff /tmp/sambamba_out.tsv (gzcat /tmp/tiny.per-base.bed.gz | psub)
```