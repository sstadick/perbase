# Notes on testing

* Check output is sorted for single chromosme: `awk 'NR != 1' output.tsv | sort -n --check`
* Check against sambamba (mate unaware, no filtering):

  ```bash
  sambamba depth base -t (nproc) -F ""  indexed.bam | awk -F '\t' 'BEGIN{OFS=FS} NR == 1 {print $0} NR != 1 {$2 += 1; $3 -= $9; print $0}' | cut -f1-9 > /tmp/sambamba_out.tsv

  ./target/release/perbase simple-depth indexed.bam | cut -f1-7,10-11 > /tmp/perbase_out.tsv

  vimdiff /tmp/perbase_out.tsv /tmp/sambamba_out.tsv
  ```

Add the `-m` flag to each to check with mate algorithms.