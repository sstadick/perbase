# CHANGELOG

## 0.8.1

- Bugfix `--keep-zeros` and `--no-merge` had an incorrect conditional resulting in all positions being filtered out.

## 0.8.0

- Added `--keep-zeros` option to prevent truncating regions that have 0 depth at the start / end of the interval
- Added `--skip-merging-intervals` to prevent the merging of overlapping intervals if a BED/VCF/BCF file is used to specify regions. This will produce duplicate entries
- [Improved](https://github.com/sstadick/perbase/issues/42) error message for invalid BED file
- [Improved](https://github.com/sstadick/perbase/issues/41) handling of 3-column BED files (they now work)
- Moved to `gzp` for multithreaded BGZF output writing when specified. Also added compression-threads and compression-level options
- Updated all deps

## 0.7.4

- Add `--bed-format` flag to `only-depth` to output a bed-like file with the depth in the 5th column and no headers.

## 0.7.3

- Update all dependencies, specifically htslib to hopefully sort out
  openssl errors in the conda build process.

## 0.7.2

- Use published version of noodles so we can publish on crates.io

## 0.7.1

- Fix regression in release profile

## 0.7.0

- Added `--min-base-quality/-Q` flag to the `base-depth` tool. When present, this flag will cause a base quality check of reach base. If the quality is less thn the specified minium quality the depth will be counted as an `N` instead of an `[A, C, T, G]`. If this flag is not set the behavior is not changed.
- Added the `--bgzip/-Z` flag to all subcommands which bgzips the output. Note that this does not index the bgzipped file at this time. This can still be done with `tabix -S 1 -s 1 -b 2 -e 2 ./output.tsv.gz`
- Changed log level of "Batch Processing ...", "Processing ...", and "Processing region" logging statements to `trace`. Set `RUST_LOG=trace` environment variable to restore previous logging verbosity.
- Added "Processing TID ..." log statement at `info` level (on by default).
- All `perbase` commands now gracefully handle broken pipes and exit 0.
