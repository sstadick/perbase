# CHANGELOG

## 1.3.0

- [Fix](https://github.com/sstadick/perbase/pull/95): typo in `count_mate_resolutions` from @camlloyd
- [chore](https://github.com/sstadick/perbase/issues/98): added CONTRIBUTING.md per @mbhall88 suggestion
- [fix](https://github.com/sstadick/perbase/issues/97): handle missing ref files gracefully
- [fix](https://github.com/sstadick/perbase/commit/76e38e6f8dc52dc3da9f87f1fa17391224ec35f8): many small typos, improved help messages, fixed conflicting short opts

## 1.2.0

- Fix: Handle BAM records with empty SEQ fields (`*` in SAM format) in `base-depth`. Previously this would cause an out-of-bounds error. Now these reads still count toward depth, and their bases are counted as `N`. ([#92](https://github.com/sstadick/perbase/pull/92) by @ghuls)
- Chore: Update deps and fixup lints

## 1.1.0

- Fix/Feat: For both MapQualBaseQualN and BaseQualMapQualN only resolve to N when the bases are ambiguous, otherwise return the consensus base.

## 1.0.0

- Feat: expand mate-fix resolution strategies
  - Adds new output columns for IUPAC bases, which are only used for some of the mate-fix resolution strats
  - By default the "original" mate fix strat will be used and remains backward compatible

## 0.10.3

- **Bugfix**: Fixed mate selection in `--mate-fix` to correctly prefer first mate over second mate when MAPQ scores are equal (issue #82)

## 0.8.3

- Updated dependencies - to fixed version of gzp

## 0.8.2

- Updated dependencies

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

## 0.10.0

- [fix](https://github.com/sstadick/perbase/pull/73) from @biermanr for [71](https://github.com/sstadick/perbase/issues/71) (--keep-zeros skip/dup loci)

## 0.10.1

- [fix](https://github.com/sstadick/perbase/issues/74) from @nkkarpov updates smartstring version wich fixes some UB in old smartstring.

## 0.10.2

- [fix](https://github.com/sstadick/perbase/pull/78) Update htslib version v0.39 by @davidecarlson
