# 0.6.4

- Added `--min-base-quality/-Q` flag to the `base-depth` tool. When present, this flag will cause a base quality check of reach base. If the quality is less thn the specified minium quality the depth will be counted as an `N` instead of an `[A, C, T, G]`. If this flag is not set the behavior is not changed.
- Added the `--bgzip/-Z` flag to all subcommands which bgzips the output. Note that this does not index the bgzipped file at this time. This can still be done with `tabix -S 1 -s 1 -b 2 -e 2 ./output.tsv.gz`
- Changed log level of "Batch Processing ...", "Processing ...", and "Processing region" logging statements to `trace`. Set `RUST_LOG=trace` environment variable to restore previous logging verbosity.
- Added "Processing TID ..." log statement at `info` level (on by default).
- All `perbase` commands now gracefully handle broken pipes and exit 0.
