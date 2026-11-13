# Notes

## Data

Data is `ERR324/ERR3240157/HG00157.final.cram` from 1KG project. https://www.internationalgenome.org/data-portal/sample/HG00157

```
mkdir data
mkdir ref
# To fetch the run
pixi run fetch-data
```

## Benchmark

```
mkdir outputs
pixi run bench
# View results in `outputs`
```

### Incremental Rust/Goose pipeline comparison

The raw measurements, commands, binary hashes, output hashes, and native-build
source hashes are in [rust-pipeline-benchmarks.json](rust-pipeline-benchmarks.json).
These are local chromosome-20 experiments, **not the paper's whole-genome results**.
The input is HG00157 chromosome 20 (64,444,167 bp; 15,887,841 mapped records),
retaining the original 3,366-reference BAM header. Empty references are still
scheduled. Each run uses `base-depth -m -t 32 -c 1000000` and writes uncompressed
TSV to the same local filesystem. Hashing is outside timing.

Earlier stages, two runs each; time and peak RSS are means:

| Stage | Rust compiler / target | Elapsed | Peak RSS |
| --- | --- | ---: | ---: |
| Original Rust v1.4.1 | 1.97.1 / default | 51.85 s | 10.39 GiB |
| Optimized Rust mate groups | 1.97.1 / default | 34.05 s | 10.65 GiB |
| Optimized Rust + Seqair | 1.98.1 / default | 27.32 s | 15.68 GiB |
| Goose v0.1.3 | native CPU | 5.82 s | 5.20 GiB |

First matched-native comparison, two sequential runs each with order reversed
on the second round. Both Rust binaries use Rust 1.98.1, the release profile
(`opt-level=3`, fat LTO, one codegen unit), and `-C target-cpu=native`.
Goose v0.1.3 uses `-O3 -march=native`.

| Stage | Elapsed runs | Mean elapsed | Mean peak RSS | Mean CPU time (user + system) |
| --- | --- | ---: | ---: | ---: |
| Optimized Rust + Seqair (`8415130`) | 31.60, 31.87 s | 31.74 s | 15.57 GiB | 193.78 s |
| Rust + Goose application pipeline | 13.96, 14.32 s | 14.15 s | 5.94 GiB | 101.25 s |
| Goose v0.1.3 | 5.82, 5.96 s | 5.90 s | 5.20 GiB | 74.34 s |

All six outputs were 3,627,671,504 bytes and had SHA-256
`b4e1e41f38b69f5ae99c2a9b5fa652773dde07df41f4467572a6a2bc8de93ed7`,
also matching every earlier stage. Compare changes within a cohort: the old
27.32-second measurement is not the native-build control.

The new Rust path reuses worker readers and allocations, groups mates with
reusable index arrays, formats TSV in counting workers, and transfers bounded
64 KiB blocks to an ordered writer. Region boundaries and counting rules are
unchanged. This controls the application architecture, **not every library
implementation**: Rust Seqair 0.10.0 retains auxiliary tags and builds primary
mate links that Goose disables. The remaining runtime gap has not been isolated
by profiling; it is not a measurement of language overhead alone.

To repeat, build the baseline from `8415130` in a separate worktree and the
pipeline checkpoint from `cedbe16`, saving each executable under a distinct name:

```sh
RUSTFLAGS='-C target-cpu=native' cargo +1.98.1 build --release --features seqair --locked
python3 paper/bench_pipeline.py /path/to/paper-chr20-full.bam \
  --baseline /path/to/baseline-native --pipeline /path/to/pipeline-native \
  --goose /path/to/goose-perbase-v013 --output /path/to/new-results-directory
```

The runner makes two passes, verifies identical output, and retains logs and
JSON while removing the large TSVs. Do not run compilation or other benchmarks
concurrently. It uses the existing OS cache, not a cold-cache protocol.

### Seqair fork checkpoint

The next stage, **optimized rust + goose optimizations + seqair fork**, uses
`sstadick/seqair` branch `perbase-lean-pileup`, commit
`471a38e735c32bc9e5a1065c45bb362ede735e89`. It starts at upstream
`seqair-v0.10.0`, not the fork's stale main branch. The fork allows callers to
omit stored auxiliary tags and skip automatic primary-mate linking; defaults
remain unchanged. Perbase enables both opt-outs because its counter does not
read tags and already resolves same-name records itself.

Fresh comparison of native builds, two sequential runs per binary, reversed
order on the second round, with the same workload and compiler flags above:

| Stage | Elapsed runs | Mean elapsed | Mean peak RSS | Mean CPU time (user + system) |
| --- | --- | ---: | ---: | ---: |
| Optimized Rust + Seqair (`8415130`) | 31.69, 32.16 s | 31.93 s | 15.68 GiB | 193.83 s |
| + Goose application pipeline (`cedbe16`) | 14.19, 14.28 s | 14.24 s | 5.94 GiB | 101.01 s |
| + Seqair fork | 13.90, 14.04 s | 13.97 s | 5.05 GiB | 97.74 s |
| Goose v0.1.3 | 5.90, 6.09 s | 6.00 s | 5.20 GiB | 73.73 s |

All eight outputs have the same length and SHA-256 listed above. The fork
reduces mean peak RSS by 15.1%; elapsed time is only 1.9% lower in this small
sample. Goose remains 2.33× faster. These two library differences account for
some memory overhead but little of the runtime gap; attributing the remainder
requires profiling. This is not a whole-genome paper rerun.

Validation: 134 Perbase tests including its doctest; 137 targeted Seqair tests
and 51 doctests (13 intentionally ignored); strict Clippy on the changed Seqair
library and test targets. Another 155 CLI cases matched both the preceding
Seqair checkpoint and HTSlib, covering all nine mate strategies, BAI/CSI,
six CRAM encodings, filters, depth caps, references, overlapping regions,
empty inputs, tiny output queues, and BGZF output.

A separate `CG`-tag CIGAR probe confirmed an existing limitation: Seqair's
record-store path processes the placeholder CIGAR rather than expanding the
tag. Both Rust Seqair checkpoints agree but differ from HTSlib on that input.
The fork changes retention only, not decoding, and does not fix this limitation.

Both repositories are local checkpoints, not published releases. For this
experiment, Cargo uses the sibling `../seqair` checkout, including its
`seqair-types` crate. It must be at the fork commit above; Cargo's path dependency
does not enforce that revision. Once the fork is published, replace these paths
with a pinned Git dependency before distributing this Perbase branch.

Build this stage with the native command above, then add
`--fork /path/to/fork-native` to `bench_pipeline.py`. Raw measurements and source
hashes are under `fork_comparison` in `rust-pipeline-benchmarks.json`.

## Generate Plots

After running benchmarks, you can generate publication-quality plots for the paper:

```
pixi run plot-benchmarks
# or
pixi run generate-plots
```

This will create benchmark comparison plots in the `outputs` directory:
- `benchmark_comparison.png` - Performance comparison plot (PNG format)
- `benchmark_comparison.pdf` - Performance comparison plot (PDF format)

The plots show the runtime performance comparison between perbase and sambamba for both standard and mate-fix modes.


## Submission

https://joss.readthedocs.io/en/latest/submitting.html
- generate metadata
- double check formatting
- fill out submission
