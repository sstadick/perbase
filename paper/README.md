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

Fresh matched-native comparison, two sequential runs each with order reversed
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
pipeline version from this tree, saving each executable under a distinct name:

```sh
RUSTFLAGS='-C target-cpu=native' cargo +1.98.1 build --release --features seqair --locked
python3 paper/bench_pipeline.py /path/to/paper-chr20-full.bam \
  --baseline /path/to/baseline-native --pipeline /path/to/pipeline-native \
  --goose /path/to/goose-perbase-v013 --output /path/to/new-results-directory
```

The runner makes two passes, verifies identical output, and retains logs and
JSON while removing the large TSVs. Do not run compilation or other benchmarks
concurrently. It uses the existing OS cache, not a cold-cache protocol.

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
