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

#### Published checkpoints and reproduction

The experiment branches are
[Perbase `optimize-base-depth-mates`](https://github.com/sstadick/perbase/tree/optimize-base-depth-mates)
and [Seqair `perbase-lean-pileup`](https://github.com/sstadick/seqair/tree/perbase-lean-pileup).
These are experimental branches, not new package releases. The related
[Perbase paper](https://joss.theoj.org/papers/10.21105/joss.09774) provides context;
the chromosome-20 measurements below are separate experiments.

| Checkpoint | Perbase commit | Rust Seqair source |
| --- | --- | --- |
| Optimized mate groups | [b1e7185](https://github.com/sstadick/perbase/commit/b1e7185) | HTSlib backend |
| Optimized Rust + Seqair | [8415130](https://github.com/sstadick/perbase/commit/8415130) | Upstream 0.10.0 |
| + Goose application pipeline | [cedbe16](https://github.com/sstadick/perbase/commit/cedbe16) | Upstream 0.10.0 |
| + Seqair fork | [70cb3a8](https://github.com/sstadick/perbase/commit/70cb3a8) | [471a38e](https://github.com/sstadick/seqair/commit/471a38e735c32bc9e5a1065c45bb362ede735e89) |
| + Remaining optimizations | [8689fe5](https://github.com/sstadick/perbase/commit/8689fe5) | [6da4be2](https://github.com/sstadick/seqair/commit/6da4be2d59f87c0e076e7cfc8b2f561fe6909396) |
| + Full-block buffering | [076471d](https://github.com/sstadick/perbase/commit/076471d) | Same fork revision |
| Shared HTSlib/Seqair pipeline | [1e68316](https://github.com/sstadick/perbase/commit/1e68316) | Same fork revision |

The final branch pins both Rust Seqair crates to the published `6da4be2` Git
revision. Earlier fork checkpoints use `../seqair` path dependencies: reproduce
them with a sibling checkout at the matching revision in the table. Checkpoint
source hashes match the recorded measured source; the final dependency pin and
documentation do not change the counting implementation.

Build the final branch with Rust 1.98.1 and Clang 18.1.3:

```sh
CC=clang CFLAGS='-O3 -march=native' CARGO_BUILD_JOBS=8 \
  RUSTFLAGS='-C target-cpu=native' \
  cargo +1.98.1 build --release --features seqair --locked
```

For the latest backend comparison, use that **same executable** with and
without `--seqair`. Set `BAM` to the chromosome-20 input, `RUST` to the executable
above, and `GOOSE` to Goose Perbase v0.1.3 built with `-O3 -march=native`.
Run on an otherwise idle machine; keep output on the same filesystem:

```sh
OUTPUT=$(mktemp -d ./comparison.XXXXXX)
for variant in htslib seqair goose goose seqair htslib; do
  case "$variant" in
    htslib) set -- "$RUST" ;;
    seqair) set -- "$RUST" ;;
    goose) set -- "$GOOSE" ;;
  esac
  set -- "$@" base-depth "$BAM" -m -t 32 -c 1000000 -o "$OUTPUT/comparison.tsv"
  if [ "$variant" = seqair ]; then set -- "$@" --seqair; fi
  printf '\n%s\n' "$variant"
  /usr/bin/time -f '%e seconds; %U user; %S system; %M KiB peak RSS' "$@" || exit
  sha256sum "$OUTPUT/comparison.tsv"
  wc -c < "$OUTPUT/comparison.tsv"
  rm "$OUTPUT/comparison.tsv"
done
rmdir "$OUTPUT"
```

Each output must have the length and SHA-256 below. Hashing is outside timing;
this is a warm-cache-oriented protocol, not a cold-cache test. Do not use a
whole-genome BAM and expect chromosome-20 timings.

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

At this historical checkpoint, Cargo uses the sibling `../seqair` checkout,
including its `seqair-types` crate. It must be at the fork commit above; Cargo's
path dependency does not enforce that revision. The final branch replaces
these paths with the pinned Git dependency described above.

Build this stage with the native command above, then add
`--fork /path/to/fork-native` to `bench_pipeline.py`. Raw measurements and source
hashes are under `fork_comparison` in `rust-pipeline-benchmarks.json`.

### Remaining Rust optimizations

The next experiment keeps the same counting rules and workload. Workers reuse
one output row rather than copying its reference name at every position, reuse
one spare TSV block per output slot, and group records using Seqair's cached
name hashes. Read-filter and depth-admission metadata now live in one worker
array rather than also being copied into the record store.

The Seqair fork adds preparation for a sequential pileup walk, without building
the unused arbitrary-window query index. Such queries remain correct through a
linear scan. It also offers lazy reference-match checking: the existing method
still computes it eagerly, while counting callers can avoid that work without
changing any alignment fields or the answer returned when explicitly queried.

Depth admission uses cached alignment ends and skips its heap when the total
number of records cannot reach the cap. **Memory tradeoff:** admission now runs
after loading a region, like Goose. A low `--max-depth` therefore no longer
bounds the retained read payload; very deep regions can use more memory even
though emitted counts remain capped by the same HTSlib admission rules.

The first attempt retained all returned output buffers. It took 14.04 s versus
13.95 s for the preceding fork checkpoint and increased mean peak RSS from
5.04 to 5.98 GiB. That was not a speedup. Limiting recycling to one spare block
per slot reduced RSS to 5.41 GiB, but still overallocated short partial blocks.
The final version also allocates those tails at their actual size, recycles
only full 64 KiB blocks, and never blocks the writer when returning them.
The two intermediate attempts are preserved under
`remaining_optimization_attempts` in `rust-pipeline-benchmarks.json`.

Final fresh comparison, with the same native builds, workload, two runs per
binary, reversed second-round order, and output validation used above:

| Stage | Elapsed runs | Mean elapsed | Mean peak RSS | Mean CPU time (user + system) |
| --- | --- | ---: | ---: | ---: |
| Optimized Rust + Seqair (`8415130`) | 31.82, 32.02 s | 31.92 s | 15.65 GiB | 193.71 s |
| + Goose application pipeline (`cedbe16`) | 14.19, 14.11 s | 14.15 s | 5.94 GiB | 101.13 s |
| + Seqair fork (`70cb3a8`) | 13.93, 13.98 s | 13.96 s | 5.05 GiB | 96.12 s |
| + Remaining optimizations (working trees) | 14.17, 14.14 s | 14.15 s | 4.94 GiB | 95.97 s |
| Goose v0.1.3 | 5.96, 5.85 s | 5.91 s | 5.19 GiB | 73.37 s |

All ten outputs were byte-identical, with the same 3,627,671,504-byte size and
SHA-256 reported above. The final pass uses about 2.1% less peak memory than
the previous fork checkpoint but takes 1.4% longer in this two-run sample:
**there is no demonstrated runtime improvement**. Goose remains 2.40× faster.
This controls the high-level pipeline, not every library implementation, and
does not isolate language overhead. This is HG00157 chromosome 20, not a
whole-genome reproduction of the paper. Final raw measurements, compiler flags,
binary hashes, and source hashes are under `remaining_optimizations_comparison`
in `rust-pipeline-benchmarks.json`.

Validation: 134 feature-enabled Perbase tests and 131 default-feature tests;
1,015 Seqair library tests (one ignored), 101 integration tests, and 51 doctests
(13 ignored). The final output-buffer change also passed the streaming
regression test and 231 three-way CLI comparisons against the preceding fork
and HTSlib, plus eight failed-sink/closed-pipe cases. Those comparisons cover all
nine mate strategies, depth caps, BAM/CRAM, filters, references, uncovered
positions, duplicate regions, output queues, and compression. The existing
`CG`-tag limitation above remains outside this equivalence claim.

Hardware-counter profiling was unavailable (`perf_event_paranoid=4`). Timing
does not establish which inner loop explains the remaining Rust/Goose gap.
These measurements were taken from working trees, subsequently preserved in
Perbase `8689fe5` and Seqair `6da4be2` as linked above. The recorded source hashes
identify the measured files. These experimental checkpoints are not releases.

To reproduce with both matching source trees, build with the native command
above and add `--extra 'rust + remaining optimizations=/path/to/new-binary'`
alongside the existing baseline, pipeline, fork, and Goose arguments.

### Correcting the output-buffer comparison

The earlier high-level pipeline comparison missed a concrete mismatch: Goose
limits each output mailbox by **bytes**, while Rust limits its queue by
**messages**. Rust was sending a full 64 KiB block followed by a tiny row tail
on most flushes. Both consumed one queue slot. With this input's 26-character
maximum reference name, the intended 90 MB budget therefore held roughly
45 MB of output. Workers blocked before completing their regions, leaving the
ordered writer waiting for serial counting tails.

Rust now retains partial blocks until more rows arrive or the region ends.
The worker count, genomic jobs, nominal budget, counting rules, and output order
are unchanged. Only the final data message in a region may be partial.

Two-run, reversed-order comparison with otherwise unchanged binaries:

| Stage | Elapsed runs | Mean elapsed | Mean peak RSS |
| --- | --- | ---: | ---: |
| Rust before buffer correction | 14.03, 14.03 s | 14.03 s | 4.94 GiB |
| Rust with full-block buffering | 6.92, 7.18 s | 7.05 s | 5.34 GiB |
| Goose v0.1.3 | 5.99, 5.80 s | 5.90 s | 5.19 GiB |

Separate instrumented runs support the diagnosis: summed worker time inside
channel sends fell from 219.34 to 0.01 seconds, and writer receive time fell
from 8.71 to 3.36 seconds. Worker durations overlap; they are not process elapsed
time. Send/receive timings include channel overhead, not only blocked time.
The instrumented runs are diagnostics, not rows in the benchmark table.

We also matched native C compilation to Goose's Clang 18.1.3 with
`CC=clang CFLAGS='-O3 -march=native'`. Rust itself remains rustc 1.98.1 with
`-C target-cpu=native`, fat LTO, and one codegen unit. This changes compilation
of native dependencies, not Rust's compiler. The earlier libdeflate objects
identified GCC 13.3.0; the new build log and object metadata confirm Clang.
Both programs use native libdeflate 1.26.

| Matched native-C comparison | Rust elapsed runs | Goose elapsed runs | Rust mean | Goose mean | Rust / Goose |
| --- | --- | --- | ---: | ---: | ---: |
| 30 counting workers (`-t 32`) | 6.93, 6.85 s | 5.88, 6.04 s | 6.89 s | 5.96 s | 1.16× |
| One counting worker (`-t 3`) | 63.31, 63.53 s | 43.99, 44.84 s | 63.42 s | 44.42 s | 1.43× |

In the same native-C cohort, the full-block Rust build with its previous C
compiler took 7.17 and 7.11 seconds. Two runs are insufficient to attribute a
small difference to compiler settings; the large, demonstrated change is the
output-buffer correction. The single-worker comparison shows that residual
work differences remain even without competition among counting workers.

The libraries are closer than a claim of fundamentally different pileup
algorithms would suggest: both use forward CIGAR cursors, stable entry order,
gap skipping, and SIMD-capable sequence expansion. Bounds checks exist in both.
But they are **not implementation-identical**:

- Rust's measured `PileupAlignment` is 44 bytes; Goose's generated-C entry stride
  is 27 bytes. Rust still materializes read summaries and library-mate metadata
  unused by Perbase. Some cached fields also avoid later record lookups, so
  compacting entries is a hypothesis to test, not a guaranteed speedup.
- Rust uses a decoded-block cache and a growing compressed-input window;
  Goose uses current-block reuse and fixed read-ahead buffers.
- Goose validates auxiliary fields and resolves long-CIGAR `CG` escapes even
  when tag retention is disabled. Rust's record-store `CG` limitation remains.
  Fewer validation checks are not an established explanation for Goose's speed.

All 16 outputs across the buffering, native-C, and single-worker experiments
have the existing byte count and SHA-256. All 134 Perbase feature-enabled tests
passed in both debug and release; the release tests used the matched native-C
build. No new Seqair library changes were made in this investigation.

Raw results, source and binary hashes, native build settings, and diagnostic
timings are under `buffer_parity` in `rust-pipeline-benchmarks.json`. Profiling
logs and instrumentation are retained under `/srv/goose-perbase/rust-parity/`.
The correction is preserved in Perbase `076471d`, linked above. The remaining
1.16× parallel gap is an implementation result on HG00157 chromosome 20,
**not a measured language
penalty or a whole-genome paper reproduction**.

### HTSlib versus Seqair with the same application pipeline

Previously, omitting `--seqair` also selected the older batching and CSV-output
pipeline. That would not isolate the backend comparison. In a build with the
`seqair` Cargo feature, both backends now share worker scheduling, genomic jobs,
TSV formatting, full-block output buffering, and ordered writing. HTSlib uses
its existing worker-local indexed readers and optimized mate counter; Seqair
uses its region store and indexed mate groups. Reading, pileup generation,
and backend-specific counting therefore remain part of the measured difference.
Builds without the Cargo feature still use the legacy HTSlib pipeline.

Fresh comparison on the same HG00157 chr20 BAM, `base-depth -m -t 32 -c 1000000`,
uncompressed TSV, two sequential runs per variant with reversed second-round
order. Both Rust variants use the **same executable**, built with rustc 1.98.1,
native CPU targeting, fat LTO, and native C dependencies built by Clang 18.1.3
with `-O3 -march=native`. Only the `--seqair` argument differs.

| Backend / implementation | Elapsed runs | Mean elapsed | Mean peak RSS | Mean CPU time (user + system) |
| --- | --- | ---: | ---: | ---: |
| Optimized Rust + HTSlib | 17.76, 17.51 s | 17.64 s | 2.51 GiB | 259.15 s |
| Optimized Rust + Seqair fork | 7.04, 7.27 s | 7.16 s | 5.38 GiB | 105.29 s |
| Goose v0.1.3 | 5.90, 5.85 s | 5.88 s | 5.20 GiB | 73.85 s |

The Seqair-backed Rust path is **2.46× faster** here and uses **2.14× the peak
memory** of the HTSlib-backed path. This is an end-to-end counting comparison,
not a decoder-only benchmark or a claim about every possible HTSlib design.
All six outputs match the existing byte count and SHA-256.

Validation: 134 feature-enabled release tests, 131 default-feature tests, and
134 three-way CLI cases comparing the previous HTSlib path, streaming HTSlib,
and streaming Seqair. Cases cover all nine mate strategies, depth caps, BAI,
CSI, CRAM, filters, reference bases, zeros, duplicate intervals, small queues,
and compressed output. The streaming test now exercises both backends,
including output failure. The known Seqair `CG` limitation remains unchanged.

Raw runs and source/build provenance are recorded under
`shared_pipeline_backend_comparison` in `rust-pipeline-benchmarks.json`.
Logs and the measured executable are under `/srv/goose-perbase/rust-backends/`.
The application changes are preserved in Perbase `1e68316`, linked above.
The final branch additionally pins the fork dependency for reproducible builds;
these experiments do not constitute a package release.

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
