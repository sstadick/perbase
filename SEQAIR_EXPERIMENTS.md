# Seqair experiments

This file tracks perbase-driven experiments against branches in `sstadick/seqair`.

Baseline for these experiments is perbase `feat/seqair-pileup-aggregation`, which pins the seqair aggregation prototype from `sstadick/seqair@6d9251c` and uses `PileupEngine::pileup_with` for non-mate `base-depth --seqair-pileup`.

Benchmark data set unless otherwise noted:

- BAM: `paper/data/HG00157.chr1_10mb.bam`
- BED: `paper/data/hg00157_chr1_10mb.bed`
- build: `cargo build --release --features seqair-pileup`
- htslib build env on this machine:
  - `SDKROOT="$(xcrun --show-sdk-path)"`
  - `BINDGEN_EXTRA_CLANG_ARGS="-isysroot $(xcrun --show-sdk-path)"`

## Current baseline: pileup accumulator

Branches / commits:

- seqair: `feat/pileup-column-aggregation` @ `6d9251c`
- perbase: `feat/seqair-pileup-aggregation` @ `6d55e59`

Effect:

- Non-mate `base-depth --seqair-pileup` uses `PileupEngine::pileup_with` and avoids materializing a public `PileupColumn` for perbase counting.
- Mate-aware `base-depth -m --seqair-pileup` still uses materialized columns.
- `only-depth --seqair` is unchanged.

Latest benchmark from PR #108:

| Mode | htslib | seqair | Output |
|---|---:|---:|---|
| `base-depth` | `5.147 ± 0.461 s` | `5.104 ± 0.166 s` | exact parity, SHA `150f6165...` |
| `base-depth -m` | `51.045 s` | `33.402 s` | known 12 sparse default `-F 0` mate-order diffs |
| `base-depth -m -F 2304` | `51.072 s` | `33.770 s` | exact parity, SHA `393a5787...` |
| `only-depth` | `1.203 ± 0.157 s` | `1.263 ± 0.125 s` | exact parity, SHA `0d98d1ab...` |
| `only-depth -x` | `859.3 ± 15.6 ms` | `1.200 ± 0.053 s` | exact parity, SHA `876b7692...` |

Takeaway: the accumulator helped non-mate `base-depth`, but only modestly. The downstream alignment pass / materialized public column was not the dominant cost.

## Experiment: reusable pileup alignment buffer

Idea: pass a reusable `Vec<PileupAlignment>` into the materialized pileup path to avoid one alignment-vector allocation per emitted column.

Branches / commits:

- seqair: `exp/pileups-reusable-buffer` @ `c1491a3`
- perbase: `exp/seqair-pileups-reusable-buffer` @ `7ed1308`

Seqair API sketch:

- `PileupColumn<'store, U, A = Vec<PileupAlignment>>` became generic over alignment storage.
- `PileupEngine::pileups_into(&mut Vec<PileupAlignment>) -> Option<PileupColumn<'_, U, &[PileupAlignment]>>` fills caller-owned storage.

Perbase integration:

- Mate-aware `base-depth -m --seqair-pileup` uses `engine.pileups_into(&mut alignments)`.
- Non-mate `base-depth --seqair-pileup` is still on the custom accumulator path.

Validation:

- seqair: `cargo check -p seqair`
- seqair: targeted `pileups_into_reuses_caller_alignment_buffer` test passed.
- perbase: `cargo check --features seqair-pileup`
- perbase: empty-SEQ seqair/htslib regression passed.
- perbase: release build passed.

Benchmark:

| Mode | htslib | seqair reusable buffer | Output |
|---|---:|---:|---|
| `base-depth -m` | `50.336 s` | `32.917 s` | known 12 sparse default `-F 0` diffs |
| `base-depth -m -F 2304` | `49.407 s` | `32.995 s` | exact parity, SHA `393a5787...` |

Comparison to baseline:

- Baseline `base-depth -m --seqair-pileup`: `33.402 s`
- Reusable-buffer `base-depth -m --seqair-pileup`: `32.917 s`

Takeaway: small improvement at best, within single-run noise. Worth keeping as an API option if the API shape is acceptable, but not a major speed lever for this dataset.

## Experiment: raw BAM record fetch visitor for only-depth

Idea: `only-depth` does not need a pileup-ready `RecordStore`; it only needs raw filter fields and CIGAR intervals. Add a seqair BAM API that visits raw fetched records after index/tid/overlap checks without decoding SEQ/QUAL or copying slabs.

Branches / commits:

- seqair: `exp/raw-record-fetch` @ `dcf26f3`
- perbase: `exp/seqair-raw-record-fetch` @ `641ab9e`

Seqair API sketch:

- `IndexedBamReader::fetch_raw_records(tid, start, end, &mut visitor)`
- Visitor receives `FilterRawFields<'_>` with parsed header fields and borrowed raw qname/CIGAR/SEQ/QUAL/AUX slices.
- No `RecordStore`, no decoded bases, no qual/aux/name slab copies.

Perbase integration:

- Non-mate `only-depth --seqair` uses raw CIGAR bytes directly.
- Non-mate `only-depth -x --seqair` uses raw record start/end fields directly.
- Mate-aware only-depth still falls back to the existing store path for now.

Validation:

- seqair: `cargo check -p seqair`
- seqair: `raw_record_fetch_matches_store_fetch_positions` passed.
- perbase: `cargo check --features seqair-pileup`
- perbase: all six `seqair_only_depth*` parity tests passed.
- perbase: release build passed.

Benchmark:

| Mode | htslib | seqair raw-record | Output |
|---|---:|---:|---|
| `only-depth` | `1.165 ± 0.224 s` | `831.1 ± 29.4 ms` | exact parity, SHA `0d98d1ab...` |
| `only-depth -x` | `856.9 ± 17.6 ms` | `854.2 ± 54.5 ms` | exact parity, SHA `876b7692...` |

Comparison to baseline seqair:

- Baseline `only-depth --seqair`: `1.263 ± 0.125 s`
- Raw-record `only-depth --seqair`: `0.831 ± 0.029 s`
- Baseline `only-depth -x --seqair`: `1.200 ± 0.053 s`
- Raw-record `only-depth -x --seqair`: `0.854 ± 0.055 s`

Takeaway: this is the clearest win so far. It turns `only-depth` from slower than htslib into faster for normal mode and roughly tied for fast mode, with exact output parity.

## Experiment: active-record CIGAR cursor

Idea: the pileup hot loop calls `CigarMapping::pos_info_at(pos)` for every active record and every column. For complex CIGARs this may scan or binary-search compact ops repeatedly. Since pileup positions are monotonic, keep a per-active-record cursor into the compact CIGAR ops.

Branches / commits:

- seqair: `exp/pileup-cigar-cursor` @ `e7c5a85`
- perbase: `exp/seqair-pileup-cigar-cursor` @ `832ca17`

Seqair implementation:

- Added `CigarMapping::pos_info_at_cursor(pos, &mut cursor)`.
- `ActiveRecord` stores a `cigar_cursor`.
- Pileup uses cursor lookup while iterating positions.

Validation:

- seqair: `cargo check -p seqair`
- seqair: cursor-vs-random-access unit test passed.
- seqair: full `cargo test -p seqair --test pileup` passed.
- perbase: `cargo check --features seqair-pileup`
- perbase: release build passed.

Benchmark:

| Mode | htslib | seqair CIGAR cursor | Output |
|---|---:|---:|---|
| `base-depth` | `5.655 ± 0.572 s` | `5.139 ± 0.038 s` | exact parity, SHA `150f6165...` |
| `base-depth -m` | `49.607 s` | `32.967 s` | known 12 sparse default `-F 0` diffs |

Comparison to baseline seqair:

- Baseline non-mate accumulator: `5.104 ± 0.166 s`
- CIGAR cursor non-mate: `5.139 ± 0.038 s`
- Baseline mate-aware: `33.402 s`
- CIGAR cursor mate-aware: `32.967 s`

Takeaway: no clear win on this short-read HG00157 slice. Likely most reads are simple linear CIGARs, where seqair was already using the fast path. This may still matter on RNA/long-read/complex-CIGAR data, but it is not a priority for this perbase benchmark.

## Experiment: cold-slab decode profile

Idea: keep using `RecordStore`, but skip cold slabs that downstream pileup consumers do not need. For base-depth, aux tags are unused; non-mate base-depth also does not need qnames.

Branches / commits:

- seqair: `exp/decode-cold-slab-profile` @ `95990c4`
- perbase: `exp/seqair-decode-profile` @ `ea8ae9d`

Seqair API sketch:

- `DecodeProfile::{FULL, PILEUP_NO_AUX, PILEUP_NO_NAMES_OR_AUX}`
- `RecordStore::push_raw_with_profile(...)`
- `IndexedBamReader::fetch_into_customized_with_profile(...)`

Important limitation:

- This prototype only skips qname and aux slabs. It does not skip bases/qualities because the current `SlimRecord` layout assumes `seq_len` indexes into dense bases/qual slabs; safely projecting those out needs an additional presence bit/sentinel design.

Perbase integration:

- BAM `base-depth --seqair-pileup` uses `PILEUP_NO_NAMES_OR_AUX`.
- BAM `base-depth -m --seqair-pileup` uses `PILEUP_NO_AUX` because mate fixing still needs qnames.
- SAM/CRAM/ref-backed `Readers` paths use the existing full decode path.

Validation:

- seqair: `cargo check -p seqair`
- seqair: `decode_profile_can_skip_qname_and_aux_slabs` passed.
- perbase: `cargo check --features seqair-pileup`
- perbase: empty-SEQ seqair/htslib regression passed.
- perbase: release build passed.

Benchmark:

| Mode | htslib | seqair decode profile | Output |
|---|---:|---:|---|
| `base-depth` | `5.839 ± 0.534 s` | `5.254 ± 0.138 s` | exact parity, SHA `150f6165...` |
| `base-depth -m` | `50.255 s` | `32.814 s` | known 12 sparse default `-F 0` diffs |

Comparison to baseline seqair:

- Baseline non-mate accumulator: `5.104 ± 0.166 s`
- Decode-profile non-mate: `5.254 ± 0.138 s`
- Baseline mate-aware: `33.402 s`
- Decode-profile mate-aware: `32.814 s`

Takeaway: skipping aux/qname alone is not a clear non-mate win and may be noise. There may be a small mate-aware improvement from skipping aux, but the larger projection win is the raw-record visitor for `only-depth`. A full SEQ/QUAL projection would need a larger `RecordStore` layout change.

## Experiment still open: mate-aware aggregation

Status: not implemented yet.

Reason to save for last:

- perbase mate fixing groups by qname and then applies order-sensitive tie-breaking.
- Current default `base-depth -m -F 0` still has 12 sparse backend-order-sensitive diffs; `-F 2304` has exact parity.
- A custom mate-aware accumulator is possible, but it likely needs either:
  - an owned per-column qname/alignment wrapper in perbase, or
  - a seqair accumulator API that can lend `AlignmentView`s during finish without materializing a public `PileupColumn`.

Expected next approach:

1. First make perbase mate-fix tie-breaking deterministic across backends.
2. Then prototype a mate-aware accumulator that groups by qname inside the accumulator and emits `PileupPosition` directly.
3. Benchmark against the materialized-column mate-aware path.

## Current ranking of wins

1. **Raw BAM record fetch visitor for only-depth**: clear win, exact parity.
2. **Reusable pileup alignment buffer**: small/noisy mate-aware improvement; API may still be useful.
3. **Cold-slab decode profile**: no clear non-mate win; maybe small mate-aware win; less important than raw projection.
4. **CIGAR cursor**: no clear win on this short-read dataset; maybe revisit for complex-CIGAR data.
5. **Mate-aware aggregation**: not yet tested; should wait until tie-breaking is deterministic or explicitly modeled.
