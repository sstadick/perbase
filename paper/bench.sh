#!/usr/bin/env bash

# Create output directory if it doesn't exist
mkdir -p ./data/outputs

# Base-depth comparison

bench_base_depth () {
	hyperfine \
		--runs 3 \
		--export-csv ./data/outputs/bench_results_base_depth.csv \
		--command-name "perbase base-depth 30X" \
		'../target/release/perbase base-depth -o ./data/outputs/perbase_base_depth_30x.tsv ./data/HG00157.final.bam' \
		--command-name "sambamba depth base 30X" \
		'sambamba depth base -t $(nproc) -F ""  -o ./data/outputs/sambamba_depth_base_30x.tsv ./data/HG00157.final.bam'
		#--command-name "bam-readcount 30X" \
		#'bam-readcount ./data/HG00157.final.bam > ./data/outputs/bam-readcount_30x.tsv'
}

bench_base_depth_with_mate_fix () {
	hyperfine \
		--runs 3 \
		--export-csv ./data/outputs/bench_results_base_depth_mate_fix.csv \
		--command-name "perbase base-depth mate-fix 30X" \
		'../target/release/perbase base-depth -m -o ./data/outputs/perbase_base_depth_mate_fix_30x.tsv ./data/HG00157.final.bam' \
		--command-name "sambamba depth base mate-fix 30X" \
		'sambamba depth base -t $(nproc)  -F "" -m -o ./data/outputs/sambamba_depth_base_mate_fix_30x.tsv ./data/HG00157.final.bam'
}

bench_base_depth_with_mate_fix_new () {
	hyperfine \
		--runs 1 \
		--export-csv ./data/outputs/bench_results_base_depth_mate_fix_new_v_old.csv \
		--command-name "perbase base-depth mate-fix 30X Original" \
		'../../temp/perbase/target/release/perbase base-depth -m --mate-resolution-strategy Original -o ./data/outputs/perbase_base_depth_mate_fix_30x.tsv ./data/HG00157.final.bam' \
		--command-name "perbase base-depth mate-fix 30X N" \
		'../../temp/perbase/target/release/perbase base-depth -m --mate-resolution-strategy N -o ./data/outputs/perbase_base_depth_mate_fix_30x.tsv ./data/HG00157.final.bam' \
		--command-name "perbase base-depth mate-fix 30X IUPAC" \
		'../../temp/perbase/target/release/perbase base-depth -m --mate-resolution-strategy IUPAC -o ./data/outputs/perbase_base_depth_mate_fix_30x.tsv ./data/HG00157.final.bam' \
		--command-name "perbase base-depth mate-fix 30X MapQualBaseQualN" \
		'../../temp/perbase/target/release/perbase base-depth -m --mate-resolution-strategy MapQualBaseQualN -o ./data/outputs/perbase_base_depth_mate_fix_30x.tsv ./data/HG00157.final.bam' \
		--command-name "perbase base-depth mate-fix 30X MapQualBaseQualIUPAC" \
		'../../temp/perbase/target/release/perbase base-depth -m --mate-resolution-strategy MapQualBaseQualIUPAC -o ./data/outputs/perbase_base_depth_mate_fix_30x.tsv ./data/HG00157.final.bam' \
		--command-name "perbase base-depth mate-fix 30X MapQualBaseQualFirstInPair" \
		'../../temp/perbase/target/release/perbase base-depth -m --mate-resolution-strategy MapQualBaseQualFirstInPair -o ./data/outputs/perbase_base_depth_mate_fix_30x.tsv ./data/HG00157.final.bam' \
		--command-name "perbase base-depth mate-fix 30X BaseQualMapQualN" \
		'../../temp/perbase/target/release/perbase base-depth -m --mate-resolution-strategy BaseQualMapQualN -o ./data/outputs/perbase_base_depth_mate_fix_30x.tsv ./data/HG00157.final.bam' \
		--command-name "perbase base-depth mate-fix 30X BaseQualMapQualIUPAC" \
		'../../temp/perbase/target/release/perbase base-depth -m --mate-resolution-strategy BaseQualMapQualIUPAC -o ./data/outputs/perbase_base_depth_mate_fix_30x.tsv ./data/HG00157.final.bam' \
		--command-name "perbase base-depth mate-fix 30X BaseQualMapQualFirstInPair" \
		'../../temp/perbase/target/release/perbase base-depth -m --mate-resolution-strategy BaseQualMapQualFirstInPair -o ./data/outputs/perbase_base_depth_mate_fix_30x.tsv ./data/HG00157.final.bam' \
		--command-name "perbase base-depth mate-fix 30X OLD" \
		'../target/release/perbase base-depth -m -o ./data/outputs/perbase_base_depth_mate_fix_30x.tsv ./data/HG00157.final.bam' \
		--command-name "sambamba depth base mate-fix 30X" \
		'sambamba depth base -t $(nproc)  -F "" -m -o ./data/outputs/sambamba_depth_base_mate_fix_30x.tsv ./data/HG00157.final.bam'
}


bench_base_depth_with_mate_fix_new
