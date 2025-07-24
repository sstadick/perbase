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


# Only-depth comparison

bench_only_depth () {
	hyperfine \
		--export-csv ./data/outputs/bench_results_only_depth.csv \
		--command-name "perbase only-depth 30X" \
		'../target/release/perbase only-depth --bgzip -o ./data/outputs/perbase_only_depth_30x.tsv.gz ./data/HG00157.final.bam' \
		--command-name "mosdepth 30X" \
		'mosdepth -t $(nproc) ./data/outputs/mosdepth_30x ./data/HG00157.final.bam'
}

#bench_base_depth_with_mate_fix
bench_only_depth
