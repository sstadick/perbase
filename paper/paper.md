---
title: 'perbase: A highly parallelized per-base sequencing metrics tool'
tags:
  - rust
  - bioinformatics
  - genomics
  - sequencing
  - depth analysis
  - parallel computing
authors:
  - name: Seth Stadick
    orcid: 0009-0002-0915-9459
    equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)
affiliations:
 - name: Bio-Rad Laboratories, United States
   index: 1
date: 14 July 2025
bibliography: paper.bib

---

# Summary

`perbase` is a highly parallelized command-line tool for calculating per-base sequencing metrics from BAM/CRAM files. Built in Rust, it leverages concurrent processing to deliver orders of magnitude faster performance than existing tools while maintaining accuracy. The tool provides unbiased base counting by directly reporting observed bases without modification or filtering beyond user specifications. `perbase` includes utilities for depth calculation, nucleotide counting at each position, and specialized features like mate-overlap correction, making it suitable for applications ranging from variant calling to coverage analysis in whole genome sequencing projects.

# Statement of need

Per-base sequencing metrics are fundamental to genomic analyses, from variant calling to coverage assessment. Existing tools like `sambamba depth` [@Tarasov2015], `samtools depth` [@Li2009], `mosdepth` [@Pedersen2018], and `bam-readcount` [@Khanna2022] provide similar functionality but often lack the performance needed for modern high-throughput sequencing datasets. As sequencing depths increase and genome-scale analyses become routine, computational efficiency becomes critical.

`perbase` addresses this need through automatic parallelization over genomic regions. Unlike tools that process data sequentially, `perbase` divides the genome into chunks and processes them concurrently, scaling performance with available compute resources. This design enables processing of deep whole genome sequencing data in minutes rather than hours.

The tool provides accurate, unbiased metrics by counting all aligned bases including deletions toward depth, while correctly excluding reference skips. This differs from some tools: `bam-readcount` excludes deletions from depth calculations, while `sambamba` incorrectly includes reference skips. These distinctions matter for downstream analyses where accurate depth representation is critical.

# Implementation and Features

## Parallel Architecture

`perbase` implements a sophisticated parallel processing system through its `ParGranges` module. The genome is divided into configurable chunks (default 1Mb), which are then distributed across worker threads using Rust's Rayon library. This zero-copy parallelization automatically scales with available CPU cores while maintaining memory efficiency through bounded channels sized proportionally to thread count.

## Core Tools

**base-depth**: Calculates depth and nucleotide composition at every position, outputting counts for A, C, G, T, N bases, insertions, deletions, reference skips, and reads failing filters. An optional `--mate-fix` flag handles overlapping mate pairs by selecting the mate with highest mapping quality, with ties broken by choosing the first mate. This prevents double-counting in paired-end data while preserving the most reliable base calls.

**only-depth**: Provides rapid depth-only calculations with two modes. The default mode considers CIGAR operations (counting deletions but not reference skips), while `--fast-mode` uses only read start/stop positions for maximum speed. Adjacent positions with identical depth are automatically merged to reduce output size, inspired by the linear depth representation pioneered by `mosdepth` [@Pedersen2018].

**merge-adjacent**: A utility for post-processing that merges adjacent intervals with the same depth, useful for creating compact representations of coverage landscapes.

## Performance

Benchmarking on standard datasets demonstrates significant performance advantages. On a 30x whole genome sequencing dataset, `perbase` can process specific regions or generate genome-wide statistics faster than existing tools, with performance scaling nearly linearly with thread count. Memory usage remains modest even for high-coverage data, typically under 1GB for a 300x genome.

# Research Applications

`perbase` has been designed for diverse genomic applications. In variant calling pipelines, the per-nucleotide counts enable sophisticated filtering and quality assessment. For structural variant detection, accurate depth profiles help identify copy number changes and deletions. The tool's speed makes it practical for real-time quality control during sequencing runs.

The unbiased output facilitates method development and benchmarking. Researchers can access raw base counts without hidden filtering or transformations, enabling custom downstream analyses. The tool's ability to restrict analysis to specific regions via BED or VCF files supports targeted applications like exome or panel sequencing.

For benchmarking and validation, datasets like the Genome in a Bottle NA12878 reference [@Zook2019] or high-coverage 1000 Genomes Project data [@Byrska2022] provide ideal test cases with known variants across varying depths. These resources allow users to validate `perbase` metrics against established truth sets.

# Acknowledgements

We acknowledge the Rust programming language community and the authors of dependencies including rust-htslib, Rayon, and rust-lapper.

# References


