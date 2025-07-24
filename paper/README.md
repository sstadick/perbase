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

