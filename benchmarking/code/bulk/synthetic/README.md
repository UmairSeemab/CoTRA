# Synthetic bulk RNA-seq benchmark code

This folder contains the final scripts used for the synthetic bulk RNA-seq
computational performance benchmark.

## Design
- 20,000 genes
- 6, 12, 24, 48, or 96 samples
- balanced two-group designs
- DESeq2 and edgeR
- five independent measured R processes per condition
- separate preliminary viability run excluded from summary statistics
- single-thread execution where applicable

## Files
- `benchmark_bulk.R`: executes one benchmark condition/run
- `run_bulk_synthetic_linux.sh`: Linux runner for the complete benchmark
- `summarize_benchmarks.R`: summarizes measured runs

The corresponding machine-readable outputs are in:
`Benchmarking_outputs/bulk/synthetic/`.
