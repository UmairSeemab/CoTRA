# CoTRA computational performance benchmarking code

This directory contains the final scripts used to generate, execute, summarize,
and visualize the computational performance benchmarks reported for CoTRA.

## Structure

```text
Benchmarking_code/
├── README.md
├── bulk/
│   ├── synthetic/
│   └── real/
├── scRNA/
│   ├── synthetic/
│   └── real/
└── figure_generation/
```

## Benchmark categories

1. Synthetic bulk RNA-seq scaling benchmark
2. Real retinal bulk RNA-seq benchmark
3. Synthetic scRNA-seq scaling benchmark
4. Real retinal scRNA-seq benchmark

All measured benchmark conditions were executed in five independent R processes.
A separate preliminary viability run was performed before the measured
replicates and was excluded from reported summary statistics.

Thread counts for OMP, OpenBLAS, MKL, VECLIB, and NUMEXPR were restricted to
one where applicable.

The matching machine-readable results are provided in `Benchmarking_outputs/`.
