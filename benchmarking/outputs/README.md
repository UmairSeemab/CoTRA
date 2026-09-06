# CoTRA computational performance benchmarking outputs

This directory contains the machine-readable outputs underlying the
computational performance benchmarking reported for CoTRA.

## Structure

```text
Benchmarking_outputs/
├── README.md
├── environment/
├── bulk/
│   ├── synthetic/
│   │   └── raw_measured/
│   └── real/
│       └── raw_measured/
├── scRNA/
│   ├── synthetic/
│   │   └── raw_measured/
│   └── real/
│       └── raw_measured/
├── step_level/
├── figure_source_data/
└── preliminary_viability_runs/
```

## Supplementary table mapping

- S1 Table: computational environment
- S2 Table: synthetic bulk RNA-seq benchmark
- S3 Table: real retinal bulk RNA-seq benchmark
- S4 Table: synthetic scRNA-seq benchmark
- S5 Table: real retinal scRNA-seq benchmark
- S6 Table: unified step-level timings

## Measurement design

Each dataset size and analysis configuration was evaluated in five independent
measured R processes. A separate preliminary viability run was performed before
the measured replicates and excluded from all reported summary statistics.

Analytical workflow runtime was measured around the relevant analysis steps.
Peak resident memory was obtained from GNU `/usr/bin/time -v` and therefore
represents the maximum memory footprint of the complete benchmark R process.

The corresponding scripts are distributed under `Benchmarking_code/`.
