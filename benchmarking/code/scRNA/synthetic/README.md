# Synthetic scRNA-seq benchmark code

This folder contains the final v3 scripts used for the synthetic scRNA-seq
computational performance benchmark.

## Design
- 20,000 genes
- 2,500, 5,000, 10,000, 25,000, or 50,000 cells
- 12 defined truth populations
- 150 designated marker genes per population
- 2,000 HVGs
- 50 PCs calculated
- PCs 1-30 used downstream
- Louvain clustering, resolution 0.5
- marker timing evaluated separately using the known truth labels
- five measured runs per dataset size

## Files
- `generate_synthetic_scRNA_v2.R`
- `benchmark_scRNA_synthetic_v3.R`
- `run_scRNA_synthetic_v3_linux.sh`
- `summarize_scRNA_synthetic_v3.R`

Only the final v3 benchmarking workflow is included here. Earlier exploratory
versions are intentionally excluded.

The corresponding machine-readable outputs are in:
`Benchmarking_outputs/scRNA/synthetic/`.
