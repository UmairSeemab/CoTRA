# Real retinal scRNA-seq benchmark code

This folder contains the final scripts used for the real retinal scRNA-seq
computational performance benchmark.

## Dataset reconstruction
The WT non-treated and vehicle-treated rd10 samples from GSE234797 were used.
Processed cell identifiers were mapped back to the original H5 matrices to
recover raw integer counts for the retained cells.

Final complete benchmark dataset:
- 18,533 genes
- 10,430 cells
- 4,478 WT
- 5,952 rd10

Nested condition-stratified subsets:
- 2,500
- 5,000
- 7,500
- 10,000
- 10,430 cells

## Files
- `download_GSE234797_processed.sh`
- `prepare_real_scRNA_exact_raw_v2.R`
- `benchmark_scRNA_real.R`
- `run_scRNA_real_exact_linux.sh`
- `summarize_scRNA_real.R`

The corresponding machine-readable outputs are in:
`Benchmarking_outputs/scRNA/real/`.
