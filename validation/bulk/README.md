# Bulk implementation-validation code

This folder reproduces the S7 CoTRA bulk implementation-concordance analysis.

## Required inputs

1. `raw_gene_counts.tsv`
2. CoTRA full DESeq2 export:
   `BulkDE_AllResults_rd10_vs_WT_DESeq2.csv`
3. CoTRA full edgeR export:
   `BulkDE_AllResults_rd10_vs_WT_edgeR.csv`

## What the script does

`validate_bulk_concordance.R`:

- reads the same retinal raw-count matrix used by CoTRA;
- identifies WT and rd10 samples;
- reruns DESeq2 directly using the current CoTRA filtering and contrast;
- reruns edgeR directly using the current CoTRA filtering and QL workflow;
- applies the same significance thresholds;
- documents the CoTRA DESeq2 `NA -> 1` p/padj post-processing;
- compares tested genes, significant genes, log2FC, significance direction,
  Jaccard overlap, and exact numeric concordance;
- repeats the scripted comparison five times;
- recreates the machine-readable S7 summary and run-level files.

## Run

```bash
chmod +x run_bulk_validation.sh
./run_bulk_validation.sh
```

Edit paths near the top of the shell script if required.

## Main outputs

- `Table_S7_Bulk_Implementation_Concordance.csv`
- `run_level_concordance.csv`
- `direct_DESeq2_results.csv`
- `direct_edgeR_results.csv`
- `validation_sessionInfo.txt`
