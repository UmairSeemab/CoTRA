#!/usr/bin/env bash
set -euo pipefail

# CoTRA bulk implementation-validation runner.
#
# Edit the three file paths below if your repository or exported CoTRA results
# are stored elsewhere.

RAW_COUNTS="/home/umair/Downloads/CoTRA/data/bulkRNA/raw_gene_counts.tsv"
COTRA_DESEQ2="${HOME}/Downloads/BulkDE_AllResults_rd10_vs_WT_DESeq2.csv"
COTRA_EDGER="${HOME}/Downloads/BulkDE_AllResults_rd10_vs_WT_edgeR.csv"

OUTDIR="/home/umair/Downloads/CoTRA/benchmark/validation_bulk_implementation"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

Rscript --vanilla validate_bulk_concordance.R \
  --input "$RAW_COUNTS" \
  --cotra-deseq2 "$COTRA_DESEQ2" \
  --cotra-edger "$COTRA_EDGER" \
  --reference WT \
  --comparison rd10 \
  --reference-regex '^WT' \
  --comparison-regex '^rd10' \
  --padj 0.05 \
  --lfc 1 \
  --repeats 5 \
  --outdir "$OUTDIR"
