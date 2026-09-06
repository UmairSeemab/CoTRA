#!/usr/bin/env bash
set -euo pipefail

# Adjust only the four CoTRA-export paths below if your browser saved them
# somewhere other than ~/Downloads.

WT_H5="/home/umair/Downloads/CoTRA/data/scRNA/GSM7474906_Wild_Type_non_treated_feature_bc_matrix.h5"
RD10_H5="/home/umair/Downloads/CoTRA/data/scRNA/GSM7474907_Rd10_Female_vehicle_feature_bc_matrix.h5"

COTRA_CLUSTERS="${HOME}/Downloads/CoTRA_cluster_assignments_20260906_161659.csv"
COTRA_MARKERS="${HOME}/Downloads/CoTRA_all_cluster_markers_20260906_161825.csv"
COTRA_SESSION="${HOME}/Downloads/CoTRA_scRNA_Clustering_session_20260906_161609.rds"

OUTDIR="/home/umair/Downloads/CoTRA/benchmark/validation_scRNA_Seurat"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

Rscript --vanilla validate_scRNA_Seurat_concordance.R \
  --wt-h5 "$WT_H5" \
  --rd10-h5 "$RD10_H5" \
  --cotra-clusters "$COTRA_CLUSTERS" \
  --cotra-markers "$COTRA_MARKERS" \
  --cotra-session "$COTRA_SESSION" \
  --outdir "$OUTDIR"
