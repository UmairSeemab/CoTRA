# scRNA-seq implementation-validation code

This folder reproduces the S8 CoTRA versus direct Seurat concordance analysis.

## Required biological test inputs

- `GSM7474906_Wild_Type_non_treated_feature_bc_matrix.h5`
- `GSM7474907_Rd10_Female_vehicle_feature_bc_matrix.h5`

## Required CoTRA validation exports

- `CoTRA_cluster_assignments_20260906_161659.csv`
- `CoTRA_all_cluster_markers_20260906_161825.csv`
- `CoTRA_scRNA_Clustering_session_20260906_161609.rds`

The large RDS is required only to recover the CoTRA Seurat object for direct HVG
and PCA concordance checks. It does not need to be distributed as an S8 output
file if the validation can be rerun from CoTRA.

## Direct Seurat settings

- LogNormalize, scale factor 10,000
- 2,000 HVGs using `vst`
- ScaleData on HVGs
- 50 PCs calculated
- PCs 1-7 used downstream
- FindNeighbors `k.param = 20`
- Louvain algorithm 1
- resolution 0.5
- random seed 1234
- FindAllMarkers with Wilcoxon
- `only.pos = TRUE`
- `min.pct = 0.10`
- `logfc.threshold = 0.25`
- `max.cells.per.ident = 500`

## Metrics

- input gene/cell equality
- cluster number
- adjusted Rand index
- normalized mutual information
- HVG intersection/Jaccard
- PCA concordance for PCs 1-7
- marker gene-cluster pair Jaccard
- significant marker pair Jaccard
- marker log2FC Pearson/Spearman correlation

## Run

```bash
chmod +x run_scRNA_Seurat_validation.sh
./run_scRNA_Seurat_validation.sh
```

Edit file paths near the top of the runner if needed.
