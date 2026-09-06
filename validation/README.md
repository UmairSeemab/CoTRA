# CoTRA implementation-validation code

This directory contains the scripts used to reproduce the implementation-
concordance analyses reported in S7 Table and S8 Table.

## Structure

```text
Implementation_validation_code/
├── README.md
├── bulk/
│   ├── README.md
│   ├── validate_bulk_concordance.R
│   └── run_bulk_validation.sh
└── scRNA/
    ├── README.md
    ├── validate_scRNA_Seurat_concordance.R
    └── run_scRNA_Seurat_validation.sh
```

## Relationship to validation outputs

The matching machine-readable result files are stored separately under:

```text
Implementation_validation_outputs/
├── bulk/
└── scRNA/
```

Recommended repository layout:

```text
validation/
├── Implementation_validation_code/
└── Implementation_validation_outputs/
```

## Scope

The validation does not attempt to revalidate the statistical methodology of
DESeq2, edgeR, or Seurat. It tests whether the CoTRA graphical implementation
reproduces the corresponding direct scripted implementations when the same
inputs and analytical parameters are used.

## Reported results

Bulk:
- exact DESeq2 and edgeR tested/significant gene sets
- Jaccard = 1.000
- directional concordance = 100%
- log2FC Pearson/Spearman = 1.000

scRNA:
- identical 15,106-cell clustering into 16 clusters
- ARI = 1.000
- NMI = 1.000
- 2,000/2,000 HVGs matched
- PCA PC1-PC7 correlation = 1.000
- 38,491 marker gene-cluster pairs matched
- marker Jaccard and log2FC correlations = 1.000
