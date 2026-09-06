# CoTRA implementation-validation outputs

This directory contains machine-readable output files underlying the
reproducibility and implementation-concordance analyses reported for CoTRA.

## Directory structure

Implementation_validation_outputs/
├── bulk/
│   ├── Table_S7_Bulk_Implementation_Concordance.csv
│   └── run_level_concordance.csv
└── scRNA/
    ├── Table_S8_scRNA_Implementation_Concordance.csv
    ├── PCA_concordance_PC1-PC7.csv
    ├── cluster_contingency_matrix.csv
    └── cluster_label_mapping.csv

## Bulk RNA-seq validation

The bulk validation compared CoTRA outputs with direct scripted executions of
DESeq2 and edgeR using the same retinal count matrix (`raw_gene_counts.tsv`),
sample grouping, contrast, filtering rules, and significance thresholds.

Key results:
- DESeq2: 19,982 tested genes in both workflows
- DESeq2: 2,518 significant genes in both workflows
- edgeR: 18,922 tested genes in both workflows
- edgeR: 2,322 significant genes in both workflows
- Significant-gene Jaccard index = 1.000 for both methods
- Direction concordance = 100%
- log2FC Pearson and Spearman correlations = 1.000

For DESeq2, 127 genes had undefined raw and adjusted P values in the direct
package output. CoTRA converts these NA values to 1 before significance
filtering. This documented post-processing did not alter the significant-gene
set.

### Files

`bulk/Table_S7_Bulk_Implementation_Concordance.csv`
Publication-level summary corresponding to S7 Table.

`bulk/run_level_concordance.csv`
Run-level concordance statistics from the direct scripted validation runs.

## scRNA-seq validation

The scRNA-seq validation compared CoTRA with a direct standalone Seurat workflow
using the same WT and rd10 10x Genomics HDF5 files and matched analysis
parameters.

Key results:
- 32,285 genes
- 15,106 cells
- 7,217 WT cells
- 7,889 rd10 cells
- 16 clusters in both workflows
- Adjusted Rand Index = 1.000
- Normalized Mutual Information = 1.000
- 2,000/2,000 highly variable genes matched
- PCA PC1-PC7 correlations = 1.000
- 38,491 marker gene-cluster pairs matched
- Marker-set Jaccard index = 1.000
- Marker avg_log2FC Pearson and Spearman correlations = 1.000

### Files

`scRNA/Table_S8_scRNA_Implementation_Concordance.csv`
Publication-level summary corresponding to S8 Table.

`scRNA/PCA_concordance_PC1-PC7.csv`
Per-component PCA concordance for PCs 1-7.

`scRNA/cluster_contingency_matrix.csv`
Cell-level contingency matrix comparing CoTRA and direct Seurat cluster labels.

`scRNA/cluster_label_mapping.csv`
Mapping between direct Seurat cluster labels and CoTRA cluster labels based on
maximum cell overlap.

## Intended use

These files can be distributed through the CoTRA GitHub repository and/or an
archival repository such as Zenodo as machine-readable supporting data for the
implementation-validation results.

Large session RDS objects and full intermediate outputs are intentionally not
included because the compact files here are sufficient to verify the reported
headline concordance statistics.
