# scRNA-seq implementation-validation outputs

- `Table_S8_scRNA_Implementation_Concordance.csv`: S8 summary table.
- `PCA_concordance_PC1-PC7.csv`: PCA score concordance.
- `cluster_contingency_matrix.csv`: cell-level cluster cross-tabulation.
- `cluster_label_mapping.csv`: direct Seurat-to-CoTRA cluster mapping.

Input test datasets:
- `GSM7474906_Wild_Type_non_treated_feature_bc_matrix.h5`
- `GSM7474907_Rd10_Female_vehicle_feature_bc_matrix.h5`

The validation compared CoTRA with a direct standalone Seurat workflow using
matched inputs and parameters.
