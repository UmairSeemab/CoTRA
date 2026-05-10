library(shiny)
library(bs4Dash)
library(shinyjs)
library(shinyWidgets)
library(Seurat)


# -----------------------------
#  Bulk module loading
# -----------------------------
source("modules/bulk/upload.R")
source("modules/bulk/mds.R")
source("modules/bulk/groups.R")
source("modules/bulk/qc_visualization.R")
source("modules/bulk/de_analysis.R")
source("modules/bulk/annotation.R")
source("modules/bulk/bulk_enrichment.R")
source("modules/bulk/report.R")

# -----------------------------
#  scRNA module loading
# -----------------------------
source("modules/sc/import.R")
source("modules/sc/qc.R")
source("modules/sc/gene_selection.R")
source("modules/sc/pca.R")
source("modules/sc/tsne_umap.R")
source("modules/sc/clustering.R")
source("modules/sc/markers.R")
source("modules/sc/annotation.R")
source("modules/sc/differential_abundance.R")
source("modules/sc/trajectory.R")
source("modules/sc/pathway_activity.R")
source("modules/sc/cell_communication.R")
#source("modules/sc/velocity.R")


options(shiny.maxRequestSize = 300 * 1024^2)

server <- function(input, output, session) {
  
  useShinyjs()
  
  # ==========================================================
  #                BULK RNA-SEQ PIPELINE
  # ==========================================================
  
  bulk_data   <- mod_bulk_upload_server("bulk_upload")
  mod_bulk_mds_server("bulk_mds", bulk_data)
  bulk_groups <- mod_bulk_groups_server("bulk_groups", bulk_data)
  mod_bulk_qc_server("bulk_qc", bulk_data, bulk_groups)
  de_results  <- mod_bulk_de_server("bulk_de", bulk_data, bulk_groups)
  bulk_annot  <- mod_bulk_annotation_server("bulk_annot", de_results)
  bulk_enrich <- mod_bulk_enrich_server("bulk_enrich", deg_data = de_results, wgcna_output = reactive({ bulk_wgcna }))
  bulk_report <- mod_bulk_report_server("bulk_report", bulk_reactive = bulk_data, groups_reactive = bulk_groups, de_results = de_results, annot_reactive = bulk_annot, enrich_state = bulk_enrich)

    
  # ==========================================================
  #             scRNA-SEQ PIPELINE (ONLY FIRST 3 MODULES)
  # ==========================================================
  
  # 1. Import
  sc_in <- mod_sc_import_server("sc_import")
  
  # 2. QC (takes Seurat from import)
  sc_qc <- mod_sc_qc_server(
    "sc_qc",
    seurat_r = sc_in$seurat
  )
  
  # 3. HVG / Gene selection
  sc_hvg <- mod_sc_gene_selection_server(
    "sc_gene_select",
    seurat_r = sc_qc$seurat,
    sc_state = cotra_state
  )
  
  # 4. PCA  (IMPORTANT: this returns a reactive, not a list)
  sc_pca <- mod_sc_pca_server(
    "sc_pca",
    seurat_r = sc_hvg$seurat,
    sc_state = cotra_state
  )
  
  # 5. t-SNE & UMAP  (IMPORTANT: pass reactive directly, NOT $seurat)
  sc_tsne_umap <- mod_sc_tsne_umap_server(
    "sc_tsne_umap",
    seurat_r = sc_pca,
    sc_state = cotra_state
  )
  
  # 6. Clustering  
  sc_cluster <- mod_sc_clustering_server(
    "sc_cluster",
    seurat_r = sc_tsne_umap$seurat, # ✔ now correct
    sc_state = cotra_state
  )
  
  # 7. Markers
  sc_markers <- mod_sc_markers_server(
    "sc_markers",
    seurat_r = sc_cluster$seurat,   # ✔ MUST COME FROM CLUSTERING MODULE
    sc_state = cotra_state
  )
  
  # Annotation
  sc_annot <- mod_sc_annotation_server(
    "sc_annot",
    seurat_r = sc_markers$seurat,
    sc_state = cotra_state
  )
  
  # 8. Differential Abundance
  sc_da <- mod_sc_differential_abundance_server(
    "sc_da",
    seurat_r = sc_annot$seurat,
    sc_state = cotra_state
  )
  
  # 8. Trajectory, independent
  sc_trajectory <- mod_sc_trajectory_server(
    "sc_trajectory",
    seurat_r = sc_markers$seurat,
    sc_state = cotra_state
  )
  
  # 9. Pathway Activity, independent from trajectory
  sc_pathway <- mod_sc_pathway_server(
    "sc_pathway",
    seurat_r = sc_markers$seurat,
    sc_state = cotra_state
  )
  
  # 10. Cell-Cell Communication
  sc_comm <- mod_sc_comm_server(
    "sc_comm",
    seurat_r = sc_markers$seurat,
    sc_state = cotra_state
  )
  
  # 10. RNA Velocity
  #sc_velocity <- mod_sc_velocity_server(
   # "sc_velocity",
   # seurat_r = sc_markers$seurat,
  #  sc_state = cotra_state
  #)
}

shinyServer(server)
