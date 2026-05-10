library(shiny)
library(bs4Dash)
library(shinyjs)
library(shinyalert)

# -------------------------
#  Load Bulk RNA modules
# -------------------------

source("modules/bulk/upload.R")
source("modules/bulk/mds.R")
source("modules/bulk/groups.R")
source("modules/bulk/qc_visualization.R")
source("modules/bulk/de_analysis.R")
source("modules/bulk/annotation.R")
source("modules/bulk/bulk_enrichment.R")
source("modules/bulk/report.R")

# -------------------------
#  Load scRNA modules
# -------------------------

source("modules/sc/import.R")
# source("modules/sc/demultiplexing.R")
source("modules/sc/qc.R")
source("modules/sc/gene_selection.R")
source("modules/sc/pca.R")
source("modules/sc/tsne_umap.R")
source("modules/sc/clustering.R")
# source("modules/sc/doublets.R")
# source("modules/sc/batch_integration.R")
# source("modules/sc/subset_recluster.R")
# source("modules/sc/annotation.R")
source("modules/sc/markers.R")
# source("modules/sc/featureplots.R")
# source("modules/sc/feature_viewer.R")
source("modules/sc/annotation.R")
source("modules/sc/differential_abundance.R")
source("modules/sc/trajectory.R")
source("modules/sc/pathway_activity.R")
source("modules/sc/cell_communication.R")
#source("modules/sc/velocity.R")
# source("modules/sc/report.R")
# source("modules/sc/export.R")

# ==========================================================
#                        UI PAGE
# ==========================================================

ui <- bs4DashPage(
  title = "CoTRA: Comprehensive Toolkit for RNA-seq Analysis",
  header = bs4DashNavbar(title = tags$span("")),
  
  sidebar = bs4DashSidebar(
    div(style="padding:5px; text-align:left;",
        tags$img(src="cotra_logo.png", style="height:80px; width:auto;")),
    
    skin = "light",
    status = "primary",
    title = "CoTRA Menu",
    width = 300,
    
    bs4SidebarMenu(
      
      bs4SidebarMenuSubItem(
        "Home",
        tabName = "home",
        icon = icon("home")
      ),
      
      # =====================================================
      #                  BULK RNA MENU
      # =====================================================
      bs4SidebarMenuItem(
        text = "Bulk RNA-seq",
        icon = icon("dna"),
        startExpanded = TRUE,
        
        bs4SidebarMenuSubItem("Upload Data",         tabName = "bulk_upload",   icon = icon("file-upload")),
        bs4SidebarMenuSubItem("Create Groups",       tabName = "bulk_groups",    icon = icon("users")),
        bs4SidebarMenuSubItem("QC / Visualization",  tabName = "bulk_qc",        icon = icon("chart-area")),
        bs4SidebarMenuSubItem("DE Analysis",         tabName = "bulk_de",        icon = icon("chart-line")),
        bs4SidebarMenuSubItem("Gene Annotation", tabName = "bulk_annot", icon = icon("book")),
        bs4SidebarMenuSubItem("Enrichment Analysis", tabName = "bulk_enrich",    icon = icon("flask")),
        bs4SidebarMenuSubItem("Report",              tabName = "bulk_report",    icon = icon("file"))
        
        
      ),

      
      # =====================================================
      #                SINGLE CELL RNA MENU
      # =====================================================
      bs4SidebarMenuItem(
        text          = "Single cell RNA",
        icon          = icon("microscope"),
        tabName       = "sc_main",
        startExpanded = TRUE,
        
        bs4SidebarMenuSubItem("Data import",  tabName = "sc_import", icon = icon("upload")),
        #bs4SidebarMenuSubItem("Demultiplexing", tabName = "sc_demux", icon = icon("vials")),
        bs4SidebarMenuSubItem("Quality Control", tabName = "sc_qc", icon = icon("heartbeat")),
        bs4SidebarMenuSubItem("Gene Selection", tabName = "sc_vargenes", icon = icon("dna")),
        bs4SidebarMenuSubItem("Variable Genes + PCA", tabName = "sc_vargene_pca", icon = icon("chart-line")),
        bs4SidebarMenuSubItem("t-SNE + UMAP ", tabName = "sc_tsne_umap", icon = icon("project-diagram")),
        bs4SidebarMenuSubItem("Clustering", tabName = "sc_cluster", icon = icon("vials")),
        #bs4SidebarMenuSubItem("Subset + Recluster", tabName = "sc_subset", icon = icon("crop")),
        bs4SidebarMenuSubItem("Markers", tabName = "sc_markers", icon = icon("bolt")),
        bs4SidebarMenuSubItem("Annotation", tabName = "sc_annot", icon = icon("tags")),
        bs4SidebarMenuSubItem("Differential Abundance", tabName = "sc_da", icon = icon("chart-bar")),
        #bs4SidebarMenuSubItem("Feature Plots", tabName = "sc_featureplots", icon = icon("palette")),
        #bs4SidebarMenuSubItem("Marker Viewer", tabName = "sc_feature", icon = icon("eye")),
        bs4SidebarMenuSubItem("Trajectory", tabName = "sc_trajectory", icon = icon("space-shuttle")),
        bs4SidebarMenuSubItem("Pathway Activity", tabName = "sc_pathway", icon = icon("fire")),
        bs4SidebarMenuSubItem("Cell Communication", tabName = "sc_comm", icon = icon("satellite-dish")),
        #bs4SidebarMenuSubItem("RNA Velocity", tabName = "sc_velocity", icon = icon("bolt")),
        bs4SidebarMenuSubItem("Report", tabName = "sc_report", icon = icon("file-alt")),
        bs4SidebarMenuSubItem("Export", tabName = "sc_export", icon = icon("download"))
      )
    )
  ),
  
  # ==========================================================
  #                        BODY
  # ==========================================================
  
  body = bs4DashBody(
    useShinyjs(),
    tags$link(rel="stylesheet", type="text/css", href="css/dark.css"),
    tags$link(rel="stylesheet", type="text/css", href="custom.css"),
    
    bs4TabItems(
      
      # HOME
      bs4TabItem(
        tabName = "home",
        fluidRow(
          bs4Card(
            title = "Welcome to CoTRA",
            status = "primary",
            width = 12,
            solidHeader = TRUE,
            h3("CoTRA: Comprehensive Toolkit for RNA-seq Analysis"),
            p("CoTRA provides complete workflows for bulk and single-cell RNA-seq analysis."),
            hr(),
            h4("Workflow Overview"),
            tags$img(src="cotra_flowchart.svg", width="100%",
                     style="border:1px solid #ccc; padding:8px; border-radius:6px;")
          )
        )
      ),

      # ---------------------- BULK RNA ----------------------
      bs4TabItem("bulk_upload",  mod_bulk_upload_ui("bulk_upload")),
      bs4TabItem("bulk_mds",     mod_bulk_mds_ui("bulk_mds")),
      bs4TabItem("bulk_groups",  mod_bulk_groups_ui("bulk_groups")),
      bs4TabItem("bulk_qc",      mod_bulk_qc_ui("bulk_qc")),
      bs4TabItem("bulk_de",      mod_bulk_de_ui("bulk_de")),
      bs4TabItem("bulk_annot",   mod_bulk_annotation_ui("bulk_annot")),
      bs4TabItem("bulk_enrich", mod_bulk_enrichment_ui("bulk_enrich")),
      bs4TabItem("bulk_report",  mod_bulk_report_ui("bulk_report")),
      
      # ---------------------- scRNA -------------------------
      # ---------------------- scRNA -------------------------
      bs4TabItem("sc_import", mod_sc_import_ui("sc_import")),
      bs4TabItem("sc_qc", mod_sc_qc_ui("sc_qc")),
      
      # Gene selection OK
      bs4TabItem("sc_vargenes", mod_sc_gene_selection_ui("sc_gene_select")),
      
      # PCA FIXED
      bs4TabItem("sc_vargene_pca", mod_sc_pca_ui("sc_pca")),
      
      # tSNE + UMAP FIXED
      bs4TabItem("sc_tsne_umap", mod_sc_tsne_umap_ui("sc_tsne_umap")),
      
      # Clustering
      bs4TabItem( tabName = "sc_cluster", mod_sc_clustering_ui("sc_cluster")),
      
      # Gene Markers
      bs4TabItem("sc_markers", mod_sc_markers_ui("sc_markers")),
      
      # Annotation
      bs4TabItem("sc_annot", mod_sc_annotation_ui("sc_annot")),
      
      # Differential Abundance
      bs4TabItem("sc_da", mod_sc_differential_abundance_ui("sc_da")),
      
      # Trajectory
      bs4TabItem("sc_trajectory", mod_sc_trajectory_ui("sc_trajectory")),
      
      #Pathway Activity
      bs4TabItem("sc_pathway", mod_sc_pathway_ui("sc_pathway")),
      
      # Cell-Cell Communication
      bs4TabItem("sc_comm", mod_sc_comm_ui("sc_comm"))
      
      # RNA velocity
      #bs4TabItem("sc_velocity", mod_sc_velocity_ui("sc_velocity"))
      # All others will be activated later one by one
    )
  ),
  
  controlbar = NULL,
  
  footer = bs4DashFooter("CoTRA: Comprehensive Toolkit for RNA-seq Analysis")
)

shinyUI(ui)
