# ===============================================================
# modules/sc/export.R
# CoTRA scRNA-seq Export Module
# - Export Seurat object (RDS)
# - Export markers table (CSV)
# - PDF summary report (UMAP + markers + QC)
# ===============================================================

library(Seurat)
library(ggplot2)
library(gridExtra)
library(cowplot)

mod_sc_export_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    bs4Card(
      title = "Export scRNA-seq Project",
      status = "primary",
      solidHeader = TRUE,
      width = 12,
      
      fluidRow(
        column(4,
               textInput(ns("export_name"), "Project Name", value = "scRNA_project")
        ),
        column(4,
               downloadButton(ns("download_rds"), "Download Seurat RDS",
                              class = "btn btn-success btn-block")
        ),
        column(4,
               downloadButton(ns("download_markers"), "Download Markers CSV",
                              class = "btn btn-info btn-block")
        )
      ),
      
      hr(),
      
      fluidRow(
        column(4,
               actionButton(ns("generate_pdf"), "Generate PDF Report",
                            class = "btn btn-warning btn-block")
        ),
        column(4,
               downloadButton(ns("download_pdf"), "Download PDF Report",
                              class = "btn btn-danger btn-block")
        )
      )
    )
  )
}



mod_sc_export_server <- function(id, seurat_in, markers_in, colors) {
  moduleServer(id, function(input, output, session) {
    
    rv <- reactiveValues(
      pdf_ready = NULL,
      pdf_file = NULL
    )
    
    # ============================================================
    # 1. Download RDS
    # ============================================================
    output$download_rds <- downloadHandler(
      filename = function() {
        paste0(input$export_name, ".rds")
      },
      content = function(file) {
        seu <- seurat_in()
        req(seu)
        saveRDS(seu, file)
      }
    )
    
    # ============================================================
    # 2. Download Marker CSV
    # ============================================================
    output$download_markers <- downloadHandler(
      filename = function() {
        paste0(input$export_name, "_markers.csv")
      },
      content = function(file) {
        mk <- markers_in()
        req(mk)
        write.csv(mk, file, row.names = FALSE)
      }
    )
    
    
    # ============================================================
    # 3. PDF Report Generator
    # ============================================================
    observeEvent(input$generate_pdf, {
      seu <- seurat_in()
      mk  <- markers_in()
      req(seu)
      
      showNotification("Generating PDF report...", type = "message")
      
      # File path inside temp directory
      pdf_path <- tempfile(fileext = ".pdf")
      rv$pdf_file <- pdf_path
      
      pdf(pdf_path, width = 10, height = 12)
      
      # Title
      plot.new()
      text(0.5, 0.9, paste("scRNA-seq Report:", input$export_name), cex = 1.8)
      text(0.5, 0.86, paste("Generated:", Sys.Date()), cex = 1.2)
      
      # --------------------------------------------------------------
      # UMAP by final annotation
      # --------------------------------------------------------------
      if (!is.null(seu$final_annotation)) {
        p1 <- DimPlot(seu, reduction="umap", group.by="final_annotation",
                      cols=colors()) + ggtitle("UMAP by Annotation")
        print(p1)
      }
      
      # --------------------------------------------------------------
      # UMAP by cluster
      # --------------------------------------------------------------
      p2 <- DimPlot(seu, reduction="umap", group.by="seurat_clusters") +
        ggtitle("UMAP by Clusters")
      print(p2)
      
      # --------------------------------------------------------------
      # QC distributions
      # --------------------------------------------------------------
      if ("nFeature_RNA" %in% colnames(seu@meta.data)) {
        p3 <- VlnPlot(seu, features = "nFeature_RNA", pt.size = 0.1) +
          ggtitle("nFeature RNA")
        print(p3)
      }
      
      if ("nCount_RNA" %in% colnames(seu@meta.data)) {
        p4 <- VlnPlot(seu, features = "nCount_RNA", pt.size = 0.1) +
          ggtitle("nCount RNA")
        print(p4)
      }
      
      if ("percent.mito" %in% colnames(seu@meta.data)) {
        p5 <- VlnPlot(seu, features = "percent.mito", pt.size = 0.1) +
          ggtitle("Percent Mito")
        print(p5)
      }
      
      # --------------------------------------------------------------
      # Top marker heatmap
      # --------------------------------------------------------------
      if (!is.null(mk)) {
        top_genes <- mk %>%
          group_by(cluster) %>%
          top_n(5, avg_log2FC)
        
        try({
          p6 <- DoHeatmap(seu, features = unique(top_genes$gene)) +
            ggtitle("Top 5 Marker Genes per Cluster")
          print(p6)
        }, silent = TRUE)
      }
      
      dev.off()
      
      rv$pdf_ready <- TRUE
      showNotification("PDF report created. Click Download.", type="message")
    })
    
    
    # ============================================================
    # 4. Download PDF
    # ============================================================
    output$download_pdf <- downloadHandler(
      filename = function() {
        paste0(input$export_name, "_report.pdf")
      },
      content = function(file) {
        req(rv$pdf_ready)
        file.copy(rv$pdf_file, file)
      }
    )
    
    return(list())
  })
}
