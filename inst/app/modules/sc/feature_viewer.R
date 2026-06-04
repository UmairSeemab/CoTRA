# ==============================================================
# modules/sc/feature_viewer.R
# CoTRA scRNA-seq: Interactive gene visualization
# Uses unified final_annotation + global color palette
# ==============================================================

library(Seurat)
library(dplyr)

mod_sc_feature_viewer_ui <- function(id) {
  ns <- NS(id)
  tagList(
    
    # --------------------- Input panel ---------------------
    bs4Card(
      title = "Interactive Marker Viewer",
      status = "primary",
      solidHeader = TRUE,
      width = 12,
      
      fluidRow(
        column(
          6,
          textInput(
            ns("genes"),
            "Gene list",
            placeholder = "RHO, NRL, PDE6B"
          )
        ),
        column(
          3,
          selectInput(
            ns("reduction"),
            "Dimensionality reduction",
            choices = c("UMAP" = "umap",
                        "t-SNE" = "tsne",
                        "PCA" = "pca"),
            selected = "umap"
          )
        ),
        column(
          3,
          selectInput(
            ns("group_by"),
            "Group for violin plot",
            choices = NULL
          )
        )
      ),
      
      fluidRow(
        column(
          3,
          actionButton(
            ns("apply"),
            "Update plots",
            class = "btn btn-primary btn-block"
          )
        )
      )
    ),
    
    # --------------------- Visualization panel ---------------------
    bs4Card(
      title = "Marker Expression",
      status = "info",
      solidHeader = TRUE,
      width = 12,
      tabsetPanel(
        tabPanel("FeaturePlot",
                 plotOutput(ns("feature_plot"), height = "650px")
        ),
        tabPanel("Violin Plot",
                 plotOutput(ns("vln"), height = "650px")
        )
      )
    )
  )
}

mod_sc_feature_viewer_server <- function(id, seurat_in, colors) {
  moduleServer(id, function(input, output, session) {
    
    rv <- reactiveValues(
      seurat = NULL,
      genes  = NULL
    )
    
    # ------------------------------------------------------------
    # Load Seurat object & update metadata choices
    # ------------------------------------------------------------
    observe({
      rv$seurat <- seurat_in()
      req(rv$seurat)
      
      meta <- rv$seurat@meta.data
      
      # Metadata options = all factor or character variables
      ok <- vapply(meta, function(x) is.factor(x) || is.character(x), logical(1))
      choices <- colnames(meta)[ok]
      
      # Fallback if none
      if (!"final_annotation" %in% choices)
        choices <- c("final_annotation", choices)
      
      updateSelectInput(session, "group_by",
                        choices = choices,
                        selected = choices[1])
    })
    
    # ------------------------------------------------------------
    # Parse gene list
    # ------------------------------------------------------------
    observeEvent(input$apply, {
      req(rv$seurat)
      
      txt <- input$genes
      if (txt == "" || is.null(txt)) {
        showNotification("Enter at least one gene name", type = "error")
        return(NULL)
      }
      
      # Clean and split text
      g <- unlist(strsplit(txt, "[,; \n\t]+"))
      g <- unique(g[nchar(g) > 0])
      
      # Keep only genes present in the dataset
      g_present <- intersect(g, rownames(rv$seurat))
      
      if (length(g_present) == 0) {
        showNotification("None of the entered genes are found in dataset", type = "error")
        return(NULL)
      }
      
      rv$genes <- g_present
      showNotification(paste("Using", length(g_present), "genes"), type = "message")
    })
    
    # ------------------------------------------------------------
    # FeaturePlot (UMAP / t-SNE / PCA)
    # ------------------------------------------------------------
    output$feature_plot <- renderPlot({
      req(rv$seurat, rv$genes)
      
      obj <- rv$seurat
      red <- input$reduction
      
      available <- Reductions(obj)
      if (!(red %in% available)) {
        # auto fallback
        if ("umap" %in% available) red <- "umap"
        else if ("tsne" %in% available) red <- "tsne"
        else red <- "pca"
      }
      
      FeaturePlot(
        obj,
        features = rv$genes,
        reduction = red,
        cols = c("lightgrey", "red"),
        order = TRUE
      )
    })
    
    # ------------------------------------------------------------
    # Violin Plot
    # ------------------------------------------------------------
    output$vln <- renderPlot({
      req(rv$seurat, rv$genes)
      
      group <- input$group_by
      
      VlnPlot(
        rv$seurat,
        features = rv$genes,
        group.by = group,
        cols = colors()
      )
    })
  })
}
