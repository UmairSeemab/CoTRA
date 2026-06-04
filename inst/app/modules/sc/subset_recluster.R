# ==========================================================
# subset_recluster.R
# CoTRA scRNA-seq Module
# Subsetting + Re-clustering
# ==========================================================
# This module lets the user:
#   • Select clusters to keep
#   • Optionally give a new subset name
#   • Create a new Seurat object subset
#   • Re-run HVG selection, PCA, UMAP, and clustering
#   • Return the new Seurat object to server.R
#
# Dependencies: Seurat, bs4Dash, shiny, shinyWidgets, plotly
# No dependency on doublet detection or annotations.
# ==========================================================

# ================= UI ======================

# modules/sc/subset_recluster.R
# Subset selected clusters or metadata groups and re run HVG, PCA, UMAP and clustering

mod_sc_subset_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    bs4Card(
      title = "Subset and Re cluster",
      status = "primary",
      solidHeader = TRUE,
      width = 12,
      collapsible = TRUE,
      
      fluidRow(
        column(
          width = 4,
          selectInput(
            ns("cluster_source"),
            "Cluster / annotation column",
            choices = c("Idents" = "idents"),
            selected = "idents"
          )
        ),
        column(
          width = 4,
          uiOutput(ns("cluster_choices_ui"))
        ),
        column(
          width = 4,
          textInput(
            ns("subset_name"),
            "Name of new subset object",
            value = "subset_reclustered"
          )
        )
      ),
      
      fluidRow(
        column(
          width = 3,
          numericInput(
            ns("nfeatures"),
            "Number of variable features",
            value = 2000,
            min = 200,
            max = 5000,
            step = 100
          )
        ),
        column(
          width = 3,
          numericInput(
            ns("ndims"),
            "Number of PCA dimensions",
            value = 30,
            min = 5,
            max = 100,
            step = 1
          )
        ),
        column(
          width = 3,
          numericInput(
            ns("resolution"),
            "Clustering resolution",
            value = 0.8,
            min = 0.1,
            max = 3,
            step = 0.1
          )
        ),
        column(
          width = 3,
          actionButton(
            ns("run_subset"),
            "Create subset and re cluster",
            class = "btn btn-primary btn-block"
          )
        )
      )
    ),
    
    bs4Card(
      title = "Subset overview and quick plots",
      status = "info",
      solidHeader = TRUE,
      width = 12,
      collapsible = TRUE,
      
      fluidRow(
        column(
          width = 4,
          verbatimTextOutput(ns("subset_summary"))
        ),
        column(
          width = 4,
          plotOutput(ns("subset_umap"), height = 400)
        ),
        column(
          width = 4,
          plotOutput(ns("subset_pca_scatter"), height = 400)
        )
      )
    )
  )
}


mod_sc_subset_server <- function(id, seurat_in) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rv <- reactiveValues(
      parent = NULL,
      subset = NULL
    )
    
    observe({
      obj <- seurat_in()
      req(obj)
      rv$parent <- obj
      
      cols <- colnames(obj@meta.data)
      factor_cols <- cols[vapply(obj@meta.data, is.factor, logical(1))]
      choices <- c("Idents" = "idents", factor_cols)
      
      updateSelectInput(
        session,
        "cluster_source",
        choices = choices,
        selected = isolate(input$cluster_source %||% "idents")
      )
    })
    
    output$cluster_choices_ui <- renderUI({
      req(rv$parent)
      obj <- rv$parent
      
      if (input$cluster_source == "idents") {
        clusters <- levels(Idents(obj))
      } else {
        col <- input$cluster_source
        if (!col %in% colnames(obj@meta.data)) {
          return(NULL)
        }
        clusters <- sort(unique(as.character(obj@meta.data[[col]])))
      }
      
      pickerInput(
        ns("clusters_to_keep"),
        "Clusters or groups to keep",
        choices = clusters,
        multiple = TRUE,
        options = list(`actions-box` = TRUE)
      )
    })
    
    observeEvent(input$run_subset, {
      req(rv$parent)
      obj <- rv$parent
      req(input$clusters_to_keep)
      subset_name <- input$subset_name
      if (is.null(subset_name) || subset_name == "") {
        subset_name <- "subset_reclustered"
      }
      
      showNotification("Creating subset and running HVG, PCA, UMAP and clustering", type = "message")
      
      if (input$cluster_source == "idents") {
        cells_keep <- WhichCells(obj, idents = input$clusters_to_keep)
      } else {
        col <- input$cluster_source
        cells_keep <- rownames(obj@meta.data)[obj@meta.data[[col]] %in% input$clusters_to_keep]
      }
      
      if (length(cells_keep) < 50) {
        showNotification("Subset has fewer than 50 cells, aborting", type = "error")
        return(NULL)
      }
      
      sub_obj <- subset(obj, cells = cells_keep)
      
      DefaultAssay(sub_obj) <- "RNA"
      
      sub_obj <- NormalizeData(
        sub_obj,
        normalization.method = "LogNormalize",
        scale.factor = 10000,
        verbose = FALSE
      )
      
      sub_obj <- FindVariableFeatures(
        sub_obj,
        selection.method = "vst",
        nfeatures = input$nfeatures,
        verbose = FALSE
      )
      
      all_feats <- VariableFeatures(sub_obj)
      sub_obj <- ScaleData(sub_obj, features = all_feats, verbose = FALSE)
      
      sub_obj <- RunPCA(
        sub_obj,
        features = all_feats,
        npcs = input$ndims,
        verbose = FALSE
      )
      
      sub_obj <- FindNeighbors(
        sub_obj,
        dims = seq_len(input$ndims),
        verbose = FALSE
      )
      
      sub_obj <- FindClusters(
        sub_obj,
        resolution = input$resolution,
        verbose = FALSE
      )
      
      sub_obj <- RunUMAP(
        sub_obj,
        dims = seq_len(input$ndims),
        reduction = "pca",
        reduction.name = "umap",
        reduction.key = "UMAP_",
        verbose = FALSE
      )
      
      sub_obj$parent_cluster_source <- input$cluster_source
      if (input$cluster_source == "idents") {
        sub_obj$parent_cluster <- as.character(Idents(obj)[Cells(sub_obj)])
      } else {
        col <- input$cluster_source
        sub_obj$parent_cluster <- as.character(obj@meta.data[Cells(sub_obj), col])
      }
      
      sub_obj@misc$subset_name <- subset_name
      
      rv$subset <- sub_obj
      
      showNotification(
        paste0(
          "Subset ", subset_name,
          " created. Cells: ", ncol(sub_obj),
          ", clusters: ", length(levels(Idents(sub_obj)))
        ),
        type = "message"
      )
    })
    
    output$subset_summary <- renderPrint({
      req(rv$subset)
      obj <- rv$subset
      cat("Subset name:", obj@misc$subset_name %||% "subset_reclustered", "\n")
      cat("Cells:", ncol(obj), "\n")
      cat("Genes:", nrow(obj), "\n")
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        cat("\nCluster sizes:\n")
        print(table(obj$seurat_clusters))
      }
      if ("parent_cluster" %in% colnames(obj@meta.data)) {
        cat("\nParent cluster mapping:\n")
        print(table(obj$parent_cluster, obj$seurat_clusters))
      }
    })
    
    output$subset_umap <- renderPlot({
      req(rv$subset)
      obj <- rv$subset
      if (!"umap" %in% Reductions(obj)) {
        plot.new()
        title("UMAP not present on subset")
        return()
      }
      DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
        ggplot2::ggtitle("Subset UMAP, new clusters")
    })
    
    output$subset_pca_scatter <- renderPlot({
      req(rv$subset)
      obj <- rv$subset
      if (!"pca" %in% Reductions(obj)) {
        plot.new()
        title("PCA not present on subset")
        return()
      }
      EmbDim <- Embeddings(obj[["pca"]])[, 1:2, drop = FALSE]
      df <- as.data.frame(EmbDim)
      df$cluster <- if ("seurat_clusters" %in% colnames(obj@meta.data))
        obj$seurat_clusters else "all"
      
      ggplot2::ggplot(
        df,
        ggplot2::aes(
          x = PC_1,
          y = PC_2,
          color = cluster
        )
      ) +
        ggplot2::geom_point(size = 0.4, alpha = 0.8) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle("Subset PCA, first two components")
    })
    
    list(
      seurat = reactive(rv$subset)
    )
  })
}
