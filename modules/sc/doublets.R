# ==========================================================
# File: doublets.R
# Module: mod_sc_doublets_ui / mod_sc_doublets_server
# Purpose:
#   Detect doublets using DoubletFinder or scDblFinder
#   Visualize doublet scores + classifications
#   Remove predicted doublets
#
# Dependencies:
#   library(Seurat)
#   library(DoubletFinder)
#   library(scDblFinder)
# ==========================================================
mod_sc_doublets_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    bs4Card(
      title = "Doublet Detection",
      status = "primary",
      solidHeader = TRUE,
      width = 12,
      collapsible = TRUE,
      
      fluidRow(
        column(
          width = 3,
          selectInput(
            ns("method"),
            "Doublet detection method",
            choices = c(
              "DoubletFinder" = "doubletfinder",
              "scDblFinder"   = "scdblfinder"
            ),
            selected = "doubletfinder"
          )
        ),
        column(
          width = 3,
          numericInput(
            ns("doublet_rate"),
            "Expected doublet rate",
            value = 0.05,
            min = 0,
            max = 0.3,
            step = 0.01
          )
        ),
        column(
          width = 3,
          numericInput(
            ns("pN"),
            "pN (DoubletFinder)",
            value = 0.25,
            min = 0.01,
            max = 1,
            step = 0.01
          )
        ),
        column(
          width = 3,
          numericInput(
            ns("pK"),
            "pK (DoubletFinder)",
            value = 0.005,
            min = 0.001,
            max = 0.2,
            step = 0.001
          )
        )
      ),
      
      fluidRow(
        column(
          width = 3,
          checkboxInput(
            ns("use_cluster_as_group"),
            "Use clusters as grouping variable",
            value = TRUE
          )
        ),
        column(
          width = 3,
          checkboxInput(
            ns("remove_doublets"),
            "Create filtered object without doublets",
            value = TRUE
          )
        ),
        column(
          width = 3,
          actionButton(
            ns("run_doublets"),
            "Run doublet detection",
            class = "btn btn-primary btn-block"
          )
        )
      )
    ),
    
    bs4Card(
      title = "Doublet summary and plots",
      status = "info",
      solidHeader = TRUE,
      width = 12,
      collapsible = TRUE,
      
      tabsetPanel(
        tabPanel(
          "Summary",
          br(),
          verbatimTextOutput(ns("summary_text")),
          DT::DTOutput(ns("table_doublets"))
        ),
        tabPanel(
          "UMAP plot",
          br(),
          plotOutput(ns("umap_doublets"), height = 450),
          br(),
          downloadButton(ns("download_umap_png"), "Download UMAP PNG")
        ),
        tabPanel(
          "QC scatter",
          br(),
          plotOutput(ns("qc_scatter"), height = 450),
          br(),
          downloadButton(ns("download_qc_png"), "Download QC PNG")
        )
      )
    )
  )
}


mod_sc_doublets_server <- function(id, seurat_in) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rv <- reactiveValues(
      seurat = NULL,
      doublets_df = NULL
    )
    
    observe({
      obj <- seurat_in()
      req(obj)
      rv$seurat <- obj
    })
    
    has_doubletfinder <- reactive({
      requireNamespace("DoubletFinder", quietly = TRUE)
    })
    
    has_scdblfinder <- reactive({
      requireNamespace("scDblFinder", quietly = TRUE) &&
        requireNamespace("SingleCellExperiment", quietly = TRUE)
    })
    
    observeEvent(input$run_doublets, {
      req(rv$seurat)
      obj <- rv$seurat
      
      if (!"pca" %in% Reductions(obj)) {
        showNotification("PCA missing. Run PCA first.", type = "error")
        return(NULL)
      }
      if (!"RNA_snn" %in% names(obj@graphs)) {
        showNotification("Neighbor graph missing. Run clustering first.", type = "error")
        return(NULL)
      }
      
      n_cells <- ncol(obj)
      n_doublets <- round(n_cells * input$doublet_rate)
      
      if (input$method == "doubletfinder") {
        
        if (!has_doubletfinder()) {
          showNotification("DoubletFinder not installed.", type = "error")
          return(NULL)
        }
        
        library(DoubletFinder)
        pcs_use <- seq_len(min(50, ncol(obj@reductions$pca@cell.embeddings)))
        
        obj <- DoubletFinder::paramSweep_v3(obj, PCs = pcs_use, sct = FALSE)
        sweep_stats <- DoubletFinder::summarizeSweep(obj, GT = FALSE)
        
        best_pK <- input$pK
        
        obj <- DoubletFinder::doubletFinder_v3(
          obj,
          PCs = pcs_use,
          pN = input$pN,
          pK = best_pK,
          nExp = n_doublets,
          reuse.pANN = FALSE,
          sct = FALSE
        )
        
        meta_names <- colnames(obj@meta.data)
        pann_col  <- meta_names[grepl("^pANN_", meta_names)][1]
        class_col <- meta_names[grepl("^DF.classifications_", meta_names)][1]
        
        obj$doublet_score  <- obj@meta.data[[pann_col]]
        obj$doublet_class  <- as.character(obj@meta.data[[class_col]])
        obj$is_doublet     <- obj$doublet_class == "Doublet"
        
      } else if (input$method == "scdblfinder") {
        
        if (!has_scdblfinder()) {
          showNotification("scDblFinder not installed.", type = "error")
          return(NULL)
        }
        
        library(scDblFinder)
        library(SingleCellExperiment)
        
        sce <- as.SingleCellExperiment(obj)
        sce <- scDblFinder::scDblFinder(sce)
        
        obj$doublet_score <- sce$scDblFinder.score
        obj$doublet_class <- as.character(sce$scDblFinder.class)
        obj$is_doublet    <- obj$doublet_class == "doublet"
      }
      
      df <- data.frame(
        cell = colnames(obj),
        cluster = if ("seurat_clusters" %in% colnames(obj@meta.data))
          obj$seurat_clusters else NA,
        doublet_score = obj$doublet_score,
        doublet_class = obj$doublet_class,
        is_doublet = obj$is_doublet,
        stringsAsFactors = FALSE
      )
      
      rv$seurat <- obj
      rv$doublets_df <- df
      
      showNotification(
        paste0("Doublet detection finished. Predicted doublets: ",
               sum(obj$is_doublet), " / ", n_cells),
        type = "message"
      )
    })
    
    seurat_filtered <- reactive({
      req(rv$seurat)
      obj <- rv$seurat
      if (!"is_doublet" %in% colnames(obj@meta.data)) return(obj)
      if (!isTRUE(input$remove_doublets)) return(obj)
      subset(obj, subset = !is_doublet)
    })
    
    output$summary_text <- renderPrint({
      req(rv$seurat)
      obj <- rv$seurat
      if (!"is_doublet" %in% colnames(obj@meta.data)) {
        cat("Run doublet detection first.")
        return()
      }
      
      n_cells <- ncol(obj)
      n_doublets <- sum(obj$is_doublet)
      rate <- n_doublets / n_cells
      
      cat("Total cells:", n_cells, "\n")
      cat("Doublets:", n_doublets, "\n")
      cat("Doublet rate:", round(rate, 4), "\n\n")
      
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        print(table(obj$seurat_clusters, obj$is_doublet))
      }
    })
    
    output$table_doublets <- DT::renderDT({
      req(rv$doublets_df)
      DT::datatable(rv$doublets_df, options = list(pageLength = 20, scrollX = TRUE))
    })
    
    output$umap_doublets <- renderPlot({
      req(rv$seurat)
      obj <- rv$seurat
      if (!"umap" %in% Reductions(obj)) {
        plot.new(); title("UMAP missing.")
        return()
      }
      if (!"doublet_class" %in% colnames(obj@meta.data)) {
        plot.new(); title("Doublet class missing.")
        return()
      }
      DimPlot(obj, reduction = "umap", group.by = "doublet_class", pt.size = 0.3)
    })
    
    output$qc_scatter <- renderPlot({
      req(rv$seurat)
      obj <- rv$seurat
      df <- obj@meta.data
      
      if (!all(c("nFeature_RNA", "nCount_RNA") %in% colnames(df))) {
        plot.new(); title("QC metrics missing.")
        return()
      }
      if (!"is_doublet" %in% colnames(df)) {
        plot.new(); title("Doublet info missing.")
        return()
      }
      
      ggplot2::ggplot(df, ggplot2::aes(
        x = nCount_RNA, y = nFeature_RNA, color = is_doublet
      )) +
        ggplot2::geom_point(size = 0.4, alpha = 0.7) +
        ggplot2::scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
        ggplot2::theme_bw()
    })
    
    output$download_umap_png <- downloadHandler(
      filename = function() paste0("doublets_umap_", Sys.Date(), ".png"),
      content = function(file) {
        req(rv$seurat)
        obj <- rv$seurat
        png(file, width = 1000, height = 800)
        print(DimPlot(obj, reduction = "umap", group.by = "doublet_class"))
        dev.off()
      }
    )
    
    output$download_qc_png <- downloadHandler(
      filename = function() paste0("doublets_qc_", Sys.Date(), ".png"),
      content = function(file) {
        req(rv$seurat)
        df <- rv$seurat@meta.data
        png(file, width = 1000, height = 800)
        print(
          ggplot2::ggplot(df, ggplot2::aes(
            x = nCount_RNA, y = nFeature_RNA, color = is_doublet
          )) +
            ggplot2::geom_point(size = 0.4, alpha = 0.7)
        )
        dev.off()
      }
    )
    
    list(
      seurat = reactive(rv$seurat),
      seurat_filtered = seurat_filtered,
      doublets = reactive(rv$doublets_df)
    )
  })
}

