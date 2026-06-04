# ==========================================================
# modules/sc/tsne_umap.R
# CoTRA scRNA-seq t-SNE + UMAP Module
#
# Features:
# - Runs UMAP and/or t-SNE from PCA embeddings
# - Uses recommended PCs from PCA module when available
# - User control for dimensions, UMAP n.neighbors/min.dist,
#   and t-SNE perplexity
# - Color by metadata, QC metrics, or later cluster labels
# - PCA existence and PC number checks
# - UMAP and t-SNE plots
# - Download plots as PDF/SVG
# - Download embedding coordinates as CSV
# - Collapsible help and interpretation sections
# - Session export to outputs/scRNA/sessioninfo/
# - Downstream-safe return object for clustering
# ==========================================================

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

mod_sc_tsne_umap_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("t-SNE and UMAP"),
    br(),
    
    bs4Dash::bs4Card(
      title = "Help: how to use t-SNE and UMAP",
      width = 12,
      solidHeader = TRUE,
      status = "info",
      collapsible = TRUE,
      collapsed = TRUE,
      icon = icon("circle-info"),
      tags$h4("Purpose"),
      tags$p("This module converts PCA structure into two-dimensional maps for visualizing cells. UMAP and t-SNE help inspect cell populations, sample effects, gradients, and potential outliers."),
      tags$h4("Recommended workflow"),
      tags$ol(
        tags$li("Run QC, gene selection, and PCA first."),
        tags$li("Start with the recommended number of PCs from the PCA module."),
        tags$li("Run UMAP first for general visualization."),
        tags$li("Use t-SNE as a complementary view if needed."),
        tags$li("Continue to clustering after embeddings are generated.")
      ),
      tags$h4("UMAP versus t-SNE"),
      tags$ul(
        tags$li(tags$b("UMAP:"), " usually preserves both local and broader structure better and is commonly used before clustering visualization."),
        tags$li(tags$b("t-SNE:"), " focuses strongly on local neighborhood structure and can separate groups visually, but distances between distant groups are harder to interpret."),
        tags$li("Neither method proves cell identity. Use marker genes and metadata for interpretation.")
      ),
      tags$h4("Important caution"),
      tags$p("UMAP and t-SNE are visualization methods. Cluster separation, shape, and distance can change with parameters. Do not treat visual distance alone as biological distance.")
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "Embedding settings",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        fluidRow(
          column(
            3,
            checkboxGroupInput(
              ns("methods_to_run"),
              "Methods to run",
              choices = c("UMAP" = "umap", "t-SNE" = "tsne"),
              selected = c("umap")
            )
          ),
          column(
            3,
            numericInput(ns("dims_start"), "First PC", value = 1, min = 1, max = 100, step = 1)
          ),
          column(
            3,
            numericInput(ns("dims_end"), "Last PC", value = 30, min = 2, max = 100, step = 1)
          ),
          column(
            3,
            br(),
            actionButton(ns("use_recommended_pcs"), "Use PCA recommendation", class = "btn-secondary")
          )
        ),
        hr(),
        h4("UMAP parameters"),
        fluidRow(
          column(3, numericInput(ns("umap_n_neighbors"), "n.neighbors", value = 30, min = 5, max = 200, step = 5)),
          column(3, numericInput(ns("umap_min_dist"), "min.dist", value = 0.3, min = 0.01, max = 0.99, step = 0.01)),
          column(3, selectInput(ns("umap_metric"), "Metric", choices = c("cosine", "euclidean", "correlation"), selected = "cosine")),
          column(3, numericInput(ns("random_seed"), "Random seed", value = 1234, min = 1, max = 999999, step = 1))
        ),
        hr(),
        h4("t-SNE parameters"),
        fluidRow(
          column(3, numericInput(ns("tsne_perplexity"), "Perplexity", value = 30, min = 5, max = 100, step = 5)),
          column(3, numericInput(ns("tsne_max_iter"), "Max iterations", value = 1000, min = 250, max = 5000, step = 250)),
          column(3, checkboxInput(ns("check_duplicates"), "Check duplicate cells", TRUE)),
          column(3, br(), actionButton(ns("run_embeddings"), "Run t-SNE/UMAP", class = "btn-primary"))
        ),
        hr(),
        fluidRow(
          column(4, downloadButton(ns("download_embedding_session"), "Save embedding session (.rds)", class = "btn-success")),
          column(8, uiOutput(ns("session_save_hint")))
        ),
        br(),
        uiOutput(ns("embedding_status_box"))
      )
    ),
    
    bs4Dash::bs4Card(
      title = "Interpretation: embedding parameters",
      width = 12,
      solidHeader = TRUE,
      status = "info",
      collapsible = TRUE,
      collapsed = TRUE,
      icon = icon("circle-info"),
      tags$h4("How to choose PCs"),
      tags$ul(
        tags$li("Start with the recommended PC number from the PCA module."),
        tags$li("Too few PCs can hide real cell populations."),
        tags$li("Too many PCs can add noise or technical variation."),
        tags$li("If UMAP is dominated by batch or QC metrics, revisit QC, PCA, or regression choices.")
      ),
      tags$h4("UMAP parameters"),
      tags$ul(
        tags$li(tags$b("n.neighbors:"), " larger values emphasize broader structure. Smaller values emphasize local structure."),
        tags$li(tags$b("min.dist:"), " lower values make clusters appear tighter. Higher values spread cells more smoothly."),
        tags$li(tags$b("metric:"), " cosine is commonly useful for PCA-based single-cell UMAP.")
      ),
      tags$h4("t-SNE perplexity"),
      tags$p("Perplexity controls the effective neighborhood size. It should be smaller than the number of cells and often works between 20 and 50 for many datasets.")
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "Embedding visualization",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        fluidRow(
          column(4, uiOutput(ns("color_by_ui"))),
          column(4, checkboxInput(ns("show_labels"), "Show labels when identity/class is selected", FALSE)),
          column(4, sliderInput(ns("point_size"), "Point size", min = 0.1, max = 3, value = 0.6, step = 0.1))
        ),
        tabsetPanel(
          tabPanel(
            "UMAP",
            br(),
            shinycssloaders::withSpinner(plotly::plotlyOutput(ns("umap_plot")), type = 4),
            br(),
            fluidRow(
              column(3, downloadButton(ns("dl_umap_pdf"), "UMAP PDF")),
              column(3, downloadButton(ns("dl_umap_svg"), "UMAP SVG")),
              column(3, downloadButton(ns("dl_umap_coords"), "UMAP coordinates CSV"))
            ),
            br(),
            bs4Dash::bs4Card(
              title = "Interpretation: UMAP plot",
              width = 12,
              solidHeader = TRUE,
              status = "info",
              collapsible = TRUE,
              collapsed = TRUE,
              icon = icon("circle-info"),
              tags$p("UMAP shows cells arranged by similarity in PCA space. Neighboring cells usually have similar expression profiles."),
              tags$ul(
                tags$li("Separated islands may represent cell types, cell states, batch effects, or technical artifacts."),
                tags$li("Gradients may represent continuous transitions such as differentiation, activation, or stress."),
                tags$li("Color by sample, condition, batch, and QC metrics before biological interpretation."),
                tags$li("Use marker genes and clustering to support cell identity assignments.")
              )
            )
          ),
          tabPanel(
            "t-SNE",
            br(),
            shinycssloaders::withSpinner(plotly::plotlyOutput(ns("tsne_plot")), type = 4),
            br(),
            fluidRow(
              column(3, downloadButton(ns("dl_tsne_pdf"), "t-SNE PDF")),
              column(3, downloadButton(ns("dl_tsne_svg"), "t-SNE SVG")),
              column(3, downloadButton(ns("dl_tsne_coords"), "t-SNE coordinates CSV"))
            ),
            br(),
            bs4Dash::bs4Card(
              title = "Interpretation: t-SNE plot",
              width = 12,
              solidHeader = TRUE,
              status = "info",
              collapsible = TRUE,
              collapsed = TRUE,
              icon = icon("circle-info"),
              tags$p("t-SNE emphasizes local neighborhoods. It can show compact groups, but global distances and group sizes are difficult to interpret."),
              tags$ul(
                tags$li("Nearby cells are likely similar."),
                tags$li("Distances between far-separated islands should not be overinterpreted."),
                tags$li("Different perplexity or random seed can change the layout."),
                tags$li("Use t-SNE as a complementary view, not as the only basis for conclusions.")
              )
            )
          ),
          tabPanel(
            "Embedding summary",
            br(),
            DT::DTOutput(ns("embedding_summary_table")),
            br(),
            downloadButton(ns("dl_embedding_summary"), "Download embedding summary CSV")
          )
        )
      )
    )
  )
}

mod_sc_tsne_umap_server <- function(id, seurat_r, sc_state = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    embedding_seu <- reactiveVal(NULL)
    embedding_ready <- reactiveVal(FALSE)
    umap_plot_gg <- reactiveVal(NULL)
    tsne_plot_gg <- reactiveVal(NULL)
    embedding_summary <- reactiveVal(NULL)
    last_session_path <- reactiveVal(NULL)
    
    get_seu <- reactive({
      if (!is.null(sc_state) && !is.null(sc_state$seurat)) return(sc_state$seurat)
      seurat_r()
    })
    
    safe_theme <- function() {
      if (exists("theme_cotra", mode = "function")) {
        theme_cotra()
      } else {
        ggplot2::theme_bw()
      }
    }
    
    has_pca <- function(seu) {
      "pca" %in% names(seu@reductions)
    }
    
    pca_npcs <- function(seu) {
      emb <- tryCatch(Seurat::Embeddings(seu, reduction = "pca"), error = function(e) NULL)
      if (is.null(emb)) return(0)
      ncol(emb)
    }
    
    get_recommended_pcs <- function(seu) {
      rec <- tryCatch(seu@misc$CoTRA_pca$recommended_pcs, error = function(e) NULL)
      if (is.null(rec) || !is.numeric(rec) || length(rec) == 0 || is.na(rec)) {
        rec <- min(30, pca_npcs(seu))
      }
      rec <- max(2, min(as.integer(rec[1]), pca_npcs(seu)))
      rec
    }
    
    observe({
      seu <- get_seu()
      if (is.null(seu)) return(NULL)
      if (!has_pca(seu)) return(NULL)
      rec <- get_recommended_pcs(seu)
      isolate({
        updateNumericInput(session, "dims_end", value = rec, max = max(2, pca_npcs(seu)))
        updateNumericInput(session, "dims_start", value = 1, max = max(1, pca_npcs(seu)))
      })
    })
    
    observeEvent(input$use_recommended_pcs, {
      seu <- get_seu()
      req(seu)
      if (!has_pca(seu)) {
        showModal(modalDialog(title = "PCA missing", "Run PCA before using recommended PCs.", easyClose = TRUE))
        return(NULL)
      }
      rec <- get_recommended_pcs(seu)
      updateNumericInput(session, "dims_start", value = 1)
      updateNumericInput(session, "dims_end", value = rec)
    })
    
    validate_dims <- function(seu) {
      if (!has_pca(seu)) {
        return(list(ok = FALSE, message = "PCA reduction was not found. Run PCA first."))
      }
      
      total_pcs <- pca_npcs(seu)
      d1 <- as.integer(input$dims_start %||% 1)
      d2 <- as.integer(input$dims_end %||% min(30, total_pcs))
      
      if (is.na(d1) || is.na(d2)) {
        return(list(ok = FALSE, message = "PC dimensions must be numeric."))
      }
      if (d1 < 1) d1 <- 1
      if (d2 > total_pcs) d2 <- total_pcs
      if (d1 > d2) {
        return(list(ok = FALSE, message = "First PC cannot be larger than last PC."))
      }
      if (length(seq(d1, d2)) < 2) {
        return(list(ok = FALSE, message = "At least two PCs are required."))
      }
      
      list(ok = TRUE, dims = seq(d1, d2), total_pcs = total_pcs)
    }
    
    color_choices <- reactive({
      seu <- embedding_seu()
      if (is.null(seu)) seu <- get_seu()
      if (is.null(seu)) return(c("None" = "none"))
      md <- seu@meta.data
      c("None" = "none", "Active identity" = "ident", colnames(md))
    })
    
    output$color_by_ui <- renderUI({
      choices <- color_choices()
      selected <- if ("seurat_clusters" %in% choices) "seurat_clusters" else if ("orig.ident" %in% choices) "orig.ident" else "none"
      selectInput(ns("color_by"), "Color cells by", choices = choices, selected = selected)
    })
    
    build_embedding_df <- function(seu, reduction = "umap") {
      emb <- tryCatch(Seurat::Embeddings(seu, reduction = reduction), error = function(e) NULL)
      if (is.null(emb) || ncol(emb) < 2) return(NULL)
      
      df <- as.data.frame(emb[, 1:2, drop = FALSE])
      colnames(df) <- c("Dim1", "Dim2")
      df$cell <- rownames(df)
      
      color_by <- input$color_by %||% "none"
      md <- seu@meta.data
      
      if (color_by == "ident") {
        df$color_value <- as.character(Seurat::Idents(seu)[rownames(df)])
        color_label <- "Identity"
      } else if (!is.null(color_by) && color_by != "none" && color_by %in% colnames(md)) {
        df$color_value <- md[rownames(df), color_by]
        color_label <- color_by
      } else {
        df$color_value <- "Cells"
        color_label <- "Cells"
      }
      
      df$color_label <- color_label
      df
    }
    
    make_embedding_plot <- function(seu, reduction = "umap") {
      df <- build_embedding_df(seu, reduction = reduction)
      if (is.null(df)) return(NULL)
      
      title <- if (reduction == "umap") "UMAP" else "t-SNE"
      xlab <- if (reduction == "umap") "UMAP_1" else "tSNE_1"
      ylab <- if (reduction == "umap") "UMAP_2" else "tSNE_2"
      ps <- input$point_size %||% 0.6
      
      p <- ggplot2::ggplot(df, ggplot2::aes(x = Dim1, y = Dim2, color = color_value, text = cell)) +
        ggplot2::geom_point(size = ps, alpha = 0.85) +
        ggplot2::labs(title = title, x = xlab, y = ylab, color = unique(df$color_label)[1]) +
        safe_theme()
      
      if (isTRUE(input$show_labels) && unique(df$color_label)[1] != "Cells") {
        label_df <- stats::aggregate(cbind(Dim1, Dim2) ~ color_value, data = df, FUN = median)
        p <- p + ggplot2::geom_text(
          data = label_df,
          ggplot2::aes(x = Dim1, y = Dim2, label = color_value),
          inherit.aes = FALSE,
          size = 3,
          fontface = "bold"
        )
      }
      
      p
    }
    
    output$embedding_status_box <- renderUI({
      seu <- embedding_seu()
      base <- get_seu()
      
      if (is.null(base)) {
        return(tags$div(style = "color:#777;", "No Seurat object available."))
      }
      if (!has_pca(base)) {
        return(tags$div(style = "color:#B00020;", "PCA is missing. Run PCA before t-SNE/UMAP."))
      }
      if (is.null(seu) || !isTRUE(embedding_ready())) {
        return(tags$div(style = "color:#777;", paste0("PCA found with ", pca_npcs(base), " PCs. Run t-SNE/UMAP to generate embeddings.")))
      }
      
      has_umap <- "umap" %in% names(seu@reductions)
      has_tsne <- "tsne" %in% names(seu@reductions)
      info <- tryCatch(seu@misc$CoTRA_tsne_umap, error = function(e) NULL)
      
      tags$div(
        tags$b("Embedding status: "), "completed", tags$br(),
        tags$b("UMAP available: "), ifelse(has_umap, "Yes", "No"), tags$br(),
        tags$b("t-SNE available: "), ifelse(has_tsne, "Yes", "No"), tags$br(),
        tags$b("PCs used: "), if (!is.null(info)) paste(info$dims, collapse = ", ") else "NA"
      )
    })
    
    output$session_save_hint <- renderUI({
      if (!isTRUE(embedding_ready()) || is.null(embedding_seu())) {
        return(tags$span(style = "color:#777;", "Run UMAP/t-SNE before saving this step."))
      }
      if (is.null(last_session_path())) {
        return(tags$span(style = "color:#777;", "Save the current embedding state as a CoTRA project RDS."))
      }
      tags$span(style = "color:#2E7D32;", paste("Last saved:", last_session_path()))
    })
    
    observeEvent(input$run_embeddings, {
      seu <- get_seu()
      req(seu)
      
      methods <- input$methods_to_run %||% character(0)
      if (length(methods) == 0) {
        showModal(modalDialog(title = "No method selected", "Select UMAP, t-SNE, or both.", easyClose = TRUE))
        return(NULL)
      }
      
      dim_check <- validate_dims(seu)
      if (!isTRUE(dim_check$ok)) {
        showModal(modalDialog(title = "PC check failed", dim_check$message, easyClose = TRUE))
        return(NULL)
      }
      
      dims <- dim_check$dims
      n_cells <- ncol(seu)
      
      withProgress(message = "Running t-SNE/UMAP", value = 0.1, {
        set.seed(as.integer(input$random_seed %||% 1234))
        
        if ("umap" %in% methods) {
          incProgress(0.35, detail = "Running UMAP")
          
          nn <- as.integer(input$umap_n_neighbors %||% 30)
          nn <- max(2, min(nn, n_cells - 1))
          mdist <- as.numeric(input$umap_min_dist %||% 0.3)
          mdist <- max(0.001, min(mdist, 0.99))
          
          seu <- tryCatch(
            Seurat::RunUMAP(
              seu,
              reduction = "pca",
              dims = dims,
              n.neighbors = nn,
              min.dist = mdist,
              metric = input$umap_metric %||% "cosine",
              reduction.name = "umap",
              reduction.key = "UMAP_",
              seed.use = as.integer(input$random_seed %||% 1234),
              verbose = FALSE
            ),
            error = function(e) {
              showModal(modalDialog(title = "RunUMAP failed", conditionMessage(e), easyClose = TRUE))
              NULL
            }
          )
          if (is.null(seu)) return(NULL)
        }
        
        if ("tsne" %in% methods) {
          incProgress(0.7, detail = "Running t-SNE")
          
          perplexity <- as.numeric(input$tsne_perplexity %||% 30)
          max_perplexity <- floor((n_cells - 1) / 3)
          perplexity <- max(2, min(perplexity, max_perplexity))
          
          seu <- tryCatch(
            Seurat::RunTSNE(
              seu,
              reduction = "pca",
              dims = dims,
              perplexity = perplexity,
              max_iter = as.integer(input$tsne_max_iter %||% 1000),
              check_duplicates = isTRUE(input$check_duplicates),
              reduction.name = "tsne",
              reduction.key = "tSNE_",
              seed.use = as.integer(input$random_seed %||% 1234),
              verbose = FALSE
            ),
            error = function(e) {
              showModal(modalDialog(title = "RunTSNE failed", conditionMessage(e), easyClose = TRUE))
              NULL
            }
          )
          if (is.null(seu)) return(NULL)
        }
        
        summary_df <- data.frame(
          method = c("UMAP", "t-SNE"),
          available = c("umap" %in% names(seu@reductions), "tsne" %in% names(seu@reductions)),
          cells = ncol(seu),
          pcs_used = paste(dims, collapse = ","),
          random_seed = as.integer(input$random_seed %||% 1234),
          stringsAsFactors = FALSE
        )
        
        seu@misc$CoTRA_tsne_umap <- list(
          methods = methods,
          dims = dims,
          umap = list(
            n.neighbors = as.integer(input$umap_n_neighbors %||% 30),
            min.dist = as.numeric(input$umap_min_dist %||% 0.3),
            metric = input$umap_metric %||% "cosine"
          ),
          tsne = list(
            perplexity = as.numeric(input$tsne_perplexity %||% 30),
            max_iter = as.integer(input$tsne_max_iter %||% 1000),
            check_duplicates = isTRUE(input$check_duplicates)
          ),
          random_seed = as.integer(input$random_seed %||% 1234),
          summary = summary_df,
          created_at = Sys.time()
        )
        
        embedding_seu(seu)
        embedding_ready(TRUE)
        embedding_summary(summary_df)
        umap_plot_gg(make_embedding_plot(seu, "umap"))
        tsne_plot_gg(make_embedding_plot(seu, "tsne"))
        
        if (!is.null(sc_state)) {
          sc_state$seurat <- seu
          sc_state$parameters$tsne_umap <- seu@misc$CoTRA_tsne_umap
        }
        
        incProgress(1)
      })
    })
    
    observe({
      seu <- embedding_seu()
      if (is.null(seu)) return(NULL)
      umap_plot_gg(make_embedding_plot(seu, "umap"))
      tsne_plot_gg(make_embedding_plot(seu, "tsne"))
    })
    
    output$umap_plot <- plotly::renderPlotly({
      p <- umap_plot_gg()
      if (is.null(p)) return(NULL)
      plotly::ggplotly(p, tooltip = c("text", "x", "y", "color"))
    })
    
    output$tsne_plot <- plotly::renderPlotly({
      p <- tsne_plot_gg()
      if (is.null(p)) return(NULL)
      plotly::ggplotly(p, tooltip = c("text", "x", "y", "color"))
    })
    
    output$embedding_summary_table <- DT::renderDT({
      df <- embedding_summary()
      if (is.null(df)) {
        df <- data.frame(message = "Run UMAP/t-SNE to create embedding summary.", stringsAsFactors = FALSE)
      }
      DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE))
    })
    
    coords_df <- function(reduction) {
      seu <- embedding_seu()
      req(seu)
      emb <- tryCatch(Seurat::Embeddings(seu, reduction = reduction), error = function(e) NULL)
      req(emb)
      df <- as.data.frame(emb)
      df$cell <- rownames(df)
      df <- df[, c("cell", setdiff(colnames(df), "cell")), drop = FALSE]
      df
    }
    
    save_plot <- function(plot_reactive, file, device = "pdf", width = 7, height = 6) {
      p <- plot_reactive()
      req(p)
      ggplot2::ggsave(file, plot = p, device = device, width = width, height = height)
    }
    
    output$dl_umap_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_UMAP_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(umap_plot_gg, file, "pdf")
    )
    
    output$dl_umap_svg <- downloadHandler(
      filename = function() paste0("CoTRA_UMAP_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(umap_plot_gg, file, "svg")
    )
    
    output$dl_tsne_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_tSNE_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(tsne_plot_gg, file, "pdf")
    )
    
    output$dl_tsne_svg <- downloadHandler(
      filename = function() paste0("CoTRA_tSNE_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(tsne_plot_gg, file, "svg")
    )
    
    output$dl_umap_coords <- downloadHandler(
      filename = function() paste0("CoTRA_UMAP_coordinates_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) utils::write.csv(coords_df("umap"), file, row.names = FALSE)
    )
    
    output$dl_tsne_coords <- downloadHandler(
      filename = function() paste0("CoTRA_tSNE_coordinates_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) utils::write.csv(coords_df("tsne"), file, row.names = FALSE)
    )
    
    output$dl_embedding_summary <- downloadHandler(
      filename = function() paste0("CoTRA_embedding_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        df <- embedding_summary()
        req(df)
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    
    build_embedding_project <- function(seu_out, out_file = NULL) {
      list(
        seurat = seu_out,
        project_meta = list(
          project_type = "CoTRA_scRNA_project",
          step = "tSNE_UMAP",
          saved_at = Sys.time(),
          cotra_module = "modules/sc/tsne_umap.R",
          output_path = out_file
        ),
        cells_meta = seu_out@meta.data,
        step = "tSNE_UMAP",
        saved_at = Sys.time(),
        genes_list = list(
          selected_genes = tryCatch(seu_out@misc$CoTRA_gene_selection$selected_genes, error = function(e) character(0)),
          hvgs = tryCatch(seu_out@misc$CoTRA_gene_selection$hvg_genes, error = function(e) character(0)),
          pca_features = tryCatch(seu_out@misc$CoTRA_pca$features, error = function(e) character(0))
        ),
        parameters = list(
          qc = tryCatch(seu_out@misc$CoTRA_qc, error = function(e) NULL),
          gene_selection = tryCatch(seu_out@misc$CoTRA_gene_selection, error = function(e) NULL),
          pca = tryCatch(seu_out@misc$CoTRA_pca, error = function(e) NULL),
          tsne_umap = tryCatch(seu_out@misc$CoTRA_tsne_umap, error = function(e) NULL)
        ),
        tables = list(
          embedding_summary = embedding_summary(),
          umap_coordinates = tryCatch(coords_df("umap"), error = function(e) NULL),
          tsne_coordinates = tryCatch(coords_df("tsne"), error = function(e) NULL)
        ),
        cotra_version = if (exists("cotra_version")) cotra_version else "unknown"
      )
    }
    
    output$download_embedding_session <- downloadHandler(
      filename = function() {
        paste0("CoTRA_scRNA_UMAP_TSNE_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
      },
      content = function(file) {
        seu_out <- embedding_seu()
        req(seu_out)
        
        out_dir <- file.path("outputs", "scRNA", "sessioninfo")
        if (!dir.exists(out_dir)) {
          dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        }
        
        out_file <- file.path(
          out_dir,
          paste0("CoTRA_scRNA_UMAP_TSNE_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
        )
        
        project_obj <- build_embedding_project(seu_out, out_file)
        
        saveRDS(project_obj, out_file)
        file.copy(out_file, file, overwrite = TRUE)
        last_session_path(out_file)
      }
    )
    
    return(list(
      seurat = reactive({
        if (!is.null(embedding_seu())) embedding_seu() else get_seu()
      }),
      embedding_ready = reactive(embedding_ready()),
      embedding_summary = reactive(embedding_summary()),
      session_path = reactive(last_session_path())
    ))
  })
}
