# ==========================================================
# modules/sc/clustering.R
# CoTRA scRNA-seq Clustering Module
#
# Features:
# - FindNeighbors using PCA dimensions
# - FindClusters with resolution control
# - Uses recommended PCs from PCA when available
# - UMAP/t-SNE plots colored by clusters
# - Cluster size table
# - Cluster composition by metadata
# - Cluster renaming
# - Cluster merging
# - Download cluster tables and plots
# - Save clustering session RDS to outputs/scRNA/sessioninfo/
# - Help and interpretation sections
# - Downstream-safe return object for markers
# ==========================================================

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

mod_sc_clustering_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("Clustering"),
    br(),
    
    bs4Dash::bs4Card(
      title = "Help: how to use clustering",
      width = 12,
      solidHeader = TRUE,
      status = "info",
      collapsible = TRUE,
      collapsed = TRUE,
      icon = icon("circle-info"),
      tags$h4("Purpose"),
      tags$p("This module groups cells with similar expression profiles using PCA space. Clusters can then be inspected with UMAP/t-SNE, marker genes, and metadata composition."),
      tags$h4("Recommended workflow"),
      tags$ol(
        tags$li("Run QC, gene selection, PCA, and UMAP/t-SNE first."),
        tags$li("Start with the recommended number of PCs from PCA."),
        tags$li("Run clustering with a moderate resolution such as 0.4 to 0.8."),
        tags$li("Inspect cluster sizes and composition by sample, condition, or batch."),
        tags$li("Rename or merge clusters only after checking marker genes and metadata."),
        tags$li("Continue to marker detection after clustering.")
      ),
      tags$h4("Main controls"),
      tags$ul(
        tags$li(tags$b("PC dimensions:"), " controls which PCA components are used to build the neighbor graph."),
        tags$li(tags$b("Resolution:"), " higher values usually create more clusters. Lower values create fewer clusters."),
        tags$li(tags$b("Algorithm:"), " controls the graph clustering method used by Seurat."),
        tags$li(tags$b("Cluster renaming:"), " lets you assign biological labels or clearer names."),
        tags$li(tags$b("Cluster merging:"), " lets you combine clusters that likely represent the same cell population.")
      ),
      tags$h4("Important caution"),
      tags$p("Clustering is not cell type annotation. Use marker genes, known biology, and metadata before naming or merging clusters.")
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "Clustering settings",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        fluidRow(
          column(3, numericInput(ns("dims_start"), "First PC", value = 1, min = 1, max = 100, step = 1)),
          column(3, numericInput(ns("dims_end"), "Last PC", value = 30, min = 2, max = 100, step = 1)),
          column(3, br(), actionButton(ns("use_recommended_pcs"), "Use PCA recommendation", class = "btn-secondary")),
          column(3, numericInput(ns("resolution"), "Resolution", value = 0.5, min = 0.01, max = 5, step = 0.05))
        ),
        fluidRow(
          column(
            3,
            selectInput(
              ns("algorithm"),
              "Clustering algorithm",
              choices = c(
                "Louvain" = 1,
                "Louvain with multilevel refinement" = 2,
                "SLM" = 3,
                "Leiden" = 4
              ),
              selected = 1
            )
          ),
          column(3, numericInput(ns("k_param"), "k.param for FindNeighbors", value = 20, min = 5, max = 100, step = 5)),
          column(3, numericInput(ns("random_seed"), "Random seed", value = 1234, min = 1, max = 999999, step = 1)),
          column(3, br(), actionButton(ns("run_clustering"), "Run clustering", class = "btn-primary"))
        ),
        hr(),
        fluidRow(
          column(4, downloadButton(ns("download_clustering_session"), "Save clustering session (.rds)", class = "btn-success")),
          column(8, uiOutput(ns("session_save_hint")))
        ),
        br(),
        uiOutput(ns("clustering_status_box"))
      )
    ),
    
    bs4Dash::bs4Card(
      title = "Interpretation: clustering settings",
      width = 12,
      solidHeader = TRUE,
      status = "info",
      collapsible = TRUE,
      collapsed = TRUE,
      icon = icon("circle-info"),
      tags$h4("Resolution"),
      tags$ul(
        tags$li("Low resolution gives fewer, broader clusters."),
        tags$li("High resolution gives more, smaller clusters."),
        tags$li("Very high resolution can split one biological cell type into many subclusters."),
        tags$li("Very low resolution can merge distinct cell types into one cluster.")
      ),
      tags$h4("PC choice"),
      tags$ul(
        tags$li("Use enough PCs to capture meaningful biological structure."),
        tags$li("Avoid adding PCs dominated by technical effects such as sequencing depth or mitochondrial percentage."),
        tags$li("Use the PCA recommendation as a starting point, then inspect marker genes.")
      ),
      tags$h4("Renaming and merging"),
      tags$p("Rename clusters when marker genes support a biological identity. Merge clusters only when they show similar markers, similar metadata composition, and no clear biological separation.")
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "Cluster visualization",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        fluidRow(
          column(4, selectInput(ns("plot_reduction"), "Reduction", choices = c("UMAP" = "umap", "t-SNE" = "tsne"), selected = "umap")),
          column(4, selectInput(ns("cluster_label_column"), "Cluster label", choices = c("seurat_clusters"), selected = "seurat_clusters")),
          column(4, sliderInput(ns("point_size"), "Point size", min = 0.1, max = 3, value = 0.6, step = 0.1))
        ),
        fluidRow(
          column(4, checkboxInput(ns("show_labels"), "Show cluster labels", TRUE)),
          column(4, checkboxInput(ns("use_custom_labels"), "Use renamed/merged labels if available", TRUE)),
          column(4, br(), actionButton(ns("refresh_cluster_plot"), "Refresh plot", class = "btn-secondary"))
        ),
        br(),
        shinycssloaders::withSpinner(plotly::plotlyOutput(ns("cluster_plot")), type = 4),
        br(),
        fluidRow(
          column(4, downloadButton(ns("dl_cluster_plot_pdf"), "Cluster plot PDF")),
          column(4, downloadButton(ns("dl_cluster_plot_svg"), "Cluster plot SVG")),
          column(4, downloadButton(ns("dl_cluster_coords"), "Cluster coordinates CSV"))
        )
      )
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "Cluster tables",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        tabsetPanel(
          tabPanel(
            "Cluster sizes",
            br(),
            DT::DTOutput(ns("cluster_size_table")),
            br(),
            downloadButton(ns("dl_cluster_sizes"), "Download cluster sizes CSV")
          ),
          tabPanel(
            "Composition by metadata",
            br(),
            fluidRow(
              column(6, uiOutput(ns("composition_metadata_ui"))),
              column(6, downloadButton(ns("dl_cluster_composition"), "Download composition CSV"))
            ),
            br(),
            DT::DTOutput(ns("cluster_composition_table")),
            br(),
            shinycssloaders::withSpinner(plotly::plotlyOutput(ns("cluster_composition_plot")), type = 4)
          ),
          tabPanel(
            "Cell assignments",
            br(),
            DT::DTOutput(ns("cluster_assignment_table")),
            br(),
            downloadButton(ns("dl_cluster_assignments"), "Download cell assignments CSV")
          )
        )
      )
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "Rename and merge clusters",
        width = 12,
        solidHeader = TRUE,
        status = "warning",
        collapsible = TRUE,
        collapsed = TRUE,
        icon = icon("pen-to-square"),
        tags$h4("Rename clusters"),
        fluidRow(
          column(4, uiOutput(ns("rename_cluster_ui"))),
          column(4, textInput(ns("new_cluster_name"), "New cluster name", value = "")),
          column(4, br(), actionButton(ns("apply_cluster_rename"), "Apply rename", class = "btn-warning"))
        ),
        hr(),
        tags$h4("Merge clusters"),
        fluidRow(
          column(4, uiOutput(ns("merge_clusters_ui"))),
          column(4, textInput(ns("merged_cluster_name"), "Merged cluster name", value = "Merged_cluster")),
          column(4, br(), actionButton(ns("apply_cluster_merge"), "Merge selected clusters", class = "btn-danger"))
        ),
        hr(),
        tags$h4("Current cluster label map"),
        DT::DTOutput(ns("cluster_label_map_table")),
        br(),
        downloadButton(ns("dl_cluster_label_map"), "Download label map CSV")
      )
    ),
    
    bs4Dash::bs4Card(
      title = "Interpretation: cluster results",
      width = 12,
      solidHeader = TRUE,
      status = "info",
      collapsible = TRUE,
      collapsed = TRUE,
      icon = icon("circle-info"),
      tags$h4("How to interpret cluster plots"),
      tags$ul(
        tags$li("Clusters separated on UMAP or t-SNE may represent different cell types, cell states, or technical effects."),
        tags$li("Cluster size alone does not prove biological importance."),
        tags$li("Very small clusters may represent rare cells, doublets, damaged cells, or overclustering."),
        tags$li("Large clusters may contain multiple biological states if resolution is too low.")
      ),
      tags$h4("How to use composition tables"),
      tags$ul(
        tags$li("A cluster dominated by one sample or batch may reflect technical variation."),
        tags$li("A cluster enriched in one condition may reflect biology, but study design and batch must be checked."),
        tags$li("Use marker genes before assigning final cell type names.")
      )
    )
  )
}

mod_sc_clustering_server <- function(id, seurat_r, sc_state = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    clustered_seu <- reactiveVal(NULL)
    clustering_ready <- reactiveVal(FALSE)
    cluster_plot_gg <- reactiveVal(NULL)
    cluster_size_df <- reactiveVal(NULL)
    cluster_assignment_df <- reactiveVal(NULL)
    cluster_label_map <- reactiveVal(data.frame(original_cluster = character(0), current_label = character(0), action = character(0), stringsAsFactors = FALSE))
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
    
    has_reduction <- function(seu, reduction) {
      reduction %in% names(seu@reductions)
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
      if (is.null(seu) || !has_pca(seu)) return(NULL)
      rec <- get_recommended_pcs(seu)
      isolate({
        updateNumericInput(session, "dims_start", value = 1, max = max(1, pca_npcs(seu)))
        updateNumericInput(session, "dims_end", value = rec, max = max(2, pca_npcs(seu)))
      })
    })
    
    observeEvent(input$use_recommended_pcs, {
      seu <- get_seu()
      req(seu)
      if (!has_pca(seu)) {
        showModal(modalDialog(title = "PCA missing", "Run PCA before clustering.", easyClose = TRUE))
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
    
    update_cluster_label_map <- function(seu) {
      if (!"seurat_clusters" %in% colnames(seu@meta.data)) return(NULL)
      cl <- sort(unique(as.character(seu@meta.data$seurat_clusters)))
      old <- cluster_label_map()
      out <- data.frame(
        original_cluster = cl,
        current_label = cl,
        action = "original",
        stringsAsFactors = FALSE
      )
      if (!is.null(old) && nrow(old) > 0) {
        for (i in seq_len(nrow(out))) {
          hit <- which(old$original_cluster == out$original_cluster[i])
          if (length(hit) > 0) {
            out$current_label[i] <- old$current_label[hit[1]]
            out$action[i] <- old$action[hit[1]]
          }
        }
      }
      cluster_label_map(out)
    }
    
    apply_label_map_to_seurat <- function(seu) {
      map <- cluster_label_map()
      if (is.null(map) || nrow(map) == 0 || !"seurat_clusters" %in% colnames(seu@meta.data)) return(seu)
      cl <- as.character(seu@meta.data$seurat_clusters)
      named <- stats::setNames(map$current_label, map$original_cluster)
      labels <- unname(named[cl])
      labels[is.na(labels)] <- cl[is.na(labels)]
      seu@meta.data$CoTRA_cluster_label <- factor(labels)
      seu
    }
    
    make_cluster_tables <- function(seu) {
      if (!"seurat_clusters" %in% colnames(seu@meta.data)) return(NULL)
      
      md <- seu@meta.data
      cluster_label <- if ("CoTRA_cluster_label" %in% colnames(md)) as.character(md$CoTRA_cluster_label) else as.character(md$seurat_clusters)
      
      size_df <- as.data.frame(table(cluster = cluster_label), stringsAsFactors = FALSE)
      colnames(size_df) <- c("cluster", "cells")
      size_df$percent <- round(size_df$cells / sum(size_df$cells) * 100, 3)
      size_df <- size_df[order(size_df$cluster), , drop = FALSE]
      
      assign_df <- data.frame(
        cell = rownames(md),
        seurat_cluster = as.character(md$seurat_clusters),
        cluster_label = cluster_label,
        stringsAsFactors = FALSE
      )
      
      cluster_size_df(size_df)
      cluster_assignment_df(assign_df)
      NULL
    }
    
    output$clustering_status_box <- renderUI({
      base <- get_seu()
      seu <- clustered_seu()
      
      if (is.null(base)) {
        return(tags$div(style = "color:#777;", "No Seurat object available."))
      }
      if (!has_pca(base)) {
        return(tags$div(style = "color:#B00020;", "PCA is missing. Run PCA before clustering."))
      }
      if (is.null(seu) || !isTRUE(clustering_ready())) {
        return(tags$div(style = "color:#777;", paste0("PCA found with ", pca_npcs(base), " PCs. Run clustering to create cluster labels.")))
      }
      
      info <- tryCatch(seu@misc$CoTRA_clustering, error = function(e) NULL)
      n_clusters <- length(unique(as.character(seu@meta.data$seurat_clusters)))
      
      tags$div(
        tags$b("Clustering status: "), "completed", tags$br(),
        tags$b("Clusters: "), n_clusters, tags$br(),
        tags$b("Resolution: "), if (!is.null(info)) info$resolution else "NA", tags$br(),
        tags$b("PCs used: "), if (!is.null(info)) paste(info$dims, collapse = ", ") else "NA"
      )
    })
    
    output$session_save_hint <- renderUI({
      if (!isTRUE(clustering_ready()) || is.null(clustered_seu())) {
        return(tags$span(style = "color:#777;", "Run clustering before saving this step."))
      }
      if (is.null(last_session_path())) {
        return(tags$span(style = "color:#777;", "Save the current clustering state as a CoTRA project RDS."))
      }
      tags$span(style = "color:#2E7D32;", paste("Last saved:", last_session_path()))
    })
    
    observeEvent(input$run_clustering, {
      seu <- get_seu()
      req(seu)
      
      dim_check <- validate_dims(seu)
      if (!isTRUE(dim_check$ok)) {
        showModal(modalDialog(title = "PC check failed", dim_check$message, easyClose = TRUE))
        return(NULL)
      }
      
      dims <- dim_check$dims
      resolution <- as.numeric(input$resolution %||% 0.5)
      k_param <- as.integer(input$k_param %||% 20)
      k_param <- max(2, min(k_param, ncol(seu) - 1))
      algorithm <- as.integer(input$algorithm %||% 1)
      
      withProgress(message = "Running clustering", value = 0.1, {
        set.seed(as.integer(input$random_seed %||% 1234))
        
        incProgress(0.35, detail = "Finding neighbors")
        seu <- tryCatch(
          Seurat::FindNeighbors(
            seu,
            reduction = "pca",
            dims = dims,
            k.param = k_param,
            verbose = FALSE
          ),
          error = function(e) {
            showModal(modalDialog(title = "FindNeighbors failed", conditionMessage(e), easyClose = TRUE))
            NULL
          }
        )
        if (is.null(seu)) return(NULL)
        
        incProgress(0.7, detail = "Finding clusters")
        seu <- tryCatch(
          Seurat::FindClusters(
            seu,
            resolution = resolution,
            algorithm = algorithm,
            random.seed = as.integer(input$random_seed %||% 1234),
            verbose = FALSE
          ),
          error = function(e) {
            showModal(modalDialog(title = "FindClusters failed", conditionMessage(e), easyClose = TRUE))
            NULL
          }
        )
        if (is.null(seu)) return(NULL)
        
        if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
          seu@meta.data$seurat_clusters <- as.character(Seurat::Idents(seu))
        }
        
        update_cluster_label_map(seu)
        seu <- apply_label_map_to_seurat(seu)
        Seurat::Idents(seu) <- seu@meta.data$CoTRA_cluster_label
        
        seu@misc$CoTRA_clustering <- list(
          dims = dims,
          resolution = resolution,
          algorithm = algorithm,
          k_param = k_param,
          random_seed = as.integer(input$random_seed %||% 1234),
          n_clusters = length(unique(as.character(seu@meta.data$seurat_clusters))),
          label_map = cluster_label_map(),
          created_at = Sys.time()
        )
        
        make_cluster_tables(seu)
        clustered_seu(seu)
        clustering_ready(TRUE)
        cluster_plot_gg(make_cluster_plot(seu))
        
        if (!is.null(sc_state)) {
          sc_state$seurat <- seu
          sc_state$parameters$clustering <- seu@misc$CoTRA_clustering
        }
        
        incProgress(1)
      })
    })
    
    output$cluster_label_column <- renderUI({
      NULL
    })
    
    observe({
      seu <- clustered_seu()
      if (is.null(seu)) return(NULL)
      choices <- c("seurat_clusters")
      if ("CoTRA_cluster_label" %in% colnames(seu@meta.data)) choices <- c(choices, "CoTRA_cluster_label")
      updateSelectInput(session, "cluster_label_column", choices = choices, selected = if ("CoTRA_cluster_label" %in% choices) "CoTRA_cluster_label" else "seurat_clusters")
    })
    
    build_reduction_df <- function(seu) {
      reduction <- input$plot_reduction %||% "umap"
      if (!has_reduction(seu, reduction)) return(NULL)
      
      emb <- tryCatch(Seurat::Embeddings(seu, reduction = reduction), error = function(e) NULL)
      if (is.null(emb) || ncol(emb) < 2) return(NULL)
      
      md <- seu@meta.data
      label_col <- input$cluster_label_column %||% "seurat_clusters"
      if (isTRUE(input$use_custom_labels) && "CoTRA_cluster_label" %in% colnames(md)) label_col <- "CoTRA_cluster_label"
      if (!label_col %in% colnames(md)) label_col <- "seurat_clusters"
      
      df <- as.data.frame(emb[, 1:2, drop = FALSE])
      colnames(df) <- c("Dim1", "Dim2")
      df$cell <- rownames(df)
      df$cluster <- as.character(md[rownames(df), label_col])
      df$label_col <- label_col
      df$reduction <- reduction
      df
    }
    
    make_cluster_plot <- function(seu) {
      df <- build_reduction_df(seu)
      if (is.null(df)) return(NULL)
      
      reduction_label <- if (unique(df$reduction)[1] == "umap") "UMAP" else "t-SNE"
      xlab <- if (unique(df$reduction)[1] == "umap") "UMAP_1" else "tSNE_1"
      ylab <- if (unique(df$reduction)[1] == "umap") "UMAP_2" else "tSNE_2"
      ps <- input$point_size %||% 0.6
      
      p <- ggplot2::ggplot(df, ggplot2::aes(x = Dim1, y = Dim2, color = cluster, text = cell)) +
        ggplot2::geom_point(size = ps, alpha = 0.85) +
        ggplot2::labs(title = paste0(reduction_label, " colored by cluster"), x = xlab, y = ylab, color = "Cluster") +
        safe_theme()
      
      if (isTRUE(input$show_labels)) {
        label_df <- stats::aggregate(cbind(Dim1, Dim2) ~ cluster, data = df, FUN = median)
        p <- p + ggplot2::geom_text(
          data = label_df,
          ggplot2::aes(x = Dim1, y = Dim2, label = cluster),
          inherit.aes = FALSE,
          size = 3,
          fontface = "bold"
        )
      }
      
      p
    }
    
    observeEvent(input$refresh_cluster_plot, {
      seu <- clustered_seu()
      req(seu)
      cluster_plot_gg(make_cluster_plot(seu))
    })
    
    observe({
      seu <- clustered_seu()
      if (is.null(seu)) return(NULL)
      cluster_plot_gg(make_cluster_plot(seu))
    })
    
    output$cluster_plot <- plotly::renderPlotly({
      p <- cluster_plot_gg()
      if (is.null(p)) return(NULL)
      plotly::ggplotly(p, tooltip = c("text", "x", "y", "color"))
    })
    
    output$cluster_size_table <- DT::renderDT({
      df <- cluster_size_df()
      if (is.null(df)) df <- data.frame(message = "Run clustering to create cluster sizes.", stringsAsFactors = FALSE)
      DT::datatable(df, options = list(pageLength = 20, scrollX = TRUE))
    })
    
    output$cluster_assignment_table <- DT::renderDT({
      df <- cluster_assignment_df()
      if (is.null(df)) df <- data.frame(message = "Run clustering to create cell assignments.", stringsAsFactors = FALSE)
      DT::datatable(df, options = list(pageLength = 20, scrollX = TRUE))
    })
    
    metadata_choices <- reactive({
      seu <- clustered_seu()
      if (is.null(seu)) return(character(0))
      md <- seu@meta.data
      cols <- colnames(md)[vapply(md, function(x) is.factor(x) || is.character(x) || is.logical(x), logical(1))]
      setdiff(cols, c("seurat_clusters", "CoTRA_cluster_label"))
    })
    
    output$composition_metadata_ui <- renderUI({
      choices <- metadata_choices()
      if (length(choices) == 0) return(tags$span("No categorical metadata columns found."))
      selectInput(ns("composition_metadata"), "Metadata column", choices = choices, selected = choices[1])
    })
    
    cluster_composition_df <- reactive({
      seu <- clustered_seu()
      req(seu)
      md_col <- input$composition_metadata
      if (is.null(md_col) || !md_col %in% colnames(seu@meta.data)) return(NULL)
      
      md <- seu@meta.data
      cluster_label <- if ("CoTRA_cluster_label" %in% colnames(md)) as.character(md$CoTRA_cluster_label) else as.character(md$seurat_clusters)
      comp <- as.data.frame(table(cluster = cluster_label, group = as.character(md[[md_col]])), stringsAsFactors = FALSE)
      colnames(comp) <- c("cluster", "group", "cells")
      totals <- stats::aggregate(cells ~ cluster, data = comp, sum)
      colnames(totals)[2] <- "cluster_total"
      comp <- merge(comp, totals, by = "cluster", all.x = TRUE)
      comp$percent_within_cluster <- round(comp$cells / comp$cluster_total * 100, 3)
      comp[order(comp$cluster, comp$group), , drop = FALSE]
    })
    
    output$cluster_composition_table <- DT::renderDT({
      df <- cluster_composition_df()
      if (is.null(df)) df <- data.frame(message = "Select metadata to view cluster composition.", stringsAsFactors = FALSE)
      DT::datatable(df, options = list(pageLength = 20, scrollX = TRUE))
    })
    
    output$cluster_composition_plot <- plotly::renderPlotly({
      df <- cluster_composition_df()
      if (is.null(df) || !"cells" %in% colnames(df)) return(NULL)
      p <- ggplot2::ggplot(df, ggplot2::aes(x = cluster, y = percent_within_cluster, fill = group)) +
        ggplot2::geom_col(position = "stack") +
        ggplot2::labs(x = "Cluster", y = "Percent within cluster", fill = input$composition_metadata) +
        safe_theme()
      plotly::ggplotly(p)
    })
    
    output$rename_cluster_ui <- renderUI({
      map <- cluster_label_map()
      if (is.null(map) || nrow(map) == 0) return(tags$span("Run clustering first."))
      selectInput(ns("cluster_to_rename"), "Cluster to rename", choices = map$current_label, selected = map$current_label[1])
    })
    
    output$merge_clusters_ui <- renderUI({
      map <- cluster_label_map()
      if (is.null(map) || nrow(map) == 0) return(tags$span("Run clustering first."))
      selectizeInput(ns("clusters_to_merge"), "Clusters to merge", choices = unique(map$current_label), selected = NULL, multiple = TRUE)
    })
    
    observeEvent(input$apply_cluster_rename, {
      seu <- clustered_seu()
      req(seu)
      old_label <- input$cluster_to_rename
      new_label <- trimws(input$new_cluster_name %||% "")
      if (is.null(old_label) || !nzchar(new_label)) {
        showModal(modalDialog(title = "Rename failed", "Select a cluster and enter a new name.", easyClose = TRUE))
        return(NULL)
      }
      
      map <- cluster_label_map()
      map$current_label[map$current_label == old_label] <- new_label
      map$action[map$current_label == new_label] <- "renamed"
      cluster_label_map(map)
      
      seu <- apply_label_map_to_seurat(seu)
      Seurat::Idents(seu) <- seu@meta.data$CoTRA_cluster_label
      seu@misc$CoTRA_clustering$label_map <- map
      
      make_cluster_tables(seu)
      clustered_seu(seu)
      cluster_plot_gg(make_cluster_plot(seu))
      
      if (!is.null(sc_state)) sc_state$seurat <- seu
    })
    
    observeEvent(input$apply_cluster_merge, {
      seu <- clustered_seu()
      req(seu)
      to_merge <- input$clusters_to_merge
      new_label <- trimws(input$merged_cluster_name %||% "")
      if (is.null(to_merge) || length(to_merge) < 2 || !nzchar(new_label)) {
        showModal(modalDialog(title = "Merge failed", "Select at least two clusters and enter a merged cluster name.", easyClose = TRUE))
        return(NULL)
      }
      
      map <- cluster_label_map()
      map$current_label[map$current_label %in% to_merge] <- new_label
      map$action[map$current_label == new_label] <- "merged"
      cluster_label_map(map)
      
      seu <- apply_label_map_to_seurat(seu)
      Seurat::Idents(seu) <- seu@meta.data$CoTRA_cluster_label
      seu@misc$CoTRA_clustering$label_map <- map
      
      make_cluster_tables(seu)
      clustered_seu(seu)
      cluster_plot_gg(make_cluster_plot(seu))
      
      if (!is.null(sc_state)) sc_state$seurat <- seu
    })
    
    output$cluster_label_map_table <- DT::renderDT({
      df <- cluster_label_map()
      if (is.null(df) || nrow(df) == 0) df <- data.frame(message = "Run clustering to create label map.", stringsAsFactors = FALSE)
      DT::datatable(df, options = list(pageLength = 20, scrollX = TRUE))
    })
    
    coords_df <- function() {
      seu <- clustered_seu()
      req(seu)
      df <- build_reduction_df(seu)
      req(df)
      df[, c("cell", "Dim1", "Dim2", "cluster", "reduction"), drop = FALSE]
    }
    
    save_plot <- function(plot_reactive, file, device = "pdf", width = 7, height = 6) {
      p <- plot_reactive()
      req(p)
      ggplot2::ggsave(file, plot = p, device = device, width = width, height = height)
    }
    
    output$dl_cluster_plot_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_cluster_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(cluster_plot_gg, file, "pdf")
    )
    
    output$dl_cluster_plot_svg <- downloadHandler(
      filename = function() paste0("CoTRA_cluster_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(cluster_plot_gg, file, "svg")
    )
    
    output$dl_cluster_coords <- downloadHandler(
      filename = function() paste0("CoTRA_cluster_coordinates_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) utils::write.csv(coords_df(), file, row.names = FALSE)
    )
    
    output$dl_cluster_sizes <- downloadHandler(
      filename = function() paste0("CoTRA_cluster_sizes_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        df <- cluster_size_df()
        req(df)
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    
    output$dl_cluster_assignments <- downloadHandler(
      filename = function() paste0("CoTRA_cluster_assignments_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        df <- cluster_assignment_df()
        req(df)
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    
    output$dl_cluster_composition <- downloadHandler(
      filename = function() paste0("CoTRA_cluster_composition_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        df <- cluster_composition_df()
        req(df)
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    
    output$dl_cluster_label_map <- downloadHandler(
      filename = function() paste0("CoTRA_cluster_label_map_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        df <- cluster_label_map()
        req(df)
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    
    build_clustering_project <- function(seu_out, out_file = NULL) {
      list(
        seurat = seu_out,
        project_meta = list(
          project_type = "CoTRA_scRNA_project",
          step = "Clustering",
          saved_at = Sys.time(),
          cotra_module = "modules/sc/clustering.R",
          output_path = out_file
        ),
        cells_meta = seu_out@meta.data,
        step = "Clustering",
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
          tsne_umap = tryCatch(seu_out@misc$CoTRA_tsne_umap, error = function(e) NULL),
          clustering = tryCatch(seu_out@misc$CoTRA_clustering, error = function(e) NULL)
        ),
        tables = list(
          cluster_sizes = cluster_size_df(),
          cluster_assignments = cluster_assignment_df(),
          cluster_label_map = cluster_label_map(),
          cluster_composition = tryCatch(cluster_composition_df(), error = function(e) NULL)
        ),
        cotra_version = if (exists("cotra_version")) cotra_version else "unknown"
      )
    }
    
    output$download_clustering_session <- downloadHandler(
      filename = function() {
        paste0("CoTRA_scRNA_Clustering_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
      },
      content = function(file) {
        seu_out <- clustered_seu()
        req(seu_out)
        
        out_dir <- file.path("outputs", "scRNA", "sessioninfo")
        if (!dir.exists(out_dir)) {
          dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        }
        
        out_file <- file.path(
          out_dir,
          paste0("CoTRA_scRNA_Clustering_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
        )
        
        project_obj <- build_clustering_project(seu_out, out_file)
        saveRDS(project_obj, out_file)
        file.copy(out_file, file, overwrite = TRUE)
        last_session_path(out_file)
      }
    )
    
    return(list(
      seurat = reactive({
        if (!is.null(clustered_seu())) clustered_seu() else get_seu()
      }),
      clustering_ready = reactive(clustering_ready()),
      cluster_sizes = reactive(cluster_size_df()),
      cluster_assignments = reactive(cluster_assignment_df()),
      cluster_label_map = reactive(cluster_label_map()),
      session_path = reactive(last_session_path())
    ))
  })
}
