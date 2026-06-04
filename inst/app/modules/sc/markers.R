# ==========================================================
# modules/sc/markers.R
# CoTRA scRNA-seq Marker Gene Module
#
# Features:
# - FindAllMarkers for all clusters
# - FindMarkers for selected cluster comparisons
# - Uses renamed or merged cluster labels when available
# - Cluster marker table
# - Top marker heatmap
# - DotPlot for top markers
# - FeaturePlot for selected markers
# - Volcano plot for pairwise cluster comparison
# - Download marker tables and plots
# - Save marker session RDS to outputs/scRNA/sessioninfo/
# - Help and interpretation sections
# - Downstream-safe return for annotation/report/enrichment
#
# Suggested improvements included:
# - Marker source selection: RNA data or SCT data
# - Pairwise comparison mode
# - Minimum cell check per cluster
# - Gene availability check for plots
# - Conservative defaults for logFC and adjusted p-value
# ==========================================================

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

mod_sc_markers_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("Marker genes"),
    br(),
    
    bs4Dash::bs4Card(
      title = "Help: how to use marker analysis",
      width = 12,
      solidHeader = TRUE,
      status = "info",
      collapsible = TRUE,
      collapsed = TRUE,
      icon = icon("circle-info"),
      tags$h4("Purpose"),
      tags$p("This module identifies genes enriched in clusters or between selected cluster groups. These markers help interpret cluster identity, biological state, and technical effects."),
      tags$h4("Recommended workflow"),
      tags$ol(
        tags$li("Run clustering first."),
        tags$li("Use renamed or merged cluster labels if you curated clusters."),
        tags$li("Run FindAllMarkers to identify markers for every cluster."),
        tags$li("Use pairwise comparison to compare specific clusters."),
        tags$li("Inspect top marker heatmap, DotPlot, FeaturePlot, and volcano plot."),
        tags$li("Use known marker databases or literature before assigning final cell type names.")
      ),
      tags$h4("Main outputs"),
      tags$ul(
        tags$li(tags$b("All markers:"), " marker genes for every cluster compared against all other cells."),
        tags$li(tags$b("Pairwise markers:"), " differential markers between two selected clusters or groups."),
        tags$li(tags$b("Heatmap:"), " expression pattern of top markers across cells."),
        tags$li(tags$b("DotPlot:"), " average expression and percent expressed across clusters."),
        tags$li(tags$b("FeaturePlot:"), " spatial expression of selected genes on UMAP or t-SNE."),
        tags$li(tags$b("Volcano plot:"), " pairwise fold-change and statistical significance.")
      ),
      tags$h4("Important caution"),
      tags$p("Marker genes suggest cluster identity. They do not prove cell identity alone. Always combine marker genes with prior biology, sample metadata, and QC checks.")
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "Marker analysis settings",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        fluidRow(
          column(3, uiOutput(ns("assay_ui"))),
          column(3, uiOutput(ns("identity_ui"))),
          column(3, selectInput(ns("test_use"), "Statistical test", choices = c("wilcox", "MAST", "LR", "DESeq2"), selected = "wilcox")),
          column(3, numericInput(ns("min_cells_per_cluster"), "Minimum cells per cluster", value = 10, min = 3, max = 100, step = 1))
        ),
        fluidRow(
          column(3, numericInput(ns("logfc_threshold"), "logFC threshold", value = 0.25, min = 0, max = 5, step = 0.05)),
          column(3, numericInput(ns("min_pct"), "Minimum pct", value = 0.1, min = 0, max = 1, step = 0.05)),
          column(3, numericInput(ns("p_adj_cutoff"), "Adjusted p-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
          column(3, checkboxInput(ns("only_pos"), "Only positive markers", TRUE))
        ),
        fluidRow(
          column(3, numericInput(ns("max_cells_per_ident"), "Max cells per identity", value = 500, min = 50, max = 10000, step = 50)),
          column(3, numericInput(ns("random_seed"), "Random seed", value = 1234, min = 1, max = 999999, step = 1)),
          column(3, checkboxInput(ns("use_latent_vars"), "Use latent vars if available", FALSE)),
          column(3, br(), actionButton(ns("run_all_markers"), "Run all markers", class = "btn-primary"))
        ),
        hr(),
        fluidRow(
          column(4, downloadButton(ns("download_marker_session"), "Save marker session (.rds)", class = "btn-success")),
          column(8, uiOutput(ns("session_save_hint")))
        ),
        br(),
        uiOutput(ns("marker_status_box"))
      )
    ),
    
    bs4Dash::bs4Card(
      title = "Interpretation: marker settings",
      width = 12,
      solidHeader = TRUE,
      status = "info",
      collapsible = TRUE,
      collapsed = TRUE,
      icon = icon("circle-info"),
      tags$h4("How to choose thresholds"),
      tags$ul(
        tags$li(tags$b("logFC threshold:"), " higher values keep stronger markers but may miss subtle markers."),
        tags$li(tags$b("minimum pct:"), " removes genes expressed in very few cells."),
        tags$li(tags$b("adjusted p-value:"), " controls statistical significance after multiple testing."),
        tags$li(tags$b("only positive markers:"), " useful for identifying genes enriched in each cluster."),
        tags$li(tags$b("max cells per identity:"), " speeds up analysis by downsampling large clusters.")
      ),
      tags$h4("Recommended defaults"),
      tags$p("Start with Wilcoxon test, logFC threshold 0.25, minimum pct 0.1, and adjusted p-value 0.05. Then adjust based on dataset size and biological question.")
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "All cluster markers",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        uiOutput(ns("all_marker_summary")),
        br(),
        DT::DTOutput(ns("all_marker_table")),
        br(),
        fluidRow(
          column(4, downloadButton(ns("dl_all_markers"), "Download all markers CSV")),
          column(4, downloadButton(ns("dl_top_markers"), "Download top markers CSV")),
          column(4, downloadButton(ns("dl_marker_summary"), "Download marker summary CSV"))
        )
      )
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "Pairwise marker comparison",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        fluidRow(
          column(4, uiOutput(ns("pairwise_ident_1_ui"))),
          column(4, uiOutput(ns("pairwise_ident_2_ui"))),
          column(4, br(), actionButton(ns("run_pairwise_markers"), "Run pairwise comparison", class = "btn-primary"))
        ),
        br(),
        DT::DTOutput(ns("pairwise_marker_table")),
        br(),
        fluidRow(
          column(6, downloadButton(ns("dl_pairwise_markers"), "Download pairwise markers CSV")),
          column(6, downloadButton(ns("dl_volcano_pdf"), "Volcano PDF"))
        )
      )
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "Marker visualizations",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        tabsetPanel(
          tabPanel(
            "Top marker heatmap",
            br(),
            fluidRow(
              column(3, numericInput(ns("top_n_heatmap"), "Top markers per cluster", value = 5, min = 1, max = 30, step = 1)),
              column(3, numericInput(ns("heatmap_max_cells"), "Max cells to plot", value = 1000, min = 100, max = 10000, step = 100)),
              column(3, actionButton(ns("make_heatmap"), "Create heatmap", class = "btn-secondary")),
              column(3, tagList(
                downloadButton(ns("dl_heatmap_pdf"), "Heatmap PDF"),
                br(), br(),
                downloadButton(ns("dl_heatmap_svg"), "Heatmap SVG")
              ))
            ),
            br(),
            shinycssloaders::withSpinner(plotly::plotlyOutput(ns("marker_heatmap"), height = "750px"), type = 4),
            br(),
            bs4Dash::bs4Card(
              title = "Interpretation: marker heatmap",
              width = 12,
              solidHeader = TRUE,
              status = "info",
              collapsible = TRUE,
              collapsed = TRUE,
              icon = icon("circle-info"),
              tags$p("The heatmap shows expression of top markers across cells. Good markers show strong expression in one cluster and lower expression elsewhere."),
              tags$ul(
                tags$li("Clear blocks suggest cluster-specific gene programs."),
                tags$li("Shared blocks suggest related clusters or possible overclustering."),
                tags$li("Stress, mitochondrial, or ribosomal genes among top markers may indicate technical effects.")
              )
            )
          ),
          tabPanel(
            "DotPlot",
            br(),
            fluidRow(
              column(4, textAreaInput(ns("dotplot_genes"), "Genes for DotPlot", value = "", rows = 5, placeholder = "Leave empty to use top markers")),
              column(3, numericInput(ns("top_n_dotplot"), "Top markers per cluster if empty", value = 5, min = 1, max = 20, step = 1)),
              column(2, br(), actionButton(ns("make_dotplot"), "Create DotPlot", class = "btn-secondary")),
              column(3, br(), tagList(
                downloadButton(ns("dl_dotplot_pdf"), "DotPlot PDF"),
                br(), br(),
                downloadButton(ns("dl_dotplot_svg"), "DotPlot SVG")
              ))
            ),
            br(),
            shinycssloaders::withSpinner(plotly::plotlyOutput(ns("marker_dotplot"), height = "650px"), type = 4),
            br(),
            bs4Dash::bs4Card(
              title = "Interpretation: DotPlot",
              width = 12,
              solidHeader = TRUE,
              status = "info",
              collapsible = TRUE,
              collapsed = TRUE,
              icon = icon("circle-info"),
              tags$p("DotPlot summarizes marker expression across clusters. Dot size shows the percentage of cells expressing the gene. Color shows average expression."),
              tags$ul(
                tags$li("Large and dark dots in one cluster suggest a strong marker."),
                tags$li("Markers expressed across many clusters are less specific."),
                tags$li("Use DotPlot to compare known marker panels across clusters.")
              )
            )
          ),
          tabPanel(
            "FeaturePlot",
            br(),
            fluidRow(
              column(4, textAreaInput(ns("featureplot_genes"), "Genes for FeaturePlot", value = "", rows = 5, placeholder = "Example:\nRHO\nPDE6B")),
              column(3, selectInput(ns("feature_reduction"), "Reduction", choices = c("UMAP" = "umap", "t-SNE" = "tsne"), selected = "umap")),
              column(2, br(), actionButton(ns("make_featureplot"), "Create FeaturePlot", class = "btn-secondary")),
              column(3, br(), tagList(
                downloadButton(ns("dl_featureplot_pdf"), "FeaturePlot PDF"),
                br(), br(),
                downloadButton(ns("dl_featureplot_svg"), "FeaturePlot SVG")
              ))
            ),
            uiOutput(ns("featureplot_missing_genes")),
            br(),
            shinycssloaders::withSpinner(plotly::plotlyOutput(ns("marker_featureplot"), height = "700px"), type = 4),
            br(),
            bs4Dash::bs4Card(
              title = "Interpretation: FeaturePlot",
              width = 12,
              solidHeader = TRUE,
              status = "info",
              collapsible = TRUE,
              collapsed = TRUE,
              icon = icon("circle-info"),
              tags$p("FeaturePlot maps gene expression onto UMAP or t-SNE. It helps show where marker genes are expressed spatially in the embedding."),
              tags$ul(
                tags$li("Cluster-localized expression supports cluster identity."),
                tags$li("Diffuse expression across many regions suggests a broad or non-specific marker."),
                tags$li("FeaturePlot should be interpreted with marker tables and DotPlot.")
              )
            )
          ),
          tabPanel(
            "Volcano plot",
            br(),
            shinycssloaders::withSpinner(plotly::plotlyOutput(ns("volcano_plot")), type = 4),
            br(),
            downloadButton(ns("dl_volcano_svg"), "Volcano SVG"),
            br(),
            bs4Dash::bs4Card(
              title = "Interpretation: volcano plot",
              width = 12,
              solidHeader = TRUE,
              status = "info",
              collapsible = TRUE,
              collapsed = TRUE,
              icon = icon("circle-info"),
              tags$p("The volcano plot shows pairwise marker strength and significance. Genes far from zero on the x-axis and high on the y-axis are stronger candidates."),
              tags$ul(
                tags$li("Positive logFC genes are enriched in the first selected group."),
                tags$li("Negative logFC genes are enriched in the second selected group."),
                tags$li("Very significant genes with small logFC may be statistically robust but biologically subtle.")
              )
            )
          )
        )
      )
    )
  )
}

mod_sc_markers_server <- function(id, seurat_r, sc_state = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    markers_seu <- reactiveVal(NULL)
    all_markers <- reactiveVal(NULL)
    pairwise_markers <- reactiveVal(NULL)
    marker_summary <- reactiveVal(NULL)
    heatmap_plot_gg <- reactiveVal(NULL)
    dotplot_gg <- reactiveVal(NULL)
    featureplot_gg <- reactiveVal(NULL)
    volcano_plot_gg <- reactiveVal(NULL)
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
    
    parse_gene_text <- function(x) {
      if (is.null(x) || !nzchar(x)) return(character(0))
      vals <- unlist(strsplit(x, "[\\s,;]+"))
      vals <- trimws(vals)
      vals <- vals[nzchar(vals)]
      unique(vals)
    }
    
    assay_choices <- reactive({
      seu <- get_seu()
      if (is.null(seu)) return(c("RNA"))
      assays <- tryCatch(Seurat::Assays(seu), error = function(e) "RNA")
      if (length(assays) == 0) assays <- "RNA"
      assays
    })
    
    output$assay_ui <- renderUI({
      choices <- assay_choices()
      selected <- if ("SCT" %in% choices) "SCT" else if ("RNA" %in% choices) "RNA" else choices[1]
      selectInput(ns("assay_use"), "Assay", choices = choices, selected = selected)
    })
    
    identity_choices <- reactive({
      seu <- get_seu()
      if (is.null(seu)) return(c("Active identity" = "ident"))
      md <- seu@meta.data
      choices <- c("Active identity" = "ident")
      if ("CoTRA_cluster_label" %in% colnames(md)) choices <- c(choices, "Renamed/merged cluster labels" = "CoTRA_cluster_label")
      if ("seurat_clusters" %in% colnames(md)) choices <- c(choices, "Seurat clusters" = "seurat_clusters")
      choices
    })
    
    output$identity_ui <- renderUI({
      choices <- identity_choices()
      selected <- if ("CoTRA_cluster_label" %in% choices) "CoTRA_cluster_label" else if ("seurat_clusters" %in% choices) "seurat_clusters" else "ident"
      selectInput(ns("identity_col"), "Identity for markers", choices = choices, selected = selected)
    })
    
    set_marker_idents <- function(seu) {
      id_col <- input$identity_col %||% "ident"
      if (id_col != "ident" && id_col %in% colnames(seu@meta.data)) {
        Seurat::Idents(seu) <- as.factor(seu@meta.data[[id_col]])
      }
      seu
    }
    
    valid_idents <- function(seu) {
      ids <- as.character(Seurat::Idents(seu))
      tab <- table(ids)
      keep <- names(tab)[tab >= as.integer(input$min_cells_per_cluster %||% 10)]
      sort(keep)
    }
    
    output$pairwise_ident_1_ui <- renderUI({
      seu <- get_seu()
      if (is.null(seu)) return(tags$span("No object available."))
      seu <- set_marker_idents(seu)
      ids <- valid_idents(seu)
      if (length(ids) == 0) return(tags$span("No identities with enough cells."))
      selectInput(ns("ident_1"), "Group 1", choices = ids, selected = ids[1])
    })
    
    output$pairwise_ident_2_ui <- renderUI({
      seu <- get_seu()
      if (is.null(seu)) return(tags$span("No object available."))
      seu <- set_marker_idents(seu)
      ids <- valid_idents(seu)
      if (length(ids) == 0) return(tags$span("No identities with enough cells."))
      selected <- if (length(ids) >= 2) ids[2] else ids[1]
      selectInput(ns("ident_2"), "Group 2", choices = ids, selected = selected)
    })
    
    output$marker_status_box <- renderUI({
      seu <- get_seu()
      if (is.null(seu)) {
        return(tags$div(style = "color:#777;", "No Seurat object available."))
      }
      if (!"seurat_clusters" %in% colnames(seu@meta.data) && !"CoTRA_cluster_label" %in% colnames(seu@meta.data)) {
        return(tags$div(style = "color:#B00020;", "No clustering labels found. Run clustering first."))
      }
      am <- all_markers()
      if (is.null(am)) {
        return(tags$div(style = "color:#777;", "Cluster labels found. Run all markers to start marker analysis."))
      }
      tags$div(
        tags$b("Marker status: "), "completed", tags$br(),
        tags$b("Total marker rows: "), nrow(am), tags$br(),
        tags$b("Clusters/groups: "), length(unique(am$cluster))
      )
    })
    
    output$session_save_hint <- renderUI({
      if (is.null(all_markers()) && is.null(pairwise_markers())) {
        return(tags$span(style = "color:#777;", "Run marker analysis before saving this step."))
      }
      if (is.null(last_session_path())) {
        return(tags$span(style = "color:#777;", "Save the current marker state as a CoTRA project RDS."))
      }
      tags$span(style = "color:#2E7D32;", paste("Last saved:", last_session_path()))
    })
    
    summarize_markers <- function(df) {
      if (is.null(df) || nrow(df) == 0) return(NULL)
      p_col <- if ("p_val_adj" %in% colnames(df)) "p_val_adj" else if ("p_val" %in% colnames(df)) "p_val" else NULL
      lfc_col <- if ("avg_log2FC" %in% colnames(df)) "avg_log2FC" else if ("avg_logFC" %in% colnames(df)) "avg_logFC" else NULL
      
      out <- as.data.frame(table(cluster = df$cluster), stringsAsFactors = FALSE)
      colnames(out) <- c("cluster", "marker_rows")
      
      if (!is.null(p_col)) {
        sig <- df[df[[p_col]] <= as.numeric(input$p_adj_cutoff %||% 0.05), , drop = FALSE]
        sig_tab <- as.data.frame(table(cluster = sig$cluster), stringsAsFactors = FALSE)
        colnames(sig_tab) <- c("cluster", "significant_markers")
        out <- merge(out, sig_tab, by = "cluster", all.x = TRUE)
        out$significant_markers[is.na(out$significant_markers)] <- 0
      }
      
      if (!is.null(lfc_col)) {
        top <- do.call(rbind, lapply(split(df, df$cluster), function(x) {
          x <- x[order(-x[[lfc_col]]), , drop = FALSE]
          data.frame(cluster = unique(x$cluster)[1], top_gene = x$gene[1], top_logFC = x[[lfc_col]][1], stringsAsFactors = FALSE)
        }))
        out <- merge(out, top, by = "cluster", all.x = TRUE)
      }
      
      out[order(out$cluster), , drop = FALSE]
    }
    
    top_markers_df <- reactive({
      df <- all_markers()
      if (is.null(df) || nrow(df) == 0) return(NULL)
      lfc_col <- if ("avg_log2FC" %in% colnames(df)) "avg_log2FC" else if ("avg_logFC" %in% colnames(df)) "avg_logFC" else NULL
      p_col <- if ("p_val_adj" %in% colnames(df)) "p_val_adj" else if ("p_val" %in% colnames(df)) "p_val" else NULL
      if (!is.null(p_col)) df <- df[df[[p_col]] <= as.numeric(input$p_adj_cutoff %||% 0.05), , drop = FALSE]
      if (!is.null(lfc_col)) df <- df[order(df$cluster, -df[[lfc_col]]), , drop = FALSE]
      do.call(rbind, lapply(split(df, df$cluster), function(x) head(x, 20)))
    })
    
    observeEvent(input$run_all_markers, {
      seu <- get_seu()
      req(seu)
      
      withProgress(message = "Finding all cluster markers", value = 0.1, {
        assay <- input$assay_use %||% Seurat::DefaultAssay(seu)
        if (!assay %in% Seurat::Assays(seu)) {
          showModal(modalDialog(title = "Assay not found", paste("Assay", assay, "was not found."), easyClose = TRUE))
          return(NULL)
        }
        Seurat::DefaultAssay(seu) <- assay
        seu <- set_marker_idents(seu)
        
        ids <- valid_idents(seu)
        if (length(ids) < 2) {
          showModal(modalDialog(title = "Not enough groups", "At least two identities with enough cells are required.", easyClose = TRUE))
          return(NULL)
        }
        
        if ((input$test_use %||% "wilcox") == "MAST" && !requireNamespace("MAST", quietly = TRUE)) {
          showModal(modalDialog(title = "MAST unavailable", "Install MAST or use Wilcoxon.", easyClose = TRUE))
          return(NULL)
        }
        
        latent_vars <- NULL
        if (isTRUE(input$use_latent_vars)) {
          candidates <- c("nCount_RNA", "percent.mito", "percent.ribo")
          latent_vars <- candidates[candidates %in% colnames(seu@meta.data)]
          if (length(latent_vars) == 0) latent_vars <- NULL
        }
        
        incProgress(0.5, detail = "Running FindAllMarkers")
        set.seed(as.integer(input$random_seed %||% 1234))
        
        markers <- tryCatch(
          Seurat::FindAllMarkers(
            seu,
            assay = assay,
            only.pos = isTRUE(input$only_pos),
            test.use = input$test_use %||% "wilcox",
            logfc.threshold = as.numeric(input$logfc_threshold %||% 0.25),
            min.pct = as.numeric(input$min_pct %||% 0.1),
            max.cells.per.ident = as.integer(input$max_cells_per_ident %||% 500),
            latent.vars = latent_vars,
            verbose = FALSE
          ),
          error = function(e) {
            showModal(modalDialog(title = "FindAllMarkers failed", conditionMessage(e), easyClose = TRUE))
            NULL
          }
        )
        
        if (is.null(markers)) return(NULL)
        if (nrow(markers) == 0) {
          showModal(modalDialog(title = "No markers found", "No markers passed the current thresholds. Try lower logFC or min.pct.", easyClose = TRUE))
        }
        
        markers$cluster <- as.character(markers$cluster)
        summary_df <- summarize_markers(markers)
        
        seu@misc$CoTRA_markers <- list(
          all_markers = markers,
          marker_summary = summary_df,
          assay = assay,
          identity_col = input$identity_col,
          test_use = input$test_use,
          only_pos = isTRUE(input$only_pos),
          logfc_threshold = as.numeric(input$logfc_threshold %||% 0.25),
          min_pct = as.numeric(input$min_pct %||% 0.1),
          p_adj_cutoff = as.numeric(input$p_adj_cutoff %||% 0.05),
          max_cells_per_ident = as.integer(input$max_cells_per_ident %||% 500),
          latent_vars = latent_vars,
          created_at = Sys.time()
        )
        
        markers_seu(seu)
        all_markers(markers)
        marker_summary(summary_df)
        
        if (!is.null(sc_state)) {
          sc_state$seurat <- seu
          sc_state$parameters$markers <- seu@misc$CoTRA_markers
        }
        
        incProgress(1)
      })
    })
    
    observeEvent(input$run_pairwise_markers, {
      seu <- get_seu()
      req(seu)
      assay <- input$assay_use %||% Seurat::DefaultAssay(seu)
      Seurat::DefaultAssay(seu) <- assay
      seu <- set_marker_idents(seu)
      
      ident1 <- input$ident_1
      ident2 <- input$ident_2
      if (is.null(ident1) || is.null(ident2) || ident1 == ident2) {
        showModal(modalDialog(title = "Invalid comparison", "Select two different groups.", easyClose = TRUE))
        return(NULL)
      }
      
      withProgress(message = "Finding pairwise markers", value = 0.2, {
        markers <- tryCatch(
          Seurat::FindMarkers(
            seu,
            assay = assay,
            ident.1 = ident1,
            ident.2 = ident2,
            test.use = input$test_use %||% "wilcox",
            logfc.threshold = as.numeric(input$logfc_threshold %||% 0.25),
            min.pct = as.numeric(input$min_pct %||% 0.1),
            max.cells.per.ident = as.integer(input$max_cells_per_ident %||% 500),
            verbose = FALSE
          ),
          error = function(e) {
            showModal(modalDialog(title = "FindMarkers failed", conditionMessage(e), easyClose = TRUE))
            NULL
          }
        )
        
        if (is.null(markers)) return(NULL)
        markers$gene <- rownames(markers)
        markers$comparison <- paste0(ident1, "_vs_", ident2)
        rownames(markers) <- NULL
        
        pairwise_markers(markers)
        if (!is.null(markers_seu())) {
          seu2 <- markers_seu()
          seu2@misc$CoTRA_markers$pairwise_markers <- markers
          markers_seu(seu2)
        }
        
        volcano_plot_gg(make_volcano_plot(markers, ident1, ident2))
        incProgress(1)
      })
    })
    
    make_volcano_plot <- function(df, ident1 = "Group1", ident2 = "Group2") {
      if (is.null(df) || nrow(df) == 0) return(NULL)
      lfc_col <- if ("avg_log2FC" %in% colnames(df)) "avg_log2FC" else if ("avg_logFC" %in% colnames(df)) "avg_logFC" else NULL
      p_col <- if ("p_val_adj" %in% colnames(df)) "p_val_adj" else if ("p_val" %in% colnames(df)) "p_val" else NULL
      if (is.null(lfc_col) || is.null(p_col)) return(NULL)
      
      df$neg_log10_padj <- -log10(pmax(df[[p_col]], .Machine$double.xmin))
      df$significant <- df[[p_col]] <= as.numeric(input$p_adj_cutoff %||% 0.05) & abs(df[[lfc_col]]) >= as.numeric(input$logfc_threshold %||% 0.25)
      
      ggplot2::ggplot(df, ggplot2::aes(x = .data[[lfc_col]], y = neg_log10_padj, color = significant, text = gene)) +
        ggplot2::geom_point(size = 1, alpha = 0.8) +
        ggplot2::geom_vline(xintercept = c(-as.numeric(input$logfc_threshold %||% 0.25), as.numeric(input$logfc_threshold %||% 0.25)), linetype = "dashed") +
        ggplot2::labs(
          title = paste0("Volcano: ", ident1, " vs ", ident2),
          x = "Average log2 fold-change",
          y = "-log10 adjusted p-value",
          color = "Significant"
        ) +
        safe_theme()
    }
    
    output$all_marker_summary <- renderUI({
      df <- all_markers()
      if (is.null(df)) return(tags$div(style = "color:#777;", "Run all markers to create marker table."))
      tags$div(
        tags$b("Marker rows: "), nrow(df), tags$br(),
        tags$b("Clusters/groups: "), length(unique(df$cluster)), tags$br(),
        tags$b("Assay used: "), tryCatch(markers_seu()@misc$CoTRA_markers$assay, error = function(e) "NA")
      )
    })
    
    output$all_marker_table <- DT::renderDT({
      df <- all_markers()
      if (is.null(df)) df <- data.frame(message = "Run all markers to create this table.", stringsAsFactors = FALSE)
      DT::datatable(df, options = list(pageLength = 25, scrollX = TRUE))
    })
    
    output$pairwise_marker_table <- DT::renderDT({
      df <- pairwise_markers()
      if (is.null(df)) df <- data.frame(message = "Run pairwise comparison to create this table.", stringsAsFactors = FALSE)
      DT::datatable(df, options = list(pageLength = 25, scrollX = TRUE))
    })
    
    output$volcano_plot <- plotly::renderPlotly({
      p <- volcano_plot_gg()
      if (is.null(p)) return(NULL)
      plotly::ggplotly(p, tooltip = c("text", "x", "y", "color"))
    })
    
    get_top_genes <- function(n_per_cluster = 5) {
      df <- all_markers()
      if (is.null(df) || nrow(df) == 0) return(character(0))
      lfc_col <- if ("avg_log2FC" %in% colnames(df)) "avg_log2FC" else if ("avg_logFC" %in% colnames(df)) "avg_logFC" else NULL
      p_col <- if ("p_val_adj" %in% colnames(df)) "p_val_adj" else if ("p_val" %in% colnames(df)) "p_val" else NULL
      if (!is.null(p_col)) df <- df[df[[p_col]] <= as.numeric(input$p_adj_cutoff %||% 0.05), , drop = FALSE]
      if (!is.null(lfc_col)) df <- df[order(df$cluster, -df[[lfc_col]]), , drop = FALSE]
      genes <- unique(unlist(lapply(split(df, df$cluster), function(x) head(x$gene, n_per_cluster))))
      genes[nzchar(genes)]
    }
    
    valid_plot_genes <- function(genes, seu) {
      genes <- unique(genes)
      found <- genes[genes %in% rownames(seu)]
      missing <- setdiff(genes, found)
      list(found = found, missing = missing)
    }
    
    observeEvent(input$make_heatmap, {
      seu <- markers_seu() %||% get_seu()
      req(seu)
      genes <- get_top_genes(as.integer(input$top_n_heatmap %||% 5))
      check <- valid_plot_genes(genes, seu)
      if (length(check$found) == 0) {
        showModal(modalDialog(title = "No heatmap genes", "No valid top marker genes found.", easyClose = TRUE))
        return(NULL)
      }
      max_cells <- as.integer(input$heatmap_max_cells %||% 1000)
      cells <- colnames(seu)
      if (length(cells) > max_cells) {
        set.seed(1234)
        cells <- sample(cells, max_cells)
      }
      p <- tryCatch(
        Seurat::DoHeatmap(seu, features = check$found, cells = cells, group.by = input$identity_col %||% "ident") + safe_theme(),
        error = function(e) {
          showModal(modalDialog(title = "Heatmap failed", conditionMessage(e), easyClose = TRUE))
          NULL
        }
      )
      heatmap_plot_gg(p)
    })
    
    output$marker_heatmap <- plotly::renderPlotly({
      p <- heatmap_plot_gg()
      if (is.null(p)) return(NULL)
      plotly::ggplotly(p)
    })
    
    observeEvent(input$make_dotplot, {
      seu <- markers_seu() %||% get_seu()
      req(seu)
      genes <- parse_gene_text(input$dotplot_genes)
      if (length(genes) == 0) genes <- get_top_genes(as.integer(input$top_n_dotplot %||% 5))
      check <- valid_plot_genes(genes, seu)
      if (length(check$found) == 0) {
        showModal(modalDialog(title = "No DotPlot genes", "No valid genes found for DotPlot.", easyClose = TRUE))
        return(NULL)
      }
      p <- tryCatch(
        {
          p <- Seurat::DotPlot(
            seu,
            features = check$found,
            group.by = input$identity_col %||% NULL
          ) + safe_theme()
          
          if ("RotatedAxis" %in% getNamespaceExports("Seurat")) {
            p <- p + Seurat::RotatedAxis()
          } else {
            p <- p + ggplot2::theme(
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
            )
          }
          
          p
          
        },
        error = function(e) {
          showModal(modalDialog(title = "DotPlot failed", conditionMessage(e), easyClose = TRUE))
          NULL
        }
      )
      dotplot_gg(p)
    })
    
    output$marker_dotplot <- plotly::renderPlotly({
      p <- dotplot_gg()
      if (is.null(p)) return(NULL)
      plotly::ggplotly(p)
    })
    
    feature_missing <- reactiveVal(character(0))
    
    observeEvent(input$make_featureplot, {
      seu <- markers_seu() %||% get_seu()
      req(seu)
      genes <- parse_gene_text(input$featureplot_genes)
      if (length(genes) == 0) genes <- head(get_top_genes(1), 6)
      check <- valid_plot_genes(genes, seu)
      feature_missing(check$missing)
      if (length(check$found) == 0) {
        showModal(modalDialog(title = "No FeaturePlot genes", "No valid genes found for FeaturePlot.", easyClose = TRUE))
        return(NULL)
      }
      reduction <- input$feature_reduction %||% "umap"
      if (!reduction %in% names(seu@reductions)) {
        showModal(modalDialog(title = "Reduction missing", paste(reduction, "was not found."), easyClose = TRUE))
        return(NULL)
      }
      p <- tryCatch(
        Seurat::FeaturePlot(seu, features = check$found, reduction = reduction, ncol = 2) + safe_theme(),
        error = function(e) {
          showModal(modalDialog(title = "FeaturePlot failed", conditionMessage(e), easyClose = TRUE))
          NULL
        }
      )
      featureplot_gg(p)
    })
    
    output$featureplot_missing_genes <- renderUI({
      miss <- feature_missing()
      if (length(miss) == 0) return(NULL)
      tags$div(style = "color:#B00020;", paste("Missing genes:", paste(miss, collapse = ", ")))
    })
    
    output$marker_featureplot <- plotly::renderPlotly({
      p <- featureplot_gg()
      if (is.null(p)) return(NULL)
      plotly::ggplotly(p)
    })
    
    save_plot <- function(plot_reactive, file, device = "pdf", width = 8, height = 7) {
      p <- plot_reactive()
      req(p)
      ggplot2::ggsave(file, plot = p, device = device, width = width, height = height)
    }
    
    output$dl_all_markers <- downloadHandler(
      filename = function() paste0("CoTRA_all_cluster_markers_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) utils::write.csv(all_markers(), file, row.names = FALSE)
    )
    
    output$dl_top_markers <- downloadHandler(
      filename = function() paste0("CoTRA_top_cluster_markers_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) utils::write.csv(top_markers_df(), file, row.names = FALSE)
    )
    
    output$dl_marker_summary <- downloadHandler(
      filename = function() paste0("CoTRA_marker_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) utils::write.csv(marker_summary(), file, row.names = FALSE)
    )
    
    output$dl_pairwise_markers <- downloadHandler(
      filename = function() paste0("CoTRA_pairwise_markers_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) utils::write.csv(pairwise_markers(), file, row.names = FALSE)
    )
    
    output$dl_heatmap_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_marker_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(heatmap_plot_gg, file, "pdf", width = 10, height = 9)
    )
    
    output$dl_heatmap_svg <- downloadHandler(
      filename = function() paste0("CoTRA_marker_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(heatmap_plot_gg, file, "svg", width = 10, height = 9)
    )
    
    output$dl_dotplot_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_marker_dotplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(dotplot_gg, file, "pdf", width = 10, height = 8)
    )
    
    output$dl_dotplot_svg <- downloadHandler(
      filename = function() paste0("CoTRA_marker_dotplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(dotplot_gg, file, "svg", width = 10, height = 8)
    )
    
    output$dl_featureplot_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_marker_featureplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(featureplot_gg, file, "pdf", width = 10, height = 8)
    )
    
    output$dl_featureplot_svg <- downloadHandler(
      filename = function() paste0("CoTRA_marker_featureplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(featureplot_gg, file, "svg", width = 10, height = 8)
    )
    
    output$dl_volcano_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_marker_volcano_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(volcano_plot_gg, file, "pdf", width = 8, height = 7)
    )
    
    output$dl_volcano_svg <- downloadHandler(
      filename = function() paste0("CoTRA_marker_volcano_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(volcano_plot_gg, file, "svg", width = 8, height = 7)
    )
    
    build_marker_project <- function(seu_out, out_file = NULL) {
      list(
        seurat = seu_out,
        project_meta = list(
          project_type = "CoTRA_scRNA_project",
          step = "Markers",
          saved_at = Sys.time(),
          cotra_module = "modules/sc/markers.R",
          output_path = out_file
        ),
        cells_meta = seu_out@meta.data,
        step = "Markers",
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
          clustering = tryCatch(seu_out@misc$CoTRA_clustering, error = function(e) NULL),
          markers = tryCatch(seu_out@misc$CoTRA_markers, error = function(e) NULL)
        ),
        tables = list(
          all_markers = all_markers(),
          pairwise_markers = pairwise_markers(),
          marker_summary = marker_summary(),
          top_markers = top_markers_df()
        ),
        cotra_version = if (exists("cotra_version")) cotra_version else "unknown"
      )
    }
    
    output$download_marker_session <- downloadHandler(
      filename = function() {
        paste0("CoTRA_scRNA_Markers_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
      },
      content = function(file) {
        seu_out <- markers_seu() %||% get_seu()
        req(seu_out)
        
        out_dir <- file.path("outputs", "scRNA", "sessioninfo")
        if (!dir.exists(out_dir)) {
          dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        }
        
        out_file <- file.path(
          out_dir,
          paste0("CoTRA_scRNA_Markers_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
        )
        
        project_obj <- build_marker_project(seu_out, out_file)
        saveRDS(project_obj, out_file)
        file.copy(out_file, file, overwrite = TRUE)
        last_session_path(out_file)
      }
    )
    
    return(list(
      seurat = reactive({
        if (!is.null(markers_seu())) markers_seu() else get_seu()
      }),
      all_markers = reactive(all_markers()),
      pairwise_markers = reactive(pairwise_markers()),
      marker_summary = reactive(marker_summary()),
      top_markers = reactive(top_markers_df()),
      session_path = reactive(last_session_path())
    ))
  })
}
