# ==========================================================
# modules/sc/pathway_activity.R
# CoTRA scRNA-seq Pathway Activity Module
#
# Methods:
# - UCell
# - AUCell
# - GSVA
# - VISION
#
# Outputs:
# - Pathway activity UMAP
# - Cluster-level pathway heatmap
# - Cluster-level dot plot
# - Pathway score table
# - Method summary
# - CSV downloads
# - PDF and SVG plot downloads
# - ZIP export
# - Session RDS save
#
# Return:
# - seurat
# - pathway_ready
# - scores
# - pathways
# - method
# - session_path
# ==========================================================

suppressPackageStartupMessages({
  library(shiny)
  library(Seurat)
  library(ggplot2)
  library(DT)
})

mod_sc_pathway_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("Pathway activity analysis"),
    
    div(
      class = "alert alert-warning",
      strong("Important: "),
      "Pathway activity scores estimate relative pathway activity per cell. They are not direct measurements of pathway activation. Interpret results with cell type, cluster composition, gene detection, and species gene symbols in mind."
    ),
    
    sidebarLayout(
      sidebarPanel(
        selectInput(
          ns("method"),
          "Scoring method",
          choices = c(
            "UCell" = "ucell",
            "AUCell" = "aucell",
            "GSVA" = "gsva",
            "VISION" = "vision"
          ),
          selected = "ucell"
        ),
        
        selectInput(
          ns("species"),
          "Species for MSigDB",
          choices = c(
            "Human" = "Homo sapiens",
            "Mouse" = "Mus musculus"
          ),
          selected = "Homo sapiens"
        ),
        
        selectInput(
          ns("pathway_source"),
          "Pathway set",
          choices = c(
            "Hallmark MSigDB" = "H",
            "C2 canonical pathways MSigDB" = "C2",
            "C5 GO biological process MSigDB" = "C5",
            "Custom GMT upload" = "custom"
          ),
          selected = "H"
        ),
        
        conditionalPanel(
          condition = sprintf("input['%s'] == 'custom'", ns("pathway_source")),
          fileInput(ns("gmt_file"), "Upload GMT file", accept = c(".gmt", "text/plain"))
        ),
        
        uiOutput(ns("reduction_ui")),
        
        uiOutput(ns("cluster_col_ui")),
        
        checkboxInput(
          ns("uppercase_genes"),
          "Convert pathway genes and object genes to uppercase for matching",
          value = FALSE
        ),
        
        numericInput(
          ns("min_genes"),
          "Minimum matched genes per pathway",
          value = 5,
          min = 2,
          max = 100,
          step = 1
        ),
        
        numericInput(
          ns("max_pathways"),
          "Maximum pathways to score",
          value = 100,
          min = 5,
          max = 1000,
          step = 5
        ),
        
        numericInput(
          ns("top_n_plot"),
          "Top pathways for heatmap and dotplot",
          value = 20,
          min = 5,
          max = 100,
          step = 5
        ),
        
        selectInput(
          ns("selected_pathway"),
          "Pathway to plot",
          choices = NULL
        ),
        
        checkboxInput(
          ns("interactive_plots"),
          "Use interactive Plotly plots when possible",
          value = FALSE
        ),
        
        actionButton(
          ns("run_pathway"),
          "Run pathway scoring",
          class = "btn-primary"
        ),
        
        br(), br(),
        
        downloadButton(ns("download_scores"), "Download pathway scores CSV"),
        br(), br(),
        downloadButton(ns("download_summary"), "Download pathway summary CSV"),
        br(), br(),
        downloadButton(ns("download_session"), "Download pathway session RDS"),
        br(), br(),
        downloadButton(ns("download_all_zip"), "Download all pathway outputs ZIP")
      ),
      
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Summary",
            br(),
            htmlOutput(ns("method_summary")),
            br(),
            DTOutput(ns("summary_table"))
          ),
          
          tabPanel(
            "Pathway UMAP",
            br(),
            uiOutput(ns("feature_plot_ui")),
            br(),
            downloadButton(ns("download_feature_pdf"), "PDF"),
            downloadButton(ns("download_feature_svg"), "SVG")
          ),
          
          tabPanel(
            "Heatmap",
            br(),
            plotOutput(ns("pathway_heatmap"), height = "650px"),
            br(),
            downloadButton(ns("download_heatmap_pdf"), "PDF"),
            downloadButton(ns("download_heatmap_svg"), "SVG")
          ),
          
          tabPanel(
            "DotPlot",
            br(),
            plotOutput(ns("pathway_dotplot"), height = "650px"),
            br(),
            downloadButton(ns("download_dotplot_pdf"), "PDF"),
            downloadButton(ns("download_dotplot_svg"), "SVG")
          ),
          
          tabPanel(
            "Scores",
            br(),
            DTOutput(ns("score_table"))
          ),
          
          tabPanel(
            "Pathways",
            br(),
            DTOutput(ns("pathway_table"))
          ),
          
          tabPanel(
            "Help",
            br(),
            h4("When to use this module"),
            p("Use pathway activity analysis when you want to summarize gene-level expression into pathway-level activity scores per cell."),
            
            h4("Recommended method"),
            tags$ul(
              tags$li("UCell is a good default for single-cell data because it is rank-based and robust to sparsity."),
              tags$li("AUCell is also rank-based and useful for gene set activity scoring."),
              tags$li("GSVA can be slower and may be less ideal for very sparse single-cell matrices."),
              tags$li("VISION is useful for signature scoring, but it can be heavier to run.")
            ),
            
            h4("Interpretation"),
            tags$table(
              class = "table table-bordered table-striped",
              tags$thead(
                tags$tr(
                  tags$th("Output"),
                  tags$th("Meaning"),
                  tags$th("Common mistake")
                )
              ),
              tags$tbody(
                tags$tr(
                  tags$td("Pathway UMAP"),
                  tags$td("Shows pathway score distribution across cells."),
                  tags$td("Calling a pathway active without checking gene overlap.")
                ),
                tags$tr(
                  tags$td("Heatmap"),
                  tags$td("Shows mean pathway score by cluster."),
                  tags$td("Ignoring cluster size and cell type composition.")
                ),
                tags$tr(
                  tags$td("DotPlot"),
                  tags$td("Shows pathway score strength and detection fraction by cluster."),
                  tags$td("Treating scores as absolute pathway activation.")
                ),
                tags$tr(
                  tags$td("Matched genes"),
                  tags$td("Number of pathway genes detected in the Seurat object."),
                  tags$td("Using human pathways on mouse-style gene symbols without conversion.")
                )
              )
            )
          )
        )
      )
    )
  )
}

mod_sc_pathway_server <- function(id, seurat_r, sc_state = NULL) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    `%||%` <- function(a, b) {
      if (!is.null(a)) a else b
    }
    
    clean_name <- function(x) {
      x <- gsub("[^A-Za-z0-9_]+", "_", x)
      x <- gsub("_+", "_", x)
      x <- gsub("^_|_$", "", x)
      x
    }
    
    get_obj <- reactive({
      obj <- seurat_r()
      validate(need(inherits(obj, "Seurat"), "Input must be a Seurat object."))
      obj
    })
    
    get_assay_data_safe <- function(obj, assay = NULL, layer = "data", slot = "data") {
      assay <- assay %||% DefaultAssay(obj)
      
      out <- tryCatch({
        GetAssayData(obj, assay = assay, layer = layer)
      }, error = function(e) NULL)
      
      if (!is.null(out) && nrow(out) > 0 && ncol(out) > 0) {
        return(out)
      }
      
      out <- tryCatch({
        suppressWarnings(GetAssayData(obj, assay = assay, slot = slot))
      }, error = function(e) NULL)
      
      if (is.null(out) || nrow(out) == 0 || ncol(out) == 0) {
        stop("Could not access normalized assay data.")
      }
      
      out
    }
    
    get_cluster_col <- function(obj) {
      if ("CoTRA_cluster_label" %in% colnames(obj@meta.data)) {
        return("CoTRA_cluster_label")
      }
      
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        return("seurat_clusters")
      }
      
      obj$CoTRA_active_ident <- as.character(Idents(obj))
      "CoTRA_active_ident"
    }
    
    output$reduction_ui <- renderUI({
      obj <- get_obj()
      reds <- names(obj@reductions)
      validate(need(length(reds) > 0, "No reductions found. Run UMAP or tSNE first."))
      
      selected <- if ("umap" %in% reds) "umap" else reds[1]
      
      selectInput(
        ns("reduction"),
        "Reduction for plotting",
        choices = reds,
        selected = selected
      )
    })
    
    output$cluster_col_ui <- renderUI({
      obj <- get_obj()
      cols <- colnames(obj@meta.data)
      
      selected <- get_cluster_col(obj)
      
      selectInput(
        ns("cluster_col"),
        "Cluster or group column",
        choices = cols,
        selected = selected
      )
    })
    
    pathway_res <- reactiveVal(NULL)
    
    read_gmt <- function(path) {
      lines <- readLines(path, warn = FALSE)
      gmt <- lapply(lines, function(z) {
        parts <- strsplit(z, "\t")[[1]]
        genes <- unique(parts[-c(1, 2)])
        genes <- genes[nzchar(genes)]
        genes
      })
      
      names(gmt) <- vapply(lines, function(z) strsplit(z, "\t")[[1]][1], character(1))
      gmt
    }
    
    load_msigdb <- function(species, category) {
      validate(need(requireNamespace("msigdbr", quietly = TRUE), "Install msigdbr first."))
      
      msig <- msigdbr::msigdbr(
        species = species,
        category = category
      )
      
      validate(need(nrow(msig) > 0, "No MSigDB pathways found for the selected species and category."))
      
      split(msig$gene_symbol, msig$gs_name)
    }
    
    prepare_pathways <- function(obj, pathways) {
      validate(need(length(pathways) > 0, "No pathways found."))
      
      object_genes <- rownames(obj)
      
      if (isTRUE(input$uppercase_genes)) {
        object_gene_map <- setNames(object_genes, toupper(object_genes))
        pathways <- lapply(pathways, toupper)
      } else {
        object_gene_map <- setNames(object_genes, object_genes)
      }
      
      pathways <- lapply(pathways, function(g) {
        g <- unique(g)
        matched <- intersect(g, names(object_gene_map))
        unique(unname(object_gene_map[matched]))
      })
      
      pathway_sizes <- vapply(pathways, length, integer(1))
      pathways <- pathways[pathway_sizes >= input$min_genes]
      
      validate(need(length(pathways) > 0, "No pathways passed the minimum matched gene filter."))
      
      pathway_sizes <- pathway_sizes[names(pathways)]
      pathways <- pathways[order(pathway_sizes, decreasing = TRUE)]
      
      if (length(pathways) > input$max_pathways) {
        pathways <- pathways[seq_len(input$max_pathways)]
      }
      
      pathways
    }
    
    check_inputs <- function(obj) {
      validate(need(input$method %in% c("ucell", "aucell", "gsva", "vision"), "Select a valid pathway method."))
      
      if (input$method == "ucell") {
        validate(need(requireNamespace("UCell", quietly = TRUE), "Install UCell first."))
      }
      
      if (input$method == "aucell") {
        validate(need(requireNamespace("AUCell", quietly = TRUE), "Install AUCell first."))
      }
      
      if (input$method == "gsva") {
        validate(need(requireNamespace("GSVA", quietly = TRUE), "Install GSVA first."))
      }
      
      if (input$method == "vision") {
        validate(need(requireNamespace("VISION", quietly = TRUE), "Install VISION first."))
      }
      
      if (isTRUE(input$interactive_plots)) {
        validate(need(requireNamespace("plotly", quietly = TRUE), "Install plotly first or disable interactive plots."))
      }
      
      validate(need(!is.null(input$reduction), "Select a reduction."))
      validate(need(input$reduction %in% names(obj@reductions), paste0("Reduction not found: ", input$reduction)))
      
      emb <- Embeddings(obj, input$reduction)
      validate(need(ncol(emb) >= 2, "Selected reduction must have at least two dimensions."))
      
      validate(need(!is.null(input$cluster_col), "Select a cluster or group column."))
      validate(need(input$cluster_col %in% colnames(obj@meta.data), "Selected cluster column was not found."))
      
      mat_ok <- tryCatch({
        mat <- get_assay_data_safe(obj)
        nrow(mat) > 0 && ncol(mat) > 0
      }, error = function(e) FALSE)
      
      validate(need(mat_ok, "No normalized data found. Run NormalizeData or SCTransform first."))
      
      TRUE
    }
    
    score_ucell <- function(mat, pathways) {
      score <- UCell::ScoreSignatures_UCell(
        matrix = mat,
        features = pathways,
        chunk.size = 500
      )
      
      as.data.frame(score)
    }
    
    score_aucell <- function(mat, pathways) {
      rankings <- AUCell::AUCell_buildRankings(mat, plotStats = FALSE, verbose = FALSE)
      auc <- AUCell::AUCell_calcAUC(pathways, rankings)
      score <- t(SummarizedExperiment::assay(auc))
      as.data.frame(score)
    }
    
    score_gsva <- function(mat, pathways) {
      mat_dense <- as.matrix(mat)
      
      score <- tryCatch({
        param <- GSVA::gsvaParam(mat_dense, pathways)
        GSVA::gsva(param)
      }, error = function(e) {
        GSVA::gsva(mat_dense, pathways, method = "gsva", verbose = FALSE)
      })
      
      as.data.frame(t(score))
    }
    
    score_vision <- function(mat, pathways) {
      mat_dense <- as.matrix(mat)
      
      vis <- VISION::Vision(mat_dense, signatures = pathways)
      vis <- VISION::analyze(vis)
      score <- VISION::getSignatureScores(vis)
      
      as.data.frame(score)
    }
    
    add_scores_to_seurat <- function(obj, scores, method) {
      scores <- as.data.frame(scores)
      
      common_cells <- intersect(rownames(scores), colnames(obj))
      validate(need(length(common_cells) > 0, "Score matrix cells do not match Seurat object cells."))
      
      scores <- scores[common_cells, , drop = FALSE]
      
      prefix <- paste0("CoTRA_Pathway_", toupper(method), "_")
      colnames(scores) <- paste0(prefix, clean_name(colnames(scores)))
      
      obj <- Seurat::AddMetaData(obj, metadata = scores)
      
      list(
        seurat = obj,
        scores = scores
      )
    }
    
    build_plot_df <- function(obj, scores) {
      emb <- as.data.frame(Embeddings(obj, input$reduction)[, 1:2, drop = FALSE])
      colnames(emb) <- c("dim1", "dim2")
      emb$cell <- rownames(emb)
      
      emb$cluster <- as.character(obj@meta.data[emb$cell, input$cluster_col])
      
      scores$cell <- rownames(scores)
      
      merge(emb, scores, by = "cell", all.x = TRUE)
    }
    
    summarize_by_cluster <- function(score_df, cluster_vec) {
      score_cols <- setdiff(colnames(score_df), "cell")
      
      tmp <- score_df
      tmp$cluster <- cluster_vec[match(rownames(score_df), names(cluster_vec))]
      
      mean_df <- aggregate(
        tmp[, score_cols, drop = FALSE],
        by = list(cluster = tmp$cluster),
        FUN = function(x) mean(x, na.rm = TRUE)
      )
      
      mean_long <- reshape(
        mean_df,
        varying = score_cols,
        v.names = "mean_score",
        timevar = "pathway",
        times = score_cols,
        direction = "long"
      )
      
      rownames(mean_long) <- NULL
      
      frac_df <- aggregate(
        tmp[, score_cols, drop = FALSE],
        by = list(cluster = tmp$cluster),
        FUN = function(x) mean(x > median(x, na.rm = TRUE), na.rm = TRUE)
      )
      
      frac_long <- reshape(
        frac_df,
        varying = score_cols,
        v.names = "fraction_high",
        timevar = "pathway",
        times = score_cols,
        direction = "long"
      )
      
      rownames(frac_long) <- NULL
      
      out <- merge(mean_long, frac_long, by = c("cluster", "pathway"), all = TRUE)
      out
    }
    
    observeEvent(input$run_pathway, {
      obj <- get_obj()
      check_inputs(obj)
      
      res <- withProgress(message = "Running pathway activity analysis", value = 0, {
        
        incProgress(0.15, detail = "Loading pathways")
        
        pathways_raw <- if (input$pathway_source == "custom") {
          validate(need(!is.null(input$gmt_file), "Upload a GMT file first."))
          read_gmt(input$gmt_file$datapath)
        } else {
          load_msigdb(input$species, input$pathway_source)
        }
        
        pathways <- prepare_pathways(obj, pathways_raw)
        
        incProgress(0.35, detail = "Preparing expression matrix")
        
        mat <- get_assay_data_safe(obj)
        mat <- mat[rownames(mat) %in% unique(unlist(pathways)), , drop = FALSE]
        
        validate(need(nrow(mat) > 0, "No pathway genes were found in the expression matrix."))
        
        incProgress(0.55, detail = "Scoring pathways")
        
        scores <- switch(
          input$method,
          ucell = score_ucell(mat, pathways),
          aucell = score_aucell(mat, pathways),
          gsva = score_gsva(mat, pathways),
          vision = score_vision(mat, pathways)
        )
        
        added <- add_scores_to_seurat(obj, scores, input$method)
        obj <- added$seurat
        scores <- added$scores
        
        pathway_summary <- data.frame(
          pathway = names(pathways),
          matched_genes = vapply(pathways, length, integer(1)),
          score_column = colnames(scores),
          method = input$method,
          source = input$pathway_source,
          species = input$species,
          stringsAsFactors = FALSE
        )
        
        cluster_vec <- as.character(obj@meta.data[[input$cluster_col]])
        names(cluster_vec) <- rownames(obj@meta.data)
        
        cluster_summary <- summarize_by_cluster(scores, cluster_vec)
        
        incProgress(0.8, detail = "Saving pathway session")
        
        obj@misc$CoTRA_pathway_activity <- list(
          method = input$method,
          species = input$species,
          pathway_source = input$pathway_source,
          reduction = input$reduction,
          cluster_col = input$cluster_col,
          pathways = pathways,
          scores = scores,
          pathway_summary = pathway_summary,
          cluster_summary = cluster_summary,
          created = Sys.time()
        )
        
        dir.create("outputs/scRNA/sessioninfo", recursive = TRUE, showWarnings = FALSE)
        
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        session_path <- file.path(
          "outputs/scRNA/sessioninfo",
          paste0("CoTRA_scRNA_PathwayActivity_session_", timestamp, ".rds")
        )
        
        saveRDS(
          list(
            seurat = obj,
            pathway_activity = obj@misc$CoTRA_pathway_activity,
            sessionInfo = sessionInfo()
          ),
          session_path
        )
        
        incProgress(1, detail = "Done")
        
        list(
          seurat = obj,
          pathway_ready = TRUE,
          method = input$method,
          species = input$species,
          pathway_source = input$pathway_source,
          reduction = input$reduction,
          cluster_col = input$cluster_col,
          pathways = pathways,
          scores = scores,
          pathway_summary = pathway_summary,
          cluster_summary = cluster_summary,
          session_path = session_path
        )
      })
      
      pathway_res(res)
      
      updateSelectInput(
        session,
        "selected_pathway",
        choices = colnames(res$scores),
        selected = colnames(res$scores)[1]
      )
    })
    
    get_res <- reactive({
      res <- pathway_res()
      validate(need(!is.null(res), "Run pathway scoring first."))
      res
    })
    
    selected_score_col <- reactive({
      res <- get_res()
      selected <- input$selected_pathway
      
      if (is.null(selected) || !selected %in% colnames(res$scores)) {
        return(colnames(res$scores)[1])
      }
      
      selected
    })
    
    top_pathways <- reactive({
      res <- get_res()
      scores <- res$scores
      
      vars <- apply(scores, 2, stats::var, na.rm = TRUE)
      names(sort(vars, decreasing = TRUE))[seq_len(min(input$top_n_plot, length(vars)))]
    })
    
    feature_plot_obj <- reactive({
      res <- get_res()
      score_col <- selected_score_col()
      
      df <- build_plot_df(res$seurat, res$scores)
      
      ggplot(df, aes(dim1, dim2, color = .data[[score_col]])) +
        geom_point(size = 0.8, alpha = 0.9) +
        scale_color_viridis_c(na.value = "grey80") +
        theme_classic(base_size = 13) +
        labs(
          title = score_col,
          x = paste0(toupper(res$reduction), "_1"),
          y = paste0(toupper(res$reduction), "_2"),
          color = "Score"
        )
    })
    
    heatmap_plot_obj <- reactive({
      res <- get_res()
      keep <- top_pathways()
      
      df <- res$cluster_summary
      df <- df[df$pathway %in% keep, , drop = FALSE]
      
      ggplot(df, aes(cluster, pathway, fill = mean_score)) +
        geom_tile(color = "white") +
        scale_fill_viridis_c() +
        theme_classic(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8)
        ) +
        labs(
          title = "Mean pathway activity by cluster",
          x = "Cluster",
          y = "Pathway",
          fill = "Mean score"
        )
    })
    
    dotplot_obj <- reactive({
      res <- get_res()
      keep <- top_pathways()
      
      df <- res$cluster_summary
      df <- df[df$pathway %in% keep, , drop = FALSE]
      
      ggplot(df, aes(cluster, pathway)) +
        geom_point(aes(size = fraction_high, color = mean_score), alpha = 0.85) +
        scale_color_viridis_c() +
        theme_classic(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8)
        ) +
        labs(
          title = "Pathway activity dot plot",
          x = "Cluster",
          y = "Pathway",
          size = "Fraction high",
          color = "Mean score"
        )
    })
    
    output$feature_plot_ui <- renderUI({
      if (isTRUE(input$interactive_plots)) {
        plotly::plotlyOutput(ns("feature_plotly"), height = "650px")
      } else {
        plotOutput(ns("pathway_featureplot"), height = "650px")
      }
    })
    
    output$pathway_featureplot <- renderPlot({
      feature_plot_obj()
    })
    
    output$pathway_heatmap <- renderPlot({
      heatmap_plot_obj()
    })
    
    output$pathway_dotplot <- renderPlot({
      dotplot_obj()
    })
    
    if (requireNamespace("plotly", quietly = TRUE)) {
      output$feature_plotly <- plotly::renderPlotly({
        plotly::ggplotly(feature_plot_obj())
      })
    }
    
    output$method_summary <- renderUI({
      res <- get_res()
      
      HTML(paste0(
        "<div class='alert alert-info'>",
        "<b>Method:</b> ", toupper(res$method), "<br>",
        "<b>Species:</b> ", res$species, "<br>",
        "<b>Pathway source:</b> ", res$pathway_source, "<br>",
        "<b>Scored pathways:</b> ", ncol(res$scores), "<br>",
        "<b>Cells:</b> ", nrow(res$scores), "<br>",
        "<b>Cluster column:</b> ", res$cluster_col, "<br>",
        "<b>Reduction:</b> ", res$reduction, "<br>",
        "<b>Session path:</b> ", res$session_path,
        "</div>"
      ))
    })
    
    output$summary_table <- renderDT({
      res <- get_res()
      
      summary_df <- data.frame(
        metric = c(
          "Method",
          "Species",
          "Pathway source",
          "Scored pathways",
          "Cells",
          "Cluster column",
          "Reduction",
          "Session path"
        ),
        value = c(
          toupper(res$method),
          res$species,
          res$pathway_source,
          ncol(res$scores),
          nrow(res$scores),
          res$cluster_col,
          res$reduction,
          res$session_path
        ),
        stringsAsFactors = FALSE
      )
      
      datatable(summary_df, options = list(dom = "t", paging = FALSE), rownames = FALSE)
    })
    
    output$score_table <- renderDT({
      res <- get_res()
      datatable(
        cbind(cell = rownames(res$scores), res$scores),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    })
    
    output$pathway_table <- renderDT({
      res <- get_res()
      datatable(
        res$pathway_summary,
        options = list(pageLength = 20, scrollX = TRUE),
        rownames = FALSE
      )
    })
    
    save_plot <- function(plot_obj, file, width = 9, height = 7) {
      ext <- tools::file_ext(file)
      
      if (tolower(ext) == "svg") {
        validate(need(requireNamespace("svglite", quietly = TRUE), "Install svglite to export SVG files."))
      }
      
      ggplot2::ggsave(
        filename = file,
        plot = plot_obj,
        width = width,
        height = height,
        units = "in",
        limitsize = FALSE
      )
    }
    
    output$download_scores <- downloadHandler(
      filename = function() paste0("CoTRA_pathway_scores_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        res <- get_res()
        write.csv(cbind(cell = rownames(res$scores), res$scores), file, row.names = FALSE)
      }
    )
    
    output$download_summary <- downloadHandler(
      filename = function() paste0("CoTRA_pathway_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        write.csv(get_res()$pathway_summary, file, row.names = FALSE)
      }
    )
    
    output$download_session <- downloadHandler(
      filename = function() basename(get_res()$session_path),
      content = function(file) {
        file.copy(get_res()$session_path, file, overwrite = TRUE)
      }
    )
    
    output$download_feature_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_pathway_umap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(feature_plot_obj(), file)
    )
    
    output$download_feature_svg <- downloadHandler(
      filename = function() paste0("CoTRA_pathway_umap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(feature_plot_obj(), file)
    )
    
    output$download_heatmap_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_pathway_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(heatmap_plot_obj(), file, width = 10, height = 8)
    )
    
    output$download_heatmap_svg <- downloadHandler(
      filename = function() paste0("CoTRA_pathway_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(heatmap_plot_obj(), file, width = 10, height = 8)
    )
    
    output$download_dotplot_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_pathway_dotplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(dotplot_obj(), file, width = 10, height = 8)
    )
    
    output$download_dotplot_svg <- downloadHandler(
      filename = function() paste0("CoTRA_pathway_dotplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(dotplot_obj(), file, width = 10, height = 8)
    )
    
    output$download_all_zip <- downloadHandler(
      filename = function() paste0("CoTRA_pathway_activity_outputs_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
      content = function(file) {
        res <- get_res()
        
        tmpdir <- tempfile("CoTRA_pathway_activity_outputs_")
        dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
        
        write.csv(cbind(cell = rownames(res$scores), res$scores), file.path(tmpdir, "CoTRA_pathway_scores.csv"), row.names = FALSE)
        write.csv(res$pathway_summary, file.path(tmpdir, "CoTRA_pathway_summary.csv"), row.names = FALSE)
        write.csv(res$cluster_summary, file.path(tmpdir, "CoTRA_pathway_cluster_summary.csv"), row.names = FALSE)
        
        file.copy(res$session_path, file.path(tmpdir, basename(res$session_path)), overwrite = TRUE)
        
        save_plot(feature_plot_obj(), file.path(tmpdir, "CoTRA_pathway_umap.pdf"))
        save_plot(heatmap_plot_obj(), file.path(tmpdir, "CoTRA_pathway_heatmap.pdf"), width = 10, height = 8)
        save_plot(dotplot_obj(), file.path(tmpdir, "CoTRA_pathway_dotplot.pdf"), width = 10, height = 8)
        
        if (requireNamespace("svglite", quietly = TRUE)) {
          save_plot(feature_plot_obj(), file.path(tmpdir, "CoTRA_pathway_umap.svg"))
          save_plot(heatmap_plot_obj(), file.path(tmpdir, "CoTRA_pathway_heatmap.svg"), width = 10, height = 8)
          save_plot(dotplot_obj(), file.path(tmpdir, "CoTRA_pathway_dotplot.svg"), width = 10, height = 8)
        }
        
        oldwd <- getwd()
        on.exit(setwd(oldwd), add = TRUE)
        setwd(tmpdir)
        
        utils::zip(zipfile = file, files = list.files(tmpdir, recursive = FALSE))
      }
    )
    
    return(list(
      seurat = reactive({
        res <- get_res()
        res$seurat
      }),
      
      pathway_ready = reactive({
        res <- get_res()
        isTRUE(res$pathway_ready)
      }),
      
      scores = reactive({
        res <- get_res()
        res$scores
      }),
      
      pathways = reactive({
        res <- get_res()
        res$pathways
      }),
      
      method = reactive({
        res <- get_res()
        res$method
      }),
      
      session_path = reactive({
        res <- get_res()
        res$session_path
      })
    ))
  })
}