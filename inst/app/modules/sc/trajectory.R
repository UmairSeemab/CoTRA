# ==========================================================
# modules/sc/trajectory.R
# CoTRA scRNA-seq Trajectory Module
#
# Backends:
# - Monocle3, default
# - Slingshot, optional
#
# Main outputs:
# - Trajectory UMAP/reduction plot with real backend trajectory overlay
# - Pseudotime UMAP/reduction plot
# - Cluster trajectory plot
# - Pseudotime distribution by cluster
# - Gene expression over pseudotime
# - Top pseudotime-associated genes table
# - Method summary
# - CSV downloads
# - PDF and SVG plot downloads
# - ZIP export
# - Session RDS save
#
# Return:
# - seurat
# - trajectory_ready
# - pseudotime_table
# - trajectory_genes
# - gene_expression
# - backend
# - trajectory_object
# - session_path
# ==========================================================

suppressPackageStartupMessages({
  library(shiny)
  library(Seurat)
  library(ggplot2)
  library(DT)
})

mod_sc_trajectory_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("Trajectory analysis"),
    
    div(
      class = "alert alert-warning",
      strong("Important: "),
      "Trajectory analysis is useful only when cells are expected to follow a continuous biological process, such as differentiation, development, activation, regeneration, degeneration, or treatment response."
    ),
    
    sidebarLayout(
      sidebarPanel(
        selectInput(
          ns("backend"),
          "Trajectory backend",
          choices = c("Monocle3" = "Monocle3", "Slingshot" = "Slingshot"),
          selected = "Monocle3"
        ),
        
        uiOutput(ns("reduction_ui")),
        
        uiOutput(ns("root_cluster_ui")),
        
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Monocle3'", ns("backend")),
          textAreaInput(
            ns("root_cells"),
            "Optional root cells",
            placeholder = "Paste cell barcodes, one per line or comma separated",
            rows = 3
          ),
          checkboxInput(
            ns("rerun_graph"),
            "Rerun Monocle3 graph learning",
            value = TRUE
          )
        ),
        
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Slingshot'", ns("backend")),
          checkboxInput(
            ns("allow_multiple_lineages"),
            "Allow multiple Slingshot lineages",
            value = TRUE
          ),
          uiOutput(ns("lineage_plot_ui"))
        ),
        
        hr(),
        
        radioButtons(
          ns("gene_test_mode"),
          "Genes to test for pseudotime association",
          choices = c(
            "Top variable genes" = "variable",
            "Selected genes only" = "selected",
            "All genes" = "all"
          ),
          selected = "variable"
        ),
        
        numericInput(
          ns("n_variable_genes"),
          "Number of variable genes to test",
          value = 2000,
          min = 100,
          max = 10000,
          step = 100
        ),
        
        textAreaInput(
          ns("genes"),
          "Genes for expression over pseudotime",
          placeholder = "Example: RHO, NRL, CRX",
          rows = 3
        ),
        
        numericInput(
          ns("top_n_genes"),
          "Number of top pseudotime genes to show",
          value = 50,
          min = 5,
          max = 500,
          step = 5
        ),
        
        checkboxInput(
          ns("store_backend_object"),
          "Store full backend object in Seurat misc",
          value = FALSE
        ),
        
        checkboxInput(
          ns("interactive_plots"),
          "Use interactive Plotly plots when possible",
          value = FALSE
        ),
        
        actionButton(
          ns("run"),
          "Run trajectory",
          class = "btn-primary"
        ),
        
        br(), br(),
        
        downloadButton(ns("download_pseudotime"), "Download pseudotime CSV"),
        br(), br(),
        downloadButton(ns("download_genes"), "Download trajectory genes CSV"),
        br(), br(),
        downloadButton(ns("download_gene_expr"), "Download gene expression CSV"),
        br(), br(),
        downloadButton(ns("download_session"), "Download trajectory session RDS"),
        br(), br(),
        downloadButton(ns("download_all_zip"), "Download all trajectory outputs ZIP")
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
            "Trajectory UMAP",
            br(),
            uiOutput(ns("trajectory_plot_ui")),
            br(),
            downloadButton(ns("download_trajectory_pdf"), "PDF"),
            downloadButton(ns("download_trajectory_svg"), "SVG")
          ),
          
          tabPanel(
            "Pseudotime UMAP",
            br(),
            uiOutput(ns("pseudotime_plot_ui")),
            br(),
            downloadButton(ns("download_pseudotime_pdf"), "PDF"),
            downloadButton(ns("download_pseudotime_svg"), "SVG")
          ),
          
          tabPanel(
            "Cluster trajectory",
            br(),
            plotOutput(ns("cluster_plot"), height = "650px"),
            br(),
            downloadButton(ns("download_cluster_pdf"), "PDF"),
            downloadButton(ns("download_cluster_svg"), "SVG")
          ),
          
          tabPanel(
            "Pseudotime distribution",
            br(),
            plotOutput(ns("distribution_plot"), height = "550px"),
            br(),
            downloadButton(ns("download_distribution_pdf"), "PDF"),
            downloadButton(ns("download_distribution_svg"), "SVG")
          ),
          
          tabPanel(
            "Gene expression",
            br(),
            uiOutput(ns("gene_expression_plot_ui")),
            br(),
            downloadButton(ns("download_gene_expression_pdf"), "PDF"),
            downloadButton(ns("download_gene_expression_svg"), "SVG")
          ),
          
          tabPanel(
            "Pseudotime genes",
            br(),
            DTOutput(ns("trajectory_gene_table"))
          ),
          
          tabPanel(
            "Help",
            br(),
            h4("When to use trajectory analysis"),
            p("Use this module when cells are expected to follow a continuous biological process."),
            tags$ul(
              tags$li("Good examples: differentiation, development, activation, regeneration, degeneration, treatment response."),
              tags$li("Poor examples: unrelated cell types, strong batch effects, discrete disease groups without transition states.")
            ),
            
            h4("Inputs"),
            tags$ul(
              tags$li("A clustered Seurat object."),
              tags$li("A reduced dimension embedding, usually UMAP."),
              tags$li("Cluster labels from CoTRA_cluster_label if available, otherwise seurat_clusters or active identity."),
              tags$li("Root cluster, and optionally root cells for Monocle3.")
            ),
            
            h4("Backends"),
            tags$ul(
              tags$li("Monocle3 learns a principal graph and orders cells in pseudotime."),
              tags$li("Slingshot infers lineage curves using reduced dimensions and cluster labels.")
            ),
            
            h4("Output interpretation"),
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
                  tags$td("Trajectory plot"),
                  tags$td("Shows cells and inferred trajectory structure."),
                  tags$td("Assuming every visible path is biologically valid.")
                ),
                tags$tr(
                  tags$td("Pseudotime"),
                  tags$td("Relative cell ordering from the selected root."),
                  tags$td("Treating pseudotime as real chronological time.")
                ),
                tags$tr(
                  tags$td("Pseudotime genes"),
                  tags$td("Genes whose expression changes with pseudotime."),
                  tags$td("Ignoring cluster composition, lineage, and batch effects.")
                ),
                tags$tr(
                  tags$td("Root cluster"),
                  tags$td("Defines the start of the trajectory."),
                  tags$td("Changing root cluster and comparing direction as if fixed.")
                )
              )
            )
          )
        )
      )
    )
  )
}

mod_sc_trajectory_server <- function(id, seurat_r, sc_state = NULL) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    `%||%` <- function(a, b) {
      if (!is.null(a)) a else b
    }
    
    clean_tokens <- function(x) {
      if (is.null(x) || !nzchar(x)) return(character(0))
      out <- unique(trimws(unlist(strsplit(x, "[,;\\n\\r\\t ]+"))))
      out[nzchar(out)]
    }
    
    clean_genes <- function(x) clean_tokens(x)
    clean_cells <- function(x) clean_tokens(x)
    
    get_obj <- reactive({
      obj <- seurat_r()
      validate(need(inherits(obj, "Seurat"), "Input must be a Seurat object."))
      
      if (!"CoTRA_cluster_label" %in% colnames(obj@meta.data) &&
          !"seurat_clusters" %in% colnames(obj@meta.data)) {
        obj$CoTRA_active_ident <- as.character(Idents(obj))
      }
      
      obj
    })
    
    get_cluster_col <- function(obj) {
      if ("CoTRA_cluster_label" %in% colnames(obj@meta.data)) {
        return("CoTRA_cluster_label")
      }
      
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        return("seurat_clusters")
      }
      
      if ("CoTRA_active_ident" %in% colnames(obj@meta.data)) {
        return("CoTRA_active_ident")
      }
      
      NULL
    }
    
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
      
      if (is.null(out)) {
        stop("Could not access assay data layer or slot: ", layer)
      }
      
      out
    }
    
    get_counts_safe <- function(obj, assay = NULL) {
      assay <- assay %||% DefaultAssay(obj)
      
      out <- tryCatch({
        GetAssayData(obj, assay = assay, layer = "counts")
      }, error = function(e) NULL)
      
      if (!is.null(out) && nrow(out) > 0 && ncol(out) > 0) {
        return(out)
      }
      
      out <- tryCatch({
        suppressWarnings(GetAssayData(obj, assay = assay, slot = "counts"))
      }, error = function(e) NULL)
      
      if (!is.null(out) && nrow(out) > 0 && ncol(out) > 0) {
        return(out)
      }
      
      get_assay_data_safe(obj, assay = assay, layer = "data", slot = "data")
    }
    
    output$reduction_ui <- renderUI({
      obj <- get_obj()
      reds <- names(obj@reductions)
      validate(need(length(reds) > 0, "No reductions found. Run PCA, UMAP, tSNE, Harmony, or another reduction first."))
      
      selected <- if ("umap" %in% reds) "umap" else reds[1]
      
      selectInput(
        ns("reduction"),
        "Reduction",
        choices = reds,
        selected = selected
      )
    })
    
    output$root_cluster_ui <- renderUI({
      obj <- get_obj()
      cluster_col <- get_cluster_col(obj)
      
      validate(need(!is.null(cluster_col), "No cluster labels found."))
      
      clusters <- sort(unique(as.character(obj@meta.data[[cluster_col]])))
      
      selectInput(
        ns("root_cluster"),
        "Root or start cluster",
        choices = clusters,
        selected = clusters[1]
      )
    })
    
    trajectory_res <- reactiveVal(NULL)
    
    output$lineage_plot_ui <- renderUI({
      res <- trajectory_res()
      if (is.null(res) || !identical(res$backend, "Slingshot")) {
        return(NULL)
      }
      
      lineages <- unique(res$pseudotime_table$lineage)
      lineages <- lineages[!is.na(lineages)]
      
      if (length(lineages) == 0) {
        return(NULL)
      }
      
      selectInput(
        ns("lineage_to_plot"),
        "Lineage to highlight",
        choices = c("All", lineages),
        selected = "All"
      )
    })
    
    check_inputs <- function(obj) {
      validate(need(input$backend %in% c("Monocle3", "Slingshot"), "Select a valid backend."))
      
      if (input$backend == "Monocle3") {
        validate(need(requireNamespace("monocle3", quietly = TRUE), "Install monocle3 first."))
        validate(need(requireNamespace("SingleCellExperiment", quietly = TRUE), "Install SingleCellExperiment first."))
        validate(need(requireNamespace("SummarizedExperiment", quietly = TRUE), "Install SummarizedExperiment first."))
        validate(need(requireNamespace("igraph", quietly = TRUE), "Install igraph first."))
      }
      
      if (input$backend == "Slingshot") {
        validate(need(requireNamespace("slingshot", quietly = TRUE), "Install slingshot first."))
        validate(need(requireNamespace("SingleCellExperiment", quietly = TRUE), "Install SingleCellExperiment first."))
      }
      
      if (isTRUE(input$interactive_plots)) {
        validate(need(requireNamespace("plotly", quietly = TRUE), "Install plotly first or disable interactive plots."))
      }
      
      validate(need(!is.null(input$reduction), "Select a reduction."))
      validate(need(input$reduction %in% names(obj@reductions), paste0("Reduction not found: ", input$reduction)))
      
      emb <- Embeddings(obj, input$reduction)
      validate(need(ncol(emb) >= 2, "Selected reduction must have at least two dimensions."))
      validate(need(nrow(emb) >= 30, "Trajectory analysis needs at least 30 cells."))
      validate(need(all(rownames(emb) %in% colnames(obj)), "Reduction cell names do not match the Seurat object."))
      
      cluster_col <- get_cluster_col(obj)
      validate(need(!is.null(cluster_col), "No clusters found. Run clustering first."))
      
      clusters <- as.character(obj@meta.data[[cluster_col]])
      validate(need(length(unique(clusters)) >= 2, "Trajectory analysis needs at least two clusters."))
      validate(need(input$root_cluster %in% unique(clusters), "Selected root cluster was not found."))
      
      root_n <- sum(clusters == input$root_cluster, na.rm = TRUE)
      validate(need(root_n >= 3, "Root cluster has fewer than 3 cells. Select another root cluster."))
      
      assay <- DefaultAssay(obj)
      normalized_ok <- tryCatch({
        mat <- get_assay_data_safe(obj, assay = assay, layer = "data", slot = "data")
        ncol(mat) > 0 && nrow(mat) > 0
      }, error = function(e) FALSE)
      
      validate(need(normalized_ok, "No normalized data found. Run NormalizeData or SCTransform first."))
      
      selected_genes <- clean_genes(input$genes)
      if (length(selected_genes) > 0) {
        missing_genes <- setdiff(selected_genes, rownames(obj))
        validate(need(length(missing_genes) == 0, paste("Missing genes:", paste(missing_genes, collapse = ", "))))
      }
      
      root_cells <- clean_cells(input$root_cells %||% "")
      if (length(root_cells) > 0) {
        missing_cells <- setdiff(root_cells, colnames(obj))
        validate(need(length(missing_cells) == 0, paste("Missing root cells:", paste(missing_cells, collapse = ", "))))
      }
      
      TRUE
    }
    
    build_plot_df <- function(obj, pseudotime_table) {
      emb <- as.data.frame(Embeddings(obj, input$reduction)[, 1:2, drop = FALSE])
      colnames(emb) <- c("dim1", "dim2")
      emb$cell <- rownames(emb)
      
      cluster_col <- get_cluster_col(obj)
      emb$cluster <- as.character(obj@meta.data[emb$cell, cluster_col])
      
      out <- merge(emb, pseudotime_table, by = c("cell", "cluster"), all.x = TRUE)
      out
    }
    
    extract_monocle3_graph_df <- function(cds) {
      out <- tryCatch({
        graph <- monocle3::principal_graph(cds)[["UMAP"]]
        coords <- monocle3::principal_graph_aux(cds)[["UMAP"]]$dp_mst
        
        if (is.null(graph) || is.null(coords)) {
          return(data.frame())
        }
        
        coords <- as.data.frame(t(coords))
        coords$node <- rownames(coords)
        colnames(coords)[1:2] <- c("dim1", "dim2")
        
        edges <- igraph::as_data_frame(graph, what = "edges")
        colnames(edges)[1:2] <- c("from", "to")
        
        edge_df <- do.call(rbind, lapply(seq_len(nrow(edges)), function(i) {
          p1 <- coords[match(edges$from[i], coords$node), , drop = FALSE]
          p2 <- coords[match(edges$to[i], coords$node), , drop = FALSE]
          
          if (nrow(p1) == 0 || nrow(p2) == 0 || any(is.na(c(p1$dim1, p1$dim2, p2$dim1, p2$dim2)))) {
            return(NULL)
          }
          
          data.frame(
            x = p1$dim1,
            y = p1$dim2,
            xend = p2$dim1,
            yend = p2$dim2,
            stringsAsFactors = FALSE
          )
        }))
        
        if (is.null(edge_df)) data.frame() else edge_df
      }, error = function(e) data.frame())
      
      out
    }
    
    extract_slingshot_curve_df <- function(sds) {
      out <- tryCatch({
        curves <- slingshot::slingCurves(sds)
        
        if (length(curves) == 0) {
          return(data.frame())
        }
        
        do.call(rbind, lapply(seq_along(curves), function(i) {
          curve <- curves[[i]]
          s <- as.data.frame(curve$s[, 1:2, drop = FALSE])
          colnames(s) <- c("dim1", "dim2")
          s$order <- seq_len(nrow(s))
          s$lineage <- paste0("Lineage", i)
          s
        }))
      }, error = function(e) data.frame())
      
      out
    }
    
    run_monocle3 <- function(obj) {
      assay <- DefaultAssay(obj)
      counts <- get_counts_safe(obj, assay = assay)
      
      cell_meta <- obj@meta.data
      gene_meta <- data.frame(
        gene_short_name = rownames(counts),
        row.names = rownames(counts),
        stringsAsFactors = FALSE
      )
      
      cds <- monocle3::new_cell_data_set(
        expression_data = counts,
        cell_metadata = cell_meta,
        gene_metadata = gene_meta
      )
      
      cds <- monocle3::preprocess_cds(
        cds,
        num_dim = min(50, max(2, ncol(cds) - 1))
      )
      
      emb <- Embeddings(obj, input$reduction)[, 1:2, drop = FALSE]
      colnames(emb) <- c("UMAP_1", "UMAP_2")
      emb <- emb[colnames(cds), , drop = FALSE]
      
      SingleCellExperiment::reducedDims(cds)$UMAP <- emb
      
      cluster_col <- get_cluster_col(obj)
      clusters <- factor(obj@meta.data[colnames(cds), cluster_col])
      names(clusters) <- colnames(cds)
      
      cds@clusters$UMAP$clusters <- clusters
      cds@clusters$UMAP$partitions <- factor(rep(1, ncol(cds)), levels = 1)
      names(cds@clusters$UMAP$partitions) <- colnames(cds)
      
      can_reuse_graph <- FALSE
      old <- obj@misc$CoTRA_trajectory %||% NULL
      
      if (!isTRUE(input$rerun_graph) &&
          !is.null(old) &&
          identical(old$backend, "Monocle3") &&
          !is.null(old$backend_object)) {
        old_graph <- tryCatch(monocle3::principal_graph(old$backend_object)[["UMAP"]], error = function(e) NULL)
        if (!is.null(old_graph)) {
          cds <- old$backend_object
          can_reuse_graph <- TRUE
        }
      }
      
      if (!can_reuse_graph) {
        cds <- tryCatch({
          monocle3::learn_graph(cds, use_partition = FALSE)
        }, error = function(e) {
          stop("Monocle3 graph learning failed: ", e$message)
        })
      }
      
      root_cells <- clean_cells(input$root_cells %||% "")
      
      if (length(root_cells) == 0) {
        root_cells <- colnames(obj)[as.character(obj@meta.data[[cluster_col]]) == input$root_cluster]
      }
      
      cds <- monocle3::order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)
      
      pt <- monocle3::pseudotime(cds)
      
      pseudotime_table <- data.frame(
        cell = names(pt),
        pseudotime = as.numeric(pt),
        mean_pseudotime = as.numeric(pt),
        lineage = "Monocle3_graph",
        cluster = as.character(obj@meta.data[names(pt), cluster_col]),
        backend = "Monocle3",
        stringsAsFactors = FALSE
      )
      
      graph_df <- extract_monocle3_graph_df(cds)
      
      list(
        pseudotime_table = pseudotime_table,
        backend_object = cds,
        graph_df = graph_df,
        curve_df = data.frame(),
        n_lineages = 1
      )
    }
    
    run_slingshot <- function(obj) {
      emb <- Embeddings(obj, input$reduction)[, 1:2, drop = FALSE]
      cluster_col <- get_cluster_col(obj)
      clusters <- as.factor(obj@meta.data[rownames(emb), cluster_col])
      
      sds <- slingshot::slingshot(
        data = emb,
        clusterLabels = clusters,
        start.clus = input$root_cluster,
        stretch = 2
      )
      
      pt_mat <- slingshot::slingPseudotime(sds)
      lineage_names <- colnames(pt_mat)
      
      if (is.null(lineage_names)) {
        lineage_names <- paste0("Lineage", seq_len(ncol(pt_mat)))
        colnames(pt_mat) <- lineage_names
      }
      
      mean_pt <- apply(pt_mat, 1, function(x) {
        if (all(is.na(x))) return(NA_real_)
        mean(x, na.rm = TRUE)
      })
      
      assigned_lineage <- apply(pt_mat, 1, function(x) {
        if (all(is.na(x))) return(NA_character_)
        lineage_names[which.min(ifelse(is.na(x), Inf, x))]
      })
      
      pseudotime_table <- data.frame(
        cell = rownames(pt_mat),
        pseudotime = as.numeric(mean_pt),
        mean_pseudotime = as.numeric(mean_pt),
        lineage = assigned_lineage,
        cluster = as.character(obj@meta.data[rownames(pt_mat), cluster_col]),
        backend = "Slingshot",
        stringsAsFactors = FALSE
      )
      
      for (ln in lineage_names) {
        pseudotime_table[[paste0("pseudotime_", make.names(ln))]] <- as.numeric(pt_mat[, ln])
      }
      
      if (!isTRUE(input$allow_multiple_lineages) && ncol(pt_mat) > 1) {
        first_lineage <- lineage_names[1]
        pseudotime_table$pseudotime <- as.numeric(pt_mat[, first_lineage])
        pseudotime_table$mean_pseudotime <- pseudotime_table$pseudotime
        pseudotime_table$lineage <- first_lineage
      }
      
      curve_df <- extract_slingshot_curve_df(sds)
      
      list(
        pseudotime_table = pseudotime_table,
        backend_object = sds,
        graph_df = data.frame(),
        curve_df = curve_df,
        n_lineages = length(lineage_names)
      )
    }
    
    get_genes_for_testing <- function(obj) {
      selected <- clean_genes(input$genes)
      mat <- get_assay_data_safe(obj, assay = DefaultAssay(obj), layer = "data", slot = "data")
      
      if (input$gene_test_mode == "selected") {
        return(intersect(selected, rownames(mat)))
      }
      
      if (input$gene_test_mode == "all") {
        return(rownames(mat))
      }
      
      variable_genes <- tryCatch(VariableFeatures(obj), error = function(e) character(0))
      variable_genes <- intersect(variable_genes, rownames(mat))
      
      if (length(variable_genes) == 0) {
        vars <- tryCatch({
          Matrix::rowMeans(mat ^ 2) - Matrix::rowMeans(mat) ^ 2
        }, error = function(e) {
          apply(as.matrix(mat), 1, stats::var)
        })
        
        variable_genes <- names(sort(vars, decreasing = TRUE))
      }
      
      head(variable_genes, input$n_variable_genes)
    }
    
    test_one_gene_lm <- function(expr, pt) {
      fit <- tryCatch(summary(lm(expr ~ pt)), error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      
      coef_table <- fit$coefficients
      if (!"pt" %in% rownames(coef_table)) return(NULL)
      
      list(
        statistic = unname(coef_table["pt", "t value"]),
        estimate = unname(coef_table["pt", "Estimate"]),
        p_value = unname(coef_table["pt", "Pr(>|t|)"]),
        method = "lm"
      )
    }
    
    test_one_gene_gam <- function(expr, pt) {
      if (!requireNamespace("mgcv", quietly = TRUE)) {
        return(test_one_gene_lm(expr, pt))
      }
      
      dat <- data.frame(expr = expr, pt = pt)
      
      fit <- tryCatch({
        mgcv::gam(expr ~ s(pt, k = min(6, max(3, length(unique(pt)) - 1))), data = dat, method = "REML")
      }, error = function(e) NULL)
      
      if (is.null(fit)) {
        return(test_one_gene_lm(expr, pt))
      }
      
      sm <- tryCatch(summary(fit), error = function(e) NULL)
      if (is.null(sm) || is.null(sm$s.table) || nrow(sm$s.table) == 0) {
        return(test_one_gene_lm(expr, pt))
      }
      
      list(
        statistic = unname(sm$s.table[1, "F"]),
        estimate = suppressWarnings(cor(expr, pt, method = "spearman")),
        p_value = unname(sm$s.table[1, "p-value"]),
        method = "gam"
      )
    }
    
    find_pseudotime_genes <- function(obj, pseudotime_table) {
      mat <- get_assay_data_safe(obj, assay = DefaultAssay(obj), layer = "data", slot = "data")
      
      cells <- pseudotime_table$cell[!is.na(pseudotime_table$pseudotime)]
      cells <- intersect(cells, colnames(mat))
      
      if (length(cells) < 10) {
        return(data.frame())
      }
      
      pt <- pseudotime_table$pseudotime[match(cells, pseudotime_table$cell)]
      keep <- !is.na(pt)
      cells <- cells[keep]
      pt <- pt[keep]
      
      if (length(unique(pt)) < 5) {
        return(data.frame())
      }
      
      genes_to_test <- get_genes_for_testing(obj)
      genes_to_test <- intersect(genes_to_test, rownames(mat))
      
      if (length(genes_to_test) == 0) {
        return(data.frame())
      }
      
      results <- vector("list", length(genes_to_test))
      
      withProgress(message = "Testing pseudotime genes", value = 0, {
        for (i in seq_along(genes_to_test)) {
          if (i %% 100 == 0) {
            incProgress(100 / length(genes_to_test), detail = paste("Tested", i, "of", length(genes_to_test), "genes"))
          }
          
          gene <- genes_to_test[i]
          expr <- as.numeric(mat[gene, cells])
          
          if (length(unique(expr)) < 3) {
            next
          }
          
          tst <- test_one_gene_gam(expr, pt)
          
          if (is.null(tst)) {
            next
          }
          
          direction <- ifelse(is.na(tst$estimate), NA_character_, ifelse(tst$estimate >= 0, "increasing", "decreasing"))
          
          results[[i]] <- data.frame(
            gene = gene,
            statistic = tst$statistic,
            estimate = tst$estimate,
            p_value = tst$p_value,
            correlation = suppressWarnings(cor(expr, pt, method = "spearman")),
            direction = direction,
            method = tst$method,
            backend = input$backend,
            stringsAsFactors = FALSE
          )
        }
      })
      
      gene_result <- do.call(rbind, results)
      
      if (is.null(gene_result) || nrow(gene_result) == 0) {
        return(data.frame())
      }
      
      gene_result$padj <- p.adjust(gene_result$p_value, method = "BH")
      gene_result <- gene_result[order(gene_result$padj, gene_result$p_value), , drop = FALSE]
      rownames(gene_result) <- NULL
      
      head(gene_result, input$top_n_genes)
    }
    
    build_gene_expr_df <- function(obj, pseudotime_table, genes) {
      if (length(genes) == 0) {
        return(data.frame())
      }
      
      mat <- get_assay_data_safe(obj, assay = DefaultAssay(obj), layer = "data", slot = "data")
      
      cells <- pseudotime_table$cell[!is.na(pseudotime_table$pseudotime)]
      cells <- intersect(cells, colnames(mat))
      
      genes <- intersect(genes, rownames(mat))
      
      if (length(cells) == 0 || length(genes) == 0) {
        return(data.frame())
      }
      
      out <- lapply(genes, function(gene) {
        data.frame(
          cell = cells,
          gene = gene,
          expression = as.numeric(mat[gene, cells]),
          pseudotime = pseudotime_table$pseudotime[match(cells, pseudotime_table$cell)],
          cluster = pseudotime_table$cluster[match(cells, pseudotime_table$cell)],
          lineage = pseudotime_table$lineage[match(cells, pseudotime_table$cell)],
          stringsAsFactors = FALSE
        )
      })
      
      do.call(rbind, out)
    }
    
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
    
    observeEvent(input$run, {
      obj <- get_obj()
      check_inputs(obj)
      
      res <- withProgress(message = "Running trajectory analysis", value = 0, {
        
        incProgress(0.15, detail = "Running trajectory backend")
        
        backend_result <- if (input$backend == "Monocle3") {
          run_monocle3(obj)
        } else {
          run_slingshot(obj)
        }
        
        pseudotime_table <- backend_result$pseudotime_table
        pseudotime_table <- pseudotime_table[order(pseudotime_table$pseudotime), , drop = FALSE]
        
        na_fraction <- mean(is.na(pseudotime_table$pseudotime))
        validate(need(na_fraction < 0.8, "More than 80 percent of cells have missing pseudotime. Check root cluster, clusters, reduction, and connected trajectory structure."))
        
        incProgress(0.45, detail = "Finding pseudotime-associated genes")
        
        trajectory_genes <- find_pseudotime_genes(obj, pseudotime_table)
        
        selected_genes <- clean_genes(input$genes)
        
        if (length(selected_genes) == 0 && nrow(trajectory_genes) > 0) {
          selected_genes <- head(trajectory_genes$gene, min(6, nrow(trajectory_genes)))
        }
        
        gene_expr <- build_gene_expr_df(obj, pseudotime_table, selected_genes)
        
        incProgress(0.75, detail = "Saving trajectory session")
        
        obj$CoTRA_pseudotime <- pseudotime_table$pseudotime[match(colnames(obj), pseudotime_table$cell)]
        
        trajectory_misc <- list(
          backend = input$backend,
          reduction = input$reduction,
          root_cluster = input$root_cluster,
          root_cells = clean_cells(input$root_cells %||% ""),
          pseudotime_table = pseudotime_table,
          trajectory_genes = trajectory_genes,
          gene_expression = gene_expr,
          graph_df = backend_result$graph_df,
          curve_df = backend_result$curve_df,
          n_lineages = backend_result$n_lineages,
          created = Sys.time()
        )
        
        if (isTRUE(input$store_backend_object)) {
          trajectory_misc$backend_object <- backend_result$backend_object
        }
        
        obj@misc$CoTRA_trajectory <- trajectory_misc
        
        dir.create("outputs/scRNA/sessioninfo", recursive = TRUE, showWarnings = FALSE)
        
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        session_path <- file.path(
          "outputs/scRNA/sessioninfo",
          paste0("CoTRA_scRNA_Trajectory_session_", timestamp, ".rds")
        )
        
        saveRDS(
          list(
            seurat = obj,
            trajectory = trajectory_misc,
            backend_object = backend_result$backend_object,
            sessionInfo = sessionInfo()
          ),
          session_path
        )
        
        incProgress(1, detail = "Done")
        
        list(
          seurat = obj,
          trajectory_ready = TRUE,
          backend = input$backend,
          reduction = input$reduction,
          root_cluster = input$root_cluster,
          pseudotime_table = pseudotime_table,
          trajectory_genes = trajectory_genes,
          gene_expression = gene_expr,
          backend_object = backend_result$backend_object,
          graph_df = backend_result$graph_df,
          curve_df = backend_result$curve_df,
          n_lineages = backend_result$n_lineages,
          session_path = session_path
        )
      })
      
      trajectory_res(res)
    })
    
    get_res <- reactive({
      res <- trajectory_res()
      validate(need(!is.null(res), "Run trajectory analysis first."))
      res
    })
    
    plot_df <- reactive({
      res <- get_res()
      build_plot_df(res$seurat, res$pseudotime_table)
    })
    
    trajectory_plot_obj <- reactive({
      res <- get_res()
      df <- plot_df()
      
      p <- ggplot(df, aes(dim1, dim2, color = cluster)) +
        geom_point(size = 0.8, alpha = 0.85) +
        theme_classic(base_size = 13) +
        labs(
          title = paste(res$backend, "trajectory"),
          x = paste0(toupper(res$reduction), "_1"),
          y = paste0(toupper(res$reduction), "_2"),
          color = "Cluster"
        )
      
      if (res$backend == "Monocle3" && nrow(res$graph_df) > 0) {
        p <- p +
          geom_segment(
            data = res$graph_df,
            aes(x = x, y = y, xend = xend, yend = yend),
            inherit.aes = FALSE,
            linewidth = 0.8,
            color = "black",
            alpha = 0.8
          )
      }
      
      if (res$backend == "Slingshot" && nrow(res$curve_df) > 0) {
        curve_df <- res$curve_df
        
        if (!is.null(input$lineage_to_plot) &&
            input$lineage_to_plot != "All" &&
            "lineage" %in% colnames(curve_df)) {
          curve_df <- curve_df[curve_df$lineage == input$lineage_to_plot, , drop = FALSE]
        }
        
        p <- p +
          geom_path(
            data = curve_df,
            aes(dim1, dim2, group = lineage),
            inherit.aes = FALSE,
            linewidth = 1.1,
            color = "black",
            alpha = 0.85
          )
      }
      
      p
    })
    
    pseudotime_plot_obj <- reactive({
      res <- get_res()
      df <- plot_df()
      
      ggplot(df, aes(dim1, dim2, color = pseudotime)) +
        geom_point(size = 0.8, alpha = 0.9) +
        scale_color_viridis_c(na.value = "grey80") +
        theme_classic(base_size = 13) +
        labs(
          title = "Pseudotime UMAP",
          x = paste0(toupper(res$reduction), "_1"),
          y = paste0(toupper(res$reduction), "_2"),
          color = "Pseudotime"
        )
    })
    
    cluster_plot_obj <- reactive({
      res <- get_res()
      df <- plot_df()
      
      centers <- aggregate(cbind(dim1, dim2) ~ cluster, data = df, FUN = median)
      pt_centers <- aggregate(pseudotime ~ cluster, data = df, FUN = median, na.rm = TRUE)
      centers$pseudotime <- pt_centers$pseudotime[match(centers$cluster, pt_centers$cluster)]
      centers <- centers[order(centers$pseudotime), , drop = FALSE]
      
      ggplot(df, aes(dim1, dim2, color = cluster)) +
        geom_point(size = 0.7, alpha = 0.65) +
        geom_path(data = centers, aes(dim1, dim2, group = 1), inherit.aes = FALSE, linewidth = 1.1) +
        geom_point(data = centers, aes(dim1, dim2), inherit.aes = FALSE, size = 4) +
        geom_text(data = centers, aes(dim1, dim2, label = cluster), inherit.aes = FALSE, vjust = -1) +
        theme_classic(base_size = 13) +
        labs(
          title = "Cluster trajectory",
          x = paste0(toupper(res$reduction), "_1"),
          y = paste0(toupper(res$reduction), "_2"),
          color = "Cluster"
        )
    })
    
    distribution_plot_obj <- reactive({
      df <- plot_df()
      
      ggplot(df, aes(cluster, pseudotime, fill = cluster)) +
        geom_violin(trim = FALSE, alpha = 0.75, na.rm = TRUE) +
        geom_boxplot(width = 0.12, outlier.size = 0.4, alpha = 0.7, na.rm = TRUE) +
        theme_classic(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
          title = "Pseudotime distribution by cluster",
          x = "Cluster",
          y = "Pseudotime"
        ) +
        guides(fill = "none")
    })
    
    gene_expression_plot_obj <- reactive({
      res <- get_res()
      gene_expr <- res$gene_expression
      
      validate(need(nrow(gene_expr) > 0, "No valid genes selected."))
      
      ggplot(gene_expr, aes(pseudotime, expression, color = cluster)) +
        geom_point(size = 0.7, alpha = 0.5, na.rm = TRUE) +
        geom_smooth(aes(group = gene), method = "loess", se = FALSE, color = "black", linewidth = 0.8, na.rm = TRUE) +
        facet_wrap(~ gene, scales = "free_y") +
        theme_classic(base_size = 13) +
        labs(
          title = "Gene expression over pseudotime",
          x = "Pseudotime",
          y = "Normalized expression",
          color = "Cluster"
        )
    })
    
    output$trajectory_plot_ui <- renderUI({
      if (isTRUE(input$interactive_plots)) {
        plotly::plotlyOutput(ns("trajectory_plotly"), height = "650px")
      } else {
        plotOutput(ns("trajectory_plot"), height = "650px")
      }
    })
    
    output$pseudotime_plot_ui <- renderUI({
      if (isTRUE(input$interactive_plots)) {
        plotly::plotlyOutput(ns("pseudotime_plotly"), height = "650px")
      } else {
        plotOutput(ns("pseudotime_plot"), height = "650px")
      }
    })
    
    output$gene_expression_plot_ui <- renderUI({
      if (isTRUE(input$interactive_plots)) {
        plotly::plotlyOutput(ns("gene_expression_plotly"), height = "650px")
      } else {
        plotOutput(ns("gene_expression_plot"), height = "650px")
      }
    })
    
    output$trajectory_plot <- renderPlot({
      trajectory_plot_obj()
    })
    
    output$pseudotime_plot <- renderPlot({
      pseudotime_plot_obj()
    })
    
    output$cluster_plot <- renderPlot({
      cluster_plot_obj()
    })
    
    output$distribution_plot <- renderPlot({
      distribution_plot_obj()
    })
    
    output$gene_expression_plot <- renderPlot({
      gene_expression_plot_obj()
    })
    
    if (requireNamespace("plotly", quietly = TRUE)) {
      output$trajectory_plotly <- plotly::renderPlotly({
        plotly::ggplotly(trajectory_plot_obj())
      })
      
      output$pseudotime_plotly <- plotly::renderPlotly({
        plotly::ggplotly(pseudotime_plot_obj())
      })
      
      output$gene_expression_plotly <- plotly::renderPlotly({
        plotly::ggplotly(gene_expression_plot_obj())
      })
    }
    
    output$method_summary <- renderUI({
      res <- get_res()
      pt <- res$pseudotime_table
      
      n_cells <- nrow(pt)
      n_pt <- sum(!is.na(pt$pseudotime))
      n_clusters <- length(unique(pt$cluster))
      n_genes <- nrow(res$trajectory_genes)
      
      HTML(paste0(
        "<div class='alert alert-info'>",
        "<b>Backend:</b> ", res$backend, "<br>",
        "<b>Reduction:</b> ", res$reduction, "<br>",
        "<b>Root/start cluster:</b> ", res$root_cluster, "<br>",
        "<b>Cells:</b> ", n_cells, "<br>",
        "<b>Cells with pseudotime:</b> ", n_pt, "<br>",
        "<b>Clusters:</b> ", n_clusters, "<br>",
        "<b>Lineages:</b> ", res$n_lineages, "<br>",
        "<b>Top pseudotime genes shown:</b> ", n_genes, "<br>",
        "<b>Session path:</b> ", res$session_path,
        "</div>"
      ))
    })
    
    output$summary_table <- renderDT({
      res <- get_res()
      pt <- res$pseudotime_table
      
      summary_df <- data.frame(
        metric = c(
          "Backend",
          "Reduction",
          "Root/start cluster",
          "Cells",
          "Cells with pseudotime",
          "Cells without pseudotime",
          "Clusters",
          "Lineages",
          "Trajectory genes shown",
          "Session path"
        ),
        value = c(
          res$backend,
          res$reduction,
          res$root_cluster,
          nrow(pt),
          sum(!is.na(pt$pseudotime)),
          sum(is.na(pt$pseudotime)),
          length(unique(pt$cluster)),
          res$n_lineages,
          nrow(res$trajectory_genes),
          res$session_path
        ),
        stringsAsFactors = FALSE
      )
      
      datatable(summary_df, options = list(dom = "t", paging = FALSE), rownames = FALSE)
    })
    
    output$trajectory_gene_table <- renderDT({
      res <- get_res()
      datatable(
        res$trajectory_genes,
        options = list(pageLength = 20, scrollX = TRUE),
        rownames = FALSE
      )
    })
    
    output$download_pseudotime <- downloadHandler(
      filename = function() paste0("CoTRA_pseudotime_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) write.csv(get_res()$pseudotime_table, file, row.names = FALSE)
    )
    
    output$download_genes <- downloadHandler(
      filename = function() paste0("CoTRA_trajectory_genes_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) write.csv(get_res()$trajectory_genes, file, row.names = FALSE)
    )
    
    output$download_gene_expr <- downloadHandler(
      filename = function() paste0("CoTRA_gene_expression_pseudotime_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) write.csv(get_res()$gene_expression, file, row.names = FALSE)
    )
    
    output$download_session <- downloadHandler(
      filename = function() basename(get_res()$session_path),
      content = function(file) file.copy(get_res()$session_path, file, overwrite = TRUE)
    )
    
    output$download_trajectory_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_trajectory_umap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(trajectory_plot_obj(), file)
    )
    
    output$download_trajectory_svg <- downloadHandler(
      filename = function() paste0("CoTRA_trajectory_umap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(trajectory_plot_obj(), file)
    )
    
    output$download_pseudotime_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_pseudotime_umap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(pseudotime_plot_obj(), file)
    )
    
    output$download_pseudotime_svg <- downloadHandler(
      filename = function() paste0("CoTRA_pseudotime_umap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(pseudotime_plot_obj(), file)
    )
    
    output$download_cluster_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_cluster_trajectory_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(cluster_plot_obj(), file)
    )
    
    output$download_cluster_svg <- downloadHandler(
      filename = function() paste0("CoTRA_cluster_trajectory_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(cluster_plot_obj(), file)
    )
    
    output$download_distribution_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_pseudotime_distribution_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(distribution_plot_obj(), file)
    )
    
    output$download_distribution_svg <- downloadHandler(
      filename = function() paste0("CoTRA_pseudotime_distribution_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(distribution_plot_obj(), file)
    )
    
    output$download_gene_expression_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_gene_expression_pseudotime_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(gene_expression_plot_obj(), file, width = 10, height = 7)
    )
    
    output$download_gene_expression_svg <- downloadHandler(
      filename = function() paste0("CoTRA_gene_expression_pseudotime_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(gene_expression_plot_obj(), file, width = 10, height = 7)
    )
    
    output$download_all_zip <- downloadHandler(
      filename = function() paste0("CoTRA_trajectory_outputs_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
      content = function(file) {
        res <- get_res()
        
        tmpdir <- tempfile("CoTRA_trajectory_outputs_")
        dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
        
        write.csv(res$pseudotime_table, file.path(tmpdir, "CoTRA_pseudotime.csv"), row.names = FALSE)
        write.csv(res$trajectory_genes, file.path(tmpdir, "CoTRA_trajectory_genes.csv"), row.names = FALSE)
        write.csv(res$gene_expression, file.path(tmpdir, "CoTRA_gene_expression_pseudotime.csv"), row.names = FALSE)
        
        file.copy(res$session_path, file.path(tmpdir, basename(res$session_path)), overwrite = TRUE)
        
        save_plot(trajectory_plot_obj(), file.path(tmpdir, "CoTRA_trajectory_umap.pdf"))
        save_plot(pseudotime_plot_obj(), file.path(tmpdir, "CoTRA_pseudotime_umap.pdf"))
        save_plot(cluster_plot_obj(), file.path(tmpdir, "CoTRA_cluster_trajectory.pdf"))
        save_plot(distribution_plot_obj(), file.path(tmpdir, "CoTRA_pseudotime_distribution.pdf"))
        save_plot(gene_expression_plot_obj(), file.path(tmpdir, "CoTRA_gene_expression_pseudotime.pdf"), width = 10, height = 7)
        
        if (requireNamespace("svglite", quietly = TRUE)) {
          save_plot(trajectory_plot_obj(), file.path(tmpdir, "CoTRA_trajectory_umap.svg"))
          save_plot(pseudotime_plot_obj(), file.path(tmpdir, "CoTRA_pseudotime_umap.svg"))
          save_plot(cluster_plot_obj(), file.path(tmpdir, "CoTRA_cluster_trajectory.svg"))
          save_plot(distribution_plot_obj(), file.path(tmpdir, "CoTRA_pseudotime_distribution.svg"))
          save_plot(gene_expression_plot_obj(), file.path(tmpdir, "CoTRA_gene_expression_pseudotime.svg"), width = 10, height = 7)
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
      
      trajectory_ready = reactive({
        res <- get_res()
        isTRUE(res$trajectory_ready)
      }),
      
      pseudotime_table = reactive({
        res <- get_res()
        res$pseudotime_table
      }),
      
      trajectory_genes = reactive({
        res <- get_res()
        res$trajectory_genes
      }),
      
      gene_expression = reactive({
        res <- get_res()
        res$gene_expression
      }),
      
      backend = reactive({
        res <- get_res()
        res$backend
      }),
      
      trajectory_object = reactive({
        res <- get_res()
        res$backend_object
      }),
      
      session_path = reactive({
        res <- get_res()
        res$session_path
      })
    ))
  })
}
