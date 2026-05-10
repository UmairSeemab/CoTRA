# ==========================================================
# modules/sc/annotation.R
# CoTRA scRNA-seq Cell Type Annotation Module
#
# Methods:
# - scType marker scoring
# - SingleR
# - CellTypist prediction import
# - Seurat label transfer
#
# Outputs:
# - Annotation UMAP
# - Confidence UMAP
# - Annotation composition plot
# - Annotation table
# - Marker score table
# - PDF and SVG plot downloads
# - ZIP export
# - Session RDS save
#
# Return:
# - seurat
# - annotation_ready
# - annotation_table
# - method
# - session_path
# ==========================================================

suppressPackageStartupMessages({
  library(shiny)
  library(Seurat)
  library(ggplot2)
  library(DT)
})

mod_sc_annotation_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("Cell type annotation"),
    
    div(
      class = "alert alert-warning",
      strong("Important: "),
      "Cell type annotation predicts biological labels from expression patterns. These labels should be checked with known marker genes, sample origin, clustering, and biological context."
    ),
    
    sidebarLayout(
      sidebarPanel(
        selectInput(
          ns("method"),
          "Annotation method",
          choices = c(
            "scType marker scoring" = "sctype",
            "SingleR" = "singler",
            "CellTypist prediction import" = "celltypist",
            "Seurat label transfer" = "seurat_transfer"
          ),
          selected = "sctype"
        ),
        
        uiOutput(ns("reduction_ui")),
        uiOutput(ns("cluster_col_ui")),
        
        conditionalPanel(
          condition = sprintf("input['%s'] == 'sctype'", ns("method")),
          textAreaInput(
            ns("marker_text"),
            "Marker list",
            value = "Rod: RHO, PDE6A, PDE6B, GNAT1, SAG\nCone: OPN1SW, OPN1MW, ARR3, GNAT2\nBipolar: VSX2, PRKCA, GRM6, CABP5\nMuller glia: RLBP1, GLUL, SOX9, APOE\nAmacrine: GAD1, GAD2, TFAP2A, SLC6A9\nGanglion: RBPMS, SNCG, THY1, POU4F1\nMicroglia: PTPRC, AIF1, CX3CR1, C1QA\nEndothelial: PECAM1, VWF, KDR",
            rows = 10,
            placeholder = "Format:\nCellType1: GENE1, GENE2, GENE3\nCellType2: GENE4, GENE5"
          )
        ),
        
        conditionalPanel(
          condition = sprintf("input['%s'] == 'singler'", ns("method")),
          selectInput(
            ns("singler_ref"),
            "SingleR reference",
            choices = c(
              "Human Primary Cell Atlas" = "hpca",
              "Blueprint Encode" = "blueprint",
              "Mouse RNA-seq" = "mouse_rnaseq",
              "ImmGen" = "immgen"
            ),
            selected = "hpca"
          )
        ),
        
        conditionalPanel(
          condition = sprintf("input['%s'] == 'celltypist'", ns("method")),
          fileInput(
            ns("celltypist_csv"),
            "Upload CellTypist predictions CSV",
            accept = c(".csv", "text/csv")
          ),
          textInput(
            ns("celltypist_cell_col"),
            "Cell barcode column",
            value = "cell"
          ),
          textInput(
            ns("celltypist_label_col"),
            "Predicted label column",
            value = "predicted_labels"
          ),
          textInput(
            ns("celltypist_conf_col"),
            "Confidence/probability column, optional",
            value = "conf_score"
          )
        ),
        
        conditionalPanel(
          condition = sprintf("input['%s'] == 'seurat_transfer'", ns("method")),
          fileInput(
            ns("reference_rds"),
            "Upload reference Seurat RDS",
            accept = c(".rds")
          ),
          textInput(
            ns("reference_label_col"),
            "Reference label column",
            value = "celltype"
          ),
          selectInput(
            ns("transfer_reduction"),
            "Transfer reduction",
            choices = c("pcaproject", "rpca", "cca"),
            selected = "pcaproject"
          )
        ),
        
        numericInput(
          ns("min_confidence"),
          "Minimum confidence for final labels",
          value = 0,
          min = 0,
          max = 1,
          step = 0.05
        ),
        
        checkboxInput(
          ns("write_labels_to_idents"),
          "Set predicted labels as active identities",
          value = FALSE
        ),
        
        checkboxInput(
          ns("interactive_plots"),
          "Use interactive Plotly plots when possible",
          value = FALSE
        ),
        
        actionButton(
          ns("run_annotation"),
          "Run annotation",
          class = "btn-primary"
        ),
        
        br(), br(),
        downloadButton(ns("download_annotations"), "Download annotation CSV"),
        br(), br(),
        downloadButton(ns("download_scores"), "Download score table CSV"),
        br(), br(),
        downloadButton(ns("download_session"), "Download annotation session RDS"),
        br(), br(),
        downloadButton(ns("download_all_zip"), "Download all annotation outputs ZIP")
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
            "Annotation UMAP",
            br(),
            uiOutput(ns("annotation_umap_ui")),
            br(),
            downloadButton(ns("download_annotation_umap_pdf"), "PDF"),
            downloadButton(ns("download_annotation_umap_svg"), "SVG")
          ),
          
          tabPanel(
            "Confidence UMAP",
            br(),
            uiOutput(ns("confidence_umap_ui")),
            br(),
            downloadButton(ns("download_confidence_umap_pdf"), "PDF"),
            downloadButton(ns("download_confidence_umap_svg"), "SVG")
          ),
          
          tabPanel(
            "Composition",
            br(),
            uiOutput(ns("composition_plot_ui")),
            br(),
            downloadButton(ns("download_composition_pdf"), "PDF"),
            downloadButton(ns("download_composition_svg"), "SVG")
          ),
          
          tabPanel(
            "Annotation table",
            br(),
            DTOutput(ns("annotation_table"))
          ),
          
          tabPanel(
            "Scores",
            br(),
            DTOutput(ns("score_table"))
          ),
          
          tabPanel(
            "Help",
            br(),
            h4("When to use this module"),
            p("Use cell type annotation after quality control, normalization, dimensionality reduction, clustering, and marker analysis."),
            
            h4("Methods"),
            tags$ul(
              tags$li("scType marker scoring uses user-provided marker genes and is transparent for retina or custom tissues."),
              tags$li("SingleR compares cells against reference datasets from celldex."),
              tags$li("CellTypist import lets you use predictions generated outside R, then stores them in CoTRA."),
              tags$li("Seurat label transfer maps labels from a reference Seurat object to the query object.")
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
                  tags$td("Annotation UMAP"),
                  tags$td("Predicted cell labels on the embedding."),
                  tags$td("Accepting labels without checking marker genes.")
                ),
                tags$tr(
                  tags$td("Confidence UMAP"),
                  tags$td("Relative confidence or score of assigned labels."),
                  tags$td("Treating confidence as biological certainty.")
                ),
                tags$tr(
                  tags$td("Composition plot"),
                  tags$td("Cell type frequency by predicted label or cluster."),
                  tags$td("Ignoring batch and sample composition.")
                ),
                tags$tr(
                  tags$td("Score table"),
                  tags$td("Method-specific scores or label probabilities."),
                  tags$td("Comparing scores across methods as if they are identical.")
                )
              )
            )
          )
        )
      )
    )
  )
}

mod_sc_annotation_server <- function(id, seurat_r, sc_state = NULL) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    annotation_res <- reactiveVal(NULL)
    
    `%||%` <- function(a, b) {
      if (!is.null(a)) a else b
    }
    
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
      if ("CoTRA_cluster_label" %in% colnames(obj@meta.data)) return("CoTRA_cluster_label")
      if ("seurat_clusters" %in% colnames(obj@meta.data)) return("seurat_clusters")
      if ("CoTRA_active_ident" %in% colnames(obj@meta.data)) return("CoTRA_active_ident")
      NULL
    }
    
    get_assay_data_safe <- function(obj, assay = NULL, layer = "data", slot = "data") {
      assay <- assay %||% DefaultAssay(obj)
      
      out <- tryCatch({
        GetAssayData(obj, assay = assay, layer = layer)
      }, error = function(e) NULL)
      
      if (!is.null(out) && nrow(out) > 0 && ncol(out) > 0) return(out)
      
      out <- tryCatch({
        suppressWarnings(GetAssayData(obj, assay = assay, slot = slot))
      }, error = function(e) NULL)
      
      if (is.null(out) || nrow(out) == 0 || ncol(out) == 0) {
        stop("Could not access normalized assay data.")
      }
      
      out
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
        "Cluster column",
        choices = cols,
        selected = selected
      )
    })
    
    parse_marker_text <- function(txt) {
      validate(need(nzchar(txt), "Marker list is empty."))
      
      lines <- unlist(strsplit(txt, "\n"))
      lines <- trimws(lines)
      lines <- lines[nzchar(lines)]
      
      out <- list()
      
      for (ln in lines) {
        parts <- strsplit(ln, ":", fixed = TRUE)[[1]]
        if (length(parts) < 2) next
        
        celltype <- trimws(parts[1])
        genes <- trimws(unlist(strsplit(parts[2], "[,; ]+")))
        genes <- unique(genes[nzchar(genes)])
        
        if (nzchar(celltype) && length(genes) > 0) {
          out[[celltype]] <- genes
        }
      }
      
      validate(need(length(out) > 0, "Could not parse marker list. Use format: CellType: GENE1, GENE2"))
      out
    }
    
    run_sctype <- function(obj) {
      mat <- get_assay_data_safe(obj)
      markers <- parse_marker_text(input$marker_text)
      
      object_genes <- rownames(mat)
      
      # Create uppercase mapping
      object_genes_upper <- toupper(object_genes)
      gene_map <- setNames(object_genes, object_genes_upper)
      
      score_mat <- sapply(names(markers), function(ct) {
        
        marker_genes <- toupper(markers[[ct]])
        
        matched <- gene_map[marker_genes]
        matched <- matched[!is.na(matched)]
        matched <- unique(matched)
        
        if (length(matched) == 0) {
          return(rep(NA_real_, ncol(mat)))
        }
        
        Matrix::colMeans(mat[matched, , drop = FALSE])
      })
      
      score_mat <- as.matrix(score_mat)
      rownames(score_mat) <- colnames(mat)
      
      validate(need(ncol(score_mat) > 0, "No marker genes matched the object."))
      
      best_idx <- apply(score_mat, 1, function(x) {
        if (all(is.na(x))) return(NA_integer_)
        which.max(x)
      })
      
      best_score <- apply(score_mat, 1, function(x) {
        if (all(is.na(x))) return(NA_real_)
        max(x, na.rm = TRUE)
      })
      
      second_score <- apply(score_mat, 1, function(x) {
        x <- x[!is.na(x)]
        if (length(x) < 2) return(0)
        sort(x, decreasing = TRUE)[2]
      })
      
      label <- rep(NA_character_, nrow(score_mat))
      ok <- !is.na(best_idx)
      label[ok] <- colnames(score_mat)[best_idx[ok]]
      
      confidence <- best_score - second_score
      confidence_scaled <- confidence
      
      if (all(is.finite(confidence_scaled)) && max(confidence_scaled, na.rm = TRUE) > min(confidence_scaled, na.rm = TRUE)) {
        confidence_scaled <- (confidence_scaled - min(confidence_scaled, na.rm = TRUE)) /
          (max(confidence_scaled, na.rm = TRUE) - min(confidence_scaled, na.rm = TRUE))
      }
      
      annotation_table <- data.frame(
        cell = rownames(score_mat),
        predicted_label = label,
        confidence = as.numeric(confidence_scaled),
        raw_score = as.numeric(best_score),
        method = "scType",
        stringsAsFactors = FALSE
      )
      
      score_table <- data.frame(cell = rownames(score_mat), score_mat, check.names = FALSE)
      
      list(
        annotation_table = annotation_table,
        score_table = score_table
      )
    }
    
    run_singler <- function(obj) {
      validate(need(requireNamespace("SingleR", quietly = TRUE), "Install SingleR first."))
      validate(need(requireNamespace("celldex", quietly = TRUE), "Install celldex first."))
      validate(need(requireNamespace("SummarizedExperiment", quietly = TRUE), "Install SummarizedExperiment first."))
      
      mat <- get_assay_data_safe(obj)
      
      ref <- switch(
        input$singler_ref,
        hpca = celldex::HumanPrimaryCellAtlasData(),
        blueprint = celldex::BlueprintEncodeData(),
        mouse_rnaseq = celldex::MouseRNAseqData(),
        immgen = celldex::ImmGenData()
      )
      
      pred <- SingleR::SingleR(
        test = mat,
        ref = ref,
        labels = ref$label.main
      )
      
      annotation_table <- data.frame(
        cell = rownames(pred),
        predicted_label = pred$labels,
        pruned_label = pred$pruned.labels,
        confidence = if ("delta.next" %in% colnames(pred)) pred$delta.next else NA_real_,
        method = "SingleR",
        stringsAsFactors = FALSE
      )
      
      if (all(is.finite(annotation_table$confidence), na.rm = TRUE)) {
        rng <- range(annotation_table$confidence, na.rm = TRUE)
        if (diff(rng) > 0) {
          annotation_table$confidence <- (annotation_table$confidence - rng[1]) / diff(rng)
        }
      }
      
      score_table <- as.data.frame(pred$scores)
      score_table$cell <- rownames(score_table)
      score_table <- score_table[, c("cell", setdiff(colnames(score_table), "cell")), drop = FALSE]
      
      list(
        annotation_table = annotation_table,
        score_table = score_table
      )
    }
    
    run_celltypist_import <- function(obj) {
      validate(need(!is.null(input$celltypist_csv), "Upload CellTypist prediction CSV first."))
      
      pred <- read.csv(input$celltypist_csv$datapath, stringsAsFactors = FALSE, check.names = FALSE)
      
      validate(need(input$celltypist_cell_col %in% colnames(pred), "Cell barcode column was not found."))
      validate(need(input$celltypist_label_col %in% colnames(pred), "Predicted label column was not found."))
      
      cell_col <- input$celltypist_cell_col
      label_col <- input$celltypist_label_col
      conf_col <- input$celltypist_conf_col
      
      common <- intersect(pred[[cell_col]], colnames(obj))
      validate(need(length(common) > 0, "No CellTypist cells matched the Seurat object."))
      
      pred <- pred[match(common, pred[[cell_col]]), , drop = FALSE]
      
      confidence <- rep(NA_real_, nrow(pred))
      if (nzchar(conf_col) && conf_col %in% colnames(pred)) {
        confidence <- suppressWarnings(as.numeric(pred[[conf_col]]))
      }
      
      annotation_table <- data.frame(
        cell = pred[[cell_col]],
        predicted_label = pred[[label_col]],
        confidence = confidence,
        method = "CellTypist_import",
        stringsAsFactors = FALSE
      )
      
      score_table <- pred
      
      list(
        annotation_table = annotation_table,
        score_table = score_table
      )
    }
    
    run_seurat_transfer <- function(obj) {
      validate(need(!is.null(input$reference_rds), "Upload a reference Seurat RDS first."))
      
      ref <- readRDS(input$reference_rds$datapath)
      validate(need(inherits(ref, "Seurat"), "Reference file must contain a Seurat object."))
      validate(need(input$reference_label_col %in% colnames(ref@meta.data), "Reference label column was not found."))
      
      common_features <- intersect(rownames(ref), rownames(obj))
      validate(need(length(common_features) >= 200, "Too few shared genes between reference and query."))
      
      anchors <- Seurat::FindTransferAnchors(
        reference = ref,
        query = obj,
        features = common_features,
        reduction = input$transfer_reduction
      )
      
      pred <- Seurat::TransferData(
        anchorset = anchors,
        refdata = ref@meta.data[[input$reference_label_col]],
        dims = 1:30
      )
      
      annotation_table <- data.frame(
        cell = rownames(pred),
        predicted_label = pred$predicted.id,
        confidence = pred$prediction.score.max,
        method = "Seurat_label_transfer",
        stringsAsFactors = FALSE
      )
      
      score_table <- data.frame(cell = rownames(pred), pred, check.names = FALSE)
      
      list(
        annotation_table = annotation_table,
        score_table = score_table
      )
    }
    
    check_inputs <- function(obj) {
      validate(need(input$method %in% c("sctype", "singler", "celltypist", "seurat_transfer"), "Select a valid annotation method."))
      
      if (isTRUE(input$interactive_plots)) {
        validate(need(requireNamespace("plotly", quietly = TRUE), "Install plotly first or disable interactive plots."))
      }
      
      validate(need(input$reduction %in% names(obj@reductions), "Selected reduction was not found."))
      validate(need(input$cluster_col %in% colnames(obj@meta.data), "Selected cluster column was not found."))
      
      mat_ok <- tryCatch({
        mat <- get_assay_data_safe(obj)
        nrow(mat) > 0 && ncol(mat) > 0
      }, error = function(e) FALSE)
      
      validate(need(mat_ok, "No normalized data found. Run NormalizeData or SCTransform first."))
      
      TRUE
    }
    
    observeEvent(input$run_annotation, {
      obj <- get_obj()
      check_inputs(obj)
      
      res <- withProgress(message = "Running cell type annotation", value = 0, {
        
        incProgress(0.35, detail = paste("Running", input$method))
        
        backend <- switch(
          input$method,
          sctype = run_sctype(obj),
          singler = run_singler(obj),
          celltypist = run_celltypist_import(obj),
          seurat_transfer = run_seurat_transfer(obj)
        )
        
        annotation_table <- backend$annotation_table
        score_table <- backend$score_table
        
        annotation_table$final_label <- annotation_table$predicted_label
        low_conf <- !is.na(annotation_table$confidence) & annotation_table$confidence < input$min_confidence
        annotation_table$final_label[low_conf] <- "Low_confidence"
        
        obj$CoTRA_celltype <- NA_character_
        obj$CoTRA_celltype[match(annotation_table$cell, colnames(obj))] <- annotation_table$final_label
        
        obj$CoTRA_celltype_method <- input$method
        
        if ("confidence" %in% colnames(annotation_table)) {
          obj$CoTRA_celltype_confidence <- NA_real_
          obj$CoTRA_celltype_confidence[match(annotation_table$cell, colnames(obj))] <- annotation_table$confidence
        }
        
        if (isTRUE(input$write_labels_to_idents)) {
          Idents(obj) <- "CoTRA_celltype"
        }
        
        cluster_vec <- as.character(obj@meta.data[[input$cluster_col]])
        names(cluster_vec) <- rownames(obj@meta.data)
        
        annotation_table$cluster <- cluster_vec[annotation_table$cell]
        
        annotation_table$cluster <- as.character(annotation_table$cluster)
        annotation_table$final_label <- as.character(annotation_table$final_label)
        
        annotation_table$cluster[is.na(annotation_table$cluster) | annotation_table$cluster == ""] <- "Unknown_cluster"
        annotation_table$final_label[is.na(annotation_table$final_label) | annotation_table$final_label == ""] <- "Unknown_celltype"
        
        composition <- aggregate(
          x = list(n = rep(1, nrow(annotation_table))),
          by = list(
            cluster = annotation_table$cluster,
            celltype = annotation_table$final_label
          ),
          FUN = sum
        )
        
        composition <- composition[order(composition$cluster, composition$celltype), , drop = FALSE]
        rownames(composition) <- NULL
        
        obj@misc$CoTRA_celltype_annotation <- list(
          method = input$method,
          reduction = input$reduction,
          cluster_col = input$cluster_col,
          annotation_table = annotation_table,
          score_table = score_table,
          composition = composition,
          created = Sys.time()
        )
        
        dir.create("outputs/scRNA/sessioninfo", recursive = TRUE, showWarnings = FALSE)
        
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        session_path <- file.path(
          "outputs/scRNA/sessioninfo",
          paste0("CoTRA_scRNA_CellTypeAnnotation_session_", timestamp, ".rds")
        )
        
        saveRDS(
          list(
            seurat = obj,
            annotation = obj@misc$CoTRA_celltype_annotation,
            sessionInfo = sessionInfo()
          ),
          session_path
        )
        
        incProgress(1, detail = "Done")
        
        list(
          seurat = obj,
          annotation_ready = TRUE,
          method = input$method,
          reduction = input$reduction,
          cluster_col = input$cluster_col,
          annotation_table = annotation_table,
          score_table = score_table,
          composition = composition,
          session_path = session_path
        )
      })
      
      annotation_res(res)
    })
    
    get_res <- reactive({
      res <- annotation_res()
      validate(need(!is.null(res), "Run cell type annotation first."))
      res
    })
    
    plot_df <- reactive({
      res <- get_res()
      obj <- res$seurat
      
      emb <- as.data.frame(Embeddings(obj, res$reduction)[, 1:2, drop = FALSE])
      colnames(emb) <- c("dim1", "dim2")
      emb$cell <- rownames(emb)
      emb$celltype <- obj@meta.data[emb$cell, "CoTRA_celltype"]
      emb$confidence <- obj@meta.data[emb$cell, "CoTRA_celltype_confidence"]
      emb$cluster <- as.character(obj@meta.data[emb$cell, res$cluster_col])
      
      emb
    })
    
    annotation_umap_obj <- reactive({
      df <- plot_df()
      
      ggplot(df, aes(dim1, dim2, color = celltype, text = paste(cell, celltype))) +
        geom_point(size = 0.8, alpha = 0.85) +
        theme_classic(base_size = 13) +
        labs(
          title = "Cell type annotation",
          x = "Dimension 1",
          y = "Dimension 2",
          color = "Cell type"
        )
    })
    
    confidence_umap_obj <- reactive({
      df <- plot_df()
      
      ggplot(df, aes(dim1, dim2, color = confidence, text = paste(cell, confidence))) +
        geom_point(size = 0.8, alpha = 0.85) +
        scale_color_viridis_c(na.value = "grey80") +
        theme_classic(base_size = 13) +
        labs(
          title = "Annotation confidence",
          x = "Dimension 1",
          y = "Dimension 2",
          color = "Confidence"
        )
    })
    
    composition_plot_obj <- reactive({
      res <- get_res()
      
      ggplot(res$composition, aes(cluster, n, fill = celltype)) +
        geom_col(position = "fill") +
        theme_classic(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
          title = "Cell type composition by cluster",
          x = "Cluster",
          y = "Fraction",
          fill = "Cell type"
        )
    })
    
    output$annotation_umap_ui <- renderUI({
      if (isTRUE(input$interactive_plots)) {
        plotly::plotlyOutput(ns("annotation_umap_plotly"), height = "650px")
      } else {
        plotOutput(ns("annotation_umap"), height = "650px")
      }
    })
    
    output$confidence_umap_ui <- renderUI({
      if (isTRUE(input$interactive_plots)) {
        plotly::plotlyOutput(ns("confidence_umap_plotly"), height = "650px")
      } else {
        plotOutput(ns("confidence_umap"), height = "650px")
      }
    })
    
    output$composition_plot_ui <- renderUI({
      if (isTRUE(input$interactive_plots)) {
        plotly::plotlyOutput(ns("composition_plotly"), height = "650px")
      } else {
        plotOutput(ns("composition_plot"), height = "650px")
      }
    })
    
    output$annotation_umap <- renderPlot({ annotation_umap_obj() })
    output$confidence_umap <- renderPlot({ confidence_umap_obj() })
    output$composition_plot <- renderPlot({ composition_plot_obj() })
    
    if (requireNamespace("plotly", quietly = TRUE)) {
      output$annotation_umap_plotly <- plotly::renderPlotly({
        plotly::ggplotly(annotation_umap_obj(), tooltip = "text")
      })
      
      output$confidence_umap_plotly <- plotly::renderPlotly({
        plotly::ggplotly(confidence_umap_obj(), tooltip = "text")
      })
      
      output$composition_plotly <- plotly::renderPlotly({
        plotly::ggplotly(composition_plot_obj())
      })
    }
    
    output$method_summary <- renderUI({
      res <- get_res()
      
      HTML(paste0(
        "<div class='alert alert-info'>",
        "<b>Method:</b> ", res$method, "<br>",
        "<b>Reduction:</b> ", res$reduction, "<br>",
        "<b>Cluster column:</b> ", res$cluster_col, "<br>",
        "<b>Annotated cells:</b> ", nrow(res$annotation_table), "<br>",
        "<b>Predicted cell types:</b> ", length(unique(res$annotation_table$final_label)), "<br>",
        "<b>Session path:</b> ", res$session_path,
        "</div>"
      ))
    })
    
    output$summary_table <- renderDT({
      res <- get_res()
      
      summary_df <- data.frame(
        metric = c(
          "Method",
          "Reduction",
          "Cluster column",
          "Annotated cells",
          "Predicted cell types",
          "Session path"
        ),
        value = c(
          res$method,
          res$reduction,
          res$cluster_col,
          nrow(res$annotation_table),
          length(unique(res$annotation_table$final_label)),
          res$session_path
        ),
        stringsAsFactors = FALSE
      )
      
      datatable(summary_df, options = list(dom = "t", paging = FALSE), rownames = FALSE)
    })
    
    output$annotation_table <- renderDT({
      res <- get_res()
      datatable(res$annotation_table, options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
    })
    
    output$score_table <- renderDT({
      res <- get_res()
      datatable(res$score_table, options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
    })
    
    save_plot <- function(plot_obj, file, width = 10, height = 8) {
      ext <- tolower(tools::file_ext(file))
      
      if (ext == "svg") {
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
    
    output$download_annotations <- downloadHandler(
      filename = function() paste0("CoTRA_celltype_annotations_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) write.csv(get_res()$annotation_table, file, row.names = FALSE)
    )
    
    output$download_scores <- downloadHandler(
      filename = function() paste0("CoTRA_celltype_scores_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) write.csv(get_res()$score_table, file, row.names = FALSE)
    )
    
    output$download_session <- downloadHandler(
      filename = function() basename(get_res()$session_path),
      content = function(file) file.copy(get_res()$session_path, file, overwrite = TRUE)
    )
    
    output$download_annotation_umap_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_celltype_umap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(annotation_umap_obj(), file)
    )
    
    output$download_annotation_umap_svg <- downloadHandler(
      filename = function() paste0("CoTRA_celltype_umap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(annotation_umap_obj(), file)
    )
    
    output$download_confidence_umap_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_celltype_confidence_umap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(confidence_umap_obj(), file)
    )
    
    output$download_confidence_umap_svg <- downloadHandler(
      filename = function() paste0("CoTRA_celltype_confidence_umap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(confidence_umap_obj(), file)
    )
    
    output$download_composition_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_celltype_composition_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(composition_plot_obj(), file)
    )
    
    output$download_composition_svg <- downloadHandler(
      filename = function() paste0("CoTRA_celltype_composition_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(composition_plot_obj(), file)
    )
    
    output$download_all_zip <- downloadHandler(
      filename = function() paste0("CoTRA_celltype_annotation_outputs_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
      content = function(file) {
        res <- get_res()
        
        tmpdir <- tempfile("CoTRA_celltype_annotation_outputs_")
        dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
        
        write.csv(res$annotation_table, file.path(tmpdir, "CoTRA_celltype_annotations.csv"), row.names = FALSE)
        write.csv(res$score_table, file.path(tmpdir, "CoTRA_celltype_scores.csv"), row.names = FALSE)
        write.csv(res$composition, file.path(tmpdir, "CoTRA_celltype_composition.csv"), row.names = FALSE)
        
        file.copy(res$session_path, file.path(tmpdir, basename(res$session_path)), overwrite = TRUE)
        
        save_plot(annotation_umap_obj(), file.path(tmpdir, "CoTRA_celltype_umap.pdf"))
        save_plot(confidence_umap_obj(), file.path(tmpdir, "CoTRA_celltype_confidence_umap.pdf"))
        save_plot(composition_plot_obj(), file.path(tmpdir, "CoTRA_celltype_composition.pdf"))
        
        if (requireNamespace("svglite", quietly = TRUE)) {
          save_plot(annotation_umap_obj(), file.path(tmpdir, "CoTRA_celltype_umap.svg"))
          save_plot(confidence_umap_obj(), file.path(tmpdir, "CoTRA_celltype_confidence_umap.svg"))
          save_plot(composition_plot_obj(), file.path(tmpdir, "CoTRA_celltype_composition.svg"))
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
      
      annotation_ready = reactive({
        res <- get_res()
        isTRUE(res$annotation_ready)
      }),
      
      annotation_table = reactive({
        res <- get_res()
        res$annotation_table
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