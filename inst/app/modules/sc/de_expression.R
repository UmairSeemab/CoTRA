# ==========================================================
# modules/sc/de_expression.R
# CoTRA scRNA-seq Differential Expression Module
#
# Supports:
# - Replicate-aware pseudobulk differential expression with DESeq2
# - Exploratory cell-level differential expression with Seurat::FindMarkers
# - Condition, sample, cluster, and cell-type metadata selection
# - Downloadable result tables and publication-ready volcano plots
#
# Required import metadata for multi-sample analysis:
# - sample_id
# - condition
# - source_file
# - orig.ident
# ==========================================================

mod_sc_de_ui <- function(id) {
  ns <- NS(id)

  tagList(
    h3("Single-cell differential expression"),
    p(
      "Compare two experimental conditions across all cells or within one selected cluster or cell type. ",
      "Cell-level Wilcoxon testing is the default exploratory method. Use pseudobulk DESeq2 only when each condition has biological replicates."
    ),

    fluidRow(
      column(
        width = 12,
        bs4Card(
          title = "Help: choosing the analysis method",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          tags$ul(
            tags$li(
              tags$b("Pseudobulk DESeq2:"),
              " aggregates raw counts by biological sample and compares conditions using sample-level replication. This is the recommended method for formal inference."
            ),
            tags$li(
              tags$b("Cell-level FindMarkers:"),
              " treats cells as observations. Use it for exploratory analysis when replicate-aware testing is not possible. P-values may be overly optimistic because cells from the same sample are not independent biological replicates."
            ),
            tags$li(
              "For cell-type-specific analysis, select CoTRA_celltype, CoTRA_cluster_label, or seurat_clusters and then choose one value."
            ),
            tags$li(
              "Positive log2 fold change means higher expression in Group 1 relative to Group 2."
            )
          )
        )
      )
    ),

    fluidRow(
      column(
        width = 4,
        bs4Card(
          title = "Analysis settings",
          status = "primary",
          solidHeader = TRUE,
          width = 12,

          selectInput(
            ns("analysis_method"),
            "Analysis method",
            choices = c(
              "Cell-level Seurat FindMarkers, Wilcoxon, exploratory" = "cell_level",
              "Pseudobulk DESeq2, replicate-aware" = "pseudobulk_deseq2"
            ),
            selected = "cell_level"
          ),

          selectInput(ns("condition_col"), "Condition column", choices = character(0)),
          selectInput(ns("sample_col"), "Biological sample column", choices = character(0)),
          selectInput(ns("subset_col"), "Restrict analysis to", choices = c("All cells" = ".none")),
          uiOutput(ns("subset_value_ui")),
          selectInput(ns("group1"), "Group 1, numerator", choices = character(0)),
          selectInput(ns("group2"), "Group 2, reference", choices = character(0)),

          conditionalPanel(
            condition = "input.analysis_method == 'pseudobulk_deseq2'",
            ns = ns,
            numericInput(ns("min_cells_sample"), "Minimum cells per sample", value = 20, min = 1, step = 1),
            numericInput(ns("min_total_count"), "Minimum total gene count", value = 10, min = 0, step = 1),
            numericInput(ns("min_samples_gene"), "Minimum samples expressing a gene", value = 2, min = 1, step = 1)
          ),

          conditionalPanel(
            condition = "input.analysis_method == 'cell_level'",
            ns = ns,
            selectInput(
              ns("cell_test"),
              "FindMarkers test",
              choices = c("Wilcoxon rank-sum" = "wilcox", "MAST" = "MAST"),
              selected = "wilcox"
            ),
            numericInput(ns("min_pct"), "Minimum expression fraction", value = 0.10, min = 0, max = 1, step = 0.01),
            numericInput(ns("logfc_threshold"), "FindMarkers log2FC threshold", value = 0.0, min = 0, step = 0.05)
          ),

          numericInput(ns("fdr_cutoff"), "Adjusted P-value cutoff", value = 0.05, min = 0.0001, max = 1, step = 0.01),
          numericInput(ns("fc_cutoff"), "Absolute log2FC cutoff", value = 0.25, min = 0, step = 0.05),
          actionButton(ns("run_de"), "Run comparison", icon = icon("play"), class = "btn-primary")
        )
      ),

      column(
        width = 8,
        bs4Card(
          title = "Comparison design",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          htmlOutput(ns("design_status")),
          DT::DTOutput(ns("design_table"))
        ),
        bs4Card(
          title = "Analysis status",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          htmlOutput(ns("result_status"))
        )
      )
    ),

    fluidRow(
      column(
        width = 12,
        bs4Card(
          title = "Differential expression results",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          tabsetPanel(
            tabPanel(
              "Result table",
              br(),
              downloadButton(ns("download_results_csv"), "Download CSV"),
              downloadButton(ns("download_results_rds"), "Download RDS"),
              br(), br(),
              DT::DTOutput(ns("de_table"))
            ),
            tabPanel(
              "Volcano plot",
              br(),
              fluidRow(
                column(2, downloadButton(ns("download_volcano_pdf"), "PDF")),
                column(2, downloadButton(ns("download_volcano_svg"), "SVG")),
                column(2, downloadButton(ns("download_volcano_png"), "PNG, 300 dpi"))
              ),
              br(),
              plotly::plotlyOutput(ns("volcano_plot"), height = "650px")
            ),
            tabPanel(
              "Sample cell counts",
              br(),
              plotOutput(ns("sample_count_plot"), height = "500px")
            )
          )
        )
      )
    )
  )
}


mod_sc_de_server <- function(id, seurat_r, sc_state = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    de_results <- reactiveVal(NULL)
    de_info <- reactiveVal(NULL)
    de_design <- reactiveVal(NULL)
    last_error <- reactiveVal(NULL)

    safe_theme <- function() {
      if (exists("theme_cotra", mode = "function")) {
        theme_cotra()
      } else {
        ggplot2::theme_bw(base_size = 12)
      }
    }

    get_seu <- reactive({
      obj <- seurat_r()
      validate(need(inherits(obj, "Seurat"), "A valid Seurat object is required."))
      validate(need(ncol(obj) > 0 && nrow(obj) > 0, "The Seurat object is empty."))
      obj
    })

    clean_character <- function(x) {
      x <- as.character(x)
      x[is.na(x)] <- ""
      trimws(x)
    }

    usable_metadata_columns <- reactive({
      md <- get_seu()@meta.data
      keep <- vapply(md, function(x) {
        ux <- unique(clean_character(x))
        ux <- ux[nzchar(ux)]
        length(ux) >= 2 && length(ux) <= 500
      }, logical(1))
      names(md)[keep]
    })

    observe({
      obj <- get_seu()
      md <- obj@meta.data
      cols <- usable_metadata_columns()

      condition_choices <- cols
      condition_selected <- if ("condition" %in% condition_choices) {
        "condition"
      } else if ("group" %in% condition_choices) {
        "group"
      } else if (length(condition_choices) > 0) {
        condition_choices[1]
      } else {
        character(0)
      }

      sample_choices <- names(md)
      sample_selected <- if ("sample_id" %in% sample_choices) {
        "sample_id"
      } else if ("orig.ident" %in% sample_choices) {
        "orig.ident"
      } else {
        sample_choices[1]
      }

      subset_priority <- c("CoTRA_celltype", "CoTRA_cluster_label", "seurat_clusters")
      subset_choices <- unique(c(subset_priority[subset_priority %in% cols], cols))
      subset_choices <- setdiff(subset_choices, c(condition_selected, sample_selected, "source_file"))
      subset_named <- c("All cells" = ".none", stats::setNames(subset_choices, subset_choices))

      updateSelectInput(session, "condition_col", choices = condition_choices, selected = condition_selected)
      updateSelectInput(session, "sample_col", choices = sample_choices, selected = sample_selected)
      updateSelectInput(session, "subset_col", choices = subset_named, selected = ".none")
    })

    observeEvent(input$condition_col, {
      req(input$condition_col)
      md <- get_seu()@meta.data
      validate(need(input$condition_col %in% colnames(md), "Select a valid condition column."))
      groups <- unique(clean_character(md[[input$condition_col]]))
      groups <- groups[nzchar(groups)]
      groups <- sort(groups)

      selected1 <- if (length(groups) >= 2) groups[2] else groups[1]
      selected2 <- if (length(groups) >= 1) groups[1] else character(0)
      updateSelectInput(session, "group1", choices = groups, selected = selected1)
      updateSelectInput(session, "group2", choices = groups, selected = selected2)
    }, ignoreInit = FALSE)

    output$subset_value_ui <- renderUI({
      req(input$subset_col)
      if (identical(input$subset_col, ".none")) return(NULL)

      md <- get_seu()@meta.data
      validate(need(input$subset_col %in% colnames(md), "Select a valid subset column."))
      values <- unique(clean_character(md[[input$subset_col]]))
      values <- sort(values[nzchar(values)])

      selectInput(
        ns("subset_value"),
        "Selected value",
        choices = values,
        selected = if (length(values) > 0) values[1] else character(0)
      )
    })

    selected_cell_metadata <- reactive({
      obj <- get_seu()
      md <- obj@meta.data
      md$cell <- rownames(md)

      req(input$condition_col, input$sample_col, input$group1, input$group2)
      validate(need(input$condition_col %in% colnames(md), "The selected condition column is missing."))
      validate(need(input$sample_col %in% colnames(md), "The selected sample column is missing."))
      validate(need(!identical(input$group1, input$group2), "Group 1 and Group 2 must be different."))

      condition_values <- clean_character(md[[input$condition_col]])
      sample_values <- clean_character(md[[input$sample_col]])
      keep <- condition_values %in% c(input$group1, input$group2) & nzchar(sample_values)

      if (!is.null(input$subset_col) && !identical(input$subset_col, ".none")) {
        req(input$subset_value)
        validate(need(input$subset_col %in% colnames(md), "The selected subset column is missing."))
        subset_values <- clean_character(md[[input$subset_col]])
        keep <- keep & subset_values == input$subset_value
      }

      out <- md[keep, , drop = FALSE]
      validate(need(nrow(out) > 0, "No cells match the selected comparison."))
      out
    })

    current_design <- reactive({
      md <- selected_cell_metadata()
      sample_values <- clean_character(md[[input$sample_col]])
      condition_values <- clean_character(md[[input$condition_col]])

      split_conditions <- split(condition_values, sample_values)
      inconsistent <- names(Filter(function(x) length(unique(x[nzchar(x)])) > 1, split_conditions))

      design <- aggregate(
        x = list(cells = rep(1L, nrow(md))),
        by = list(sample_id = sample_values, condition = condition_values),
        FUN = sum
      )
      design <- design[order(design$condition, design$sample_id), , drop = FALSE]
      attr(design, "inconsistent_samples") <- inconsistent
      design
    })

    output$design_status <- renderUI({
      design <- current_design()
      inconsistent <- attr(design, "inconsistent_samples")
      group_counts <- table(design$condition)

      text <- paste0(
        "<b>Selected cells:</b> ", sum(design$cells), "<br>",
        "<b>Biological samples:</b> ", nrow(design), "<br>",
        "<b>Samples per condition:</b> ",
        paste(paste(names(group_counts), as.integer(group_counts), sep = " = "), collapse = "; "), "<br>"
      )

      if (length(inconsistent) > 0) {
        text <- paste0(
          text,
          "<span style='color:#b30000;'><b>Invalid sample mapping:</b> ",
          paste(inconsistent, collapse = ", "),
          " map to more than one condition.</span>"
        )
      }

      HTML(text)
    })

    output$design_table <- DT::renderDT({
      design <- current_design()
      DT::datatable(
        design,
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE, dom = "tip")
      )
    })

    get_assay_matrix <- function(obj, assay = "RNA", layer = "counts") {
      mat <- tryCatch(
        SeuratObject::LayerData(object = obj, assay = assay, layer = layer),
        error = function(e) NULL
      )

      if (!is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0) return(mat)

      mat <- tryCatch(
        Seurat::GetAssayData(object = obj, assay = assay, slot = layer),
        error = function(e) NULL
      )

      if (!is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0) return(mat)
      NULL
    }

    prepare_selected_cells <- function() {
      md <- selected_cell_metadata()
      cells <- rownames(md)
      list(metadata = md, cells = cells)
    }

    run_pseudobulk_deseq2 <- function(obj, selected) {
      if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("The DESeq2 package is required for pseudobulk analysis.")
      }
      if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("The Matrix package is required for pseudobulk aggregation.")
      }

      md <- selected$metadata
      cells <- selected$cells
      counts <- get_assay_matrix(obj, assay = "RNA", layer = "counts")
      if (is.null(counts)) {
        stop("A non-empty RNA counts layer is required for pseudobulk analysis.")
      }
      counts <- counts[, cells, drop = FALSE]

      sample_vec <- clean_character(md[[input$sample_col]])
      condition_vec <- clean_character(md[[input$condition_col]])

      condition_by_sample <- split(condition_vec, sample_vec)
      inconsistent <- names(Filter(function(x) length(unique(x[nzchar(x)])) > 1, condition_by_sample))
      if (length(inconsistent) > 0) {
        stop(
          "Each biological sample must map to exactly one condition. Check: ",
          paste(inconsistent, collapse = ", "), "."
        )
      }

      sample_levels <- unique(sample_vec)
      cell_counts <- table(factor(sample_vec, levels = sample_levels))
      keep_samples <- names(cell_counts)[cell_counts >= input$min_cells_sample]
      if (length(keep_samples) == 0) {
        stop("No samples contain the requested minimum number of cells.")
      }

      keep_cells <- sample_vec %in% keep_samples
      counts <- counts[, keep_cells, drop = FALSE]
      md <- md[keep_cells, , drop = FALSE]
      sample_vec <- clean_character(md[[input$sample_col]])
      condition_vec <- clean_character(md[[input$condition_col]])
      sample_levels <- unique(sample_vec)

      sample_condition <- vapply(sample_levels, function(s) {
        unique(condition_vec[sample_vec == s])[1]
      }, character(1))

      samples_per_group <- table(sample_condition)
      missing_groups <- setdiff(c(input$group1, input$group2), names(samples_per_group))
      if (length(missing_groups) > 0) {
        stop("No eligible samples remain for: ", paste(missing_groups, collapse = ", "), ".")
      }
      if (any(samples_per_group[c(input$group1, input$group2)] < 2)) {
        stop(
          "Pseudobulk DESeq2 requires at least two biological samples in each condition after filtering. ",
          "Use additional replicates or choose the exploratory cell-level method."
        )
      }

      membership <- Matrix::sparseMatrix(
        i = seq_along(sample_vec),
        j = match(sample_vec, sample_levels),
        x = 1,
        dims = c(length(sample_vec), length(sample_levels)),
        dimnames = list(colnames(counts), sample_levels)
      )
      pseudobulk_counts <- counts %*% membership

      genes_keep <- Matrix::rowSums(pseudobulk_counts) >= input$min_total_count &
        Matrix::rowSums(pseudobulk_counts > 0) >= input$min_samples_gene
      pseudobulk_counts <- pseudobulk_counts[genes_keep, , drop = FALSE]
      if (nrow(pseudobulk_counts) < 2) {
        stop("Too few genes remain after pseudobulk count filtering.")
      }

      sample_data <- data.frame(
        sample_id = sample_levels,
        condition = factor(sample_condition, levels = c(input$group2, input$group1)),
        cells = as.integer(table(factor(sample_vec, levels = sample_levels))),
        row.names = sample_levels,
        stringsAsFactors = FALSE
      )

      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = round(as.matrix(pseudobulk_counts)),
        colData = sample_data,
        design = ~ condition
      )
      dds <- DESeq2::DESeq(dds, quiet = TRUE)
      res <- DESeq2::results(
        dds,
        contrast = c("condition", input$group1, input$group2),
        alpha = input$fdr_cutoff,
        independentFiltering = TRUE
      )

      result <- as.data.frame(res)
      result$gene <- rownames(result)
      rownames(result) <- NULL

      normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
      group1_samples <- rownames(sample_data)[sample_data$condition == input$group1]
      group2_samples <- rownames(sample_data)[sample_data$condition == input$group2]
      result$mean_group1 <- rowMeans(normalized_counts[result$gene, group1_samples, drop = FALSE])
      result$mean_group2 <- rowMeans(normalized_counts[result$gene, group2_samples, drop = FALSE])

      names(result)[names(result) == "log2FoldChange"] <- "log2FC"
      names(result)[names(result) == "pvalue"] <- "p_value"
      names(result)[names(result) == "padj"] <- "p_adj"

      result$direction <- "Not significant"
      significant <- !is.na(result$p_adj) & result$p_adj <= input$fdr_cutoff &
        !is.na(result$log2FC) & abs(result$log2FC) >= input$fc_cutoff
      result$direction[significant & result$log2FC > 0] <- paste0("Higher in ", input$group1)
      result$direction[significant & result$log2FC < 0] <- paste0("Higher in ", input$group2)

      result <- result[order(result$p_adj, -abs(result$log2FC), na.last = TRUE), , drop = FALSE]
      result <- result[, c(
        "gene", "baseMean", "mean_group1", "mean_group2", "log2FC",
        "lfcSE", "stat", "p_value", "p_adj", "direction"
      )]

      design <- data.frame(
        sample_id = rownames(sample_data),
        condition = as.character(sample_data$condition),
        cells = sample_data$cells,
        included = TRUE,
        stringsAsFactors = FALSE
      )

      list(
        results = result,
        design = design,
        method = "Pseudobulk DESeq2",
        tested_genes = nrow(result),
        selected_cells = ncol(counts),
        samples = nrow(sample_data),
        warning = NULL
      )
    }

    run_cell_level <- function(obj, selected) {
      if (identical(input$cell_test, "MAST") && !requireNamespace("MAST", quietly = TRUE)) {
        stop("The MAST package is required for the selected FindMarkers test.")
      }

      cells <- selected$cells
      md <- selected$metadata
      obj_sub <- subset(obj, cells = cells)
      Seurat::DefaultAssay(obj_sub) <- "RNA"

      data_layer <- get_assay_matrix(obj_sub, assay = "RNA", layer = "data")
      if (is.null(data_layer)) {
        obj_sub <- Seurat::NormalizeData(obj_sub, assay = "RNA", verbose = FALSE)
      }

      condition_values <- clean_character(md[colnames(obj_sub), input$condition_col])
      obj_sub@meta.data$CoTRA_DE_condition <- condition_values
      Seurat::Idents(obj_sub) <- obj_sub@meta.data$CoTRA_DE_condition

      result <- Seurat::FindMarkers(
        object = obj_sub,
        ident.1 = input$group1,
        ident.2 = input$group2,
        assay = "RNA",
        test.use = input$cell_test,
        min.pct = input$min_pct,
        logfc.threshold = input$logfc_threshold,
        only.pos = FALSE,
        verbose = FALSE
      )

      if (is.null(result) || nrow(result) == 0) {
        stop("FindMarkers returned no genes for the selected comparison.")
      }

      result$gene <- rownames(result)
      rownames(result) <- NULL
      if ("avg_log2FC" %in% names(result)) {
        result$log2FC <- result$avg_log2FC
      } else if ("avg_logFC" %in% names(result)) {
        result$log2FC <- result$avg_logFC
      } else {
        stop("FindMarkers output does not contain a recognized log fold-change column.")
      }
      result$p_value <- if ("p_val" %in% names(result)) result$p_val else NA_real_
      result$p_adj <- if ("p_val_adj" %in% names(result)) result$p_val_adj else stats::p.adjust(result$p_value, method = "BH")

      result$direction <- "Not significant"
      significant <- !is.na(result$p_adj) & result$p_adj <= input$fdr_cutoff &
        !is.na(result$log2FC) & abs(result$log2FC) >= input$fc_cutoff
      result$direction[significant & result$log2FC > 0] <- paste0("Higher in ", input$group1)
      result$direction[significant & result$log2FC < 0] <- paste0("Higher in ", input$group2)

      keep_cols <- c("gene", "log2FC", "pct.1", "pct.2", "p_value", "p_adj", "direction")
      keep_cols <- keep_cols[keep_cols %in% colnames(result)]
      result <- result[order(result$p_adj, -abs(result$log2FC), na.last = TRUE), keep_cols, drop = FALSE]

      design <- current_design()
      list(
        results = result,
        design = design,
        method = paste0("Cell-level Seurat FindMarkers, ", input$cell_test),
        tested_genes = nrow(result),
        selected_cells = nrow(md),
        samples = length(unique(clean_character(md[[input$sample_col]]))),
        warning = "Exploratory result: cells from the same biological sample are not independent replicates."
      )
    }

    observeEvent(input$run_de, {
      last_error(NULL)
      de_results(NULL)
      de_info(NULL)
      de_design(NULL)

      outcome <- tryCatch({
        obj <- get_seu()
        selected <- prepare_selected_cells()

        if (identical(input$analysis_method, "pseudobulk_deseq2")) {
          run_pseudobulk_deseq2(obj, selected)
        } else {
          run_cell_level(obj, selected)
        }
      }, error = function(e) e)

      if (inherits(outcome, "error")) {
        last_error(conditionMessage(outcome))
        showModal(modalDialog(
          title = "Differential expression failed",
          conditionMessage(outcome),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
        return(NULL)
      }

      de_results(outcome$results)
      de_design(outcome$design)
      info <- list(
        method = outcome$method,
        condition_column = input$condition_col,
        sample_column = input$sample_col,
        subset_column = input$subset_col,
        subset_value = if (!is.null(input$subset_col) && input$subset_col != ".none") input$subset_value else "All cells",
        group1 = input$group1,
        group2 = input$group2,
        tested_genes = outcome$tested_genes,
        selected_cells = outcome$selected_cells,
        samples = outcome$samples,
        fdr_cutoff = input$fdr_cutoff,
        log2fc_cutoff = input$fc_cutoff,
        warning = outcome$warning,
        run_at = Sys.time()
      )
      de_info(info)

      if (!is.null(sc_state)) {
        try({
          parameters <- sc_state$parameters
          if (is.null(parameters)) parameters <- list()
          parameters$de_expression <- info
          sc_state$parameters <- parameters

          tables <- sc_state$tables
          if (is.null(tables)) tables <- list()
          tables$de_expression <- outcome$results
          sc_state$tables <- tables
        }, silent = TRUE)
      }
    })

    output$result_status <- renderUI({
      if (!is.null(last_error())) {
        return(HTML(paste0("<span style='color:#b30000;'><b>Error:</b> ", last_error(), "</span>")))
      }

      info <- de_info()
      result <- de_results()
      if (is.null(info) || is.null(result)) {
        return(HTML("<span style='color:#777;'>Configure the comparison and click Run comparison.</span>"))
      }

      significant <- sum(!is.na(result$p_adj) & result$p_adj <= info$fdr_cutoff & abs(result$log2FC) >= info$log2fc_cutoff)
      up1 <- sum(result$direction == paste0("Higher in ", info$group1), na.rm = TRUE)
      up2 <- sum(result$direction == paste0("Higher in ", info$group2), na.rm = TRUE)

      warning_text <- if (!is.null(info$warning) && nzchar(info$warning)) {
        paste0("<br><span style='color:#a35a00;'><b>Interpretation:</b> ", info$warning, "</span>")
      } else {
        ""
      }

      HTML(paste0(
        "<b>Method:</b> ", info$method, "<br>",
        "<b>Contrast:</b> ", info$group1, " versus ", info$group2, "<br>",
        "<b>Subset:</b> ", info$subset_value, "<br>",
        "<b>Selected cells:</b> ", info$selected_cells, "<br>",
        "<b>Biological samples:</b> ", info$samples, "<br>",
        "<b>Genes tested:</b> ", info$tested_genes, "<br>",
        "<b>Significant genes:</b> ", significant, " (", up1, " higher in ", info$group1,
        "; ", up2, " higher in ", info$group2, ")",
        warning_text
      ))
    })

    output$de_table <- DT::renderDT({
      result <- de_results()
      validate(need(!is.null(result), "Run differential expression first."))
      DT::datatable(
        result,
        rownames = FALSE,
        filter = "top",
        options = list(pageLength = 20, scrollX = TRUE)
      )
    })

    volcano_plot_data <- reactive({
      result <- de_results()
      info <- de_info()
      validate(need(!is.null(result) && !is.null(info), "Run differential expression first."))

      plot_df <- result
      plot_df$plot_p <- plot_df$p_adj
      plot_df$plot_p[is.na(plot_df$plot_p)] <- 1
      plot_df$plot_p <- pmax(plot_df$plot_p, .Machine$double.xmin)
      plot_df$minus_log10_padj <- -log10(plot_df$plot_p)
      plot_df$label <- ""

      significant_idx <- which(
        !is.na(plot_df$p_adj) &
          plot_df$p_adj <= info$fdr_cutoff &
          abs(plot_df$log2FC) >= info$log2fc_cutoff
      )
      if (length(significant_idx) > 0) {
        label_idx <- head(
          significant_idx[order(plot_df$p_adj[significant_idx], -abs(plot_df$log2FC[significant_idx]))],
          12
        )
        plot_df$label[label_idx] <- plot_df$gene[label_idx]
      }

      plot_df$hover_text <- paste0(
        "Gene: ", plot_df$gene,
        "<br>log2FC: ", signif(plot_df$log2FC, 4),
        "<br>Adjusted P-value: ", format(plot_df$p_adj, digits = 4, scientific = TRUE),
        "<br>Direction: ", plot_df$direction
      )
      plot_df
    })

    volcano_plot_gg <- reactive({
      plot_df <- volcano_plot_data()
      info <- de_info()

      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
          x = log2FC,
          y = minus_log10_padj,
          colour = direction,
          text = hover_text
        )
      ) +
        ggplot2::geom_point(alpha = 0.70, size = 1.6, na.rm = TRUE) +
        ggplot2::geom_vline(
          xintercept = c(-info$log2fc_cutoff, info$log2fc_cutoff),
          linetype = 2,
          inherit.aes = FALSE
        ) +
        ggplot2::geom_hline(
          yintercept = -log10(info$fdr_cutoff),
          linetype = 2,
          inherit.aes = FALSE
        ) +
        ggplot2::labs(
          title = paste0(info$group1, " versus ", info$group2),
          subtitle = info$method,
          x = paste0("log2 fold change, positive = higher in ", info$group1),
          y = "-log10 adjusted P-value",
          colour = "Result"
        ) +
        safe_theme()

      if (requireNamespace("ggrepel", quietly = TRUE) && any(nzchar(plot_df$label))) {
        p <- p + ggrepel::geom_text_repel(
          ggplot2::aes(label = label),
          max.overlaps = 20,
          size = 3,
          show.legend = FALSE
        )
      }
      p
    })

    output$volcano_plot <- plotly::renderPlotly({
      validate(need(requireNamespace("plotly", quietly = TRUE), "Install the plotly package to display the interactive volcano plot."))
      plotly::ggplotly(
        volcano_plot_gg(),
        tooltip = "text",
        dynamicTicks = TRUE
      ) |>
        plotly::layout(
          hovermode = "closest",
          legend = list(orientation = "h", x = 0, y = -0.18),
          margin = list(b = 100)
        ) |>
        plotly::config(
          displaylogo = FALSE,
          responsive = TRUE,
          toImageButtonOptions = list(
            format = "png",
            filename = "CoTRA_scRNA_DE_volcano",
            scale = 2
          )
        )
    })

    sample_count_plot_gg <- reactive({
      design <- de_design()
      if (is.null(design)) design <- current_design()
      validate(need(!is.null(design) && nrow(design) > 0, "No sample design is available."))

      ggplot2::ggplot(
        design,
        ggplot2::aes(x = sample_id, y = cells, fill = condition)
      ) +
        ggplot2::geom_col() +
        ggplot2::labs(
          title = "Cells contributing to the comparison",
          x = "Biological sample",
          y = "Number of cells",
          fill = "Condition"
        ) +
        safe_theme() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    })

    output$sample_count_plot <- renderPlot({
      sample_count_plot_gg()
    }, res = 120)

    save_plot <- function(plot_obj, file, format, width = 8, height = 7) {
      if (identical(format, "pdf")) {
        grDevices::pdf(file, width = width, height = height, onefile = TRUE)
        on.exit(grDevices::dev.off(), add = TRUE)
        print(plot_obj)
      } else if (identical(format, "svg")) {
        if (requireNamespace("svglite", quietly = TRUE)) {
          svglite::svglite(file, width = width, height = height)
        } else {
          grDevices::svg(file, width = width, height = height)
        }
        on.exit(grDevices::dev.off(), add = TRUE)
        print(plot_obj)
      } else if (identical(format, "png")) {
        grDevices::png(file, width = width, height = height, units = "in", res = 300)
        on.exit(grDevices::dev.off(), add = TRUE)
        print(plot_obj)
      } else {
        stop("Unsupported plot format: ", format)
      }
    }

    timestamp_name <- function(prefix, extension) {
      paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", extension)
    }

    output$download_results_csv <- downloadHandler(
      filename = function() timestamp_name("CoTRA_scRNA_differential_expression", "csv"),
      content = function(file) {
        result <- de_results()
        req(result)
        utils::write.csv(result, file, row.names = FALSE)
      }
    )

    output$download_results_rds <- downloadHandler(
      filename = function() timestamp_name("CoTRA_scRNA_differential_expression", "rds"),
      content = function(file) {
        result <- de_results()
        info <- de_info()
        design <- de_design()
        req(result, info)
        saveRDS(list(results = result, design = design, analysis = info), file)
      }
    )

    output$download_volcano_pdf <- downloadHandler(
      filename = function() timestamp_name("CoTRA_scRNA_DE_volcano", "pdf"),
      content = function(file) save_plot(volcano_plot_gg(), file, "pdf")
    )

    output$download_volcano_svg <- downloadHandler(
      filename = function() timestamp_name("CoTRA_scRNA_DE_volcano", "svg"),
      content = function(file) save_plot(volcano_plot_gg(), file, "svg")
    )

    output$download_volcano_png <- downloadHandler(
      filename = function() timestamp_name("CoTRA_scRNA_DE_volcano", "png"),
      content = function(file) save_plot(volcano_plot_gg(), file, "png")
    )

    return(list(
      seurat = reactive(get_seu()),
      results = reactive(de_results()),
      design = reactive(de_design()),
      analysis_info = reactive(de_info()),
      volcano_plot = reactive(volcano_plot_gg()),
      volcano_data = reactive(volcano_plot_data())
    ))
  })
}
