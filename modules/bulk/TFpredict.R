# ==========================================================
# modules/bulk/TFpredict.R
# CoTRA Bulk RNA-seq Transcription Factor Prediction Module
#
# Input:
# - de_results: output reactive from mod_bulk_de_server()
# Optional:
# - bulk_reactive: uploaded or normalized expression matrix for gene count plots
# - groups_reactive: groups from mod_bulk_groups_server()
#
# Main features:
# - Uses DE results from CoTRA de_analysis.R
# - Supports Human, Mouse, Rat
# - Predicts enriched TF motifs using RcisTarget::primo() if available
# - Falls back to motif enrichment table display only when primo output exists
# - Maps SYMBOL or ENSEMBL gene IDs to ENTREZID
# - Runs separately for Up, Down, and All significant DE genes
# - Shows predicted TF motif table and TFs present in the DE table
# - Plots selected TF expression across groups when expression matrix is available
# - Downloads CSV/XLSX tables
# - Downloads gene expression plot as SVG/PDF
# ==========================================================

mod_bulk_tfpredict_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    shinyjs::useShinyjs(),
    
    actionButton(ns("help_toggle"), icon("question-circle"), class = "btn btn-info btn-sm"),
    
    shinyjs::hidden(
      div(
        id = ns("help_box"),
        style = "border:1px solid #d9d9d9; padding:12px; margin-top:10px; margin-bottom:10px; border-radius:6px; background:#f7f7f7;",
        h4("Step-by-Step Guide"),
        
        h5("1. Run differential expression first"),
        p("This module uses the DE table from CoTRA. The table should contain gene, padj, and lfc or log2FoldChange columns."),
        
        h5("2. Select organism"),
        p("Choose Human, Mouse, or Rat. Gene IDs are mapped to SYMBOL and ENTREZID using org.*.eg.db packages."),
        
        h5("3. Set DE thresholds"),
        p("TF prediction is performed separately for up-regulated, down-regulated, and all significant genes using the selected padj and log2 fold-change cutoffs."),
        
        h5("4. Run TF prediction"),
        p("Click Run TF Prediction. The module uses primo() to detect over-represented transcription factor motifs from Entrez IDs."),
        
        h5("5. Explore results"),
        p("The summary table shows how many predicted TFs are present in the DE table. Click any cell to view detailed predicted TF motifs and matching TF genes."),
        
        h5("6. Plot TF expression"),
        p("Select a TF gene from the matching TF table to view expression across sample groups when expression data and groups are available."),
        
        h5("7. Download outputs"),
        p("Download predicted TF motif tables, matching TF gene tables, and gene expression plots for reporting.")
      )
    ),
    
    br(),
    
    fluidRow(
      box(
        width = 12,
        title = "TF prediction settings",
        solidHeader = TRUE,
        status = "primary",
        collapsible = TRUE,
        
        fluidRow(
          column(3, selectInput(ns("organism"), "Organism", choices = c("Human" = "human", "Mouse" = "mouse", "Rat" = "rat"), selected = "mouse")),
          column(3, selectInput(ns("gene_col"), "Gene ID column", choices = NULL)),
          column(3, numericInput(ns("padj_cutoff"), "Padj cutoff", value = 0.05, min = 0, max = 1, step = 0.005)),
          column(3, numericInput(ns("logfc_cutoff"), "Abs log2FC cutoff", value = 1, min = 0, step = 0.1))
        ),
        
        fluidRow(
          column(3, checkboxInput(ns("show_lfc"), "Show log2 fold-change columns", value = TRUE)),
          column(3, selectInput(ns("download_type"), "Download type", choices = c("CSV" = "csv", "Excel" = "xlsx"), selected = "csv")),
          column(3, actionButton(ns("run_tf"), "Run TF Prediction", class = "btn btn-primary")),
          column(3, downloadButton(ns("download_summary"), "Download Summary", class = "btn btn-success"))
        ),
        
        br(),
        verbatimTextOutput(ns("status"))
      )
    ),
    
    fluidRow(
      box(
        width = 12,
        title = "Predicted TF summary",
        solidHeader = TRUE,
        status = "primary",
        DT::DTOutput(ns("tf_summary"))
      )
    ),
    
    fluidRow(
      box(
        width = 12,
        title = "Detailed predicted TF motifs",
        solidHeader = TRUE,
        status = "primary",
        fluidRow(
          column(6, downloadButton(ns("download_predicted_motifs"), "Download Predicted Motifs")),
          column(6, downloadButton(ns("download_matching_tfs"), "Download Matching TF Genes"))
        ),
        br(),
        textOutput(ns("selected_text")),
        br(),
        DT::DTOutput(ns("predicted_motif_table")),
        br(),
        DT::DTOutput(ns("matching_tf_table"))
      )
    ),
    
    fluidRow(
      box(
        width = 12,
        title = "Selected TF expression",
        solidHeader = TRUE,
        status = "primary",
        fluidRow(
          column(3, radioButtons(ns("log_scale"), "Use log scale", choices = c("Yes", "No"), selected = "Yes", inline = TRUE)),
          column(3, downloadButton(ns("download_gene_svg"), "Download SVG")),
          column(3, downloadButton(ns("download_gene_pdf"), "Download PDF"))
        ),
        br(),
        plotOutput(ns("gene_plot"), height = "500px")
      )
    )
  )
}

mod_bulk_tfpredict_server <- function(id,
                                      de_results,
                                      bulk_reactive = reactive(NULL),
                                      groups_reactive = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    
    rv <- reactiveValues(
      tf_result = NULL,
      selected_set = NULL,
      selected_col = NULL,
      selected_gene = NULL
    )
    
    observeEvent(input$help_toggle, {
      shinyjs::toggle(id = "help_box", anim = TRUE)
    })
    
    clean_ids <- function(x) {
      x <- trimws(as.character(x))
      sub("\\.\\d+$", "", x)
    }
    
    first_non_na <- function(x) {
      x <- as.character(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (length(x) == 0) return(NA_character_)
      x[1]
    }
    
    safe_filename <- function(x) {
      x <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
      gsub("_+", "_", x)
    }
    
    save_table_file <- function(df, file, type, sheet_name = "Sheet1") {
      if (identical(type, "xlsx")) {
        wb <- openxlsx::createWorkbook()
        openxlsx::addWorksheet(wb, sheet_name)
        openxlsx::writeData(wb, sheet = 1, x = df, colNames = TRUE, rowNames = FALSE)
        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      } else {
        utils::write.csv(df, file, row.names = FALSE)
      }
    }
    
    orgdb_obj <- reactive({
      switch(
        input$organism,
        "mouse" = org.Mm.eg.db::org.Mm.eg.db,
        "human" = org.Hs.eg.db::org.Hs.eg.db,
        "rat" = org.Rn.eg.db::org.Rn.eg.db,
        org.Mm.eg.db::org.Mm.eg.db
      )
    })
    
    primo_org <- reactive({
      switch(input$organism, "mouse" = "Mm", "human" = "Hs", "rat" = "Rn", "Mm")
    })
    
    guess_id_type <- function(ids) {
      ids <- clean_ids(ids)
      patt <- switch(
        input$organism,
        "mouse" = "^ENSMUSG[0-9]+$",
        "human" = "^ENSG[0-9]+$",
        "rat" = "^ENSRNOG[0-9]+$",
        "^ENS"
      )
      frac_ens <- mean(grepl(patt, ids), na.rm = TRUE)
      if (!is.na(frac_ens) && frac_ens > 0.5) "ENSEMBL" else "SYMBOL"
    }
    
    de_df <- reactive({
      df <- de_results()
      req(df)
      validate(need(is.data.frame(df), "DE result must be a data frame."))
      validate(need(nrow(df) > 0, "DE result is empty."))
      df
    })
    
    observeEvent(de_df(), {
      cols <- colnames(de_df())
      gene_default <- if ("gene" %in% cols) "gene" else cols[1]
      updateSelectInput(session, "gene_col", choices = cols, selected = gene_default)
    }, ignoreInit = FALSE)
    
    lfc_col <- reactive({
      df <- de_df()
      if ("lfc" %in% colnames(df)) return("lfc")
      if ("log2FoldChange" %in% colnames(df)) return("log2FoldChange")
      validate(need(FALSE, "No log fold-change column found. Expected lfc or log2FoldChange."))
    })
    
    map_de_ids <- reactive({
      df <- de_df()
      ids <- clean_ids(df[[input$gene_col]])
      validate(need(length(ids) > 0, "No gene IDs available."))
      
      keytype <- guess_id_type(ids)
      db <- orgdb_obj()
      
      mp <- suppressWarnings(
        AnnotationDbi::select(
          db,
          keys = unique(ids),
          keytype = keytype,
          columns = c("SYMBOL", "ENTREZID", "ENSEMBL")
        )
      )
      
      validate(need(!is.null(mp) && nrow(mp) > 0, "No gene IDs mapped."))
      
      mp <- as.data.frame(mp, stringsAsFactors = FALSE)
      colnames(mp)[colnames(mp) == keytype] <- "input_id"
      mp$input_id <- clean_ids(mp$input_id)
      
      if (!"SYMBOL" %in% colnames(mp)) mp$SYMBOL <- NA_character_
      if (!"ENTREZID" %in% colnames(mp)) mp$ENTREZID <- NA_character_
      if (!"ENSEMBL" %in% colnames(mp)) mp$ENSEMBL <- NA_character_
      
      mp <- mp[!is.na(mp$input_id) & nzchar(mp$input_id), , drop = FALSE]
      
      mp <- mp %>%
        dplyr::group_by(.data$input_id) %>%
        dplyr::summarise(
          SYMBOL = first_non_na(.data$SYMBOL),
          ENTREZID = first_non_na(.data$ENTREZID),
          ENSEMBL = first_non_na(.data$ENSEMBL),
          .groups = "drop"
        )
      
      mp
    })
    
    de_mapped <- reactive({
      df <- de_df()
      df$input_id <- clean_ids(df[[input$gene_col]])
      out <- dplyr::left_join(df, map_de_ids(), by = "input_id")
      out$report_lfc <- suppressWarnings(as.numeric(out[[lfc_col()]]))
      out$report_padj <- suppressWarnings(as.numeric(out$padj))
      out
    })
    
    split_gene_sets <- reactive({
      df <- de_mapped()
      validate(need("padj" %in% colnames(df), "No padj column found in DE table."))
      
      base <- !is.na(df$report_padj) &
        df$report_padj <= input$padj_cutoff &
        !is.na(df$report_lfc) &
        abs(df$report_lfc) >= input$logfc_cutoff
      
      up <- base & df$report_lfc >= input$logfc_cutoff
      down <- base & df$report_lfc <= -input$logfc_cutoff
      
      list(
        Up = df[up, , drop = FALSE],
        Down = df[down, , drop = FALSE],
        All = df[base, , drop = FALSE]
      )
    })
    
    parse_primo_tf_symbols <- function(tf_table) {
      if (is.null(tf_table) || nrow(as.data.frame(tf_table)) == 0) return(character(0))
      tf_table <- as.data.frame(tf_table, stringsAsFactors = FALSE)
      
      motif_col <- if (ncol(tf_table) >= 4) 4 else ncol(tf_table)
      raw <- as.character(tf_table[[motif_col]])
      raw <- raw[!is.na(raw) & nzchar(raw)]
      if (length(raw) == 0) return(character(0))
      
      x <- unlist(strsplit(raw, "_|::|;|,|/|\\s+"))
      x <- trimws(x)
      x <- x[nzchar(x)]
      x <- gsub("[^A-Za-z0-9.-]", "", x)
      x <- x[nzchar(x)]
      
      if (identical(input$organism, "human")) {
        x <- toupper(x)
      } else {
        x <- stringr::str_to_title(tolower(x))
      }
      
      unique(x)
    }
    
    run_primo_safe <- function(entrez_ids) {
      entrez_ids <- unique(as.character(entrez_ids))
      entrez_ids <- entrez_ids[!is.na(entrez_ids) & nzchar(entrez_ids)]
      if (length(entrez_ids) < 3) {
        return(list(table = data.frame(), tfs = character(0), message = "Too few Entrez IDs."))
      }
      
      if (!exists("primo", mode = "function")) {
        return(list(table = data.frame(), tfs = character(0), message = "Function primo() is not available. Load the package or source file that provides primo()."))
      }
      
      res <- tryCatch(
        primo(entrez_ids, inputType = "entrezID", org = primo_org()),
        error = function(e) e
      )
      
      if (inherits(res, "error")) {
        return(list(table = data.frame(), tfs = character(0), message = conditionMessage(res)))
      }
      
      tf_table <- tryCatch(as.data.frame(res$overRepresented, stringsAsFactors = FALSE), error = function(e) data.frame())
      tfs <- parse_primo_tf_symbols(tf_table)
      
      list(table = tf_table, tfs = tfs, message = "OK")
    }
    
    predicted_tf <- eventReactive(input$run_tf, {
      sets <- split_gene_sets()
      de_all <- de_mapped()
      
      out <- list()
      summary <- data.frame(
        Gene_Set = character(0),
        Input_DE_genes = integer(0),
        Entrez_mapped = integer(0),
        Predicted_motif_rows = integer(0),
        Predicted_TF_symbols = integer(0),
        TFs_present_in_DE_table = integer(0),
        Status = character(0),
        stringsAsFactors = FALSE
      )
      
      for (set_name in names(sets)) {
        df <- sets[[set_name]]
        entrez <- unique(as.character(df$ENTREZID))
        entrez <- entrez[!is.na(entrez) & nzchar(entrez)]
        
        pr <- run_primo_safe(entrez)
        predicted_symbols <- pr$tfs
        
        matching <- de_all[!is.na(de_all$SYMBOL) & de_all$SYMBOL %in% predicted_symbols, , drop = FALSE]
        matching <- matching[!duplicated(matching$SYMBOL), , drop = FALSE]
        
        motif_table <- as.data.frame(pr$table, stringsAsFactors = FALSE)
        if (nrow(motif_table) > 0) {
          colnames(motif_table) <- make.names(colnames(motif_table), unique = TRUE)
          motif_table$Predicted_TF_symbols_parsed <- paste(predicted_symbols, collapse = ";")
        }
        
        out[[set_name]] <- list(
          input_de = df,
          entrez = entrez,
          motif_table = motif_table,
          predicted_symbols = predicted_symbols,
          matching_tfs = matching,
          status = pr$message
        )
        
        summary <- rbind(
          summary,
          data.frame(
            Gene_Set = set_name,
            Input_DE_genes = nrow(df),
            Entrez_mapped = length(entrez),
            Predicted_motif_rows = nrow(motif_table),
            Predicted_TF_symbols = length(predicted_symbols),
            TFs_present_in_DE_table = nrow(matching),
            Status = pr$message,
            stringsAsFactors = FALSE
          )
        )
      }
      
      list(summary = summary, result = out)
    })
    
    output$status <- renderPrint({
      df <- tryCatch(de_mapped(), error = function(e) NULL)
      if (is.null(df)) {
        cat("No DE table available. Run DE analysis first.\n")
      } else {
        cat("DE rows:", nrow(df), "\n")
        cat("Mapped SYMBOL:", sum(!is.na(df$SYMBOL) & nzchar(df$SYMBOL)), "\n")
        cat("Mapped ENTREZID:", sum(!is.na(df$ENTREZID) & nzchar(df$ENTREZID)), "\n")
        cat("Selected organism:", input$organism, "\n")
        cat("primo() available:", exists("primo", mode = "function"), "\n")
        if (input$run_tf > 0 && !is.null(rv$tf_result)) cat("TF prediction completed.\n")
      }
    })
    
    observeEvent(input$run_tf, {
      rv$tf_result <- predicted_tf()
    })
    
    output$tf_summary <- DT::renderDT({
      req(rv$tf_result)
      DT::datatable(
        rv$tf_result$summary,
        selection = list(mode = "single", target = "cell"),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    })
    
    selected_gene_set <- reactive({
      req(rv$tf_result)
      selected <- input$tf_summary_cells_selected
      if (is.null(selected) || length(selected) == 0) return(rv$tf_result$summary$Gene_Set[1])
      row <- selected[1]
      rv$tf_result$summary$Gene_Set[row]
    })
    
    selected_result <- reactive({
      req(rv$tf_result)
      set_name <- selected_gene_set()
      rv$tf_result$result[[set_name]]
    })
    
    output$selected_text <- renderText({
      req(selected_result())
      res <- selected_result()
      paste0(
        selected_gene_set(), ": ",
        nrow(res$matching_tfs), " TF genes from predicted motifs are present in the DE table. ",
        "Status: ", res$status
      )
    })
    
    output$predicted_motif_table <- DT::renderDT({
      res <- selected_result()
      df <- res$motif_table
      validate(need(!is.null(df) && nrow(df) > 0, "No predicted motif table available."))
      DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
    })
    
    matching_tf_table <- reactive({
      res <- selected_result()
      df <- res$matching_tfs
      if (is.null(df) || nrow(df) == 0) return(data.frame())
      
      keep_cols <- colnames(df)
      if (!isTRUE(input$show_lfc)) {
        keep_cols <- keep_cols[!grepl("lfc|log2FoldChange", keep_cols, ignore.case = TRUE)]
      }
      
      df <- df[, keep_cols, drop = FALSE]
      df <- df[order(df$report_padj, na.last = TRUE), , drop = FALSE]
      df
    })
    
    output$matching_tf_table <- DT::renderDT({
      df <- matching_tf_table()
      validate(need(nrow(df) > 0, "No predicted TF genes are present in the DE table."))
      DT::datatable(
        df,
        selection = list(mode = "single", target = "row"),
        options = list(pageLength = 15, scrollX = TRUE),
        rownames = FALSE
      )
    })
    
    observeEvent(input$matching_tf_table_rows_selected, {
      df <- matching_tf_table()
      idx <- input$matching_tf_table_rows_selected
      if (!is.null(idx) && length(idx) > 0 && nrow(df) >= idx[1]) {
        if ("SYMBOL" %in% colnames(df) && !is.na(df$SYMBOL[idx[1]]) && nzchar(df$SYMBOL[idx[1]])) {
          rv$selected_gene <- df$SYMBOL[idx[1]]
        } else {
          rv$selected_gene <- as.character(df[[input$gene_col]][idx[1]])
        }
      }
    })
    
    get_matrix <- function() {
      x <- tryCatch({
        if (is.function(bulk_reactive)) bulk_reactive() else bulk_reactive
      }, error = function(e) NULL)
      if (is.null(x)) return(NULL)
      
      if (is.matrix(x) || is.data.frame(x)) return(as.data.frame(x, check.names = FALSE))
      
      if (is.list(x)) {
        candidates <- c("final", "norm", "normalized", "normalized_counts", "counts", "raw", "data", "matrix", "expr", "expression")
        for (nm in candidates) {
          if (!is.null(x[[nm]])) {
            y <- tryCatch(if (is.function(x[[nm]])) x[[nm]]() else x[[nm]], error = function(e) NULL)
            if (is.matrix(y) || is.data.frame(y)) return(as.data.frame(y, check.names = FALSE))
          }
        }
      }
      
      NULL
    }
    
    clean_matrix <- function(mat) {
      if (is.null(mat)) return(NULL)
      df <- as.data.frame(mat, check.names = FALSE)
      if (ncol(df) > 1) {
        first_col <- df[[1]]
        if (!is.numeric(first_col) && anyDuplicated(as.character(first_col)) == 0) {
          rownames(df) <- as.character(first_col)
          df <- df[, -1, drop = FALSE]
        }
      }
      mat2 <- as.matrix(df)
      suppressWarnings(storage.mode(mat2) <- "numeric")
      mat2[is.na(mat2)] <- 0
      mat2
    }
    
    get_group_df <- function() {
      x <- tryCatch({
        if (is.function(groups_reactive)) groups_reactive() else groups_reactive
      }, error = function(e) NULL)
      if (is.null(x)) return(NULL)
      if (is.data.frame(x)) return(x)
      if (is.list(x)) {
        candidates <- c("groups", "sample_info", "metadata", "meta", "coldata", "group_table")
        for (nm in candidates) {
          if (!is.null(x[[nm]])) {
            y <- tryCatch(if (is.function(x[[nm]])) x[[nm]]() else x[[nm]], error = function(e) NULL)
            if (is.data.frame(y)) return(y)
          }
        }
      }
      NULL
    }
    
    sample_groups <- function(samples) {
      out <- data.frame(Sample = samples, Group = "Unassigned", stringsAsFactors = FALSE)
      gdf <- get_group_df()
      if (is.null(gdf) || nrow(gdf) == 0) return(out)
      
      if (all(c("Group", "Samples") %in% colnames(gdf))) {
        for (i in seq_len(nrow(gdf))) {
          smp <- trimws(unlist(strsplit(as.character(gdf$Samples[i]), ";|,|\\s+")))
          smp <- smp[nzchar(smp)]
          out$Group[out$Sample %in% smp] <- as.character(gdf$Group[i])
        }
        return(out)
      }
      
      sample_col <- grep("sample|Sample|sample_id|SampleID", colnames(gdf), value = TRUE)[1]
      group_col <- grep("group|Group|condition|Condition", colnames(gdf), value = TRUE)[1]
      if (!is.na(sample_col) && !is.na(group_col)) {
        idx <- match(out$Sample, as.character(gdf[[sample_col]]))
        out$Group[!is.na(idx)] <- as.character(gdf[[group_col]])[idx[!is.na(idx)]]
      }
      out
    }
    
    gene_plot_obj <- reactive({
      req(rv$selected_gene)
      mat <- clean_matrix(get_matrix())
      validate(need(!is.null(mat), "No expression matrix available for plotting."))
      
      gene <- rv$selected_gene
      rn <- rownames(mat)
      idx <- which(rn == gene)
      if (length(idx) == 0) {
        idx <- which(toupper(rn) == toupper(gene))
      }
      validate(need(length(idx) > 0, paste("Selected gene not found in expression matrix:", gene)))
      
      counts <- as.numeric(mat[idx[1], ])
      df <- data.frame(Sample = colnames(mat), Count = counts, stringsAsFactors = FALSE)
      df <- dplyr::left_join(df, sample_groups(df$Sample), by = "Sample")
      
      p <- ggplot2::ggplot(df, ggplot2::aes(x = Group, y = Count, fill = Group)) +
        ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.15, height = 0), size = 2) +
        ggplot2::theme_bw() +
        ggplot2::labs(title = paste0("Gene expression: ", gene), x = "Group", y = "Expression") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      
      if (identical(input$log_scale, "Yes")) {
        p <- p + ggplot2::scale_y_continuous(trans = "log2")
      }
      
      p
    })
    
    output$gene_plot <- renderPlot({
      req(rv$selected_gene)
      print(gene_plot_obj())
    })
    
    output$download_gene_svg <- downloadHandler(
      filename = function() paste0("TF_expression_", safe_filename(rv$selected_gene), ".svg"),
      content = function(file) {
        svglite::svglite(file, width = 8, height = 6)
        print(gene_plot_obj())
        grDevices::dev.off()
      }
    )
    
    output$download_gene_pdf <- downloadHandler(
      filename = function() paste0("TF_expression_", safe_filename(rv$selected_gene), ".pdf"),
      content = function(file) {
        grDevices::pdf(file, width = 8, height = 6)
        print(gene_plot_obj())
        grDevices::dev.off()
      }
    )
    
    output$download_summary <- downloadHandler(
      filename = function() if (input$download_type == "xlsx") "TF_prediction_summary.xlsx" else "TF_prediction_summary.csv",
      content = function(file) {
        req(rv$tf_result)
        save_table_file(rv$tf_result$summary, file, input$download_type, "TF_summary")
      }
    )
    
    output$download_predicted_motifs <- downloadHandler(
      filename = function() {
        suffix <- if (input$download_type == "xlsx") "xlsx" else "csv"
        paste0("TF_predicted_motifs_", safe_filename(selected_gene_set()), ".", suffix)
      },
      content = function(file) {
        res <- selected_result()
        df <- res$motif_table
        save_table_file(df, file, input$download_type, "Predicted_motifs")
      }
    )
    
    output$download_matching_tfs <- downloadHandler(
      filename = function() {
        suffix <- if (input$download_type == "xlsx") "xlsx" else "csv"
        paste0("TF_matching_genes_", safe_filename(selected_gene_set()), ".", suffix)
      },
      content = function(file) {
        save_table_file(matching_tf_table(), file, input$download_type, "Matching_TFs")
      }
    )
    
    return(list(
      summary = reactive({
        req(rv$tf_result)
        rv$tf_result$summary
      }),
      result = reactive({
        req(rv$tf_result)
        rv$tf_result$result
      }),
      selected_gene = reactive(rv$selected_gene)
    ))
  })
}

mod_bulk_TFpredict_ui <- function(id) mod_bulk_tfpredict_ui(id)
mod_bulk_TFpredict_server <- function(id, de_results, bulk_reactive = reactive(NULL), groups_reactive = reactive(NULL)) {
  mod_bulk_tfpredict_server(id, de_results, bulk_reactive, groups_reactive)
}
