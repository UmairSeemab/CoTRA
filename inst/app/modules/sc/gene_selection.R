# ==========================================================
# modules/sc/gene_selection.R
# CoTRA scRNA-seq Gene Selection Module
#
# Features:
# - HVG selection using FindVariableFeatures or SCTransform
# - Custom pasted gene list
# - Uploaded gene list from CSV/TXT/TSV
# - Combine or replace selected genes
# - Seurat v4/v5 safe assay access
# - Missing gene checks
# - Gene summary table
# - Collapsible help and interpretation sections
# - Session export to outputs/scRNA/sessioninfo/
# - Downstream-safe return object for PCA and later modules
# ==========================================================

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

mod_sc_gene_selection_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("Gene selection"),
    br(),
    
    bs4Dash::bs4Card(
      title = "Help: how to use gene selection",
      width = 12,
      solidHeader = TRUE,
      status = "info",
      collapsible = TRUE,
      collapsed = TRUE,
      icon = icon("circle-info"),
      tags$h4("Purpose"),
      tags$p("This module selects genes for PCA and downstream analysis. You can use highly variable genes, custom genes, uploaded gene lists, or a combination of these."),
      tags$h4("Available gene selection modes"),
      tags$ul(
        tags$li(tags$b("Highly variable genes:"), " recommended default for unsupervised PCA, UMAP, clustering, and marker discovery."),
        tags$li(tags$b("Pasted custom genes:"), " useful when you want to focus on known markers, pathways, or a curated biological signature."),
        tags$li(tags$b("Uploaded gene list:"), " useful for reusable gene panels from CSV, TXT, or TSV files."),
        tags$li(tags$b("Combine selected genes:"), " keeps HVGs and adds valid custom or uploaded genes."),
        tags$li(tags$b("Replace with custom/uploaded genes:"), " uses only the valid custom or uploaded genes for downstream PCA.")
      ),
      tags$h4("Recommended workflow"),
      tags$ol(
        tags$li("Run QC first."),
        tags$li("Use HVG mode for the first unsupervised analysis."),
        tags$li("Add custom or uploaded genes only when you have a clear biological reason."),
        tags$li("Check the missing genes table before continuing."),
        tags$li("Continue to PCA using the final selected genes.")
      ),
      tags$h4("Important caution"),
      tags$p("Using only a small curated gene set can bias PCA, UMAP, clustering, and biological interpretation. For unbiased discovery, start with HVGs.")
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "Gene selection settings",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        
        fluidRow(
          column(
            4,
            selectInput(
              ns("selection_mode"),
              "Gene selection mode",
              choices = c(
                "HVGs only" = "hvg_only",
                "Custom genes only" = "custom_only",
                "Uploaded genes only" = "uploaded_only",
                "HVGs plus custom genes" = "hvg_custom",
                "HVGs plus uploaded genes" = "hvg_uploaded",
                "HVGs plus custom and uploaded genes" = "hvg_custom_uploaded"
              ),
              selected = "hvg_only"
            )
          ),
          column(
            4,
            selectInput(
              ns("assay_use"),
              "Assay for gene selection",
              choices = c("RNA"),
              selected = "RNA"
            )
          ),
          column(
            4,
            selectInput(
              ns("hvg_mode"),
              "HVG method",
              choices = c(
                "FindVariableFeatures" = "fvf",
                "SCTransform" = "sct"
              ),
              selected = "fvf"
            )
          )
        ),
        
        fluidRow(
          column(
            4,
            numericInput(
              ns("nfeatures"),
              "Number of HVGs",
              value = 2000,
              min = 100,
              max = 10000,
              step = 100
            )
          ),
          column(
            4,
            selectInput(
              ns("gene_case_mode"),
              "Gene name matching",
              choices = c(
                "Exact match" = "exact",
                "Case-insensitive match" = "case_insensitive"
              ),
              selected = "exact"
            )
          ),
          column(
            4,
            br(),
            actionButton(ns("run_gene_selection"), "Run gene selection", class = "btn-primary")
          )
        ),
        
        hr(),
        h4("Optional filters for final selected genes"),
        fluidRow(
          column(4, checkboxInput(ns("filter_mito"), "Remove mitochondrial genes", FALSE)),
          column(4, checkboxInput(ns("filter_ribo"), "Remove ribosomal genes", FALSE)),
          column(4, textInput(ns("filter_regex"), "Remove genes matching regex", value = ""))
        ),
        
        hr(),
        h4("Custom gene list"),
        fluidRow(
          column(
            6,
            textAreaInput(
              ns("custom_genes"),
              "Paste genes, one per line or separated by comma/space",
              value = "",
              rows = 8,
              placeholder = "Example:\nRHO\nPDE6B\nGNAT1"
            )
          ),
          column(
            6,
            fileInput(
              ns("gene_file"),
              "Upload gene list, optional",
              accept = c(".csv", ".txt", ".tsv")
            ),
            helpText("The first column will be used as gene symbols unless a column named gene, genes, symbol, or feature exists."),
            checkboxInput(ns("use_uploaded_header"), "Uploaded file has header", TRUE)
          )
        ),
        hr(),
        fluidRow(
          column(4, downloadButton(ns("download_gene_selection_session"), "Save gene selection session (.rds)", class = "btn-success")),
          column(8, uiOutput(ns("session_save_hint")))
        )
      )
    ),
    
    bs4Dash::bs4Card(
      title = "Interpretation: gene selection choices",
      width = 12,
      solidHeader = TRUE,
      status = "info",
      collapsible = TRUE,
      collapsed = TRUE,
      icon = icon("circle-info"),
      tags$h4("How to interpret HVGs"),
      tags$p("Highly variable genes are genes with stronger cell-to-cell variation than expected from their average expression. They are useful because they capture the major structure of the dataset."),
      tags$h4("When custom genes are useful"),
      tags$ul(
        tags$li("Use marker panels when you want to inspect known cell types."),
        tags$li("Use pathway genes when you want to focus on a specific biological process."),
        tags$li("Use disease genes when you want to explore disease-associated expression patterns.")
      ),
      tags$h4("Main caution"),
      tags$p("Custom gene lists can miss unexpected biology. They can also force PCA or clustering to reflect your prior assumptions. For exploratory analysis, use HVGs first.")
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "Selected gene summary",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        uiOutput(ns("gene_summary_box")),
        br(),
        DT::DTOutput(ns("gene_summary_table")),
        br(),
        downloadButton(ns("dl_gene_summary"), "Download selected genes CSV")
      )
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "Missing genes",
        width = 12,
        solidHeader = TRUE,
        status = "warning",
        collapsible = TRUE,
        collapsed = TRUE,
        icon = icon("triangle-exclamation"),
        uiOutput(ns("missing_gene_message")),
        DT::DTOutput(ns("missing_gene_table")),
        br(),
        downloadButton(ns("dl_missing_genes"), "Download missing genes CSV")
      )
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "HVG plot",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        shinycssloaders::withSpinner(plotly::plotlyOutput(ns("hvg_plot")), type = 4),
        br(),
        fluidRow(
          column(6, downloadButton(ns("dl_hvg_pdf"), "HVG plot PDF")),
          column(6, downloadButton(ns("dl_hvg_svg"), "HVG plot SVG"))
        ),
        br(),
        bs4Dash::bs4Card(
          title = "Interpretation: HVG plot",
          width = 12,
          solidHeader = TRUE,
          status = "info",
          collapsible = TRUE,
          collapsed = TRUE,
          icon = icon("circle-info"),
          tags$p("The HVG plot shows average expression and dispersion. Labeled genes are among the strongest selected variable genes."),
          tags$ul(
            tags$li("Genes with high standardized variance are more informative for PCA."),
            tags$li("Very high HVGs can include technical genes, mitochondrial genes, ribosomal genes, or stress genes."),
            tags$li("Filtering mitochondrial or ribosomal genes may help if they dominate the HVG list."),
            tags$li("A strong HVG list does not prove biological relevance. It only shows variation across cells.")
          )
        )
      )
    )
  )
}

mod_sc_gene_selection_server <- function(id, seurat_r, sc_state = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    hvg_list <- reactiveVal(character(0))
    selected_genes <- reactiveVal(character(0))
    hvg_stats <- reactiveVal(NULL)
    gene_summary <- reactiveVal(NULL)
    missing_genes <- reactiveVal(NULL)
    seu_obj <- reactiveVal(NULL)
    hvg_plot_gg <- reactiveVal(NULL)
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
    
    assay_choices <- reactive({
      seu <- get_seu()
      if (is.null(seu)) return(c("RNA"))
      assays <- tryCatch(Seurat::Assays(seu), error = function(e) "RNA")
      if (length(assays) == 0) assays <- "RNA"
      assays
    })
    
    observe({
      choices <- assay_choices()
      updateSelectInput(session, "assay_use", choices = choices, selected = if ("RNA" %in% choices) "RNA" else choices[1])
    })
    
    assay_exists <- function(seu, assay) {
      isTRUE(assay %in% Seurat::Assays(seu))
    }
    
    get_assay_matrix_safe <- function(seu, assay = "RNA", preferred = c("data", "counts", "scale.data")) {
      preferred <- unique(preferred)
      
      for (layer_name in preferred) {
        mat <- tryCatch(
          Seurat::GetAssayData(seu, assay = assay, layer = layer_name),
          error = function(e) NULL
        )
        if (!is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0) {
          return(list(matrix = mat, source = paste0("layer:", layer_name)))
        }
      }
      
      for (slot_name in preferred) {
        mat <- tryCatch(
          Seurat::GetAssayData(seu, assay = assay, slot = slot_name),
          error = function(e) NULL
        )
        if (!is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0) {
          return(list(matrix = mat, source = paste0("slot:", slot_name)))
        }
      }
      
      list(matrix = NULL, source = NA_character_)
    }
    
    assay_has_data <- function(seu, assay) {
      res <- get_assay_matrix_safe(seu, assay = assay, preferred = c("data"))
      !is.null(res$matrix)
    }
    
    parse_gene_text <- function(x) {
      if (is.null(x) || !nzchar(x)) return(character(0))
      vals <- unlist(strsplit(x, "[\\s,;]+"))
      vals <- trimws(vals)
      vals <- vals[nzchar(vals)]
      unique(vals)
    }
    
    read_uploaded_genes <- reactive({
      file <- input$gene_file
      if (is.null(file)) return(character(0))
      
      ext <- tolower(tools::file_ext(file$name))
      sep <- if (ext == "tsv" || ext == "txt") "\t" else ","
      
      df <- tryCatch(
        utils::read.table(
          file$datapath,
          header = isTRUE(input$use_uploaded_header),
          sep = sep,
          quote = "\"",
          comment.char = "",
          stringsAsFactors = FALSE,
          check.names = FALSE
        ),
        error = function(e) NULL
      )
      
      if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) {
        return(character(0))
      }
      
      cn <- tolower(colnames(df))
      preferred <- which(cn %in% c("gene", "genes", "symbol", "symbols", "feature", "features", "gene_symbol"))
      use_col <- if (length(preferred) > 0) preferred[1] else 1
      
      vals <- as.character(df[[use_col]])
      vals <- trimws(vals)
      vals <- vals[nzchar(vals) & !is.na(vals)]
      unique(vals)
    })
    
    match_genes_to_assay <- function(genes, available_genes, mode = "exact") {
      genes <- unique(trimws(genes))
      genes <- genes[nzchar(genes)]
      
      if (length(genes) == 0) {
        return(list(found = character(0), missing = character(0), map = data.frame()))
      }
      
      if (mode == "case_insensitive") {
        available_upper <- toupper(available_genes)
        gene_upper <- toupper(genes)
        idx <- match(gene_upper, available_upper)
        found <- available_genes[idx[!is.na(idx)]]
        missing <- genes[is.na(idx)]
        map <- data.frame(
          requested_gene = genes[!is.na(idx)],
          matched_gene = found,
          stringsAsFactors = FALSE
        )
      } else {
        found <- genes[genes %in% available_genes]
        missing <- genes[!genes %in% available_genes]
        map <- data.frame(
          requested_gene = found,
          matched_gene = found,
          stringsAsFactors = FALSE
        )
      }
      
      list(found = unique(found), missing = unique(missing), map = map)
    }
    
    apply_gene_filters <- function(genes) {
      genes <- unique(genes)
      removed_reason <- data.frame(gene = character(0), reason = character(0), stringsAsFactors = FALSE)
      
      if (isTRUE(input$filter_mito)) {
        rem <- genes[grepl("^MT-|^mt-", genes)]
        if (length(rem) > 0) removed_reason <- rbind(removed_reason, data.frame(gene = rem, reason = "mitochondrial", stringsAsFactors = FALSE))
        genes <- genes[!grepl("^MT-|^mt-", genes)]
      }
      
      if (isTRUE(input$filter_ribo)) {
        rem <- genes[grepl("^RP[SL]|^Rp[sl]", genes)]
        if (length(rem) > 0) removed_reason <- rbind(removed_reason, data.frame(gene = rem, reason = "ribosomal", stringsAsFactors = FALSE))
        genes <- genes[!grepl("^RP[SL]|^Rp[sl]", genes)]
      }
      
      if (!is.null(input$filter_regex) && nzchar(input$filter_regex)) {
        valid_regex <- tryCatch({
          grepl(input$filter_regex, genes, ignore.case = TRUE)
        }, error = function(e) NULL)
        
        if (is.null(valid_regex)) {
          showModal(modalDialog(
            title = "Invalid regex",
            "The regex filter is invalid. It was ignored.",
            easyClose = TRUE
          ))
        } else {
          rem <- genes[valid_regex]
          if (length(rem) > 0) removed_reason <- rbind(removed_reason, data.frame(gene = rem, reason = "regex", stringsAsFactors = FALSE))
          genes <- genes[!valid_regex]
        }
      }
      
      list(genes = unique(genes), removed = removed_reason)
    }
    
    summarize_selected_genes <- function(seu, genes, assay, source_vector, matrix_source) {
      mat_res <- get_assay_matrix_safe(seu, assay = assay, preferred = c("data", "counts"))
      mat <- mat_res$matrix
      
      if (is.null(mat)) {
        return(data.frame(
          gene = genes,
          source = source_vector[genes],
          mean_expression = NA_real_,
          pct_cells_expressed = NA_real_,
          assay = assay,
          matrix_source = matrix_source,
          stringsAsFactors = FALSE
        ))
      }
      
      genes_present <- genes[genes %in% rownames(mat)]
      if (length(genes_present) == 0) {
        return(data.frame())
      }
      
      sub <- mat[genes_present, , drop = FALSE]
      data.frame(
        gene = genes_present,
        source = source_vector[genes_present],
        mean_expression = Matrix::rowMeans(sub),
        pct_cells_expressed = Matrix::rowMeans(sub > 0) * 100,
        assay = assay,
        matrix_source = mat_res$source,
        stringsAsFactors = FALSE
      )
    }
    
    make_missing_df <- function(genes, source_name) {
      genes <- unique(genes)
      genes <- genes[nzchar(genes)]
      if (length(genes) == 0) {
        return(data.frame(
          requested_gene = character(0),
          source = character(0),
          stringsAsFactors = FALSE
        ))
      }
      data.frame(
        requested_gene = genes,
        source = rep(source_name, length(genes)),
        stringsAsFactors = FALSE
      )
    }
    
    build_gene_selection_project <- function(seu_out, out_file = NULL) {
      list(
        seurat = seu_out,
        project_meta = list(
          project_type = "CoTRA_scRNA_project",
          step = "Gene_Selection",
          saved_at = Sys.time(),
          cotra_module = "modules/sc/gene_selection.R",
          output_path = out_file
        ),
        cells_meta = seu_out@meta.data,
        step = "Gene_Selection",
        saved_at = Sys.time(),
        genes_list = list(
          selected_genes = selected_genes(),
          hvgs = hvg_list()
        ),
        parameters = list(
          gene_selection = tryCatch(seu_out@misc$CoTRA_gene_selection, error = function(e) NULL),
          qc = tryCatch(seu_out@misc$CoTRA_qc, error = function(e) NULL)
        ),
        tables = list(
          gene_summary = gene_summary(),
          missing_genes = missing_genes()
        ),
        cotra_version = if (exists("cotra_version")) cotra_version else "unknown"
      )
    }
    
    output$session_save_hint <- renderUI({
      if (length(selected_genes()) == 0) {
        return(tags$span(style = "color:#777;", "Run gene selection before saving this step."))
      }
      if (is.null(last_session_path())) {
        return(tags$span(style = "color:#777;", "Save the current gene selection state as a CoTRA project RDS."))
      }
      tags$span(style = "color:#2E7D32;", paste("Last saved:", last_session_path()))
    })
    
    output$gene_summary_box <- renderUI({
      genes <- selected_genes()
      miss <- missing_genes()
      if (length(genes) == 0) {
        return(tags$div(style = "color:#777;", "No genes selected yet. Choose a mode and click Run gene selection."))
      }
      
      n_missing <- if (is.null(miss)) 0 else nrow(miss)
      tags$div(
        tags$b("Final selected genes: "), length(genes), tags$br(),
        tags$b("Missing requested genes: "), n_missing, tags$br(),
        tags$b("Downstream status: "), ifelse(length(genes) >= 50, "Ready for PCA", "Too few genes for robust PCA")
      )
    })
    
    observeEvent(input$run_gene_selection, {
      seu <- get_seu()
      req(seu)
      
      withProgress(message = "Running gene selection", value = 0.1, {
        assay <- input$assay_use %||% "RNA"
        if (!assay_exists(seu, assay)) {
          showModal(modalDialog(
            title = "Assay not found",
            paste("Assay", assay, "was not found in the Seurat object."),
            easyClose = TRUE
          ))
          return(NULL)
        }
        
        available_genes <- rownames(seu[[assay]])
        if (length(available_genes) == 0) {
          showModal(modalDialog(title = "No genes found", "The selected assay has no feature names.", easyClose = TRUE))
          return(NULL)
        }
        
        incProgress(0.2, detail = "Checking assay data")
        
        hvg <- character(0)
        hvg_plot <- NULL
        mode <- input$selection_mode %||% "hvg_only"
        needs_hvg <- mode %in% c("hvg_only", "hvg_custom", "hvg_uploaded", "hvg_custom_uploaded")
        
        if (needs_hvg) {
          if ((input$hvg_mode %||% "fvf") == "fvf") {
            if (!assay_has_data(seu, assay)) {
              seu <- Seurat::NormalizeData(
                seu,
                assay = assay,
                normalization.method = "LogNormalize",
                scale.factor = 10000,
                verbose = FALSE
              )
            }
            
            incProgress(0.35, detail = "Finding highly variable genes")
            
            seu <- Seurat::FindVariableFeatures(
              seu,
              assay = assay,
              selection.method = "vst",
              nfeatures = input$nfeatures %||% 2000,
              verbose = FALSE
            )
            
            hvg <- Seurat::VariableFeatures(seu, assay = assay)
            hvg_plot <- tryCatch({
              p <- Seurat::VariableFeaturePlot(seu, assay = assay)
              lab <- head(hvg, 15)
              Seurat::LabelPoints(p, points = lab, repel = TRUE) + safe_theme()
            }, error = function(e) NULL)
          } else {
            incProgress(0.35, detail = "Running SCTransform")
            
            if (!requireNamespace("sctransform", quietly = TRUE)) {
              showModal(modalDialog(
                title = "SCTransform unavailable",
                "Install the sctransform package or use FindVariableFeatures.",
                easyClose = TRUE
              ))
              return(NULL)
            }
            
            seu <- Seurat::SCTransform(
              seu,
              assay = assay,
              verbose = FALSE,
              return.only.var.genes = FALSE
            )
            
            assay <- "SCT"
            available_genes <- rownames(seu[[assay]])
            hvg <- Seurat::VariableFeatures(seu, assay = assay)
            hvg_plot <- tryCatch({
              p <- Seurat::VariableFeaturePlot(seu, assay = assay)
              lab <- head(hvg, 15)
              Seurat::LabelPoints(p, points = lab, repel = TRUE) + safe_theme()
            }, error = function(e) NULL)
          }
        }
        
        incProgress(0.55, detail = "Checking custom and uploaded genes")
        
        custom_requested <- parse_gene_text(input$custom_genes)
        uploaded_requested <- read_uploaded_genes()
        
        custom_match <- match_genes_to_assay(custom_requested, available_genes, input$gene_case_mode)
        uploaded_match <- match_genes_to_assay(uploaded_requested, available_genes, input$gene_case_mode)
        
        final <- character(0)
        source_vector <- character(0)
        
        add_source <- function(genes, source_name) {
          if (length(genes) == 0) return(NULL)
          for (g in genes) {
            if (!g %in% names(source_vector)) {
              source_vector[g] <<- source_name
            } else if (!grepl(source_name, source_vector[g], fixed = TRUE)) {
              source_vector[g] <<- paste(source_vector[g], source_name, sep = "+")
            }
          }
          NULL
        }
        
        if (mode %in% c("hvg_only", "hvg_custom", "hvg_uploaded", "hvg_custom_uploaded")) {
          final <- c(final, hvg)
          add_source(hvg, "HVG")
        }
        if (mode %in% c("custom_only", "hvg_custom", "hvg_custom_uploaded")) {
          final <- c(final, custom_match$found)
          add_source(custom_match$found, "custom")
        }
        if (mode %in% c("uploaded_only", "hvg_uploaded", "hvg_custom_uploaded")) {
          final <- c(final, uploaded_match$found)
          add_source(uploaded_match$found, "uploaded")
        }
        
        final <- unique(final)
        final <- final[final %in% available_genes]
        
        filtered <- apply_gene_filters(final)
        final <- filtered$genes
        
        if (length(final) == 0) {
          showModal(modalDialog(
            title = "No genes selected",
            "No valid genes remained after matching and filtering.",
            easyClose = TRUE
          ))
          return(NULL)
        }
        
        if (length(final) < 50) {
          showModal(modalDialog(
            title = "Small gene set",
            paste("Only", length(final), "genes were selected. This may be too few for robust PCA."),
            easyClose = TRUE
          ))
        }
        
        incProgress(0.75, detail = "Creating summary")
        
        source_vector <- source_vector[final]
        source_vector[is.na(source_vector)] <- "selected"
        
        summary_df <- summarize_selected_genes(seu, final, assay, source_vector, NA_character_)
        summary_df <- summary_df[order(summary_df$source, -summary_df$mean_expression), , drop = FALSE]
        
        miss_df <- rbind(
          make_missing_df(custom_match$missing, "custom"),
          make_missing_df(uploaded_match$missing, "uploaded")
        )
        
        if (!is.null(filtered$removed) && nrow(filtered$removed) > 0) {
          miss_df <- rbind(
            miss_df,
            data.frame(
              requested_gene = filtered$removed$gene,
              source = paste0("removed_by_", filtered$removed$reason),
              stringsAsFactors = FALSE
            )
          )
        }
        
        rownames(miss_df) <- NULL
        
        Seurat::VariableFeatures(seu, assay = assay) <- final
        seu@misc$CoTRA_gene_selection <- list(
          selected_genes = final,
          hvg_genes = hvg,
          custom_requested = custom_requested,
          uploaded_requested = uploaded_requested,
          missing_genes = miss_df,
          selection_mode = mode,
          assay = assay,
          hvg_method = input$hvg_mode,
          nfeatures = input$nfeatures,
          filters = list(
            filter_mito = isTRUE(input$filter_mito),
            filter_ribo = isTRUE(input$filter_ribo),
            filter_regex = input$filter_regex
          ),
          created_at = Sys.time()
        )
        
        hvg_list(hvg)
        selected_genes(final)
        hvg_stats(summary_df)
        gene_summary(summary_df)
        missing_genes(miss_df)
        seu_obj(seu)
        hvg_plot_gg(hvg_plot)
        
        if (!is.null(sc_state)) {
          sc_state$seurat <- seu
          sc_state$genes_list[["selected_genes"]] <- final
          sc_state$genes_list[["HVG"]] <- hvg
          sc_state$parameters$gene_selection <- seu@misc$CoTRA_gene_selection
        }
        
        incProgress(1)
      })
    })
    
    output$hvg_plot <- plotly::renderPlotly({
      p <- hvg_plot_gg()
      if (is.null(p)) return(NULL)
      plotly::ggplotly(p)
    })
    
    output$gene_summary_table <- DT::renderDT({
      df <- gene_summary()
      req(df)
      DT::datatable(df, options = list(pageLength = 20, scrollX = TRUE))
    })
    
    output$missing_gene_message <- renderUI({
      df <- missing_genes()
      if (is.null(df)) {
        return(tags$p("Run gene selection to check custom or uploaded genes."))
      }
      if (nrow(df) == 0) {
        return(tags$p(style = "color:#2E7D32;", "All requested custom or uploaded genes were found. No missing genes detected."))
      }
      tags$p(style = "color:#B00020;", "Some requested genes were not found in the selected assay or were removed by filters.")
    })
    
    output$missing_gene_table <- DT::renderDT({
      df <- missing_genes()
      if (is.null(df) || nrow(df) == 0) {
        df <- data.frame(message = "No missing genes detected.", stringsAsFactors = FALSE)
      }
      DT::datatable(df, options = list(pageLength = 20, scrollX = TRUE))
    })
    
    output$dl_gene_summary <- downloadHandler(
      filename = function() paste0("CoTRA_selected_genes_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        df <- gene_summary()
        req(df)
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    
    output$dl_missing_genes <- downloadHandler(
      filename = function() paste0("CoTRA_missing_genes_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        df <- missing_genes()
        if (is.null(df)) df <- data.frame()
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    
    output$dl_hvg_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_HVG_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) {
        p <- hvg_plot_gg()
        req(p)
        ggplot2::ggsave(file, plot = p, width = 7, height = 5)
      }
    )
    
    output$dl_hvg_svg <- downloadHandler(
      filename = function() paste0("CoTRA_HVG_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) {
        p <- hvg_plot_gg()
        req(p)
        ggplot2::ggsave(file, plot = p, device = "svg", width = 7, height = 5)
      }
    )
    
    output$download_gene_selection_session <- downloadHandler(
      filename = function() {
        paste0("CoTRA_scRNA_Gene_Selection_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
      },
      content = function(file) {
        seu_out <- if (!is.null(seu_obj())) seu_obj() else get_seu()
        req(seu_out)
        
        out_dir <- file.path("outputs", "scRNA", "sessioninfo")
        if (!dir.exists(out_dir)) {
          dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        }
        
        out_file <- file.path(
          out_dir,
          paste0("CoTRA_scRNA_Gene_Selection_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
        )
        
        project_obj <- build_gene_selection_project(seu_out, out_file)
        
        saveRDS(project_obj, out_file)
        file.copy(out_file, file, overwrite = TRUE)
        last_session_path(out_file)
      }
    )
    
    return(list(
      seurat = reactive({
        if (!is.null(seu_obj())) seu_obj() else get_seu()
      }),
      selected_genes = reactive(selected_genes()),
      hvgs = reactive(hvg_list()),
      gene_summary = reactive(gene_summary()),
      missing_genes = reactive(missing_genes()),
      hvg_stats = reactive(hvg_stats()),
      ready_for_pca = reactive(length(selected_genes()) >= 50),
      session_path = reactive(last_session_path())
    ))
  })
}
