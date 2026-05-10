# CoTRA Bulk Enrichment Module
# Integrated enrichment module with:
# - ORA for Hallmark, GO BP, GO MF, GO CC, KEGG, Reactome
# - Direction-aware ORA for up-regulated and down-regulated genes
# - KEGG Pathview pathway visualization
# - Network plots using cnetplot and emapplot
# - CSV/XLSX table downloads
# - PDF and SVG plot downloads

mod_bulk_enrich_ui <- function(id) {
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
        p("This module requires a differential expression table from the bulk RNA-seq workflow. The table should contain a gene identifier column, padj, log2 fold change, and preferably a ranking statistic such as stat or lfc."),
        
        h5("2. Check organism and gene ID mapping"),
        p("The module auto-detects organism from gene identifiers when possible. It supports human, mouse, and rat. You can apply the detected organism manually. Gene identifiers are mapped to SYMBOL, ENTREZID, and ENSEMBL using org.*.eg.db packages."),
        
        h5("3. Select ranking and filtering columns"),
        p("Choose the gene ID column and ranking column. GSEA uses the ranking column. ORA uses significant genes based on adjusted p value and log2 fold change cutoffs."),
        
        h5("4. Run ORA"),
        p("The ORA tab runs over-representation analysis for Hallmark, GO, KEGG, or Reactome. Only one enrichment family is selected at a time. If GO is selected, you can choose GO BP, GO MF, GO CC, or all GO ontologies."),
        
        h5("5. Inspect KEGG pathways"),
        p("The KEGG Pathview tab uses selected KEGG pathway IDs and fold-change values to draw pathway diagrams. Select a KEGG row from the ORA table before drawing Pathview."),
        
        h5("6. Explore network plots"),
        p("The Network plots tab supports cnetplot and emapplot for the selected ORA result."),
        
        h5("7. Download results"),
        p("All result tables can be saved as CSV or Excel. Main plots can be saved as PDF and SVG for publication or reporting.")
      )
    ),
    
    fluidRow(
      box(
        width = 12,
        title = "Enrichment analysis settings",
        solidHeader = TRUE,
        status = "primary",
        collapsible = TRUE,
        
        fluidRow(
          column(3, selectInput(ns("gene_col"), "Gene ID column", choices = NULL)),
          column(3, selectInput(ns("rank_col"), "Ranking column", choices = NULL)),
          column(3, selectInput(ns("organism"), "Organism", c("Human", "Mouse", "Rat"), selected = "Human")),
          column(3, textInput(ns("detected_org"), "Auto-detected organism", value = "", placeholder = "Not detected yet"))
        ),
        
        fluidRow(
          column(3, actionButton(ns("apply_detected_org"), "Use detected organism")),
          column(3, numericInput(ns("padj_cutoff"), "Padj cutoff for ORA", value = 0.05, min = 0, max = 1)),
          column(3, numericInput(ns("logfc_cutoff"), "Abs log2FC cutoff for ORA", value = 1, min = 0)),
          column(3, checkboxInput(ns("split_direction"), "Run ORA separately for up and down genes", value = TRUE))
        ),
        
        fluidRow(
          column(12, tags$small("ORA uses significant DE genes based on padj and log2 fold-change cutoffs."))
        ),
        
        fluidRow(
          column(12, verbatimTextOutput(ns("mapping_qc")))
        )
      )
    ),
    
    fluidRow(
      box(
        width = 12,
        title = "Enrichment results",
        solidHeader = TRUE,
        status = "primary",
        
        tabsetPanel(
          id = ns("enrich_tabs"),
          
          tabPanel(
            "ORA",
            br(),
            fluidRow(
              column(
                3,
                selectInput(
                  ns("ora_type"),
                  "Select enrichment analysis",
                  choices = c(
                    "Hallmark, SYMBOL" = "hallmark",
                    "GO, SYMBOL" = "go",
                    "KEGG, ENTREZID" = "kegg",
                    "Reactome, ENTREZID" = "reactome"
                  ),
                  selected = "hallmark"
                )
              ),
              column(2, uiOutput(ns("ora_go_selector"))),
              column(2, numericInput(ns("ora_show_n"), "Top categories", value = 10, min = 5, max = 50)),
              column(2, actionButton(ns("run_ora"), "Run ORA")),
              column(3, selectInput(ns("ora_download_type"), "Download type", c("CSV" = "csv", "Excel" = "xlsx"), selected = "csv"))
            ),
            br(),
            fluidRow(
              column(3, downloadButton(ns("dl_ora_table"), "Download Table")),
              column(3, downloadButton(ns("dl_ora_selected"), "Download Selected")),
              column(3, downloadButton(ns("dl_ora_pdf"), "PDF")),
              column(3, downloadButton(ns("dl_ora_svg"), "SVG"))
            ),
            br(),
            shinycssloaders::withSpinner(DT::DTOutput(ns("ora_table"))),
            br(),
            shinycssloaders::withSpinner(DT::DTOutput(ns("ora_selected_table"))),
            br(),
            shinycssloaders::withSpinner(plotly::plotlyOutput(ns("ora_plot"), height = "650px"))
          ),
          
          tabPanel(
            "KEGG Pathview",
            br(),
            fluidRow(
              column(3, selectInput(ns("pathview_source"), "Fold-change source", c("DE ranked table" = "de"), selected = "de")),
              column(3, numericInput(ns("pathview_limit"), "Fold-change color limit", value = 2, min = 0.1)),
              column(3, actionButton(ns("run_pathview"), "Draw selected KEGG pathway")),
              column(3, downloadButton(ns("dl_pathview_png"), "Download PNG"))
            ),
            br(),
            tags$small("Select a KEGG row in the ORA table before drawing the pathway."),
            br(), br(),
            imageOutput(ns("pathview_image"), height = "700px")
          ),
          
          tabPanel(
            "Network plots",
            br(),
            fluidRow(
              column(3, selectInput(ns("nr_type"), "Plot type", choices = c("Cnetplot" = "cnet", "Emapplot" = "emap"), selected = "cnet")),
              column(3, numericInput(ns("nr_show"), "Top categories", value = 10, min = 5)),
              column(3, actionButton(ns("draw_nr"), "Draw network"))
            ),
            br(),
            fluidRow(
              column(6, downloadButton(ns("dl_nr_pdf"), "PDF")),
              column(6, downloadButton(ns("dl_nr_svg"), "SVG"))
            ),
            br(),
            plotOutput(ns("nr_plot"), height = "600px")
          )
        )
      )
    )
  )
}

mod_bulk_enrichment_ui <- function(id) mod_bulk_enrich_ui(id)

mod_bulk_enrich_server <- function(id, deg_data, wgcna_output = NULL, dds.fc = NULL, anova_table = NULL) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
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
    
    wrap_text <- function(x, width = 35) {
      vapply(as.character(x), function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
    }
    
    wrap_text_html <- function(x, width = 35) {
      vapply(as.character(x), function(z) paste(strwrap(z, width = width), collapse = "<br>"), character(1))
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
    
    detect_organism_from_ids <- function(ids) {
      ids <- clean_ids(ids)
      ids <- ids[!is.na(ids) & nzchar(ids)]
      if (length(ids) == 0) return(NA_character_)
      samp <- head(ids, 500)
      
      mouse_ens <- mean(grepl("^ENSMUSG[0-9]+$", samp), na.rm = TRUE)
      human_ens <- mean(grepl("^ENSG[0-9]+$", samp), na.rm = TRUE)
      rat_ens <- mean(grepl("^ENSRNOG[0-9]+$", samp), na.rm = TRUE)
      
      if (!is.na(mouse_ens) && mouse_ens > 0.5) return("Mouse")
      if (!is.na(human_ens) && human_ens > 0.5) return("Human")
      if (!is.na(rat_ens) && rat_ens > 0.5) return("Rat")
      
      mouse_hits <- sum(samp %in% c("A1bg", "A2m", "Aaas", "Aak1", "Aanat", "Crx", "Rho", "Gnat1"), na.rm = TRUE)
      human_hits <- sum(samp %in% c("A1BG", "A2M", "ACTB", "XIST", "RPS4Y1", "KDM5D", "DDX3Y"), na.rm = TRUE)
      
      if (mouse_hits > human_hits && mouse_hits > 0) return("Mouse")
      if (human_hits > mouse_hits && human_hits > 0) return("Human")
      
      upper_frac <- mean(grepl("^[A-Z0-9\\.-]+$", samp), na.rm = TRUE)
      title_frac <- mean(grepl("^[A-Z][a-z0-9\\.-]+$", samp), na.rm = TRUE)
      
      if (!is.na(upper_frac) && upper_frac > 0.8) return("Human")
      if (!is.na(title_frac) && title_frac > 0.5) return("Mouse")
      NA_character_
    }
    
    guess_id_type <- function(ids, organism) {
      ids <- clean_ids(ids)
      patt <- switch(organism, "Mouse" = "^ENSMUSG[0-9]+$", "Human" = "^ENSG[0-9]+$", "Rat" = "^ENSRNOG[0-9]+$", "^ENS")
      frac_ens <- mean(grepl(patt, ids), na.rm = TRUE)
      if (!is.na(frac_ens) && frac_ens > 0.5) "ENSEMBL" else "SYMBOL"
    }
    
    looks_like_symbol <- function(ids, organism) {
      identical(guess_id_type(ids, organism), "SYMBOL")
    }
    
    orgdb_obj <- reactive({
      switch(input$organism, "Mouse" = org.Mm.eg.db::org.Mm.eg.db, "Rat" = org.Rn.eg.db::org.Rn.eg.db, org.Hs.eg.db::org.Hs.eg.db)
    })
    
    species_msigdb <- reactive({
      switch(input$organism, "Mouse" = "Mus musculus", "Rat" = "Rattus norvegicus", "Homo sapiens")
    })
    
    msigdb_db_species <- reactive({
      if (input$organism == "Mouse") "MM" else "HS"
    })
    
    kegg_code <- reactive({
      switch(input$organism, "Mouse" = "mmu", "Rat" = "rno", "hsa")
    })
    
    reactome_species <- reactive({
      switch(input$organism, "Mouse" = "mouse", "Rat" = "rat", "human")
    })
    
    deg_df <- reactive({
      df <- deg_data()
      req(df)
      validate(need(is.data.frame(df), "DE result must be a data frame"))
      validate(need(nrow(df) > 0, "DE result is empty"))
      df
    })
    
    observeEvent(deg_df(), {
      df <- deg_df()
      cols <- colnames(df)
      
      rank_default <- if ("stat" %in% cols) {
        "stat"
      } else if ("lfc" %in% cols) {
        "lfc"
      } else if ("log2FoldChange" %in% cols) {
        "log2FoldChange"
      } else {
        num_cols <- cols[vapply(df, is.numeric, logical(1))]
        if (length(num_cols) > 0) num_cols[1] else cols[1]
      }
      
      gene_default <- if ("gene" %in% cols) "gene" else cols[1]
      
      updateSelectInput(session, "gene_col", choices = cols, selected = gene_default)
      updateSelectInput(session, "rank_col", choices = cols, selected = rank_default)
      
      det_org <- detect_organism_from_ids(df[[gene_default]])
      updateTextInput(session, "detected_org", value = ifelse(is.na(det_org), "", det_org))
    }, ignoreInit = FALSE)
    
    observeEvent(input$apply_detected_org, {
      req(nzchar(input$detected_org))
      if (input$detected_org %in% c("Human", "Mouse", "Rat")) {
        updateSelectInput(session, "organism", selected = input$detected_org)
      }
    })
    
    map_ids <- reactive({
      df <- deg_df()
      ids <- clean_ids(df[[input$gene_col]])
      validate(need(length(ids) > 0, "No gene IDs available"))
      
      id_type <- guess_id_type(ids, input$organism)
      db <- orgdb_obj()
      
      map_df <- suppressWarnings(
        AnnotationDbi::select(db, keys = unique(ids), keytype = id_type, columns = c("SYMBOL", "ENTREZID", "ENSEMBL"))
      )
      
      validate(need(!is.null(map_df) && nrow(map_df) > 0, "No gene IDs mapped"))
      map_df <- as.data.frame(map_df, stringsAsFactors = FALSE)
      
      colnames(map_df)[colnames(map_df) == id_type] <- "input_id"
      map_df$input_id <- clean_ids(map_df$input_id)
      
      if (!"SYMBOL" %in% colnames(map_df)) map_df$SYMBOL <- NA_character_
      if (!"ENTREZID" %in% colnames(map_df)) map_df$ENTREZID <- NA_character_
      if (!"ENSEMBL" %in% colnames(map_df)) map_df$ENSEMBL <- NA_character_
      
      map_df <- map_df[!is.na(map_df$input_id) & nzchar(map_df$input_id), , drop = FALSE]
      
      map_df <- map_df %>%
        dplyr::group_by(.data$input_id) %>%
        dplyr::summarise(
          SYMBOL = first_non_na(.data$SYMBOL),
          ENTREZID = first_non_na(.data$ENTREZID),
          ENSEMBL = first_non_na(.data$ENSEMBL),
          .groups = "drop"
        )
      
      map_df
    })
    
    deg_mapped <- reactive({
      df <- deg_df()
      df$input_id <- clean_ids(df[[input$gene_col]])
      dplyr::left_join(df, map_ids(), by = "input_id")
    })
    
    lfc_col <- reactive({
      df <- deg_mapped()
      if ("lfc" %in% colnames(df)) return("lfc")
      if ("log2FoldChange" %in% colnames(df)) return("log2FoldChange")
      validate(need(FALSE, "No log fold change column found. Expected lfc or log2FoldChange."))
    })
    
    output$mapping_qc <- renderPrint({
      df <- deg_mapped()
      cat("Total DE rows:", nrow(df), "\n")
      cat("Mapped SYMBOL:", sum(!is.na(df$SYMBOL) & nzchar(df$SYMBOL)), "\n")
      cat("Mapped ENTREZID:", sum(!is.na(df$ENTREZID) & nzchar(df$ENTREZID)), "\n")
      cat("Detected organism:", ifelse(nzchar(input$detected_org), input$detected_org, "Not detected"), "\n")
      cat("Selected organism:", input$organism, "\n")
      if ("padj" %in% colnames(df)) cat("Genes passing padj cutoff:", sum(!is.na(df$padj) & df$padj <= input$padj_cutoff), "\n")
    })
    
    de_gene_sets <- reactive({
      df <- deg_mapped()
      validate(need("padj" %in% colnames(df), "No padj column found"))
      lc <- lfc_col()
      df$lfc_for_filter <- suppressWarnings(as.numeric(df[[lc]]))
      
      base <- !is.na(df$padj) & df$padj <= input$padj_cutoff & !is.na(df$lfc_for_filter) & abs(df$lfc_for_filter) >= input$logfc_cutoff
      up <- base & df$lfc_for_filter >= input$logfc_cutoff
      down <- base & df$lfc_for_filter <= -input$logfc_cutoff
      
      get_symbols <- function(x) {
        raw_ids <- clean_ids(x[[input$gene_col]])
        if (looks_like_symbol(raw_ids, input$organism)) {
          genes <- unique(as.character(raw_ids))
        } else {
          genes <- unique(as.character(x$SYMBOL))
        }
        genes[!is.na(genes) & nzchar(genes)]
      }
      
      get_entrez <- function(x) {
        genes <- unique(as.character(x$ENTREZID))
        genes[!is.na(genes) & nzchar(genes)]
      }
      
      out <- list()
      if (isTRUE(input$split_direction)) {
        out$Up <- list(symbol = get_symbols(df[up, , drop = FALSE]), entrez = get_entrez(df[up, , drop = FALSE]))
        out$Down <- list(symbol = get_symbols(df[down, , drop = FALSE]), entrez = get_entrez(df[down, , drop = FALSE]))
      } else {
        out$Significant <- list(symbol = get_symbols(df[base, , drop = FALSE]), entrez = get_entrez(df[base, , drop = FALSE]))
      }
      out
    })
    
    normalise_msigdb_collection <- function(collection_value) {
      switch(
        as.character(collection_value),
        "H" = c("H", "hallmark"),
        "C2" = c("C2", "curated"),
        "C5" = c("C5", "ontology"),
        "hallmark" = c("hallmark", "H"),
        "curated" = c("curated", "C2"),
        "ontology" = c("ontology", "C5"),
        c(collection_value)
      )
    }
    
    get_msigdb <- function(collection_value) {
      fmls <- names(formals(msigdbr::msigdbr))
      collection_candidates <- normalise_msigdb_collection(collection_value)
      last_error <- NULL
      
      for (coll in collection_candidates) {
        args_list <- list(
          list(db_species = msigdb_db_species(), species = species_msigdb(), collection = coll),
          list(db_species = msigdb_db_species(), species = species_msigdb(), category = coll),
          list(species = species_msigdb(), collection = coll),
          list(species = species_msigdb(), category = coll)
        )
        
        for (args in args_list) {
          args <- args[names(args) %in% fmls]
          res <- tryCatch(
            do.call(msigdbr::msigdbr, args),
            error = function(e) {
              last_error <<- conditionMessage(e)
              NULL
            }
          )
          
          if (!is.null(res) && is.data.frame(res) && nrow(res) > 0) {
            return(res)
          }
        }
      }
      
      available <- tryCatch(
        paste(capture.output(print(msigdbr::msigdbr_collections())), collapse = "\n"),
        error = function(e) "Could not read msigdbr_collections()."
      )
      
      stop(paste0("Could not load MSigDB collection: ", collection_value, "\nLast msigdbr error: ", last_error, "\nAvailable collections:\n", available))
    }
    
    hallmark_term2gene <- reactive({
      msig <- get_msigdb("H")
      as.data.frame(msig[, c("gs_name", "gene_symbol")])
    })
    
    output$ora_go_selector <- renderUI({
      if (identical(input$ora_type, "go")) {
        selectInput(
          ns("ora_go_type"),
          "GO ontology",
          choices = c("GO BP" = "gobp", "GO MF" = "gomf", "GO CC" = "gocc", "All GO" = "go_all"),
          selected = "gobp"
        )
      }
    })
    
    expand_enrichment_choice <- function(choice, go_choice = NULL) {
      if (identical(choice, "go")) {
        if (is.null(go_choice) || !nzchar(go_choice)) return("gobp")
        if (identical(go_choice, "go_all")) return(c("gobp", "gomf", "gocc"))
        return(go_choice)
      }
      choice
    }
    
    run_single_ora <- function(symbols, entrez, type_choice, go_choice = NULL) {
      types <- expand_enrichment_choice(type_choice, go_choice)
      out <- list()
      db <- orgdb_obj()
      
      if ("hallmark" %in% types && length(symbols) > 0) {
        out$Hallmark <- clusterProfiler::enricher(
          gene = symbols,
          TERM2GENE = hallmark_term2gene(),
          pAdjustMethod = "BH",
          pvalueCutoff = input$padj_cutoff,
          qvalueCutoff = input$padj_cutoff
        )
      }
      
      if ("gobp" %in% types && length(symbols) > 0) {
        out$GO_BP <- clusterProfiler::enrichGO(gene = symbols, OrgDb = db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = input$padj_cutoff, qvalueCutoff = input$padj_cutoff, readable = TRUE)
      }
      
      if ("gomf" %in% types && length(symbols) > 0) {
        out$GO_MF <- clusterProfiler::enrichGO(gene = symbols, OrgDb = db, keyType = "SYMBOL", ont = "MF", pAdjustMethod = "BH", pvalueCutoff = input$padj_cutoff, qvalueCutoff = input$padj_cutoff, readable = TRUE)
      }
      
      if ("gocc" %in% types && length(symbols) > 0) {
        out$GO_CC <- clusterProfiler::enrichGO(gene = symbols, OrgDb = db, keyType = "SYMBOL", ont = "CC", pAdjustMethod = "BH", pvalueCutoff = input$padj_cutoff, qvalueCutoff = input$padj_cutoff, readable = TRUE)
      }
      
      if ("kegg" %in% types && length(entrez) > 0) {
        out$KEGG <- clusterProfiler::enrichKEGG(gene = entrez, organism = kegg_code(), pvalueCutoff = input$padj_cutoff)
      }
      
      if ("reactome" %in% types && length(entrez) > 0) {
        out$Reactome <- ReactomePA::enrichPathway(gene = entrez, organism = reactome_species(), readable = TRUE, pvalueCutoff = input$padj_cutoff)
      }
      
      out
    }
    
    flatten_enrich_list <- function(res_nested) {
      all <- list()
      for (set_name in names(res_nested)) {
        inner <- res_nested[[set_name]]
        for (type_name in names(inner)) {
          obj <- inner[[type_name]]
          if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
            d <- as.data.frame(obj)
            d$Gene_Set <- set_name
            d$Type <- type_name
            d <- d[, c("Gene_Set", "Type", setdiff(colnames(d), c("Gene_Set", "Type"))), drop = FALSE]
            all[[length(all) + 1]] <- d
          }
        }
      }
      if (length(all) == 0) return(NULL)
      do.call(rbind, all)
    }
    
    ora_res <- eventReactive(input$run_ora, {
      sets <- de_gene_sets()
      out <- list()
      for (nm in names(sets)) {
        out[[nm]] <- run_single_ora(sets[[nm]]$symbol, sets[[nm]]$entrez, input$ora_type, input$ora_go_type)
      }
      out
    })
    
    ora_df <- reactive({
      df <- flatten_enrich_list(ora_res())
      validate(need(!is.null(df) && nrow(df) > 0, "No ORA enrichment found"))
      df
    })
    
    selected_ora_row <- reactive({
      df <- ora_df()
      idx <- input$ora_table_rows_selected
      if (is.null(idx) || length(idx) == 0) return(df[1, , drop = FALSE])
      df[idx[1], , drop = FALSE]
    })
    
    selected_ora_obj <- reactive({
      sel <- selected_ora_row()
      res <- ora_res()
      set_name <- sel$Gene_Set[1]
      type_name <- sel$Type[1]
      obj <- res[[set_name]][[type_name]]
      validate(need(!is.null(obj), "Selected enrichment object not available"))
      obj
    })
    
    ora_plot_obj <- reactive({
      obj <- selected_ora_obj()
      validate(need(nrow(as.data.frame(obj)) > 0, "No ORA terms available for plotting"))
      
      obj2 <- obj
      obj2@result$Description <- wrap_text_html(obj2@result$Description, width = 35)
      
      enrichplot::dotplot(obj2, showCategory = input$ora_show_n) +
        ggplot2::labs(size = "No. of Genes") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.y = ggplot2::element_text(size = 9, lineheight = 0.9),
          axis.title.y = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(5.5, 20, 5.5, 120, "pt")
        )
    })
    
    output$ora_table <- DT::renderDT({
      req(input$run_ora > 0)
      DT::datatable(ora_df(), selection = list(mode = "single", target = "row"), options = list(pageLength = 20, scrollX = TRUE))
    })
    
    output$ora_selected_table <- DT::renderDT({
      req(input$run_ora > 0)
      DT::datatable(selected_ora_row(), options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
    })
    
    output$ora_plot <- plotly::renderPlotly({
      req(input$run_ora > 0)
      plotly::ggplotly(ora_plot_obj(), tooltip = c("x", "y", "colour"))
    })
    
    output$dl_ora_table <- downloadHandler(
      filename = function() if (input$ora_download_type == "xlsx") "ORA_results.xlsx" else "ORA_results.csv",
      content = function(file) save_table_file(ora_df(), file, input$ora_download_type, "ORA_results")
    )
    
    output$dl_ora_selected <- downloadHandler(
      filename = function() if (input$ora_download_type == "xlsx") "ORA_selected.xlsx" else "ORA_selected.csv",
      content = function(file) save_table_file(selected_ora_row(), file, input$ora_download_type, "ORA_selected")
    )
    
    output$dl_ora_pdf <- downloadHandler(
      filename = function() "ORA_plot.pdf",
      content = function(file) { grDevices::pdf(file, width = 9, height = 7); print(ora_plot_obj()); grDevices::dev.off() }
    )
    
    output$dl_ora_svg <- downloadHandler(
      filename = function() "ORA_plot.svg",
      content = function(file) { svglite::svglite(file, width = 9, height = 7); print(ora_plot_obj()); grDevices::dev.off() }
    )
    
    selected_kegg_id <- reactive({
      sel <- tryCatch(selected_ora_row(), error = function(e) NULL)
      validate(need(!is.null(sel) && "Type" %in% colnames(sel) && sel$Type[1] == "KEGG", "Select a KEGG row from the ORA table first"))
      sel$ID[1]
    })
    
    pathview_gene_vector <- reactive({
      df <- deg_mapped()
      lc <- lfc_col()
      vals <- suppressWarnings(as.numeric(df[[lc]]))
      entrez <- as.character(df$ENTREZID)
      keep <- !is.na(vals) & !is.na(entrez) & nzchar(entrez)
      vals <- vals[keep]
      entrez <- entrez[keep]
      tmp <- data.frame(entrez = entrez, val = vals, stringsAsFactors = FALSE)
      tmp <- tmp[order(abs(tmp$val), decreasing = TRUE), , drop = FALSE]
      tmp <- tmp[!duplicated(tmp$entrez), , drop = FALSE]
      v <- tmp$val
      names(v) <- tmp$entrez
      validate(need(length(v) > 0, "No Entrez fold-change vector available for Pathview"))
      v
    })
    
    pathview_file <- eventReactive(input$run_pathview, {
      pathway_id <- selected_kegg_id()
      species <- kegg_code()
      gene_vec <- pathview_gene_vector()
      old_files <- list.files(pattern = paste0(pathway_id, ".*\\.png$"), full.names = TRUE)
      if (length(old_files) > 0) unlink(old_files)
      pathview::pathview(gene.data = gene_vec, pathway.id = pathway_id, species = species, limit = list(gene = as.numeric(input$pathview_limit), cpd = 1))
      files <- list.files(pattern = paste0(pathway_id, ".*pathview\\.png$"), full.names = TRUE)
      if (length(files) == 0) files <- list.files(pattern = paste0(pathway_id, ".*\\.png$"), full.names = TRUE)
      validate(need(length(files) > 0, "Pathview image was not created"))
      files[1]
    })
    
    output$pathview_image <- renderImage({
      file <- pathview_file()
      list(src = file, contentType = "image/png", alt = "KEGG Pathview")
    }, deleteFile = FALSE)
    
    output$dl_pathview_png <- downloadHandler(
      filename = function() paste0("KEGG_pathview_", selected_kegg_id(), ".png"),
      content = function(file) file.copy(pathview_file(), file, overwrite = TRUE)
    )
    
    nr_plot_obj <- eventReactive(input$draw_nr, {
      obj <- selected_ora_obj()
      validate(need(!is.null(obj) && nrow(as.data.frame(obj)) > 0, "No ORA enrichment object available"))
      
      if (input$nr_type == "cnet") {
        p <- enrichplot::cnetplot(obj, showCategory = input$nr_show) +
          ggplot2::labs(size = "No. of Genes")
      } else {
        obj_sim <- enrichplot::pairwise_termsim(obj)
        p <- enrichplot::emapplot(obj_sim, showCategory = input$nr_show) +
          ggplot2::labs(size = "No. of Genes")
      }
      
      p
    })
    
    output$nr_plot <- renderPlot({
      req(input$draw_nr > 0)
      p <- tryCatch(nr_plot_obj(), error = function(e) NULL)
      validate(need(!is.null(p), "Run ORA, select a result, then draw the network plot."))
      print(p)
    })
    
    output$dl_nr_pdf <- downloadHandler(
      filename = function() paste0("network_plot_", input$nr_type, ".pdf"),
      content = function(file) { grDevices::pdf(file, width = 9, height = 7); print(nr_plot_obj()); grDevices::dev.off() }
    )
    
    output$dl_nr_svg <- downloadHandler(
      filename = function() paste0("network_plot_", input$nr_type, ".svg"),
      content = function(file) { svglite::svglite(file, width = 9, height = 7); print(nr_plot_obj()); grDevices::dev.off() }
    )
    
    return(list(
      ora = ora_res,
      ora_table = ora_df
    ))
  })
}

mod_bulk_enrichment_server <- function(id, deg_data, wgcna_output = NULL, dds.fc = NULL, anova_table = NULL) {
  mod_bulk_enrich_server(id, deg_data, wgcna_output, dds.fc, anova_table)
}
