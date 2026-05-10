# ==========================================================
# modules/bulk/annotation.R
# Gene annotation module
# Uses DE results from de_analysis.R as input
# Supports Human, Mouse, Rat
# Uses EnsDb packages for offline genomic annotation
# User can choose all tested genes or differentially expressed genes
# ==========================================================

mod_bulk_annotation_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    shinyjs::useShinyjs(),
    
    actionButton(ns("help_toggle"), icon("question-circle"), class = "btn btn-info btn-sm"),
    
    shinyjs::hidden(
      div(
        id = ns("help_box"),
        style = "border:1px solid #d9d9d9; padding:12px; margin-top:10px; border-radius:6px; background:#f7f7f7;",
        
        h4("Step-by-Step Guide"),
        
        h5("1. Run differential expression first"),
        p("This module uses the output of the DE analysis step. You can annotate either all tested genes or only differentially expressed genes."),
        
        h5("2. Select species"),
        p("Choose the organism for annotation. Mouse is selected by default."),
        
        h5("3. Choose gene set"),
        p("Select whether to annotate all tested genes or only differentially expressed genes. All tested genes is selected by default."),
        
        h5("4. Run annotation"),
        p("Click Run Annotation to map genes to gene symbol, gene name, Ensembl gene ID, chromosome, genomic start and end positions, strand, biotype, and gene length using offline EnsDb annotation packages."),
        
        h5("5. Explore results"),
        p("Use the annotation table and plots to inspect gene biotype, chromosome distribution, gene length, and genomic location."),
        
        h5("6. Download results"),
        p("Download the annotated gene table as CSV for downstream analysis or reporting.")
      )
    ),
    
    br(),
    
    fluidRow(
      column(
        3,
        selectInput(
          ns("species"),
          "Select organism",
          choices = c(
            "Mouse" = "mouse",
            "Human" = "human",
            "Rat" = "rat"
          ),
          selected = "mouse"
        )
      ),
      column(
        3,
        selectInput(
          ns("gene_set"),
          "Genes to annotate",
          choices = c(
            "All tested genes" = "all",
            "Differentially expressed genes" = "sig"
          ),
          selected = "all"
        )
      ),
      column(
        3,
        actionButton(
          ns("run_annot"),
          "Run Annotation",
          class = "btn btn-primary"
        )
      ),
      column(
        3,
        downloadButton(
          ns("download_table"),
          "Download CSV",
          class = "btn btn-success"
        )
      )
    ),
    
    br(),
    verbatimTextOutput(ns("status")),
    hr(),
    
    DT::DTOutput(ns("annot_table")),
    hr(),
    
    tabsetPanel(
      tabPanel(
        "Gene Biotype",
        br(),
        downloadButton(ns("download_type_pdf"), "Download PDF"),
        downloadButton(ns("download_type_svg"), "Download SVG"),
        br(), br(),
        plotly::plotlyOutput(ns("plot_type"), height = "520px")
      ),
      tabPanel(
        "Chromosome Distribution",
        br(),
        downloadButton(ns("download_chr_pdf"), "Download PDF"),
        downloadButton(ns("download_chr_svg"), "Download SVG"),
        br(), br(),
        plotly::plotlyOutput(ns("plot_chr"), height = "520px")
      ),
      tabPanel(
        "Gene Length",
        br(),
        downloadButton(ns("download_len_pdf"), "Download PDF"),
        downloadButton(ns("download_len_svg"), "Download SVG"),
        br(), br(),
        plotly::plotlyOutput(ns("plot_len"), height = "520px")
      ),
      tabPanel(
        "Genomic Location",
        br(),
        downloadButton(ns("download_ideo_pdf"), "Download PDF"),
        downloadButton(ns("download_ideo_svg"), "Download SVG"),
        br(), br(),
        plotly::plotlyOutput(ns("plot_ideo"), height = "520px")
      )
    )
  )
}

mod_bulk_annotation_server <- function(id, de_results) {
  moduleServer(id, function(input, output, session) {
    
    rv <- reactiveValues(annot = NULL)
    
    observeEvent(input$help_toggle, {
      shinyjs::toggle(id = "help_box", anim = TRUE)
    })
    
    get_ensdb <- function(species) {
      switch(
        species,
        "mouse" = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
        "human" = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
        "rat"   = EnsDb.Rnorvegicus.v79::EnsDb.Rnorvegicus.v79,
        EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79
      )
    }
    
    get_chr_levels <- function(species) {
      switch(
        species,
        "mouse" = c(as.character(1:19), "X", "Y", "MT"),
        "human" = c(as.character(1:22), "X", "Y", "MT"),
        "rat"   = c(as.character(1:20), "X", "Y", "MT"),
        c(as.character(1:19), "X", "Y", "MT")
      )
    }
    
    first_non_na <- function(x) {
      x <- x[!is.na(x) & nzchar(as.character(x))]
      if (length(x) == 0) return(NA_character_)
      as.character(x[1])
    }
    
    clean_gene_ids <- function(x) {
      x <- trimws(as.character(x))
      x <- sub("\\.\\d+$", "", x)
      x
    }
    
    guess_gene_id_type <- function(genes, species) {
      genes2 <- clean_gene_ids(genes)
      
      ens_pattern <- switch(
        species,
        "mouse" = "^ENSMUSG[0-9]+$",
        "human" = "^ENSG[0-9]+$",
        "rat"   = "^ENSRNOG[0-9]+$",
        "^ENS"
      )
      
      frac_ens <- mean(grepl(ens_pattern, genes2), na.rm = TRUE)
      if (!is.na(frac_ens) && frac_ens > 0.5) return("ensembl")
      "symbol"
    }
    
    filter_de_input <- function(de_df, gene_set) {
      if (is.null(de_df) || nrow(de_df) == 0) return(NULL)
      if (!"gene" %in% colnames(de_df)) return(NULL)
      
      if (identical(gene_set, "sig")) {
        if (!"significant" %in% colnames(de_df)) return(NULL)
        de_df <- de_df %>% dplyr::filter(significant %in% TRUE)
      }
      
      if (nrow(de_df) == 0) return(NULL)
      de_df
    }
    
    annotate_genes_ensdb <- function(genes, species) {
      edb <- get_ensdb(species)
      input_genes <- clean_gene_ids(genes)
      id_type <- guess_gene_id_type(input_genes, species)
      
      gr <- tryCatch({
        if (identical(id_type, "ensembl")) {
          ensembldb::genes(
            edb,
            filter = AnnotationFilter::GeneIdFilter(input_genes)
          )
        } else {
          ensembldb::genes(
            edb,
            filter = AnnotationFilter::GeneNameFilter(input_genes)
          )
        }
      }, error = function(e) NULL)
      
      if (is.null(gr) || length(gr) == 0) {
        return(data.frame(
          input_gene = input_genes,
          gene_symbol = NA_character_,
          gene_name = NA_character_,
          ensembl_gene_id = NA_character_,
          chromosome = NA_character_,
          start_position = NA_real_,
          end_position = NA_real_,
          strand = NA_character_,
          gene_biotype = NA_character_,
          gene_length = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
      
      meta <- as.data.frame(S4Vectors::mcols(gr), stringsAsFactors = FALSE)
      
      ann <- data.frame(
        gene_symbol = if ("gene_name" %in% colnames(meta)) as.character(meta$gene_name) else NA_character_,
        ensembl_gene_id = if ("gene_id" %in% colnames(meta)) as.character(meta$gene_id) else NA_character_,
        chromosome = as.character(GenomeInfoDb::seqnames(gr)),
        start_position = as.numeric(IRanges::start(gr)),
        end_position = as.numeric(IRanges::end(gr)),
        strand = as.character(BiocGenerics::strand(gr)),
        gene_biotype = if ("gene_biotype" %in% colnames(meta)) as.character(meta$gene_biotype) else NA_character_,
        stringsAsFactors = FALSE
      )
      
      if (identical(id_type, "ensembl")) {
        ann$input_gene <- ann$ensembl_gene_id
      } else {
        ann$input_gene <- ann$gene_symbol
      }
      
      ann$input_gene <- clean_gene_ids(ann$input_gene)
      ann$gene_symbol <- as.character(ann$gene_symbol)
      ann$gene_name <- as.character(ann$gene_symbol)
      
      ann <- ann %>%
        dplyr::group_by(input_gene) %>%
        dplyr::summarise(
          gene_symbol = first_non_na(gene_symbol),
          gene_name = first_non_na(gene_name),
          ensembl_gene_id = first_non_na(ensembl_gene_id),
          chromosome = first_non_na(chromosome),
          start_position = suppressWarnings(as.numeric(first_non_na(as.character(start_position)))),
          end_position = suppressWarnings(as.numeric(first_non_na(as.character(end_position)))),
          strand = first_non_na(strand),
          gene_biotype = first_non_na(gene_biotype),
          .groups = "drop"
        )
      
      ann$gene_length <- ifelse(
        !is.na(ann$start_position) & !is.na(ann$end_position),
        ann$end_position - ann$start_position + 1,
        NA_real_
      )
      
      missing_genes <- setdiff(input_genes, ann$input_gene)
      if (length(missing_genes) > 0) {
        ann_missing <- data.frame(
          input_gene = missing_genes,
          gene_symbol = NA_character_,
          gene_name = NA_character_,
          ensembl_gene_id = NA_character_,
          chromosome = NA_character_,
          start_position = NA_real_,
          end_position = NA_real_,
          strand = NA_character_,
          gene_biotype = NA_character_,
          gene_length = NA_real_,
          stringsAsFactors = FALSE
        )
        ann <- dplyr::bind_rows(ann, ann_missing)
      }
      
      ann
    }
    
    build_annotation <- function(de_df, species, gene_set) {
      de_df <- filter_de_input(de_df, gene_set)
      if (is.null(de_df) || nrow(de_df) == 0) return(NULL)
      
      de_df <- de_df %>%
        dplyr::mutate(
          gene = as.character(gene),
          gene_clean = clean_gene_ids(gene)
        )
      
      genes <- unique(de_df$gene_clean)
      genes <- genes[!is.na(genes) & nzchar(genes)]
      if (length(genes) == 0) return(NULL)
      
      ann <- annotate_genes_ensdb(genes, species)
      
      out <- de_df %>%
        dplyr::left_join(ann, by = c("gene_clean" = "input_gene"))
      
      chr_levels <- get_chr_levels(species)
      out$chromosome <- gsub("^chr", "", as.character(out$chromosome), ignore.case = TRUE)
      out$chromosome[out$chromosome %in% c("M", "MT")] <- "MT"
      out$chromosome[!(out$chromosome %in% chr_levels)] <- NA
      out$chromosome <- factor(out$chromosome, levels = chr_levels)
      
      out$annotation_status <- ifelse(
        !is.na(out$ensembl_gene_id) | !is.na(out$gene_symbol),
        "Mapped",
        "Unmapped"
      )
      
      out
    }
    
    observeEvent(input$run_annot, {
      output$status <- renderText("Running annotation...")
      
      de_df <- try(de_results(), silent = TRUE)
      
      if (inherits(de_df, "try-error") || is.null(de_df) || nrow(de_df) == 0) {
        rv$annot <- NULL
        output$status <- renderText("No DE results available. Run DE analysis first.")
        return(NULL)
      }
      
      if (!"gene" %in% colnames(de_df)) {
        rv$annot <- NULL
        output$status <- renderText("DE results do not contain a gene column.")
        return(NULL)
      }
      
      if (identical(input$gene_set, "sig")) {
        if (!"significant" %in% colnames(de_df)) {
          rv$annot <- NULL
          output$status <- renderText("DE results do not contain a significant column.")
          return(NULL)
        }
        
        sig_n <- sum(de_df$significant %in% TRUE, na.rm = TRUE)
        if (sig_n == 0) {
          rv$annot <- NULL
          output$status <- renderText("No significant DE genes available for annotation under the current DE thresholds.")
          return(NULL)
        }
      }
      
      ann <- tryCatch(
        build_annotation(de_df, input$species, input$gene_set),
        error = function(e) e
      )
      
      if (inherits(ann, "error")) {
        rv$annot <- NULL
        output$status <- renderText(paste("Annotation failed:", ann$message))
        return(NULL)
      }
      
      if (is.null(ann) || nrow(ann) == 0) {
        rv$annot <- NULL
        output$status <- renderText("No annotation results returned.")
        return(NULL)
      }
      
      rv$annot <- ann
      
      mapped_n <- sum(rv$annot$annotation_status == "Mapped", na.rm = TRUE)
      unique_input_n <- length(unique(rv$annot$gene))
      input_label <- if (identical(input$gene_set, "sig")) "differentially expressed genes" else "tested genes"
      
      output$status <- renderText(
        paste0(
          "Annotation complete. ",
          unique_input_n, " ", input_label, " processed, ",
          mapped_n, " rows mapped."
        )
      )
    })
    
    output$annot_table <- DT::renderDT({
      req(rv$annot)
      DT::datatable(
        rv$annot,
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 15)
      )
    })
    
    plot_data_type <- reactive({
      req(rv$annot)
      rv$annot %>%
        dplyr::filter(!is.na(gene_biotype), nzchar(gene_biotype)) %>%
        dplyr::count(gene_biotype, name = "n") %>%
        dplyr::arrange(dplyr::desc(n))
    })
    
    plot_data_chr <- reactive({
      req(rv$annot)
      rv$annot %>%
        dplyr::filter(!is.na(chromosome)) %>%
        dplyr::count(chromosome, name = "n")
    })
    
    plot_data_len <- reactive({
      req(rv$annot)
      rv$annot %>%
        dplyr::filter(!is.na(gene_length), is.finite(gene_length), gene_length > 0)
    })
    
    plot_data_ideo <- reactive({
      req(rv$annot)
      rv$annot %>%
        dplyr::filter(!is.na(chromosome), !is.na(start_position))
    })
    
    output$plot_type <- plotly::renderPlotly({
      df <- plot_data_type()
      validate(need(nrow(df) > 0, "No gene biotype data available."))
      
      plotly::plot_ly(
        df,
        labels = ~gene_biotype,
        values = ~n,
        type = "pie",
        textinfo = "label+percent"
      ) %>%
        plotly::layout(title = "Gene biotype distribution")
    })
    
    output$plot_chr <- plotly::renderPlotly({
      df <- plot_data_chr()
      validate(need(nrow(df) > 0, "No chromosome data available."))
      
      plotly::plot_ly(
        df,
        x = ~chromosome,
        y = ~n,
        type = "bar",
        text = ~n,
        hoverinfo = "text",
        hovertext = ~paste0("Chromosome: ", chromosome, "<br>Genes: ", n)
      ) %>%
        plotly::layout(
          title = "Chromosome distribution",
          xaxis = list(title = "Chromosome"),
          yaxis = list(title = "Gene count")
        )
    })
    
    output$plot_len <- plotly::renderPlotly({
      df <- plot_data_len()
      validate(need(nrow(df) > 0, "No gene length data available."))
      
      plotly::plot_ly(
        df,
        x = ~gene_length,
        type = "histogram"
      ) %>%
        plotly::layout(
          title = "Gene length distribution",
          xaxis = list(title = "Gene length"),
          yaxis = list(title = "Count")
        )
    })
    
    output$plot_ideo <- plotly::renderPlotly({
      df <- plot_data_ideo()
      validate(need(nrow(df) > 0, "No genomic location data available."))
      
      plotly::plot_ly(
        df,
        x = ~chromosome,
        y = ~start_position,
        type = "scatter",
        mode = "markers",
        text = ~paste0(
          "Input gene: ", gene,
          "<br>Gene symbol: ", gene_symbol,
          "<br>Gene name: ", gene_name,
          "<br>Ensembl: ", ensembl_gene_id,
          "<br>Chromosome: ", chromosome,
          "<br>Start: ", start_position,
          "<br>End: ", end_position,
          "<br>Strand: ", strand,
          "<br>Biotype: ", gene_biotype,
          "<br>Gene length: ", gene_length
        ),
        hoverinfo = "text",
        marker = list(size = 6)
      ) %>%
        plotly::layout(
          title = "Genomic locations of annotated genes",
          xaxis = list(title = "Chromosome"),
          yaxis = list(title = "Start position")
        )
    })
    plot_type_static <- function() {
      df <- plot_data_type()
      validate(need(nrow(df) > 0, "No gene biotype data available."))
      
      ggplot2::ggplot(df, ggplot2::aes(x = "", y = n, fill = gene_biotype)) +
        ggplot2::geom_col(width = 1) +
        ggplot2::coord_polar(theta = "y") +
        ggplot2::theme_void() +
        ggplot2::labs(title = "Gene biotype distribution", fill = "Gene biotype")
    }
    
    plot_chr_static <- function() {
      df <- plot_data_chr()
      validate(need(nrow(df) > 0, "No chromosome data available."))
      
      ggplot2::ggplot(df, ggplot2::aes(x = chromosome, y = n)) +
        ggplot2::geom_col() +
        ggplot2::theme_bw() +
        ggplot2::labs(
          title = "Chromosome distribution",
          x = "Chromosome",
          y = "Gene count"
        )
    }
    
    plot_len_static <- function() {
      df <- plot_data_len()
      validate(need(nrow(df) > 0, "No gene length data available."))
      
      ggplot2::ggplot(df, ggplot2::aes(x = gene_length)) +
        ggplot2::geom_histogram(bins = 40) +
        ggplot2::theme_bw() +
        ggplot2::labs(
          title = "Gene length distribution",
          x = "Gene length",
          y = "Count"
        )
    }
    
    plot_ideo_static <- function() {
      df <- plot_data_ideo()
      validate(need(nrow(df) > 0, "No genomic location data available."))
      
      ggplot2::ggplot(df, ggplot2::aes(x = chromosome, y = start_position)) +
        ggplot2::geom_point(size = 1.8) +
        ggplot2::theme_bw() +
        ggplot2::labs(
          title = "Genomic locations of annotated genes",
          x = "Chromosome",
          y = "Start position"
        )
    }
    
    save_static_plot <- function(plot_fun, file, type, width = 9, height = 6) {
      p <- plot_fun()
      
      if (identical(type, "pdf")) {
        grDevices::pdf(file, width = width, height = height)
        print(p)
        grDevices::dev.off()
      } else {
        svglite::svglite(file, width = width, height = height)
        print(p)
        grDevices::dev.off()
      }
    }
    
    output$download_type_pdf <- downloadHandler(
      filename = function() "gene_biotype_distribution.pdf",
      content = function(file) save_static_plot(plot_type_static, file, "pdf")
    )
    
    output$download_type_svg <- downloadHandler(
      filename = function() "gene_biotype_distribution.svg",
      content = function(file) save_static_plot(plot_type_static, file, "svg")
    )
    
    output$download_chr_pdf <- downloadHandler(
      filename = function() "chromosome_distribution.pdf",
      content = function(file) save_static_plot(plot_chr_static, file, "pdf")
    )
    
    output$download_chr_svg <- downloadHandler(
      filename = function() "chromosome_distribution.svg",
      content = function(file) save_static_plot(plot_chr_static, file, "svg")
    )
    
    output$download_len_pdf <- downloadHandler(
      filename = function() "gene_length_distribution.pdf",
      content = function(file) save_static_plot(plot_len_static, file, "pdf")
    )
    
    output$download_len_svg <- downloadHandler(
      filename = function() "gene_length_distribution.svg",
      content = function(file) save_static_plot(plot_len_static, file, "svg")
    )
    
    output$download_ideo_pdf <- downloadHandler(
      filename = function() "genomic_location.pdf",
      content = function(file) save_static_plot(plot_ideo_static, file, "pdf")
    )
    
    output$download_ideo_svg <- downloadHandler(
      filename = function() "genomic_location.svg",
      content = function(file) save_static_plot(plot_ideo_static, file, "svg")
    )
    output$download_table <- downloadHandler(
      filename = function() {
        paste0("gene_annotation_", input$gene_set, "_", input$species, ".csv")
      },
      content = function(file) {
        write.csv(rv$annot, file, row.names = FALSE)
      }
    )
    
    return(reactive(rv$annot))
  })
}