# ==========================================================
# modules/bulk/report.R
# CoTRA Bulk RNA-seq Report Module
#
# Generates:
# - HTML report
# - PDF report using pagedown::chrome_print, no LaTeX needed
# - SVG figures packed into a ZIP file
# ==========================================================

mod_bulk_report_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    shinyjs::useShinyjs(),
    
    actionButton(ns("help_toggle"), icon("question-circle"), class = "btn btn-info btn-sm"),
    
    shinyjs::hidden(
      div(
        id = ns("help_box"),
        style = "border:1px solid #d9d9d9; padding:12px; margin-top:10px; margin-bottom:10px; border-radius:6px; background:#f7f7f7;",
        
        h4("Step-by-Step Guide"),
        
        h5("1. Complete upstream analysis"),
        p("Upload the count matrix, create groups, run differential expression, and optionally run annotation and enrichment."),
        
        h5("2. Select report sections"),
        p("Choose which analysis sections to include. Sections with missing data are skipped automatically."),
        
        h5("3. Generate report"),
        p("Click Generate Report. CoTRA creates an HTML report, converts it to PDF without LaTeX, and saves generated figures as SVG."),
        
        h5("4. Download outputs"),
        p("Download the PDF report, HTML report, and ZIP file containing SVG figures.")
      )
    ),
    
    br(),
    
    fluidRow(
      box(
        width = 12,
        title = "Report settings",
        solidHeader = TRUE,
        status = "primary",
        collapsible = TRUE,
        
        fluidRow(
          column(4, textInput(ns("project_title"), "Report title", value = "CoTRA Bulk RNA-seq Analysis Report")),
          column(4, textInput(ns("author"), "Author", value = "CoTRA user")),
          column(4, textInput(ns("project_name"), "Project name", value = "Bulk RNA-seq project"))
        ),
        
        fluidRow(
          column(3, numericInput(ns("top_genes_heatmap"), "Top genes for heatmap", value = 50, min = 10, max = 300, step = 10)),
          column(3, numericInput(ns("top_enrich_terms"), "Top enrichment terms", value = 20, min = 5, max = 50, step = 5)),
          column(3, numericInput(ns("padj_cutoff"), "Padj cutoff for report", value = 0.05, min = 0, max = 1, step = 0.005)),
          column(3, numericInput(ns("logfc_cutoff"), "Abs log2FC cutoff for report", value = 1, min = 0, max = 10, step = 0.1))
        ),
        
        fluidRow(
          column(
            12,
            checkboxGroupInput(
              ns("sections"),
              "Sections to include",
              choices = c(
                "Dataset summary" = "dataset",
                "Group summary" = "groups",
                "Quality control figures" = "qc",
                "Differential expression" = "de",
                "Gene annotation" = "annotation",
                "Enrichment analysis" = "enrichment"
              ),
              selected = c("dataset", "groups", "qc", "de", "annotation", "enrichment"),
              inline = TRUE
            )
          )
        ),
        
        fluidRow(
          column(3, actionButton(ns("generate_report"), "Generate Report", class = "btn btn-primary")),
          column(3, downloadButton(ns("download_pdf"), "Download PDF", class = "btn btn-success")),
          column(3, downloadButton(ns("download_svg_zip"), "Download SVG Figures ZIP", class = "btn btn-warning")),
          column(3, downloadButton(ns("download_html"), "Download HTML", class = "btn btn-info"))
        ),
        
        br(),
        verbatimTextOutput(ns("status"))
      )
    ),
    
  )
}

mod_bulk_report_server <- function(id,
                                   bulk_reactive,
                                   de_results,
                                   annot_reactive = reactive(NULL),
                                   enrich_state = NULL,
                                   groups_reactive = NULL) {
  moduleServer(id, function(input, output, session) {
    
    rv <- reactiveValues(
      report_dir = NULL,
      report_html = NULL,
      report_pdf = NULL,
      svg_zip = NULL,
      status = "No report generated yet."
    )
    
    observeEvent(input$help_toggle, {
      shinyjs::toggle(id = "help_box", anim = TRUE)
    })
    
    output$status <- renderText({
      rv$status
    })
    
    safe_name <- function(x) {
      x <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
      x <- gsub("_+", "_", x)
      x
    }
    
    html_clean <- function(x) {
      x <- as.character(x)
      x <- gsub("&", "&amp;", x, fixed = TRUE)
      x <- gsub("<", "&lt;", x, fixed = TRUE)
      x <- gsub(">", "&gt;", x, fixed = TRUE)
      x
    }
    
    wrap_text <- function(x, width = 40) {
      vapply(
        as.character(x),
        function(z) paste(strwrap(z, width = width), collapse = "
"),
        character(1)
      )
    }
    
    resolve_value <- function(x) {
      if (is.function(x)) return(x())
      x
    }
    
    get_matrix <- function() {
      x <- tryCatch(resolve_value(bulk_reactive), error = function(e) NULL)
      if (is.null(x)) return(NULL)
      
      if (is.matrix(x) || is.data.frame(x)) {
        return(as.data.frame(x, check.names = FALSE))
      }
      
      if (is.list(x)) {
        candidates <- c("final", "norm", "normalized", "normalized_counts", "counts", "raw", "data", "matrix", "expr", "expression")
        for (nm in candidates) {
          if (!is.null(x[[nm]])) {
            y <- tryCatch(resolve_value(x[[nm]]), error = function(e) NULL)
            if (is.matrix(y) || is.data.frame(y)) {
              return(as.data.frame(y, check.names = FALSE))
            }
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
      mat2 <- mat2[rowSums(is.na(mat2)) < ncol(mat2), , drop = FALSE]
      mat2[is.na(mat2)] <- 0
      mat2
    }
    
    get_groups_df <- function() {
      if (is.null(groups_reactive)) return(NULL)
      
      x <- tryCatch(resolve_value(groups_reactive), error = function(e) NULL)
      if (is.null(x)) return(NULL)
      
      if (is.data.frame(x)) return(x)
      
      if (is.list(x)) {
        candidates <- c("groups", "sample_info", "metadata", "meta", "coldata", "group_table")
        for (nm in candidates) {
          if (!is.null(x[[nm]])) {
            y <- tryCatch(resolve_value(x[[nm]]), error = function(e) NULL)
            if (is.data.frame(y)) return(y)
          }
        }
      }
      
      NULL
    }
    
    sample_group_map <- function(samples) {
      out <- data.frame(Sample = samples, Group = "Unassigned", stringsAsFactors = FALSE)
      gdf <- get_groups_df()
      
      if (is.null(gdf) || nrow(gdf) == 0) return(out)
      
      if (all(c("Group", "Samples") %in% colnames(gdf))) {
        for (i in seq_len(nrow(gdf))) {
          smp <- trimws(unlist(strsplit(as.character(gdf$Samples[i]), ";|,|[[:space:]]+")))
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
    
    get_de <- function() {
      df <- tryCatch(de_results(), error = function(e) NULL)
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
      df
    }
    
    get_annot <- function() {
      df <- tryCatch(annot_reactive(), error = function(e) NULL)
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
      df
    }
    
    get_enrich <- function() {
      if (is.null(enrich_state)) return(NULL)
      
      df <- tryCatch({
        if (is.list(enrich_state) && !is.null(enrich_state$ora_table)) {
          enrich_state$ora_table()
        } else if (is.function(enrich_state)) {
          enrich_state()
        } else {
          NULL
        }
      }, error = function(e) NULL)
      
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
      df
    }
    
    get_lfc_col <- function(df) {
      if (is.null(df)) return(NULL)
      if ("lfc" %in% colnames(df)) return("lfc")
      if ("log2FoldChange" %in% colnames(df)) return("log2FoldChange")
      NULL
    }
    
    significant_de <- function(df) {
      if (is.null(df)) return(NULL)
      
      lfc_col <- get_lfc_col(df)
      if (is.null(lfc_col) || !"padj" %in% colnames(df)) return(NULL)
      
      df$report_lfc <- suppressWarnings(as.numeric(df[[lfc_col]]))
      df$report_padj <- suppressWarnings(as.numeric(df$padj))
      
      df[
        !is.na(df$report_padj) &
          df$report_padj <= input$padj_cutoff &
          !is.na(df$report_lfc) &
          abs(df$report_lfc) >= input$logfc_cutoff,
        ,
        drop = FALSE
      ]
    }
    
    write_table <- function(df, file) {
      if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
        utils::write.csv(df, file, row.names = FALSE)
        return(file)
      }
      NULL
    }
    
    save_plot_all <- function(p, base_file, width = 9, height = 6) {
      pdf_file <- paste0(base_file, ".pdf")
      svg_file <- paste0(base_file, ".svg")
      png_file <- paste0(base_file, ".png")
      
      grDevices::pdf(pdf_file, width = width, height = height)
      print(p)
      grDevices::dev.off()
      
      svglite::svglite(svg_file, width = width, height = height)
      print(p)
      grDevices::dev.off()
      
      ggplot2::ggsave(png_file, plot = p, width = width, height = height, dpi = 300)
      
      list(pdf = pdf_file, svg = svg_file, png = png_file)
    }
    
    make_qc_plots <- function(mat, fig_dir) {
      out <- list()
      if (is.null(mat) || nrow(mat) < 2 || ncol(mat) < 2) return(out)
      
      mat_log <- log2(mat + 1)
      vars <- apply(mat_log, 1, stats::var, na.rm = TRUE)
      if (length(vars) < 2) return(out)
      
      keep <- order(vars, decreasing = TRUE)[seq_len(min(500, length(vars)))]
      mat_top <- mat_log[keep, , drop = FALSE]
      
      pca <- tryCatch(stats::prcomp(t(mat_top), scale. = TRUE), error = function(e) NULL)
      if (!is.null(pca) && ncol(pca$x) >= 2) {
        pcs <- as.data.frame(pca$x[, 1:2, drop = FALSE])
        pcs$Sample <- rownames(pcs)
        pcs <- dplyr::left_join(pcs, sample_group_map(pcs$Sample), by = "Sample")
        ve <- 100 * summary(pca)$importance[2, 1:2]
        
        p <- ggplot2::ggplot(pcs, ggplot2::aes(x = PC1, y = PC2, color = Group, label = Sample)) +
          ggplot2::geom_point(size = 3) +
          ggplot2::geom_text(vjust = -0.8, size = 3) +
          ggplot2::theme_bw() +
          ggplot2::labs(
            title = "PCA of top variable genes",
            x = paste0("PC1 (", round(ve[1], 1), "%)"),
            y = paste0("PC2 (", round(ve[2], 1), "%)")
          )
        out$qc_pca <- save_plot_all(p, file.path(fig_dir, "qc_pca"), 8, 6)
      }
      
      cor_mat <- tryCatch(stats::cor(mat_top), error = function(e) NULL)
      if (!is.null(cor_mat)) {
        cor_df <- as.data.frame(as.table(cor_mat))
        colnames(cor_df) <- c("Sample1", "Sample2", "Correlation")
        
        p <- ggplot2::ggplot(cor_df, ggplot2::aes(Sample1, Sample2, fill = Correlation)) +
          ggplot2::geom_tile() +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
          ggplot2::labs(title = "Sample correlation heatmap", x = "", y = "")
        out$qc_correlation <- save_plot_all(p, file.path(fig_dir, "qc_sample_correlation"), 8, 7)
      }
      
      hc <- tryCatch(stats::hclust(stats::dist(t(mat_top))), error = function(e) NULL)
      if (!is.null(hc)) {
        dend_pdf <- file.path(fig_dir, "qc_sample_clustering.pdf")
        dend_svg <- file.path(fig_dir, "qc_sample_clustering.svg")
        dend_png <- file.path(fig_dir, "qc_sample_clustering.png")
        
        grDevices::pdf(dend_pdf, width = 9, height = 6)
        plot(hc, main = "Sample hierarchical clustering", xlab = "", sub = "")
        grDevices::dev.off()
        
        svglite::svglite(dend_svg, width = 9, height = 6)
        plot(hc, main = "Sample hierarchical clustering", xlab = "", sub = "")
        grDevices::dev.off()
        
        grDevices::png(dend_png, width = 2400, height = 1600, res = 300)
        plot(hc, main = "Sample hierarchical clustering", xlab = "", sub = "")
        grDevices::dev.off()
        
        out$qc_clustering <- list(pdf = dend_pdf, svg = dend_svg, png = dend_png)
      }
      
      out
    }
    
    make_de_plots <- function(de_df, mat, fig_dir) {
      out <- list()
      if (is.null(de_df) || nrow(de_df) == 0) return(out)
      
      lfc_col <- get_lfc_col(de_df)
      if (is.null(lfc_col) || !"padj" %in% colnames(de_df)) return(out)
      
      df <- de_df
      df$lfc_report <- suppressWarnings(as.numeric(df[[lfc_col]]))
      df$padj_report <- suppressWarnings(as.numeric(df$padj))
      df$neglog10padj <- -log10(df$padj_report + 1e-300)
      df$Significance <- ifelse(
        !is.na(df$padj_report) &
          df$padj_report <= input$padj_cutoff &
          !is.na(df$lfc_report) &
          abs(df$lfc_report) >= input$logfc_cutoff,
        "Significant",
        "Not significant"
      )
      
      p <- ggplot2::ggplot(df, ggplot2::aes(x = lfc_report, y = neglog10padj, color = Significance)) +
        ggplot2::geom_point(alpha = 0.7, size = 1.5) +
        ggplot2::theme_bw() +
        ggplot2::labs(title = "Volcano plot", x = "log2 fold change", y = "-log10 adjusted p value")
      out$de_volcano <- save_plot_all(p, file.path(fig_dir, "de_volcano"), 8, 6)
      
      if ("baseMean" %in% colnames(df)) {
        df$baseMean_report <- suppressWarnings(as.numeric(df$baseMean))
        p <- ggplot2::ggplot(df, ggplot2::aes(x = log10(baseMean_report + 1), y = lfc_report, color = Significance)) +
          ggplot2::geom_point(alpha = 0.7, size = 1.5) +
          ggplot2::theme_bw() +
          ggplot2::labs(title = "MA plot", x = "log10 mean expression", y = "log2 fold change")
        out$de_ma <- save_plot_all(p, file.path(fig_dir, "de_ma"), 8, 6)
      }
      
      sig <- significant_de(de_df)
      if (!is.null(sig) && nrow(sig) > 0 && !is.null(mat) && "gene" %in% colnames(sig)) {
        genes <- intersect(as.character(sig$gene), rownames(mat))
        if (length(genes) > 0) {
          genes <- genes[seq_len(min(length(genes), input$top_genes_heatmap))]
        }
        
        if (length(genes) >= 2) {
          sub <- log2(mat[genes, , drop = FALSE] + 1)
          sub <- t(scale(t(sub)))
          sub[is.na(sub)] <- 0
          heat_df <- as.data.frame(as.table(sub))
          colnames(heat_df) <- c("Gene", "Sample", "Z")
          
          p <- ggplot2::ggplot(heat_df, ggplot2::aes(x = Sample, y = Gene, fill = Z)) +
            ggplot2::geom_tile() +
            ggplot2::theme_bw() +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
              axis.text.y = ggplot2::element_text(size = 6)
            ) +
            ggplot2::labs(title = paste0("Heatmap of top significant genes, n = ", length(genes)), x = "", y = "")
          out$de_heatmap <- save_plot_all(p, file.path(fig_dir, "de_heatmap"), 10, 10)
        }
      }
      
      out
    }
    
    make_annotation_plots <- function(ann_df, fig_dir) {
      out <- list()
      if (is.null(ann_df) || nrow(ann_df) == 0) return(out)
      
      if ("gene_biotype" %in% colnames(ann_df)) {
        df <- ann_df %>%
          dplyr::filter(!is.na(gene_biotype), nzchar(gene_biotype)) %>%
          dplyr::count(gene_biotype, name = "n") %>%
          dplyr::arrange(dplyr::desc(n))
        
        if (nrow(df) > 0) {
          p <- ggplot2::ggplot(df, ggplot2::aes(x = reorder(gene_biotype, n), y = n)) +
            ggplot2::geom_col() +
            ggplot2::coord_flip() +
            ggplot2::theme_bw() +
            ggplot2::labs(title = "Gene biotype distribution", x = "Gene biotype", y = "Gene count")
          out$annotation_biotype <- save_plot_all(p, file.path(fig_dir, "annotation_gene_biotype"), 9, 6)
        }
      }
      
      if ("chromosome" %in% colnames(ann_df)) {
        df <- ann_df %>%
          dplyr::filter(!is.na(chromosome)) %>%
          dplyr::count(chromosome, name = "n")
        
        if (nrow(df) > 0) {
          p <- ggplot2::ggplot(df, ggplot2::aes(x = chromosome, y = n)) +
            ggplot2::geom_col() +
            ggplot2::theme_bw() +
            ggplot2::labs(title = "Chromosome distribution", x = "Chromosome", y = "Gene count")
          out$annotation_chromosome <- save_plot_all(p, file.path(fig_dir, "annotation_chromosome_distribution"), 9, 6)
        }
      }
      
      if ("gene_length" %in% colnames(ann_df)) {
        df <- ann_df %>%
          dplyr::filter(!is.na(gene_length), is.finite(gene_length), gene_length > 0)
        
        if (nrow(df) > 0) {
          p <- ggplot2::ggplot(df, ggplot2::aes(x = gene_length)) +
            ggplot2::geom_histogram(bins = 40) +
            ggplot2::theme_bw() +
            ggplot2::labs(title = "Gene length distribution", x = "Gene length", y = "Count")
          out$annotation_gene_length <- save_plot_all(p, file.path(fig_dir, "annotation_gene_length"), 9, 6)
        }
      }
      
      if (all(c("chromosome", "start_position") %in% colnames(ann_df))) {
        df <- ann_df %>%
          dplyr::filter(!is.na(chromosome), !is.na(start_position))
        
        if (nrow(df) > 0) {
          p <- ggplot2::ggplot(df, ggplot2::aes(x = chromosome, y = start_position)) +
            ggplot2::geom_point(size = 1.5, alpha = 0.8) +
            ggplot2::theme_bw() +
            ggplot2::labs(title = "Genomic location of annotated genes", x = "Chromosome", y = "Start position")
          out$annotation_genomic_location <- save_plot_all(p, file.path(fig_dir, "annotation_genomic_location"), 9, 6)
        }
      }
      
      out
    }
    
    make_enrichment_plots <- function(enrich_df, fig_dir) {
      out <- list()
      if (is.null(enrich_df) || nrow(enrich_df) == 0 || !"Description" %in% colnames(enrich_df)) return(out)
      
      df <- enrich_df
      if (!"p.adjust" %in% colnames(df)) df$p.adjust <- NA_real_
      if (!"Count" %in% colnames(df)) df$Count <- NA_real_
      
      df$p.adjust <- suppressWarnings(as.numeric(df$p.adjust))
      df$Count <- suppressWarnings(as.numeric(df$Count))
      df <- df[order(df$p.adjust, decreasing = FALSE), , drop = FALSE]
      df <- head(df, input$top_enrich_terms)
      df$DescriptionWrapped <- factor(wrap_text(df$Description, 45), levels = rev(wrap_text(df$Description, 45)))
      
      p <- ggplot2::ggplot(df, ggplot2::aes(x = Count, y = DescriptionWrapped, size = Count, color = p.adjust)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::labs(title = "Top enrichment terms", x = "Gene count", y = "", size = "No. of Genes", color = "Adjusted p")
      out$enrichment_dotplot <- save_plot_all(p, file.path(fig_dir, "enrichment_dotplot"), 10, 8)
      
      out
    }
    
    img_tag <- function(fig, caption) {
      if (is.null(fig) || is.null(fig$png) || !file.exists(fig$png)) return("")
      
      src <- file.path("figures", basename(fig$png))
      src <- gsub("\\", "/", src, fixed = TRUE)
      
      paste0(
        "<figure>",
        "<img src='", src, "' style='max-width:100%; height:auto;'>",
        "<figcaption>", html_clean(caption), "</figcaption>",
        "</figure>"
      )
    }
    
    table_html <- function(df, max_rows = 20) {
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
        return("<p>No table available.</p>")
      }
      
      df <- utils::head(df, max_rows)
      df[] <- lapply(df, function(x) html_clean(as.character(x)))
      paste(capture.output(print(knitr::kable(df, format = "html"))), collapse = "
")
    }
    
    build_html <- function(html_file, data, figs) {
      de_sig <- data$de_sig
      enrich <- data$enrichment
      ann <- data$annotation
      mat <- data$mat
      groups <- data$groups
      sections <- input$sections
      
      body <- c(
        "<!doctype html><html><head><meta charset='utf-8'>",
        "<title>", html_clean(input$project_title), "</title>",
        "<style>",
        "body{font-family:Arial, sans-serif; margin:40px; line-height:1.45; color:#222;}",
        "h1,h2,h3{color:#1f4e79;}",
        "table{border-collapse:collapse; width:100%; font-size:11px; margin-bottom:20px;}",
        "th,td{border:1px solid #ddd; padding:5px; text-align:left;}",
        "th{background:#f2f2f2;}",
        "figure{page-break-inside:avoid; margin:20px 0;}",
        "figcaption{font-size:12px; color:#555; margin-top:4px;}",
        ".note{background:#f7f7f7; padding:10px; border-left:4px solid #1f4e79;}",
        "</style></head><body>",
        "<h1>", html_clean(input$project_title), "</h1>",
        "<p><b>Project:</b> ", html_clean(input$project_name), "<br>",
        "<b>Author:</b> ", html_clean(input$author), "<br>",
        "<b>Generated:</b> ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "</p>",
        "<div class='note'>Padj cutoff: ", input$padj_cutoff, ", Abs log2FC cutoff: ", input$logfc_cutoff, "</div>"
      )
      
      if ("dataset" %in% sections) {
        body <- c(body, "<h2>Dataset summary</h2>")
        if (!is.null(mat)) {
          body <- c(
            body,
            "<ul>",
            "<li>Genes: ", nrow(mat), "</li>",
            "<li>Samples: ", ncol(mat), "</li>",
            "</ul>"
          )
        } else {
          body <- c(body, "<p>No count matrix available.</p>")
        }
      }
      
      if ("groups" %in% sections) {
        body <- c(body, "<h2>Group summary</h2>", table_html(groups, 50))
      }
      
      if ("qc" %in% sections) {
        body <- c(
          body,
          "<h2>Quality control</h2>",
          img_tag(figs$qc_pca, "PCA of top variable genes."),
          img_tag(figs$qc_correlation, "Sample correlation heatmap."),
          img_tag(figs$qc_clustering, "Sample hierarchical clustering.")
        )
      }
      
      if ("de" %in% sections) {
        body <- c(
          body,
          "<h2>Differential expression</h2>",
          "<ul>",
          "<li>Total DE rows: ", ifelse(is.null(data$de), "Not available", nrow(data$de)), "</li>",
          "<li>Significant genes under report thresholds: ", ifelse(is.null(de_sig), "Not available", nrow(de_sig)), "</li>",
          "</ul>",
          img_tag(figs$de_volcano, "Volcano plot."),
          img_tag(figs$de_ma, "MA plot."),
          img_tag(figs$de_heatmap, "Heatmap of top significant genes."),
          "<h3>Top significant genes</h3>",
          table_html(de_sig, 20)
        )
      }
      
      if ("annotation" %in% sections) {
        body <- c(
          body,
          "<h2>Gene annotation</h2>",
          "<p>Annotated rows: ", ifelse(is.null(ann), "Not available", nrow(ann)), "</p>",
          img_tag(figs$annotation_biotype, "Gene biotype distribution."),
          img_tag(figs$annotation_chromosome, "Chromosome distribution."),
          img_tag(figs$annotation_gene_length, "Gene length distribution."),
          img_tag(figs$annotation_genomic_location, "Genomic location of annotated genes."),
          "<h3>Annotation table preview</h3>",
          table_html(ann, 20)
        )
      }
      
      if ("enrichment" %in% sections) {
        body <- c(
          body,
          "<h2>Enrichment analysis</h2>",
          "<p>Enrichment rows: ", ifelse(is.null(enrich), "Not available", nrow(enrich)), "</p>",
          img_tag(figs$enrichment_dotplot, "Top enrichment terms."),
          "<h3>Enrichment table preview</h3>",
          table_html(enrich, 20)
        )
      }
      
      body <- c(
        body,
        "<h2>Session information</h2><pre>",
        html_clean(paste(capture.output(sessionInfo()), collapse = "
")),
        "</pre>",
        "</body></html>"
      )
      
      writeLines(body, html_file)
      html_file
    }
    
    generate_report_files <- function() {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      report_prefix <- paste0("CoTRA_bulk_report_", timestamp)
      report_dir <- file.path("outputs", report_prefix)
      if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
      
      fig_dir <- file.path(report_dir, "figures")
      table_dir <- file.path(report_dir, "tables")
      dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
      
      mat <- clean_matrix(get_matrix())
      groups <- get_groups_df()
      de <- get_de()
      de_sig <- significant_de(de)
      ann <- get_annot()
      enrich <- get_enrich()
      
      write_table(groups, file.path(table_dir, "groups.csv"))
      write_table(de, file.path(table_dir, "de_all.csv"))
      write_table(de_sig, file.path(table_dir, "de_significant.csv"))
      write_table(ann, file.path(table_dir, "annotation.csv"))
      write_table(enrich, file.path(table_dir, "enrichment.csv"))
      
      figs <- list()
      if ("qc" %in% input$sections) figs <- c(figs, make_qc_plots(mat, fig_dir))
      if ("de" %in% input$sections) figs <- c(figs, make_de_plots(de, mat, fig_dir))
      if ("annotation" %in% input$sections) figs <- c(figs, make_annotation_plots(ann, fig_dir))
      if ("enrichment" %in% input$sections) figs <- c(figs, make_enrichment_plots(enrich, fig_dir))
      
      data <- list(
        mat = mat,
        groups = groups,
        de = de,
        de_sig = de_sig,
        annotation = ann,
        enrichment = enrich
      )
      
      html_file <- file.path(report_dir, paste0(report_prefix, ".html"))
      pdf_file <- file.path(report_dir, paste0(report_prefix, ".pdf"))
      build_html(html_file, data, figs)
      
      pdf_res <- tryCatch(
        pagedown::chrome_print(input = html_file, output = pdf_file),
        error = function(e) e
      )
      
      if (inherits(pdf_res, "error")) {
        stop(paste("PDF export failed. Install pagedown and webshot2, and ensure Chrome or Chromium is available. Error:", pdf_res$message))
      }
      
      svg_files <- unlist(lapply(figs, function(x) x$svg), use.names = FALSE)
      svg_files <- svg_files[file.exists(svg_files)]
      zip_file <- file.path(report_dir, paste0(report_prefix, "_svg_figures.zip"))
      
      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(fig_dir)
      
      if (length(svg_files) > 0) {
        utils::zip(zipfile = zip_file, files = basename(svg_files))
      } else {
        writeLines("No SVG figures were generated.", "no_figures.txt")
        utils::zip(zipfile = zip_file, files = "no_figures.txt")
      }
      
      list(
        report_dir = report_dir,
        html = html_file,
        pdf = pdf_file,
        zip = zip_file,
        n_figures = length(svg_files)
      )
    }
    
    observeEvent(input$generate_report, {
      rv$status <- "Generating report."
      
      res <- tryCatch(generate_report_files(), error = function(e) e)
      
      if (inherits(res, "error")) {
        rv$report_dir <- NULL
        rv$report_html <- NULL
        rv$report_pdf <- NULL
        rv$svg_zip <- NULL
        rv$status <- paste("Report generation failed:", res$message)
        return(NULL)
      }
      
      rv$report_dir <- res$report_dir
      rv$report_html <- res$html
      rv$report_pdf <- res$pdf
      rv$svg_zip <- res$zip
      rv$status <- paste0("Report generated successfully. SVG figures: ", res$n_figures)
    })
    
    output$download_pdf <- downloadHandler(
      filename = function() {
        req(rv$report_pdf)
        basename(rv$report_pdf)
      },
      content = function(file) {
        req(rv$report_pdf)
        file.copy(rv$report_pdf, file, overwrite = TRUE)
      }
    )
    
    output$download_html <- downloadHandler(
      filename = function() {
        req(rv$report_html)
        basename(rv$report_html)
      },
      content = function(file) {
        req(rv$report_html)
        file.copy(rv$report_html, file, overwrite = TRUE)
      }
    )
    
    output$download_svg_zip <- downloadHandler(
      filename = function() {
        req(rv$svg_zip)
        basename(rv$svg_zip)
      },
      content = function(file) {
        req(rv$svg_zip)
        file.copy(rv$svg_zip, file, overwrite = TRUE)
      },
      contentType = "application/zip"
    )
    
    return(list(
      report_html = reactive(rv$report_html),
      report_pdf = reactive(rv$report_pdf),
      svg_zip = reactive(rv$svg_zip),
      report_dir = reactive(rv$report_dir)
    ))
  })
}
