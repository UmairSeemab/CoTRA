# ==========================================================
# modules/sc/report.R
# CoTRA scRNA-seq Report Module
#
# Generates:
# - HTML report
# - PDF report using pagedown, no LaTeX needed
# - ZIP export
#
# Fix included:
# - Uses relative figure paths in Rmd/HTML:
#   figures/pca.png
#   figures/umap.png
#   figures/qc_violin.png
#
# Return:
# - report_ready
# - report_path
# - report_pdf
# ==========================================================

suppressPackageStartupMessages({
  library(shiny)
  library(Seurat)
  library(ggplot2)
  library(DT)
})

mod_sc_report_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("scRNA-seq Report"),
    
    div(
      class = "alert alert-info",
      strong("Purpose: "),
      "This module generates an HTML and PDF summary report of completed single-cell RNA-seq analysis steps in CoTRA."
    ),
    
    sidebarLayout(
      sidebarPanel(
        textInput(
          ns("report_title"),
          "Report title",
          value = "CoTRA scRNA-seq Analysis Report"
        ),
        
        textInput(
          ns("project_name"),
          "Project name",
          value = "CoTRA scRNA-seq project"
        ),
        
        textAreaInput(
          ns("report_notes"),
          "User notes",
          value = "",
          rows = 5,
          placeholder = "Optional notes about samples, organism, experimental design, or analysis choices."
        ),
        
        checkboxGroupInput(
          ns("include_sections"),
          "Sections to include",
          choices = c(
            "Dataset summary" = "dataset",
            "Quality control" = "qc",
            "Gene selection" = "gene_selection",
            "PCA" = "pca",
            "t-SNE / UMAP" = "tsne_umap",
            "Clustering" = "clustering",
            "Markers" = "markers",
            "Annotation" = "annotation",
            "Differential expression" = "differential_expression",
            "Trajectory" = "trajectory",
            "Pathway activity" = "pathway",
            "Cell communication" = "communication",
            "Session information" = "session"
          ),
          selected = c(
            "dataset",
            "qc",
            "gene_selection",
            "pca",
            "tsne_umap",
            "clustering",
            "markers",
            "annotation",
            "differential_expression",
            "trajectory",
            "pathway",
            "communication",
            "session"
          )
        ),
        
        checkboxInput(
          ns("include_figures"),
          "Include available figures",
          value = TRUE
        ),
        
        checkboxInput(
          ns("include_tables"),
          "Include available tables",
          value = TRUE
        ),
        
        actionButton(
          ns("generate_report"),
          "Generate report",
          class = "btn-primary"
        ),
        
        br(), br(),
        downloadButton(ns("download_html"), "Download HTML report"),
        br(), br(),
        downloadButton(ns("download_pdf"), "Download PDF report"),
        br(), br(),
        downloadButton(ns("download_zip"), "Download report ZIP")
      ),
      
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Summary",
            br(),
            htmlOutput(ns("report_summary")),
            br(),
            DTOutput(ns("module_status_table"))
          ),
          
          tabPanel(
            "Report preview",
            br(),
            uiOutput(ns("report_preview"))
          ),
          
          tabPanel(
            "Help",
            br(),
            h4("When to use this module"),
            p("Use this module after completing one or more scRNA-seq analysis steps. The report includes only available results."),
            
            h4("Recommended workflow"),
            tags$ol(
              tags$li("Import data"),
              tags$li("Run quality control"),
              tags$li("Select variable genes"),
              tags$li("Run PCA"),
              tags$li("Run UMAP or t-SNE"),
              tags$li("Run clustering"),
              tags$li("Find markers"),
              tags$li("Run annotation"),
              tags$li("Optionally compare conditions using differential expression"),
              tags$li("Optionally run trajectory, pathway activity, and cell communication"),
              tags$li("Generate report")
            ),
            
            h4("Output files"),
            tags$ul(
              tags$li("HTML report for browser viewing."),
              tags$li("PDF report for sharing and archiving."),
              tags$li("ZIP file containing report outputs and figures.")
            ),
            
            h4("PDF export"),
            p("PDF export requires pagedown and Chrome or Chromium.")
          )
        )
      )
    )
  )
}

mod_sc_report_server <- function(
    id,
    seurat_r,
    sc_state = NULL,
    sc_markers = NULL,
    sc_annot = NULL,
    sc_de = NULL,
    sc_trajectory = NULL,
    sc_pathway = NULL,
    sc_comm = NULL,
    output_dir = NULL
) {
  moduleServer(id, function(input, output, session) {
    
    report_html_path <- reactiveVal(NULL)
    report_pdf_path  <- reactiveVal(NULL)
    report_dir_path  <- reactiveVal(NULL)
    report_ready     <- reactiveVal(FALSE)
    pdf_error_msg    <- reactiveVal(NULL)
    
    effective_output_dir <- reactive({
      if (!is.null(output_dir) && nzchar(as.character(output_dir)[1])) {
        return(cotra_ensure_output_dir(output_dir))
      }
      cotra_module_output_dir("scRNA", "Reports")
    })
    
    get_obj <- reactive({
      obj <- seurat_r()
      validate(need(inherits(obj, "Seurat"), "Input must be a Seurat object."))
      obj
    })
    
    safe_value <- function(expr, default = NULL) {
      tryCatch(expr, error = function(e) default)
    }
    
    get_module_status <- reactive({
      obj <- safe_value(get_obj(), NULL)
      
      md_cols <- if (!is.null(obj)) colnames(obj@meta.data) else character(0)
      reductions <- if (!is.null(obj)) names(obj@reductions) else character(0)
      de_ready <- !is.null(safe_value(sc_de$results(), NULL))
      
      data.frame(
        Module = c(
          "Import",
          "Quality Control",
          "Gene Selection",
          "PCA",
          "t-SNE / UMAP",
          "Clustering",
          "Markers",
          "Annotation",
          "Differential Expression",
          "Trajectory",
          "Pathway Activity",
          "Cell Communication"
        ),
        Status = c(
          if (!is.null(obj)) "Available" else "Missing",
          if (any(c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.mito") %in% md_cols)) "Available" else "Missing",
          if (!is.null(obj) && length(VariableFeatures(obj)) > 0) "Available" else "Missing",
          if ("pca" %in% reductions) "Available" else "Missing",
          if (any(c("umap", "tsne") %in% reductions)) "Available" else "Missing",
          if (any(c("seurat_clusters", "CoTRA_cluster_label") %in% md_cols)) "Available" else "Missing",
          if (!is.null(sc_markers)) "Available if completed" else "Not connected",
          if ("CoTRA_celltype" %in% md_cols) "Available" else "Missing",
          if (isTRUE(de_ready)) "Available" else if (!is.null(sc_de)) "Not completed" else "Not connected",
          if (!is.null(sc_trajectory)) "Available if completed" else "Not connected",
          if (!is.null(sc_pathway)) "Available if completed" else "Not connected",
          if (!is.null(sc_comm)) "Available if completed" else "Not connected"
        ),
        stringsAsFactors = FALSE
      )
    })
    
    output$module_status_table <- renderDT({
      datatable(
        get_module_status(),
        options = list(dom = "t", paging = FALSE),
        rownames = FALSE
      )
    })
    
    output$report_summary <- renderUI({
      if (!isTRUE(report_ready())) {
        return(HTML(
          "<div class='alert alert-warning'>No report generated yet. Click <b>Generate report</b>.</div>"
        ))
      }
      
      pdf_text <- if (!is.null(report_pdf_path()) && file.exists(report_pdf_path())) {
        report_pdf_path()
      } else {
        paste0("PDF not generated", if (!is.null(pdf_error_msg())) paste0("<br><b>PDF error:</b> ", pdf_error_msg()) else "")
      }
      
      HTML(paste0(
        "<div class='alert alert-success'>",
        "<b>Report generated successfully.</b><br>",
        "<b>HTML:</b> ", report_html_path(), "<br>",
        "<b>PDF:</b> ", pdf_text,
        "</div>"
      ))
    })
    
    output$report_preview <- renderUI({
      req(report_html_path())
      
      tags$iframe(
        src = report_html_path(),
        width = "100%",
        height = "800px",
        style = "border:1px solid #ccc;"
      )
    })
    
    observeEvent(input$generate_report, {
      obj <- get_obj()
      
      report_timestamp <- cotra_timestamp()
      run_dir <- file.path(effective_output_dir(), paste0("CoTRA_scRNA_Report_", report_timestamp))
      fig_dir <- file.path(run_dir, "figures")
      
      dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
      
      html_file <- file.path(run_dir, cotra_file_name("CoTRA_scRNA_Report", "html", timestamp = report_timestamp))
      pdf_file  <- file.path(run_dir, cotra_file_name("CoTRA_scRNA_Report", "pdf", timestamp = report_timestamp))
      rmd_file  <- file.path(run_dir, cotra_file_name("CoTRA_scRNA_Report_Template", "Rmd", timestamp = report_timestamp))
      
      withProgress(message = "Generating scRNA report", value = 0, {
        
        incProgress(0.2, detail = "Collecting report data")
        
        report_data <- collect_sc_report_data(
          obj = obj,
          sc_markers = sc_markers,
          sc_annot = sc_annot,
          sc_de = sc_de,
          sc_trajectory = sc_trajectory,
          sc_pathway = sc_pathway,
          sc_comm = sc_comm
        )

        if (!is.null(report_data$de_results) && nrow(report_data$de_results) > 0) {
          utils::write.csv(
            report_data$de_results,
            file.path(run_dir, "scRNA_differential_expression_results.csv"),
            row.names = FALSE
          )
        }
        
        incProgress(0.45, detail = "Saving figures")
        
        figure_files <- create_sc_report_figures(
          obj = obj,
          fig_dir = fig_dir,
          include_figures = isTRUE(input$include_figures),
          sc_de = sc_de
        )
        
        incProgress(0.6, detail = "Writing report template")
        
        write_sc_report_rmd(
          file = rmd_file,
          title = input$report_title,
          project_name = input$project_name,
          notes = input$report_notes,
          sections = input$include_sections,
          report_data = report_data,
          figure_files = figure_files,
          include_tables = isTRUE(input$include_tables)
        )
        
        incProgress(0.75, detail = "Rendering HTML report")
        
        validate(
          need(requireNamespace("rmarkdown", quietly = TRUE), "Install rmarkdown to generate reports.")
        )
        
        rmarkdown::render(
          input = rmd_file,
          output_file = basename(html_file),
          output_dir = run_dir,
          quiet = TRUE,
          envir = new.env(parent = globalenv())
        )
        
        incProgress(0.9, detail = "Rendering PDF report")
        
        pdf_error_msg(NULL)
        
        if (requireNamespace("pagedown", quietly = TRUE)) {
          tryCatch({
            pagedown::chrome_print(
              input = normalizePath(html_file),
              output = normalizePath(pdf_file, mustWork = FALSE)
            )
          }, error = function(e) {
            pdf_error_msg(e$message)
            pdf_file <<- NULL
          })
        } else {
          pdf_error_msg("Package pagedown is not installed.")
          pdf_file <- NULL
        }
        
        report_html_path(html_file)
        report_pdf_path(pdf_file)
        report_dir_path(run_dir)
        report_ready(TRUE)
        
        if (!is.null(sc_state) && is.environment(sc_state)) {
          sc_state$sc_report <- list(
            html = html_file,
            pdf = pdf_file,
            folder = run_dir,
            created = Sys.time()
          )
        }
        
        incProgress(1, detail = "Done")
      })
    })
    
    output$download_html <- downloadHandler(
      filename = function() {
        req(report_html_path())
        basename(report_html_path())
      },
      content = function(file) {
        req(report_html_path())
        file.copy(report_html_path(), file, overwrite = TRUE)
      }
    )
    
    output$download_pdf <- downloadHandler(
      filename = function() {
        req(report_pdf_path())
        basename(report_pdf_path())
      },
      content = function(file) {
        req(report_pdf_path())
        file.copy(report_pdf_path(), file, overwrite = TRUE)
      }
    )
    
    output$download_zip <- downloadHandler(
      filename = function() {
        cotra_file_name("CoTRA_scRNA_Report_AllFiles", "zip")
      },
      content = function(file) {
        req(report_dir_path())
        
        oldwd <- getwd()
        on.exit(setwd(oldwd), add = TRUE)
        
        setwd(report_dir_path())
        utils::zip(
          zipfile = file,
          files = list.files(report_dir_path(), recursive = TRUE)
        )
      }
    )
    
    return(list(
      report_ready = reactive({
        report_ready()
      }),
      report_path = reactive({
        report_html_path()
      }),
      report_pdf = reactive({
        report_pdf_path()
      })
    ))
  })
}

# ==========================================================
# Helper functions
# ==========================================================

collect_sc_report_data <- function(
    obj,
    sc_markers = NULL,
    sc_annot = NULL,
    sc_de = NULL,
    sc_trajectory = NULL,
    sc_pathway = NULL,
    sc_comm = NULL
) {
  md <- obj@meta.data
  
  marker_table <- tryCatch(sc_markers$markers_table(), error = function(e) NULL)
  annotation_table <- tryCatch(sc_annot$annotation_table(), error = function(e) NULL)
  de_results <- tryCatch(sc_de$results(), error = function(e) NULL)
  de_info <- tryCatch(sc_de$analysis_info(), error = function(e) NULL)
  de_design <- tryCatch(sc_de$design(), error = function(e) NULL)
  trajectory_path <- tryCatch(sc_trajectory$session_path(), error = function(e) NULL)
  pathway_ready <- tryCatch(sc_pathway$pathway_ready(), error = function(e) NULL)
  communication_ready <- tryCatch(sc_comm$comm_ready(), error = function(e) NULL)
  
  list(
    n_cells = ncol(obj),
    n_genes = nrow(obj),
    assay = DefaultAssay(obj),
    reductions = names(obj@reductions),
    metadata_columns = colnames(md),
    variable_features = length(VariableFeatures(obj)),
    clusters = if ("seurat_clusters" %in% colnames(md)) length(unique(md$seurat_clusters)) else NA,
    celltypes = if ("CoTRA_celltype" %in% colnames(md)) length(unique(md$CoTRA_celltype)) else NA,
    marker_table = marker_table,
    annotation_table = annotation_table,
    de_results = de_results,
    de_info = de_info,
    de_design = de_design,
    trajectory_path = trajectory_path,
    pathway_ready = pathway_ready,
    communication_ready = communication_ready
  )
}

create_sc_report_figures <- function(obj, fig_dir, include_figures = TRUE, sc_de = NULL) {
  if (!isTRUE(include_figures)) {
    return(list())
  }
  
  figs <- list()
  md <- obj@meta.data
  
  if (all(c("nFeature_RNA", "nCount_RNA") %in% colnames(md))) {
    features <- c("nFeature_RNA", "nCount_RNA")
    
    if ("percent.mt" %in% colnames(md)) {
      features <- c(features, "percent.mt")
    } else if ("percent.mito" %in% colnames(md)) {
      features <- c(features, "percent.mito")
    }
    
    p <- Seurat::VlnPlot(
      obj,
      features = features,
      ncol = min(length(features), 3),
      pt.size = 0
    ) + ggplot2::ggtitle("Quality control metrics")
    
    f <- file.path(fig_dir, cotra_file_name("scRNA_QC_ViolinPlot", "png"))
    ggplot2::ggsave(f, p, width = 10, height = 5, dpi = 150)
    figs$qc <- f
  }
  
  if ("pca" %in% names(obj@reductions)) {
    p <- Seurat::DimPlot(obj, reduction = "pca") +
      ggplot2::ggtitle("PCA")
    
    f <- file.path(fig_dir, cotra_file_name("scRNA_PCA_DimPlot", "png"))
    ggplot2::ggsave(f, p, width = 8, height = 6, dpi = 150)
    figs$pca <- f
  }
  
  if ("umap" %in% names(obj@reductions)) {
    group_col <- if ("CoTRA_celltype" %in% colnames(md)) {
      "CoTRA_celltype"
    } else if ("CoTRA_cluster_label" %in% colnames(md)) {
      "CoTRA_cluster_label"
    } else if ("seurat_clusters" %in% colnames(md)) {
      "seurat_clusters"
    } else {
      NULL
    }
    
    p <- if (!is.null(group_col)) {
      Seurat::DimPlot(obj, reduction = "umap", group.by = group_col, label = TRUE) +
        ggplot2::ggtitle(paste("UMAP colored by", group_col))
    } else {
      Seurat::DimPlot(obj, reduction = "umap") +
        ggplot2::ggtitle("UMAP")
    }
    
    f <- file.path(fig_dir, cotra_file_name("scRNA_UMAP_DimPlot", "png"))
    ggplot2::ggsave(f, p, width = 8, height = 6, dpi = 150)
    figs$umap <- f
  }
  
  if ("tsne" %in% names(obj@reductions)) {
    p <- Seurat::DimPlot(obj, reduction = "tsne") +
      ggplot2::ggtitle("t-SNE")
    
    f <- file.path(fig_dir, cotra_file_name("scRNA_tSNE_DimPlot", "png"))
    ggplot2::ggsave(f, p, width = 8, height = 6, dpi = 150)
    figs$tsne <- f
  }

  de_results <- tryCatch(sc_de$results(), error = function(e) NULL)
  de_info <- tryCatch(sc_de$analysis_info(), error = function(e) NULL)
  
  if (!is.null(de_results) && nrow(de_results) > 0 && !is.null(de_info)) {
    plot_df <- de_results
    plot_df$plot_p <- plot_df$p_adj
    plot_df$plot_p[is.na(plot_df$plot_p)] <- 1
    plot_df$plot_p <- pmax(plot_df$plot_p, .Machine$double.xmin)
    plot_df$minus_log10_padj <- -log10(plot_df$plot_p)
    plot_df$hover_text <- paste0(
      "Gene: ", plot_df$gene,
      "<br>log2FC: ", signif(plot_df$log2FC, 4),
      "<br>Adjusted P-value: ", format(plot_df$p_adj, digits = 4, scientific = TRUE),
      "<br>Direction: ", plot_df$direction
    )
    
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = log2FC, y = minus_log10_padj, colour = direction, text = hover_text)
    ) +
      ggplot2::geom_point(alpha = 0.70, size = 1.5, na.rm = TRUE) +
      ggplot2::geom_vline(
        xintercept = c(-de_info$log2fc_cutoff, de_info$log2fc_cutoff),
        linetype = 2,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(de_info$fdr_cutoff),
        linetype = 2,
        inherit.aes = FALSE
      ) +
      ggplot2::labs(
        title = paste0(de_info$group1, " versus ", de_info$group2),
        subtitle = de_info$method,
        x = paste0("log2 fold change, positive = higher in ", de_info$group1),
        y = "-log10 adjusted P-value",
        colour = "Result"
      ) +
      ggplot2::theme_bw(base_size = 12)
    
    f <- file.path(fig_dir, cotra_file_name("scRNA_DE_VolcanoPlot", "png"))
    ggplot2::ggsave(f, p, width = 8, height = 7, dpi = 300)
    figs$de_volcano <- f

    if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
      widget <- plotly::ggplotly(p, tooltip = "text", dynamicTicks = TRUE) |>
        plotly::layout(
          hovermode = "closest",
          legend = list(orientation = "h", x = 0, y = -0.18),
          margin = list(b = 100)
        ) |>
        plotly::config(displaylogo = FALSE, responsive = TRUE)

      f_html <- file.path(fig_dir, cotra_file_name("scRNA_DE_InteractiveVolcanoPlot", "html"))
      widget_saved <- tryCatch({
        htmlwidgets::saveWidget(widget, f_html, selfcontained = TRUE)
        file.exists(f_html)
      }, error = function(e) FALSE)
      if (isTRUE(widget_saved)) figs$de_volcano_interactive <- f_html
    }
  }
  
  figs
}

write_sc_report_rmd <- function(
    file,
    title,
    project_name,
    notes,
    sections,
    report_data,
    figure_files,
    include_tables = TRUE
) {
  q <- function(x) {
    if (is.null(x) || length(x) == 0) return("")
    paste(x, collapse = ", ")
  }
  
  report_img <- function(x) {
    if (is.null(x) || !file.exists(x)) return(NULL)
    paste0("![](figures/", basename(x), ")")
  }

  report_iframe <- function(x) {
    if (is.null(x) || !file.exists(x)) return(NULL)
    paste0(
      '<iframe src="figures/', basename(x),
      '" width="100%" height="680" style="border:1px solid #ddd;"></iframe>'
    )
  }
  
  esc <- function(x) {
    x <- as.character(x)
    x <- gsub("\\\\", "/", x)
    x <- gsub('"', "'", x)
    x
  }
  
  markdown_table <- function(df, max_rows = 20) {
    if (is.null(df) || nrow(df) == 0) return(character(0))
    preferred <- c("gene", "log2FC", "pct.1", "pct.2", "p_value", "p_adj", "direction")
    cols <- preferred[preferred %in% colnames(df)]
    if (length(cols) == 0) cols <- head(colnames(df), 7)
    tab <- head(df[, cols, drop = FALSE], max_rows)

    format_cell <- function(x) {
      if (length(x) == 0 || is.na(x)) return("")
      if (is.numeric(x)) {
        if (abs(x) > 0 && abs(x) < 0.001) return(format(x, scientific = TRUE, digits = 3))
        return(format(round(x, 4), trim = TRUE, scientific = FALSE))
      }
      gsub("\\|", "\\\\|", as.character(x))
    }

    header <- paste0("| ", paste(cols, collapse = " | "), " |")
    separator <- paste0("| ", paste(rep("---", length(cols)), collapse = " | "), " |")
    rows <- apply(tab, 1, function(row) {
      paste0("| ", paste(vapply(row, format_cell, character(1)), collapse = " | "), " |")
    })
    c(header, separator, rows)
  }
  
  lines <- c(
    "---",
    paste0("title: \"", esc(title), "\""),
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_depth: 3",
    "    theme: flatly",
    "    self_contained: false",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "```",
    "",
    paste0("# ", esc(project_name)),
    "",
    paste0("Generated: ", Sys.time()),
    "",
    "## User notes",
    "",
    ifelse(nzchar(notes), esc(notes), "No user notes provided."),
    ""
  )
  
  if ("dataset" %in% sections) {
    lines <- c(
      lines,
      "## Dataset summary",
      "",
      paste0("- Cells: ", report_data$n_cells),
      paste0("- Genes: ", report_data$n_genes),
      paste0("- Default assay: ", report_data$assay),
      paste0("- Reductions: ", q(report_data$reductions)),
      paste0("- Metadata columns: ", q(report_data$metadata_columns)),
      ""
    )
  }
  
  if ("qc" %in% sections) {
    img <- report_img(figure_files$qc)
    
    lines <- c(
      lines,
      "## Quality control",
      "",
      "QC metrics summarize library size, detected genes, and optional mitochondrial content.",
      "",
      if (!is.null(img)) img else "No QC figure available.",
      ""
    )
  }
  
  if ("gene_selection" %in% sections) {
    lines <- c(
      lines,
      "## Gene selection",
      "",
      paste0("Number of variable genes: ", report_data$variable_features),
      ""
    )
  }
  
  if ("pca" %in% sections) {
    img <- report_img(figure_files$pca)
    
    lines <- c(
      lines,
      "## PCA",
      "",
      if (!is.null(img)) img else "PCA figure not available.",
      ""
    )
  }
  
  if ("tsne_umap" %in% sections) {
    umap_img <- report_img(figure_files$umap)
    tsne_img <- report_img(figure_files$tsne)
    
    lines <- c(
      lines,
      "## t-SNE / UMAP",
      "",
      if (!is.null(umap_img)) umap_img else "UMAP figure not available.",
      "",
      if (!is.null(tsne_img)) tsne_img else "",
      ""
    )
  }
  
  if ("clustering" %in% sections) {
    lines <- c(
      lines,
      "## Clustering",
      "",
      paste0("Detected clusters: ", report_data$clusters),
      ""
    )
  }
  
  if ("markers" %in% sections) {
    lines <- c(
      lines,
      "## Marker genes",
      "",
      if (!is.null(report_data$marker_table)) {
        "Marker table was available from the marker module."
      } else {
        "Marker table was not available."
      },
      ""
    )
  }
  
  if ("annotation" %in% sections) {
    lines <- c(
      lines,
      "## Cell type annotation",
      "",
      paste0("Detected cell types: ", report_data$celltypes),
      "",
      if (!is.null(report_data$annotation_table)) {
        "Annotation table was available from the annotation module."
      } else {
        "Annotation table was not available."
      },
      ""
    )
  }
  
  if ("differential_expression" %in% sections) {
    de_img <- report_img(figure_files$de_volcano)
    de_interactive <- report_iframe(figure_files$de_volcano_interactive)
    de_info <- report_data$de_info
    de_results <- report_data$de_results

    if (!is.null(de_info) && !is.null(de_results) && nrow(de_results) > 0) {
      significant <- sum(
        !is.na(de_results$p_adj) &
          de_results$p_adj <= de_info$fdr_cutoff &
          abs(de_results$log2FC) >= de_info$log2fc_cutoff,
        na.rm = TRUE
      )
      higher_group1 <- sum(de_results$direction == paste0("Higher in ", de_info$group1), na.rm = TRUE)
      higher_group2 <- sum(de_results$direction == paste0("Higher in ", de_info$group2), na.rm = TRUE)

      lines <- c(
        lines,
        "## Differential expression",
        "",
        paste0("- Method: ", esc(de_info$method)),
        paste0("- Contrast: ", esc(de_info$group1), " versus ", esc(de_info$group2)),
        paste0("- Subset: ", esc(de_info$subset_value)),
        paste0("- Cells included: ", de_info$selected_cells),
        paste0("- Biological samples represented: ", de_info$samples),
        paste0("- Genes tested: ", de_info$tested_genes),
        paste0("- Significant genes: ", significant),
        paste0("- Higher in ", esc(de_info$group1), ": ", higher_group1),
        paste0("- Higher in ", esc(de_info$group2), ": ", higher_group2),
        "",
        if (!is.null(de_info$warning) && nzchar(de_info$warning)) paste0("Interpretation: ", esc(de_info$warning)) else "",
        "",
        if (!is.null(de_interactive)) "### Interactive volcano plot" else "",
        if (!is.null(de_interactive)) de_interactive else "",
        "",
        "### Static volcano plot",
        "",
        if (!is.null(de_img)) de_img else "Differential-expression volcano figure not available.",
        ""
      )

      if (isTRUE(include_tables)) {
        lines <- c(
          lines,
          "### Top differential-expression results",
          "",
          markdown_table(de_results, max_rows = 20),
          ""
        )
      }
    } else {
      lines <- c(
        lines,
        "## Differential expression",
        "",
        "Differential-expression results were not available. Run the Differential Expression module before generating the report.",
        ""
      )
    }
  }
  
  if ("trajectory" %in% sections) {
    lines <- c(
      lines,
      "## Trajectory analysis",
      "",
      if (!is.null(report_data$trajectory_path)) {
        paste0("Trajectory session path: ", report_data$trajectory_path)
      } else {
        "Trajectory results were not available."
      },
      ""
    )
  }
  
  if ("pathway" %in% sections) {
    lines <- c(
      lines,
      "## Pathway activity",
      "",
      if (isTRUE(report_data$pathway_ready)) {
        "Pathway activity results were available."
      } else {
        "Pathway activity results were not available."
      },
      ""
    )
  }
  
  if ("communication" %in% sections) {
    lines <- c(
      lines,
      "## Cell communication",
      "",
      if (isTRUE(report_data$communication_ready)) {
        "Cell communication results were available."
      } else {
        "Cell communication results were not available."
      },
      ""
    )
  }
  
  if ("session" %in% sections) {
    lines <- c(
      lines,
      "## Session information",
      "",
      "```{r}",
      "sessionInfo()",
      "```",
      ""
    )
  }
  
  writeLines(lines, file)
}