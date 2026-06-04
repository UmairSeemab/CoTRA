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
    sc_trajectory = NULL,
    sc_pathway = NULL,
    sc_comm = NULL,
    output_dir = "outputs/scRNA/report"
) {
  moduleServer(id, function(input, output, session) {
    
    report_html_path <- reactiveVal(NULL)
    report_pdf_path  <- reactiveVal(NULL)
    report_dir_path  <- reactiveVal(NULL)
    report_ready     <- reactiveVal(FALSE)
    pdf_error_msg    <- reactiveVal(NULL)
    
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
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
      
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      run_dir <- file.path(output_dir, paste0("CoTRA_scRNA_report_", timestamp))
      fig_dir <- file.path(run_dir, "figures")
      
      dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
      
      html_file <- file.path(run_dir, paste0("CoTRA_scRNA_report_", timestamp, ".html"))
      pdf_file  <- file.path(run_dir, paste0("CoTRA_scRNA_report_", timestamp, ".pdf"))
      rmd_file  <- file.path(run_dir, "CoTRA_scRNA_report.Rmd")
      
      withProgress(message = "Generating scRNA report", value = 0, {
        
        incProgress(0.2, detail = "Saving figures")
        
        figure_files <- create_sc_report_figures(
          obj = obj,
          fig_dir = fig_dir,
          include_figures = isTRUE(input$include_figures)
        )
        
        incProgress(0.45, detail = "Collecting report data")
        
        report_data <- collect_sc_report_data(
          obj = obj,
          sc_markers = sc_markers,
          sc_annot = sc_annot,
          sc_trajectory = sc_trajectory,
          sc_pathway = sc_pathway,
          sc_comm = sc_comm
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
        paste0("CoTRA_scRNA_report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
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
    sc_trajectory = NULL,
    sc_pathway = NULL,
    sc_comm = NULL
) {
  md <- obj@meta.data
  
  marker_table <- tryCatch(sc_markers$markers_table(), error = function(e) NULL)
  annotation_table <- tryCatch(sc_annot$annotation_table(), error = function(e) NULL)
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
    trajectory_path = trajectory_path,
    pathway_ready = pathway_ready,
    communication_ready = communication_ready
  )
}

create_sc_report_figures <- function(obj, fig_dir, include_figures = TRUE) {
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
    
    f <- file.path(fig_dir, "qc_violin.png")
    ggplot2::ggsave(f, p, width = 10, height = 5, dpi = 150)
    figs$qc <- f
  }
  
  if ("pca" %in% names(obj@reductions)) {
    p <- Seurat::DimPlot(obj, reduction = "pca") +
      ggplot2::ggtitle("PCA")
    
    f <- file.path(fig_dir, "pca.png")
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
    
    f <- file.path(fig_dir, "umap.png")
    ggplot2::ggsave(f, p, width = 8, height = 6, dpi = 150)
    figs$umap <- f
  }
  
  if ("tsne" %in% names(obj@reductions)) {
    p <- Seurat::DimPlot(obj, reduction = "tsne") +
      ggplot2::ggtitle("t-SNE")
    
    f <- file.path(fig_dir, "tsne.png")
    ggplot2::ggsave(f, p, width = 8, height = 6, dpi = 150)
    figs$tsne <- f
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
  
  esc <- function(x) {
    x <- as.character(x)
    x <- gsub("\\\\", "/", x)
    x <- gsub('"', "'", x)
    x
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