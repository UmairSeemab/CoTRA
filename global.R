# ======================================================
#                GLOBAL LIBRARIES
# ======================================================
library(shinyFiles)
library(shiny)
library(shinydashboard)
library(bs4Dash)
library(shinyWidgets)
library(shinyjs)
library(plotly)
library(DT)
library(dplyr)
library(Seurat)
library(Matrix)
library(data.table)
library(ggplot2)
library(shinyalert)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(enrichplot)
library(knitr)
library(kableExtra)
library(scales)
library(stringr)
library(shinycssloaders)
library(svglite)
library(htmlwidgets)
library(clusterProfiler)
library(fgsea)
library(msigdbr)
library(ReactomePA)
library(pathview)
library(rmarkdown)
library(SingleCellExperiment)


# ======================================================
#                GLOBAL OPTIONS
# ======================================================

options(shiny.maxRequestSize = 1000 * 1024^2)
options(stringsAsFactors = FALSE)

# ======================================================
#                PATHS
# ======================================================

cotra_path <- getwd()
ref_path   <- file.path(cotra_path, "References")

# ======================================================
#                HELP FILE
# ======================================================

help_file <- file.path(ref_path, "help_infos.csv")

help_infos <- if (file.exists(help_file)) {
  read.csv(help_file, header = TRUE, stringsAsFactors = FALSE)
} else {
  data.frame(
    key   = character(0),
    title = character(0),
    value = character(0),
    stringsAsFactors = FALSE
  )
}

# ======================================================
#         CELL CYCLE REFERENCES (IF AVAILABLE)
# ======================================================

cc_mouse_file <- file.path(ref_path, "mouse_cell_cycle_genes.rds")
cc_human_file <- file.path(ref_path, "human_cell_cycle_genes.rds")

cc.genes.mouse <- if (file.exists(cc_mouse_file)) readRDS(cc_mouse_file) else NULL
cc.genes.human <- if (file.exists(cc_human_file)) readRDS(cc_human_file) else NULL

# ======================================================
#                GENE CONVERSION HELPERS
# ======================================================

convertHumanGeneList <- function(x) {
  hfile <- file.path(ref_path, "human_db_biomart")
  mfile <- file.path(ref_path, "mouse_db_biomart")
  
  if (!file.exists(hfile) || !file.exists(mfile)) return(x)
  
  human <- readRDS(hfile)
  mouse <- readRDS(mfile)
  
  g <- biomaRt::getLDS(
    attributes   = "hgnc_symbol",
    filters      = "hgnc_symbol",
    values       = x,
    mart         = human,
    attributesL  = "mgi_symbol",
    martL        = mouse,
    uniqueRows   = TRUE
  )
  
  unique(g[, 2])
}

convertMouseGeneList <- function(x) {
  hfile <- file.path(ref_path, "human_db_biomart")
  mfile <- file.path(ref_path, "mouse_db_biomart")
  
  if (!file.exists(hfile) || !file.exists(mfile)) return(x)
  
  human <- readRDS(hfile)
  mouse <- readRDS(mfile)
  
  g <- biomaRt::getLDS(
    attributes   = "mgi_symbol",
    filters      = "mgi_symbol",
    values       = x,
    mart         = mouse,
    attributesL  = "hgnc_symbol",
    martL        = human,
    uniqueRows   = TRUE
  )
  
  unique(g[, 2])
}

# ======================================================
#                 GLOBAL SAVE FUNCTION
# ======================================================

save_cotra_project <- function(path, seurat_obj, project_meta = NULL, cells_meta = NULL) {
  saveRDS(
    list(
      seurat       = seurat_obj,
      project_meta = project_meta,
      cells_meta   = cells_meta
    ),
    path
  )
}

# ======================================================
#                  GLOBAL PLOT THEME
# ======================================================

theme_cotra <- function() {
  theme_minimal() +
    theme(
      text             = element_text(size = 14),
      plot.title       = element_text(face = "bold", size = 16),
      panel.grid.minor = element_blank()
    )
}

# ======================================================
#                 OUTPUT DIRECTORY CHECK
# ======================================================

if (!dir.exists("outputs")) dir.create("outputs")

# ======================================================
#               REACTIVE STATE STORE
# ======================================================

cotra_state <- reactiveValues(
  seurat = NULL,
  seurat_list = list(),
  seurat_selected = NULL,
  parameters = list(),
  genes_list = list(),
  output = "outputs"
)

# ======================================================
#        UNIVERSAL DOWNLOAD HANDLER FOR PLOTS
# ======================================================

downPlotUI <- function(id, label = "Download") {
  ns <- NS(id)
  
  tagList(
    br(),
    fluidRow(
      column(3, downloadButton(ns("pdf"),  paste(label, "PDF"),  class = "btn-primary")),
      column(3, downloadButton(ns("svg"),  paste(label, "SVG"),  class = "btn-primary")),
      column(3, downloadButton(ns("png"),  paste(label, "PNG"),  class = "btn-primary")),
      column(3, downloadButton(ns("html"), paste(label, "HTML"), class = "btn-primary"))
    )
  )
}

downPlotServer <- function(id, plot_obj, out_file = "plot") {
  
  moduleServer(id, function(input, output, session) {
    
    get_plot <- reactive({
      p <- plot_obj()
      if ("plotly" %in% class(p)) {
        return(list(type = "plotly", obj = p))
      }
      list(type = "ggplot", obj = p)
    })
    
    # PDF
    output$pdf <- downloadHandler(
      filename = function() paste0(out_file, ".pdf"),
      content  = function(file) {
        p <- get_plot()
        if (p$type == "ggplot") {
          pdf(file, width = 8, height = 6)
          print(p$obj)
          dev.off()
        } else {
          orca::orca(p$obj, file = file)
        }
      }
    )
    
    # SVG
    output$svg <- downloadHandler(
      filename = function() paste0(out_file, ".svg"),
      content  = function(file) {
        p <- get_plot()
        if (p$type == "ggplot") {
          svglite::svglite(file, width = 8, height = 6)
          print(p$obj)
          dev.off()
        } else {
          orca::orca(p$obj, file = file)
        }
      }
    )
    
    # PNG
    output$png <- downloadHandler(
      filename = function() paste0(out_file, ".png"),
      content  = function(file) {
        p <- get_plot()
        if (p$type == "ggplot") {
          png(file, width = 2000, height = 1500, res = 300)
          print(p$obj)
          dev.off()
        } else {
          orca::orca(p$obj, file = file)
        }
      }
    )
    
    # HTML (Interactive)
    output$html <- downloadHandler(
      filename = function() paste0(out_file, ".html"),
      content  = function(file) {
        p <- get_plot()
        if (p$type == "plotly") {
          htmlwidgets::saveWidget(as_widget(p$obj), file)
        } else {
          htmlwidgets::saveWidget(plotly::ggplotly(p$obj), file)
        }
      }
    )
  })
}

# ======================================================
#                 AUTOLOAD MODULES
# ======================================================

mod_files <- list.files(
  "modules",
  pattern     = "\\.(R|r)$",
  full.names  = TRUE,
  recursive   = TRUE
)

for (f in mod_files) {
  try(source(f), silent = TRUE)
}
