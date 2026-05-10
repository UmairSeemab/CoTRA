# ==========================================================
# modules/sc/import.R
# CoTRA scRNA-seq Import Module
#
# Production-ready import module for:
# - CellRanger matrix folders/files
# - Seurat .rds objects from Seurat v4 or Seurat v5
# - CoTRA scRNA project .rds files
#
# Features:
# - Automatic CellRanger v2/v3/v4/v5-style file detection
# - Automatic Seurat object version/layer check
# - RNA count layer validation for downstream compatibility
# - Human/mouse/other species detection
# - Basic QC metric calculation
# - Optional HTO and ADT import
# - Optional metadata import
# - Built-in user help section
# ==========================================================

mod_sc_import_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("Single-cell data import and project setup"),
    br(),
    
    fluidRow(
      column(
        width = 12,
        bs4Card(
          title = "Help: supported input formats",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          tags$div(
            tags$h4("What you can upload"),
            tags$ul(
              tags$li(tags$b("CellRanger output:"), " upload the three matrix files: matrix.mtx or matrix.mtx.gz, barcodes.tsv or barcodes.tsv.gz, and features.tsv/genes.tsv or their gzipped versions."),
              tags$li(tags$b("Seurat object:"), " upload a .rds file containing a valid Seurat object."),
              tags$li(tags$b("CoTRA project:"), " upload a .rds file saved by CoTRA. It should contain a Seurat object and optional project/cell metadata.")
            ),
            tags$h4("CellRanger version detection"),
            tags$p("CoTRA checks file names and feature files to infer the likely CellRanger format."),
            tags$ul(
              tags$li(tags$b("CellRanger v2:"), " commonly uses genes.tsv, barcodes.tsv, and matrix.mtx."),
              tags$li(tags$b("CellRanger v3 or newer:"), " commonly uses features.tsv.gz, barcodes.tsv.gz, and matrix.mtx.gz."),
              tags$li(tags$b("CellRanger v4/v5:"), " uses the same matrix exchange structure as v3 for filtered_feature_bc_matrix. CoTRA treats these as suitable if RNA features are found.")
            ),
            tags$h4("Seurat v4/v5 compatibility"),
            tags$ul(
              tags$li("CoTRA checks whether the uploaded object is a valid Seurat object."),
              tags$li("For Seurat v4, it checks the RNA assay counts slot."),
              tags$li("For Seurat v5, it checks RNA assay layers and looks for a counts layer."),
              tags$li("If a counts layer is missing, the object may still load, but downstream QC and normalization may be unreliable.")
            ),
            tags$h4("Species and QC detection"),
            tags$ul(
              tags$li("Human is detected mainly from MT- and RP genes."),
              tags$li("Mouse is detected mainly from mt- and Rp genes."),
              tags$li("If species cannot be detected, choose Human, Mouse, or Other from Species override before normalization.")
            ),
            tags$h4("Recommended input"),
            tags$p("For a new analysis, upload raw CellRanger filtered_feature_bc_matrix files. For continuing an old analysis, upload a Seurat .rds or CoTRA project .rds.")
          )
        )
      )
    ),
    
    fluidRow(
      column(
        width = 4,
        bs4Card(
          title = "Input type",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          radioButtons(
            ns("type_import"),
            label = "Select input type",
            choices = c(
              "CellRanger output" = "cellranger",
              "Seurat R object (.rds)" = "rds",
              "CoTRA project (.rds)" = "cotra"
            ),
            selected = "cellranger"
          )
        )
      ),
      column(
        width = 8,
        bs4Card(
          title = "Import status",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          htmlOutput(ns("import_status"))
        )
      )
    ),
    
    conditionalPanel(
      condition = "input.type_import == 'cellranger'",
      ns = ns,
      fluidRow(
        column(
          width = 6,
          bs4Card(
            title = "CellRanger RNA files",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            helpText("Upload exactly three files from filtered_feature_bc_matrix or raw_feature_bc_matrix."),
            fileInput(
              ns("crfiles"),
              "Upload RNA matrix files",
              multiple = TRUE,
              accept = c(".mtx", ".gz", ".tsv")
            )
          )
        ),
        column(
          width = 6,
          bs4Card(
            title = "HTO / ADT files, optional",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            helpText("Upload separate CellRanger-style files only if HTO or ADT were generated separately."),
            fileInput(
              ns("htofiles"),
              "HTO files, optional",
              multiple = TRUE,
              accept = c(".mtx", ".gz", ".tsv")
            ),
            fileInput(
              ns("adtfiles"),
              "ADT files, optional",
              multiple = TRUE,
              accept = c(".mtx", ".gz", ".tsv")
            )
          )
        )
      )
    ),
    
    conditionalPanel(
      condition = "input.type_import == 'rds'",
      ns = ns,
      fluidRow(
        column(
          width = 6,
          bs4Card(
            title = "Seurat R object",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            fileInput(
              ns("rdsfile"),
              "Upload Seurat .rds file",
              multiple = FALSE,
              accept = ".rds"
            )
          )
        )
      )
    ),
    
    conditionalPanel(
      condition = "input.type_import == 'cotra'",
      ns = ns,
      fluidRow(
        column(
          width = 6,
          bs4Card(
            title = "CoTRA project",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            fileInput(
              ns("cotrafile"),
              "Upload CoTRA project .rds",
              multiple = FALSE,
              accept = ".rds"
            )
          )
        )
      )
    ),
    
    br(),
    fluidRow(column(width = 12, uiOutput(ns("qc_panel")))),
    br(),
    fluidRow(column(width = 12, uiOutput(ns("project_panel"))))
  )
}


mod_sc_import_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rv <- reactiveValues(
      seu = NULL,
      input_check = NULL,
      cellranger_version = NULL,
      seurat_version = NULL,
      seurat_suitable = FALSE,
      hto = FALSE,
      adt = FALSE,
      project_meta = NULL,
      cells_meta = NULL,
      imported = FALSE,
      normalized = FALSE,
      messages = character(0)
    )
    
    # -------------------- helpers --------------------
    
    add_message <- function(...) {
      msg <- paste0(...)
      rv$messages <- unique(c(rv$messages, msg))
    }
    
    safe_theme <- function() {
      if (exists("theme_cotra", mode = "function")) {
        theme_cotra()
      } else {
        ggplot2::theme_bw()
      }
    }
    
    cotra_interpretation_box <- function(title, ...) {
      bs4Card(
        title = paste("Interpretation:", title),
        status = "info",
        solidHeader = TRUE,
        width = 12,
        collapsible = TRUE,
        collapsed = TRUE,
        icon = icon("circle-info"),
        ...
      )
    }
    
    show_import_error <- function(title, message) {
      showModal(modalDialog(title = title, message, easyClose = TRUE))
    }
    
    normalize_uploaded_file_names <- function(files_df) {
      if (is.null(files_df) || nrow(files_df) == 0) return(NULL)
      
      tmp_dir <- tempfile("cotra_upload_")
      dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
      
      out <- files_df
      for (i in seq_len(nrow(files_df))) {
        clean_name <- basename(files_df$name[i])
        clean_name <- gsub("\\\\", "/", clean_name)
        clean_name <- basename(clean_name)
        target <- file.path(tmp_dir, clean_name)
        ok <- file.copy(files_df$datapath[i], target, overwrite = TRUE)
        if (!isTRUE(ok)) stop("Could not copy uploaded file: ", files_df$name[i])
        out$datapath[i] <- target
        out$name[i] <- clean_name
      }
      out$upload_dir <- tmp_dir
      out
    }
    
    detect_cellranger_files <- function(files_df) {
      if (is.null(files_df) || nrow(files_df) == 0) {
        return(list(ok = FALSE, message = "No files were uploaded."))
      }
      
      fn <- tolower(basename(files_df$name))
      has_matrix <- any(fn %in% c("matrix.mtx", "matrix.mtx.gz"))
      has_barcodes <- any(fn %in% c("barcodes.tsv", "barcodes.tsv.gz"))
      has_features <- any(fn %in% c("features.tsv", "features.tsv.gz"))
      has_genes <- any(fn %in% c("genes.tsv", "genes.tsv.gz"))
      
      if (!has_matrix || !has_barcodes || (!has_features && !has_genes)) {
        return(list(
          ok = FALSE,
          message = paste(
            "Upload must contain matrix.mtx or matrix.mtx.gz, barcodes.tsv or barcodes.tsv.gz, and features.tsv/genes.tsv or gzipped equivalent."
          )
        ))
      }
      
      version <- "CellRanger v3 or newer"
      feature_file <- fn[fn %in% c("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz")][1]
      
      if (grepl("^genes.tsv", feature_file)) {
        version <- "CellRanger v2-like"
      } else if (grepl("^features.tsv", feature_file)) {
        version <- "CellRanger v3/v4/v5-like"
      }
      
      list(ok = TRUE, version = version, message = paste("Detected", version))
    }
    
    read_cellranger_triple <- function(files_df) {
      files_df <- normalize_uploaded_file_names(files_df)
      check <- detect_cellranger_files(files_df)
      if (!isTRUE(check$ok)) return(list(ok = FALSE, message = check$message))
      
      data_dir <- unique(dirname(files_df$datapath))[1]
      mat <- tryCatch(
        Seurat::Read10X(data.dir = data_dir),
        error = function(e) e
      )
      
      if (inherits(mat, "error")) {
        return(list(ok = FALSE, message = paste("Read10X failed:", conditionMessage(mat))))
      }
      
      list(ok = TRUE, mat = mat, version = check$version, message = check$message)
    }
    
    choose_rna_matrix <- function(mat) {
      if (!is.list(mat)) return(mat)
      
      possible_names <- c(
        "Gene Expression",
        "gene expression",
        "RNA",
        "RNA Gene Expression"
      )
      
      hit <- intersect(possible_names, names(mat))
      if (length(hit) > 0) return(mat[[hit[1]]])
      
      mat[[1]]
    }
    
    detect_species <- function(seu) {
      feats <- rownames(seu)
      if (length(feats) == 0) return("other")
      
      human_mito <- sum(grepl("^MT-", feats), na.rm = TRUE)
      mouse_mito <- sum(grepl("^mt-", feats), na.rm = TRUE)
      human_ribo <- sum(grepl("^RP[LS]", feats), na.rm = TRUE)
      mouse_ribo <- sum(grepl("^Rp[ls]", feats), na.rm = TRUE)
      
      human_score <- human_mito + human_ribo
      mouse_score <- mouse_mito + mouse_ribo
      
      if (human_score > mouse_score && human_score > 0) return("human")
      if (mouse_score > human_score && mouse_score > 0) return("mouse")
      "other"
    }
    
    set_species_and_qc <- function(seu, species_override = NULL) {
      sp <- species_override
      if (is.null(sp) || sp == "auto" || !nzchar(sp)) sp <- detect_species(seu)
      
      seu@meta.data$species <- sp
      
      if (sp == "mouse") {
        seu@meta.data$percent.mito <- Seurat::PercentageFeatureSet(seu, assay = "RNA", pattern = "^mt-")
        seu@meta.data$percent.ribo <- Seurat::PercentageFeatureSet(seu, assay = "RNA", pattern = "^Rp[ls]")
      } else if (sp == "human") {
        seu@meta.data$percent.mito <- Seurat::PercentageFeatureSet(seu, assay = "RNA", pattern = "^MT-")
        seu@meta.data$percent.ribo <- Seurat::PercentageFeatureSet(seu, assay = "RNA", pattern = "^RP[LS]")
      } else {
        seu@meta.data$percent.mito <- NA_real_
        seu@meta.data$percent.ribo <- NA_real_
      }
      
      seu
    }
    
    get_seurat_major_version <- function() {
      tryCatch(
        as.integer(strsplit(as.character(utils::packageVersion("SeuratObject")), "\\.")[[1]][1]),
        error = function(e) NA_integer_
      )
    }
    
    assay_has_counts <- function(seu, assay = "RNA") {
      if (!inherits(seu, "Seurat")) return(FALSE)
      if (!assay %in% Seurat::Assays(seu)) return(FALSE)
      
      ok <- tryCatch({
        mat <- Seurat::GetAssayData(seu, assay = assay, slot = "counts")
        !is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0
      }, error = function(e) FALSE)
      
      if (isTRUE(ok)) return(TRUE)
      
      ok_layer <- tryCatch({
        assay_obj <- seu[[assay]]
        layer_names <- SeuratObject::Layers(assay_obj)
        any(layer_names %in% c("counts", "counts.1"))
      }, error = function(e) FALSE)
      
      isTRUE(ok_layer)
    }
    
    get_assay_layers_text <- function(seu, assay = "RNA") {
      out <- tryCatch({
        paste(SeuratObject::Layers(seu[[assay]]), collapse = ", ")
      }, error = function(e) "not available")
      if (!nzchar(out)) out <- "not available"
      out
    }
    
    inspect_seurat_object <- function(obj) {
      if (!inherits(obj, "Seurat")) {
        return(list(ok = FALSE, suitable = FALSE, message = "Uploaded RDS is not a Seurat object."))
      }
      
      assays <- Seurat::Assays(obj)
      if (!"RNA" %in% assays) {
        return(list(
          ok = FALSE,
          suitable = FALSE,
          message = paste0("No RNA assay found. Available assays: ", paste(assays, collapse = ", "))
        ))
      }
      
      major <- get_seurat_major_version()
      layers <- get_assay_layers_text(obj, "RNA")
      has_counts <- assay_has_counts(obj, "RNA")
      
      version_text <- if (!is.na(major) && major >= 5) {
        "SeuratObject v5 environment"
      } else if (!is.na(major)) {
        "SeuratObject v4 environment"
      } else {
        "SeuratObject version unknown"
      }
      
      suitable <- has_counts && ncol(obj) > 0 && nrow(obj) > 0
      
      message <- paste0(
        version_text,
        "; assays: ", paste(assays, collapse = ", "),
        "; RNA layers/slots: ", layers,
        "; counts available: ", ifelse(has_counts, "yes", "no"),
        "; cells: ", ncol(obj),
        "; features: ", nrow(obj),
        "."
      )
      
      list(
        ok = TRUE,
        suitable = suitable,
        version = version_text,
        layers = layers,
        has_counts = has_counts,
        message = message
      )
    }
    
    safe_read_rds <- function(path) {
      tryCatch(readRDS(path), error = function(e) e)
    }
    
    qc_metrics <- reactive({
      req(rv$seu)
      seu <- rv$seu
      md <- seu@meta.data
      
      data.frame(
        Metric = c(
          "Cells",
          "Features",
          "Assays",
          "Detected species",
          "CellRanger format",
          "Seurat check",
          "Suitable for downstream",
          "Median nCount_RNA",
          "Median nFeature_RNA",
          "Median percent.mito",
          "Median percent.ribo"
        ),
        Value = c(
          ncol(seu),
          nrow(seu),
          paste(Seurat::Assays(seu), collapse = ", "),
          if ("species" %in% colnames(md)) as.character(unique(md$species))[1] else "NA",
          if (!is.null(rv$cellranger_version)) rv$cellranger_version else "NA",
          if (!is.null(rv$seurat_version)) rv$seurat_version else "Created in current session",
          ifelse(isTRUE(rv$seurat_suitable), "Yes", "Check warning"),
          if ("nCount_RNA" %in% colnames(md)) stats::median(md$nCount_RNA, na.rm = TRUE) else NA,
          if ("nFeature_RNA" %in% colnames(md)) stats::median(md$nFeature_RNA, na.rm = TRUE) else NA,
          if ("percent.mito" %in% colnames(md)) stats::median(md$percent.mito, na.rm = TRUE) else NA,
          if ("percent.ribo" %in% colnames(md)) stats::median(md$percent.ribo, na.rm = TRUE) else NA
        ),
        stringsAsFactors = FALSE
      )
    })
    
    # -------------------- status --------------------
    
    output$import_status <- renderUI({
      if (!isTRUE(rv$imported) || is.null(rv$seu)) {
        return(HTML("<span style='color:#777;'>No object imported yet.</span>"))
      }
      
      status <- paste0(
        "<b>Imported:</b> ", ncol(rv$seu), " cells and ", nrow(rv$seu), " features.<br>",
        "<b>Assays:</b> ", paste(Seurat::Assays(rv$seu), collapse = ", "), "<br>",
        "<b>Species:</b> ", unique(rv$seu@meta.data$species)[1], "<br>"
      )
      
      if (!is.null(rv$cellranger_version)) {
        status <- paste0(status, "<b>CellRanger check:</b> ", rv$cellranger_version, "<br>")
      }
      
      if (!is.null(rv$seurat_version)) {
        status <- paste0(status, "<b>Seurat check:</b> ", rv$seurat_version, "<br>")
      }
      
      status <- paste0(
        status,
        "<b>Downstream suitability:</b> ", ifelse(isTRUE(rv$seurat_suitable), "Suitable", "Check warnings"), "<br>",
        "<b>Normalized:</b> ", ifelse(isTRUE(rv$normalized), "Yes", "No"), "<br>"
      )
      
      if (length(rv$messages) > 0) {
        status <- paste0(status, "<br><b>Messages:</b><ul>", paste0("<li>", rv$messages, "</li>", collapse = ""), "</ul>")
      }
      
      HTML(status)
    })
    
    # -------------------- CellRanger import --------------------
    
    observeEvent(input$crfiles, {
      req(input$crfiles)
      
      withProgress(message = "Importing CellRanger RNA files", value = 0, {
        incProgress(0.2)
        
        res <- read_cellranger_triple(input$crfiles)
        if (!isTRUE(res$ok)) {
          show_import_error("CellRanger file check failed", res$message)
          rv$seu <- NULL
          rv$imported <- FALSE
          return(NULL)
        }
        
        incProgress(0.4, detail = res$message)
        
        mat_rna <- choose_rna_matrix(res$mat)
        if (is.null(mat_rna) || nrow(mat_rna) == 0 || ncol(mat_rna) == 0) {
          show_import_error("RNA matrix error", "Could not identify a non-empty RNA expression matrix.")
          rv$seu <- NULL
          rv$imported <- FALSE
          return(NULL)
        }
        
        seu <- Seurat::CreateSeuratObject(
          counts = mat_rna,
          assay = "RNA",
          project = "CoTRA_scRNA"
        )
        
        incProgress(0.7, detail = "Calculating QC metrics")
        
        seu <- set_species_and_qc(seu, species_override = "auto")
        
        rv$seu <- seu
        rv$input_check <- res$message
        rv$cellranger_version <- res$version
        rv$seurat_version <- "Created from CellRanger matrix"
        rv$seurat_suitable <- TRUE
        rv$imported <- TRUE
        rv$normalized <- FALSE
        rv$hto <- FALSE
        rv$adt <- FALSE
        rv$messages <- character(0)
        
        if (detect_species(seu) == "other") {
          add_message("Species could not be confidently detected. Use Species override before normalization.")
        }
        
        incProgress(1)
      })
    })
    
    # -------------------- HTO / ADT import --------------------
    
    add_modality_assay <- function(files_df, assay_name, normalize_clr = FALSE) {
      req(rv$seu, files_df)
      
      res <- read_cellranger_triple(files_df)
      if (!isTRUE(res$ok)) {
        show_import_error(paste(assay_name, "file check failed"), res$message)
        return(NULL)
      }
      
      mat <- choose_rna_matrix(res$mat)
      if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
        show_import_error(paste(assay_name, "matrix error"), "Could not identify a non-empty matrix.")
        return(NULL)
      }
      
      common <- intersect(colnames(mat), colnames(rv$seu))
      if (length(common) < 10) {
        show_import_error(
          paste(assay_name, "barcode mismatch"),
          paste("Only", length(common), "overlapping barcodes were found between RNA and", assay_name, ".")
        )
        return(NULL)
      }
      
      seu <- subset(rv$seu, cells = common)
      mat <- mat[, common, drop = FALSE]
      seu[[assay_name]] <- Seurat::CreateAssayObject(counts = mat)
      
      if (isTRUE(normalize_clr)) {
        seu <- Seurat::NormalizeData(seu, assay = assay_name, normalization.method = "CLR", margin = 2)
      }
      
      rv$seu <- seu
      TRUE
    }
    
    observeEvent(input$htofiles, {
      ok <- add_modality_assay(input$htofiles, "HTO", normalize_clr = FALSE)
      if (isTRUE(ok)) {
        rv$hto <- TRUE
        add_message("HTO assay added. RNA object was subset to barcodes shared with HTO.")
      }
    })
    
    observeEvent(input$adtfiles, {
      ok <- add_modality_assay(input$adtfiles, "ADT", normalize_clr = TRUE)
      if (isTRUE(ok)) {
        rv$adt <- TRUE
        add_message("ADT assay added and CLR-normalized. RNA object was subset to barcodes shared with ADT.")
      }
    })
    
    # -------------------- Seurat RDS import --------------------
    
    observeEvent(input$rdsfile, {
      req(input$rdsfile)
      
      withProgress(message = "Loading Seurat object", value = 0, {
        incProgress(0.3)
        obj <- safe_read_rds(input$rdsfile$datapath)
        
        if (inherits(obj, "error")) {
          show_import_error("RDS read error", conditionMessage(obj))
          rv$seu <- NULL
          rv$imported <- FALSE
          return(NULL)
        }
        
        check <- inspect_seurat_object(obj)
        if (!isTRUE(check$ok)) {
          show_import_error("Invalid Seurat object", check$message)
          rv$seu <- NULL
          rv$imported <- FALSE
          return(NULL)
        }
        
        incProgress(0.6, detail = "Checking Seurat assay structure")
        
        obj <- set_species_and_qc(obj, species_override = "auto")
        
        rv$seu <- obj
        rv$input_check <- check$message
        rv$cellranger_version <- NULL
        rv$seurat_version <- check$version
        rv$seurat_suitable <- isTRUE(check$suitable)
        rv$imported <- TRUE
        rv$normalized <- FALSE
        rv$hto <- "HTO" %in% Seurat::Assays(obj)
        rv$adt <- "ADT" %in% Seurat::Assays(obj)
        rv$messages <- character(0)
        
        add_message(check$message)
        if (!isTRUE(check$suitable)) {
          add_message("Warning: RNA counts were not clearly detected. Downstream QC or normalization may fail.")
        }
        if (detect_species(obj) == "other") {
          add_message("Species could not be confidently detected. Use Species override before normalization.")
        }
        
        incProgress(1)
      })
    })
    
    # -------------------- CoTRA project import --------------------
    
    observeEvent(input$cotrafile, {
      req(input$cotrafile)
      
      withProgress(message = "Loading CoTRA project", value = 0, {
        incProgress(0.3)
        obj <- safe_read_rds(input$cotrafile$datapath)
        
        if (inherits(obj, "error")) {
          show_import_error("RDS read error", conditionMessage(obj))
          rv$seu <- NULL
          rv$imported <- FALSE
          return(NULL)
        }
        
        if (inherits(obj, "Seurat")) {
          seu <- obj
          rv$project_meta <- NULL
          rv$cells_meta <- NULL
        } else if (is.list(obj) && "seurat" %in% names(obj) && inherits(obj$seurat, "Seurat")) {
          seu <- obj$seurat
          rv$project_meta <- obj$project_meta
          rv$cells_meta <- obj$cells_meta
        } else {
          show_import_error("Invalid CoTRA project", "The uploaded file is not a CoTRA project and does not contain a valid Seurat object.")
          rv$seu <- NULL
          rv$imported <- FALSE
          return(NULL)
        }
        
        check <- inspect_seurat_object(seu)
        if (!isTRUE(check$ok)) {
          show_import_error("Invalid Seurat object inside project", check$message)
          rv$seu <- NULL
          rv$imported <- FALSE
          return(NULL)
        }
        
        incProgress(0.6, detail = "Checking project object")
        
        seu <- set_species_and_qc(seu, species_override = "auto")
        
        rv$seu <- seu
        rv$input_check <- check$message
        rv$cellranger_version <- NULL
        rv$seurat_version <- check$version
        rv$seurat_suitable <- isTRUE(check$suitable)
        rv$imported <- TRUE
        rv$normalized <- FALSE
        rv$hto <- "HTO" %in% Seurat::Assays(seu)
        rv$adt <- "ADT" %in% Seurat::Assays(seu)
        rv$messages <- character(0)
        
        add_message(check$message)
        if (!isTRUE(check$suitable)) {
          add_message("Warning: RNA counts were not clearly detected. Downstream QC or normalization may fail.")
        }
        if (detect_species(seu) == "other") {
          add_message("Species could not be confidently detected. Use Species override before normalization.")
        }
        
        incProgress(1)
      })
    })
    
    # -------------------- QC panel --------------------
    
    output$qc_panel <- renderUI({
      if (!isTRUE(rv$imported) || is.null(rv$seu)) return(NULL)
      
      tagList(
        bs4Card(
          title = "QC summary after import",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          fluidRow(
            column(width = 4, tableOutput(ns("qc_table"))),
            column(width = 4, shinycssloaders::withSpinner(plotOutput(ns("qc_vln_counts"), height = "250px"))),
            column(width = 4, shinycssloaders::withSpinner(plotOutput(ns("qc_vln_features"), height = "250px")))
          ),
          br(),
          fluidRow(
            column(width = 4, shinycssloaders::withSpinner(plotOutput(ns("qc_vln_mito"), height = "250px"))),
            column(width = 4, shinycssloaders::withSpinner(plotOutput(ns("qc_vln_ribo"), height = "250px")))
          )
        ),
        cotra_interpretation_box(
          "QC summary and import checks",
          tags$h5("How to read this section"),
          tags$p("The QC summary confirms whether the uploaded data were imported correctly and whether the object is suitable for downstream scRNA-seq analysis."),
          tags$ul(
            tags$li(tags$b("Cells:"), " total number of barcodes or cells in the imported object."),
            tags$li(tags$b("Features:"), " total number of genes or measured features."),
            tags$li(tags$b("Assays:"), " available data layers such as RNA, ADT, or HTO."),
            tags$li(tags$b("Detected species:"), " inferred from mitochondrial and ribosomal gene naming."),
            tags$li(tags$b("Suitable for downstream:"), " shows whether RNA counts were detected for QC, normalization, and clustering.")
          ),
          tags$h5("How to interpret the QC plots"),
          tags$ul(
            tags$li(tags$b("nCount_RNA:"), " total transcript counts per cell. Very high values may indicate doublets. Very low values may indicate empty droplets or low-quality cells."),
            tags$li(tags$b("nFeature_RNA:"), " number of detected genes per cell. Low values often indicate poor-quality cells. Very high values may indicate doublets."),
            tags$li(tags$b("percent.mito:"), " fraction of mitochondrial reads. High values may indicate stressed or damaged cells."),
            tags$li(tags$b("percent.ribo:"), " fraction of ribosomal reads. Strong shifts may reflect technical or biological differences.")
          ),
          tags$h5("Recommended action"),
          tags$p("Use these plots to choose filtering thresholds in the QC module. Do not apply fixed thresholds blindly. Compare the distribution of your own dataset, sample type, and expected cell population."),
          tags$h5("Caution"),
          tags$p("This section gives technical interpretation only. It does not assign cell types or explain biology. Biological interpretation should use clustering, marker genes, metadata, and enrichment results.")
        )
      )
    })
    
    output$qc_table <- renderTable({
      qc_metrics()
    }, digits = 3)
    
    output$qc_vln_counts <- renderPlot({
      req(rv$seu)
      validate(need("nCount_RNA" %in% colnames(rv$seu@meta.data), "nCount_RNA is missing."))
      Seurat::VlnPlot(rv$seu, features = "nCount_RNA", pt.size = 0.1) + safe_theme()
    })
    
    output$qc_vln_features <- renderPlot({
      req(rv$seu)
      validate(need("nFeature_RNA" %in% colnames(rv$seu@meta.data), "nFeature_RNA is missing."))
      Seurat::VlnPlot(rv$seu, features = "nFeature_RNA", pt.size = 0.1) + safe_theme()
    })
    
    output$qc_vln_mito <- renderPlot({
      req(rv$seu)
      validate(need("percent.mito" %in% colnames(rv$seu@meta.data), "percent.mito is missing."))
      Seurat::VlnPlot(rv$seu, features = "percent.mito", pt.size = 0.1) + safe_theme()
    })
    
    output$qc_vln_ribo <- renderPlot({
      req(rv$seu)
      validate(need("percent.ribo" %in% colnames(rv$seu@meta.data), "percent.ribo is missing."))
      Seurat::VlnPlot(rv$seu, features = "percent.ribo", pt.size = 0.1) + safe_theme()
    })
    
    # -------------------- project setup panel --------------------
    
    output$project_panel <- renderUI({
      if (!isTRUE(rv$imported) || is.null(rv$seu)) return(NULL)
      
      tagList(
        bs4Card(
          title = "Project setup and export",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          fluidRow(
            column(
              width = 4,
              textInput(
                ns("proj_name"),
                "Project name",
                value = if (!is.null(rv$seu@project.name) && nzchar(rv$seu@project.name)) rv$seu@project.name else "CoTRA_scRNA_project"
              ),
              selectInput(
                ns("species_override"),
                "Species override",
                choices = c("Auto detect" = "auto", "Human" = "human", "Mouse" = "mouse", "Other" = "other"),
                selected = "auto"
              ),
              selectInput(
                ns("norm_method"),
                "Normalization method",
                choices = c("LogNormalize", "SCTransform"),
                selected = "LogNormalize"
              ),
              actionButton(ns("run_norm"), "Apply normalization", class = "btn-primary")
            ),
            column(
              width = 4,
              h4("Metadata, optional"),
              helpText("Cell metadata CSV. First column must contain cell barcodes matching the Seurat object."),
              fileInput(ns("cell_meta_file"), "Cell metadata", accept = c(".csv", ".txt")),
              br(),
              helpText("Project-level metadata CSV."),
              fileInput(ns("proj_meta_file"), "Project metadata", accept = c(".csv", ".txt"))
            ),
            column(
              width = 4,
              h4("Export"),
              helpText("Save a reusable CoTRA scRNA project."),
              br(),
              downloadButton(ns("save_project"), "Download CoTRA project (.rds)")
            )
          )
        ),
        cotra_interpretation_box(
          "project setup and normalization",
          tags$h5("Species override"),
          tags$p("Use Auto detect first. Change this only if mitochondrial or ribosomal gene names were not detected correctly."),
          tags$ul(
            tags$li(tags$b("Human:"), " uses MT- and RP gene patterns for QC metrics."),
            tags$li(tags$b("Mouse:"), " uses mt- and Rp gene patterns for QC metrics."),
            tags$li(tags$b("Other:"), " keeps mitochondrial and ribosomal QC metrics unavailable unless custom species support is added.")
          ),
          tags$h5("Normalization choice"),
          tags$ul(
            tags$li(tags$b("LogNormalize:"), " a standard, fast option suitable for most initial scRNA-seq workflows."),
            tags$li(tags$b("SCTransform:"), " useful when stronger technical normalization is needed, especially for datasets with variable sequencing depth."),
            tags$li("Use the same normalization method consistently before PCA, UMAP, clustering, marker detection, and integration steps.")
          ),
          tags$h5("Metadata upload"),
          tags$p("Cell metadata should contain cell barcodes in the first column. Only matching barcodes will be added to the Seurat object."),
          tags$h5("Recommended action"),
          tags$p("After import, check QC plots first. Then set species correctly, add metadata if available, normalize the data, and continue to variable gene selection and PCA.")
        )
      )
    })
    
    # -------------------- normalization --------------------
    
    observeEvent(input$run_norm, {
      req(rv$seu)
      
      withProgress(message = "Running normalization", value = 0, {
        incProgress(0.2)
        
        seu <- rv$seu
        seu <- set_species_and_qc(seu, species_override = input$species_override)
        
        incProgress(0.5, detail = input$norm_method)
        
        if (input$norm_method == "LogNormalize") {
          seu <- Seurat::NormalizeData(
            seu,
            assay = "RNA",
            normalization.method = "LogNormalize",
            scale.factor = 10000,
            verbose = FALSE
          )
        } else {
          if (!requireNamespace("sctransform", quietly = TRUE)) {
            show_import_error("SCTransform unavailable", "Install the sctransform package or choose LogNormalize.")
            return(NULL)
          }
          seu <- Seurat::SCTransform(
            seu,
            assay = "RNA",
            verbose = FALSE,
            return.only.var.genes = FALSE
          )
        }
        
        seu@meta.data$Project <- input$proj_name
        rv$seu <- seu
        rv$normalized <- TRUE
        
        incProgress(1)
      })
    })
    
    # -------------------- metadata import --------------------
    
    observeEvent(input$cell_meta_file, {
      req(rv$seu, input$cell_meta_file)
      
      df <- tryCatch(
        read.csv(input$cell_meta_file$datapath, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE),
        error = function(e) e
      )
      
      if (inherits(df, "error")) {
        show_import_error("Cell metadata error", conditionMessage(df))
        return(NULL)
      }
      
      common <- intersect(rownames(df), colnames(rv$seu))
      if (length(common) == 0) {
        show_import_error("Cell metadata mismatch", "No overlapping cell IDs were found between metadata and Seurat object.")
        return(NULL)
      }
      
      seu <- rv$seu
      for (col in colnames(df)) {
        cname <- col
        if (cname %in% colnames(seu@meta.data)) cname <- paste0(cname, "_meta")
        seu@meta.data[[cname]] <- NA
        seu@meta.data[common, cname] <- df[common, col]
      }
      
      rv$seu <- seu
      rv$cells_meta <- df
      add_message(paste("Cell metadata added for", length(common), "cells."))
    })
    
    observeEvent(input$proj_meta_file, {
      req(input$proj_meta_file)
      
      df <- tryCatch(
        read.csv(input$proj_meta_file$datapath, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
        error = function(e) e
      )
      
      if (inherits(df, "error")) {
        show_import_error("Project metadata error", conditionMessage(df))
        return(NULL)
      }
      
      rv$project_meta <- df
      add_message("Project metadata added.")
    })
    
    # -------------------- export CoTRA project --------------------
    
    output$save_project <- downloadHandler(
      filename = function() {
        nm <- if (!is.null(input$proj_name) && nzchar(input$proj_name)) input$proj_name else "CoTRA_scRNA_project"
        nm <- gsub("[^A-Za-z0-9_\\-]", "_", nm)
        paste0(nm, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_CoTRA_scRNA.rds")
      },
      content = function(file) {
        req(rv$seu)
        
        if (exists("save_cotra_project", mode = "function")) {
          save_cotra_project(
            path = file,
            seurat_obj = rv$seu,
            project_meta = rv$project_meta,
            cells_meta = rv$cells_meta
          )
        } else {
          saveRDS(
            list(
              seurat = rv$seu,
              project_meta = rv$project_meta,
              cells_meta = rv$cells_meta,
              import_check = rv$input_check,
              cellranger_version = rv$cellranger_version,
              seurat_version = rv$seurat_version,
              seurat_suitable = rv$seurat_suitable,
              saved_at = Sys.time()
            ),
            file
          )
        }
      }
    )
    
    # -------------------- return for downstream modules --------------------
    
    return(
      list(
        seurat = reactive(rv$seu),
        project_meta = reactive(rv$project_meta),
        cells_meta = reactive(rv$cells_meta),
        imported = reactive(rv$imported),
        normalized = reactive(rv$normalized),
        input_check = reactive(rv$input_check),
        seurat_suitable = reactive(rv$seurat_suitable)
      )
    )
  })
}
