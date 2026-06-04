cotra_cran_packages <- function() {
  c(
    "shiny", "shinyFiles", "shinydashboard", "bs4Dash", "shinyWidgets",
    "shinyjs", "plotly", "DT", "dplyr", "Matrix", "data.table",
    "ggplot2", "shinyalert", "knitr", "kableExtra", "scales",
    "stringr", "shinycssloaders", "svglite", "htmlwidgets",
    "msigdbr", "rmarkdown", "gridExtra", "cowplot", "sctransform",
    "corrplot", "igraph", "mgcv", "pagedown", "webshot2", "pheatmap",
    "RColorBrewer", "readr", "tidyr", "openxlsx", "remotes"
  )
}

cotra_bioc_packages <- function() {
  c(
    "Seurat", "SeuratObject", "clusterProfiler", "org.Hs.eg.db",
    "org.Mm.eg.db", "org.Rn.eg.db", "enrichplot", "fgsea", "ReactomePA",
    "pathview", "SingleCellExperiment", "scDblFinder", "MAST",
    "SummarizedExperiment", "slingshot", "UCell", "AUCell", "GSVA",
    "SingleR", "celldex", "limma", "edgeR", "DESeq2", "AnnotationDbi",
    "biomaRt", "GenomeInfoDb", "IRanges", "S4Vectors", "BiocGenerics",
    "ensembldb", "EnsDb.Hsapiens.v86", "EnsDb.Mmusculus.v79",
    "EnsDb.Rnorvegicus.v79"
  )
}

cotra_github_packages <- function() {
  c(
    DoubletFinder = "chris-mcginnis-ucsf/DoubletFinder",
    monocle3 = "cole-trapnell-lab/monocle3",
    CellChat = "jinworks/CellChat"
  )
}

#' Check CoTRA dependencies
#'
#' @param quiet Logical. Suppress messages.
#' @export
check_cotra_dependencies <- function(quiet = FALSE) {
  pkgs <- unique(c(cotra_cran_packages(), cotra_bioc_packages(), names(cotra_github_packages())))
  installed <- vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  missing <- names(installed)[!installed]

  result <- list(
    ok = length(missing) == 0,
    installed = names(installed)[installed],
    missing = missing
  )

  if (!quiet) {
    if (result$ok) {
      message("All CoTRA dependencies are installed.")
    } else {
      message("Missing CoTRA dependencies: ", paste(missing, collapse = ", "))
    }
  }

  result
}

#' Install CoTRA dependencies
#'
#' Installs CRAN, Bioconductor, and selected GitHub dependencies required by the CoTRA app.
#'
#' @param ask Logical. Ask before installing missing packages.
#' @param update Logical. Passed to BiocManager::install().
#' @export
install_cotra_dependencies <- function(ask = interactive(), update = FALSE) {
  options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = max(1000, getOption("timeout")))

  cran <- cotra_cran_packages()
  missing_cran <- cran[!vapply(cran, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

  if (length(missing_cran) > 0) {
    if (ask && !utils::askYesNo(paste("Install missing CRAN packages?", paste(missing_cran, collapse = ", ")))) {
      stop("CRAN dependency installation cancelled.", call. = FALSE)
    }
    utils::install.packages(missing_cran, dependencies = TRUE)
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    utils::install.packages("BiocManager")
  }

  bioc <- cotra_bioc_packages()
  missing_bioc <- bioc[!vapply(bioc, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

  if (length(missing_bioc) > 0) {
    if (ask && !utils::askYesNo(paste("Install missing Bioconductor packages?", paste(missing_bioc, collapse = ", ")))) {
      stop("Bioconductor dependency installation cancelled.", call. = FALSE)
    }
    BiocManager::install(missing_bioc, ask = FALSE, update = update)
  }

  if (!requireNamespace("remotes", quietly = TRUE)) {
    utils::install.packages("remotes")
  }

  github <- cotra_github_packages()
  missing_github <- names(github)[!vapply(names(github), requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

  if (length(missing_github) > 0) {
    if (ask && !utils::askYesNo(paste("Install missing GitHub packages?", paste(missing_github, collapse = ", ")))) {
      stop("GitHub dependency installation cancelled.", call. = FALSE)
    }
    for (pkg in missing_github) {
      remotes::install_github(github[[pkg]], dependencies = TRUE, upgrade = "never")
    }
  }

  check_cotra_dependencies(quiet = FALSE)
  invisible(TRUE)
}
