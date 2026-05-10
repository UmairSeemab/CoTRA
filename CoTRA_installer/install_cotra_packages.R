# ==========================================================
# CoTRA package installation script
# Skips packages already installed
# Minimizes user input
# ==========================================================

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  ask = FALSE,
  timeout = 1000
)

cran_pkgs <- c(
  "pagedown", "webshot2",
  "BiocManager",
  "remotes", "devtools",
  "svglite", "plotly", "dplyr",
  "msigdbr",
  "NMF",
  "Rcpp", "RcppArmadillo", "Matrix", "igraph", "h5",
  "shiny", "bs4Dash", "shinyjs", "shinyWidgets",
  "shinyalert", "shinyFiles", "shinydashboard",
  "DT", "data.table", "ggplot2",
  "knitr", "kableExtra", "scales", "stringr",
  "shinycssloaders", "htmlwidgets",
  "rmarkdown", "lme4"
)

bioc_pkgs <- c(
  "DirichletMultinomial",
  "TFBSTools",
  "ReactomePA",
  "enrichplot",
  "AnnotationDbi",
  "org.Mm.eg.db", "org.Hs.eg.db", "org.Rn.eg.db",
  "ensembldb",
  "EnsDb.Mmusculus.v79",
  "EnsDb.Hsapiens.v86",
  "EnsDb.Rnorvegicus.v79",
  "BiocGenerics",
  "DelayedArray",
  "DelayedMatrixStats",
  "limma",
  "S4Vectors",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "batchelor",
  "HDF5Array",
  "ggrastr",
  "slingshot",
  "UCell",
  "AUCell",
  "GSVA",
  "VISION",
  "pcaMethods",
  "clusterProfiler",
  "fgsea",
  "pathview"
)

github_pkgs <- list(
  azimuth = "satijalab/azimuth",
  BPCells = "bnprks/BPCells/r",
  monocle3 = "cole-trapnell-lab/monocle3",
  circlize = "jokergoo/circlize",
  ComplexHeatmap = "jokergoo/ComplexHeatmap",
  CellChat = "jinworks/CellChat",
  nichenetr = "saeyslab/nichenetr",
  liana = "saezlab/liana"
)

is_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

install_cran_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, is_installed, logical(1))]
  if (length(missing) > 0) {
    install.packages(missing, dependencies = TRUE)
  } else {
    message("All CRAN packages already installed.")
  }
}

install_bioc_if_missing <- function(pkgs) {
  if (!is_installed("BiocManager")) {
    install.packages("BiocManager")
  }
  
  missing <- pkgs[!vapply(pkgs, is_installed, logical(1))]
  if (length(missing) > 0) {
    BiocManager::install(
      missing,
      ask = FALSE,
      update = FALSE,
      dependencies = TRUE
    )
  } else {
    message("All Bioconductor packages already installed.")
  }
}

install_github_if_missing <- function(pkg, repo) {
  if (!is_installed(pkg)) {
    remotes::install_github(
      repo,
      dependencies = TRUE,
      upgrade = "never",
      build_vignettes = FALSE
    )
  } else {
    message(pkg, " already installed.")
  }
}

message("Installing missing CRAN packages...")
install_cran_if_missing(cran_pkgs)

message("Installing missing Bioconductor packages...")
install_bioc_if_missing(bioc_pkgs)

message("Installing missing GitHub packages...")
for (pkg in names(github_pkgs)) {
  install_github_if_missing(pkg, github_pkgs[[pkg]])
}

message("Checking important packages...")

check_pkgs <- unique(c(
  cran_pkgs,
  bioc_pkgs,
  names(github_pkgs)
))

status <- data.frame(
  package = check_pkgs,
  installed = vapply(check_pkgs, is_installed, logical(1))
)

print(status)

failed <- status$package[!status$installed]

if (length(failed) > 0) {
  message("These packages still failed or need manual checking:")
  print(failed)
} else {
  message("All requested packages installed successfully.")
}

message("Done.")
