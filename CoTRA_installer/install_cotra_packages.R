# ==========================================================
# CoTRA package installation script
# Compatible with R 4.6 and Bioconductor 3.23
# ==========================================================

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  ask = FALSE,
  timeout = 1000,
  Ncpus = max(1, parallel::detectCores() - 1)
)

message("R version: ", getRversion())

is_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

install_cran_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, is_installed, logical(1))]
  
  if (length(missing) > 0) {
    message("Installing CRAN packages:")
    print(missing)
    
    install.packages(
      missing,
      dependencies = TRUE,
      type = "source"
    )
  } else {
    message("All CRAN packages already installed.")
  }
}

setup_bioc <- function() {
  if (!is_installed("BiocManager")) {
    install.packages("BiocManager")
  }
  
  r_ver <- as.character(getRversion())
  
  if (startsWith(r_ver, "4.6")) {
    BiocManager::install(version = "3.23", ask = FALSE, update = TRUE)
  } else {
    message("Using default Bioconductor version for R ", r_ver)
  }
  
  message("Bioconductor version: ", BiocManager::version())
}

install_bioc_if_missing <- function(pkgs) {
  setup_bioc()
  
  missing <- pkgs[!vapply(pkgs, is_installed, logical(1))]
  
  if (length(missing) > 0) {
    message("Installing Bioconductor packages:")
    print(missing)
    
    BiocManager::install(
      missing,
      ask = FALSE,
      update = TRUE,
      dependencies = TRUE
    )
  } else {
    message("All Bioconductor packages already installed.")
  }
}

install_github_if_missing <- function(pkg, repo) {
  if (!is_installed("remotes")) {
    install.packages("remotes")
  }
  
  if (!is_installed(pkg)) {
    message("Installing GitHub package: ", pkg)
    
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

cran_core_pkgs <- c(
  "BiocManager",
  "remotes",
  "Rcpp",
  "RcppArmadillo",
  "Matrix",
  "cli",
  "fs",
  "cpp11",
  "data.table",
  "dplyr",
  "ggplot2",
  "plotly",
  "DT",
  "shiny",
  "bs4Dash",
  "shinyjs",
  "shinyWidgets",
  "shinyalert",
  "shinyFiles",
  "shinydashboard",
  "shinycssloaders",
  "htmlwidgets",
  "knitr",
  "rmarkdown",
  "scales",
  "stringr",
  "igraph",
  "openxlsx",
  "msigdbr",
  "lme4"
)

cran_optional_pkgs <- c(
  "svglite",
  "ragg",
  "textshaping",
  "pagedown",
  "webshot2",
  "kableExtra",
  "devtools",
  "pkgdown",
  "roxygen2",
  "NMF"
)

bioc_pkgs <- c(
  "AnnotationDbi",
  "BiocGenerics",
  "Biobase",
  "S4Vectors",
  "IRanges",
  "DelayedArray",
  "DelayedMatrixStats",
  "SummarizedExperiment",
  "SingleCellExperiment",
  "HDF5Array",
  "rhdf5",
  "rhdf5filters",
  "Rhdf5lib",
  "limma",
  "edgeR",
  "DESeq2",
  "batchelor",
  "ggrastr",
  "slingshot",
  "UCell",
  "AUCell",
  "GSEABase",
  "GSVA",
  "pcaMethods",
  "clusterProfiler",
  "ReactomePA",
  "enrichplot",
  "pathview",
  "fgsea",
  "DirichletMultinomial",
  "TFBSTools",
  "ensembldb",
  "org.Mm.eg.db",
  "org.Hs.eg.db",
  "org.Rn.eg.db",
  "EnsDb.Mmusculus.v79",
  "EnsDb.Hsapiens.v86",
  "EnsDb.Rnorvegicus.v79",
  "ComplexHeatmap",
  "circlize"
)

github_pkgs <- list(
  azimuth = "satijalab/azimuth",
  BPCells = "bnprks/BPCells/r",
  monocle3 = "cole-trapnell-lab/monocle3",
  CellChat = "jinworks/CellChat",
  nichenetr = "saeyslab/nichenetr",
  liana = "saezlab/liana"
)

message("Installing missing core CRAN packages...")
install_cran_if_missing(cran_core_pkgs)

message("Installing missing Bioconductor packages...")
install_bioc_if_missing(bioc_pkgs)

message("Installing missing optional CRAN packages...")
install_cran_if_missing(cran_optional_pkgs)

message("Installing missing GitHub packages...")
for (pkg in names(github_pkgs)) {
  install_github_if_missing(pkg, github_pkgs[[pkg]])
}

message("Checking installed packages...")

check_pkgs <- unique(c(
  cran_core_pkgs,
  cran_optional_pkgs,
  bioc_pkgs,
  names(github_pkgs)
))

check_pkgs <- setdiff(check_pkgs, "h5")

status <- data.frame(
  package = check_pkgs,
  installed = vapply(check_pkgs, is_installed, logical(1)),
  stringsAsFactors = FALSE
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