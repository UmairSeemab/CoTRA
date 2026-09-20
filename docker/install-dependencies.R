options(timeout = 2000)

options(
  repos = c(
    CRAN = Sys.getenv(
      "CRAN",
      unset = "https://cloud.r-project.org"
    )
  )
)

options(
  Ncpus = max(
    1L,
    parallel::detectCores(logical = TRUE) - 1L
  )
)

# ------------------------------------------------------------
# Bootstrap installers
# ------------------------------------------------------------

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Reuse CoTRA's own dependency declarations.
source(
  "/opt/CoTRA/R/dependencies.R",
  local = TRUE
)

runtime_dependencies <- cotra_runtime_dependency_types()


# ------------------------------------------------------------
# CRAN dependencies
# ------------------------------------------------------------

cran <- unique(
  c(
    cotra_cran_packages(),
    "hdf5r"
  )
)

missing_cran <- cran[
  !vapply(
    cran,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_cran) > 0) {
  message(
    "Installing CRAN packages: ",
    paste(missing_cran, collapse = ", ")
  )

  install.packages(
    missing_cran,
    dependencies = runtime_dependencies
  )
}


# ------------------------------------------------------------
# Bioconductor dependencies
# ------------------------------------------------------------

bioc <- unique(
  cotra_bioc_packages()
)

missing_bioc <- bioc[
  !vapply(
    bioc,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_bioc) > 0) {
  message(
    "Installing Bioconductor packages: ",
    paste(missing_bioc, collapse = ", ")
  )

  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}


# ------------------------------------------------------------
# Explicit CellChat Bioconductor dependency verification
# ------------------------------------------------------------

cellchat_bioc <- c(
  "ComplexHeatmap",
  "BiocNeighbors",
  "BiocGenerics"
)

missing_cellchat_bioc <- cellchat_bioc[
  !vapply(
    cellchat_bioc,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_cellchat_bioc) > 0) {
  message(
    "Installing CellChat Bioconductor requirements: ",
    paste(missing_cellchat_bioc, collapse = ", ")
  )

  BiocManager::install(
    missing_cellchat_bioc,
    ask = FALSE,
    update = FALSE
  )
}


# ------------------------------------------------------------
# GitHub dependencies
# ------------------------------------------------------------

github <- cotra_github_packages()

for (pkg in names(github)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(
      "Installing GitHub package ",
      pkg,
      " from ",
      github[[pkg]]
    )

    remotes::install_github(
      github[[pkg]],
      dependencies = runtime_dependencies,
      upgrade = "never"
    )
  }
}


# ------------------------------------------------------------
# Final verification
# ------------------------------------------------------------

all_required <- unique(
  c(
    cran,
    bioc,
    names(github)
  )
)

missing <- all_required[
  !vapply(
    all_required,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing) > 0) {
  stop(
    "Container dependency installation is incomplete. Missing: ",
    paste(missing, collapse = ", ")
  )
}

message(
  "All declared CoTRA container dependencies are installed."
)
