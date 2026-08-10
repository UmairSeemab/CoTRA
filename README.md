# Note before installation
Make sure you have latest version of R installed before installing CoTRA

# CoTRA

CoTRA, Comprehensive Toolbox for RNA Sequencing Data Analysis, is a Shiny-based graphical interface for bulk RNA-seq and single-cell RNA-seq workflows.

## Install from GitHub

```r
install.packages("remotes")
remotes::install_github("UmairSeemab/CoTRA", dependencies = TRUE)
```

## Install required analysis dependencies

```r
library(CoTRA)
CoTRA::install_cotra_dependencies()
```

## Run CoTRA

```r
library(CoTRA)
CoTRA::runCoTRA()
```

`runCoTRA()` creates a temporary writable copy of the Shiny app. This avoids writing output files inside the installed R package library.

## Platform notes

Windows, Ubuntu, and macOS users can install the package with the same R commands above.

External command-line tools such as FastQC, STAR, HTSeq, MultiQC, and Chrome or Chromium for PDF reports must be installed separately when those workflows are used.

## Bioconductor preparation

This package keeps reusable R functions under `R/` and the Shiny app under `inst/app/`. This layout is compatible with later Bioconductor preparation. Before Bioconductor submission, run:

```r
devtools::check()
BiocCheck::BiocCheck()
```
## Cloud computer installation (CSC)

```r
# ============================================================
# CoTRA installation on CSC Roihu
# r-env/452: R 4.5.2 + Bioconductor 3.22
# ============================================================

libpath <- "/projappl/project_2007629/CoTRA_Rlibs_452"

dir.create(
  libpath,
  recursive = TRUE,
  showWarnings = FALSE
)

.libPaths(
  c(
    libpath,
    .libPaths()
  )
)

options(
  repos = c(
    CRAN = "https://cloud.r-project.org"
  ),
  timeout = 2000
)

cat("R version:\n")
print(R.version.string)

# ------------------------------------------------------------
# BiocManager
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    lib = libpath
  )
}

BiocManager::install(
  version = "3.22",
  ask = FALSE,
  update = FALSE
)

cat("\nBioconductor version:\n")
print(BiocManager::version())

# ------------------------------------------------------------
# remotes
# ------------------------------------------------------------

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages(
    "remotes",
    lib = libpath
  )
}

# ------------------------------------------------------------
# Install CoTRA
# ------------------------------------------------------------

remotes::install_github(
  "UmairSeemab/CoTRA",
  dependencies = c(
    "Depends",
    "Imports",
    "LinkingTo"
  ),
  upgrade = "never",
  force = TRUE,
  lib = libpath
)

# ------------------------------------------------------------
# Test installation
# ------------------------------------------------------------

library(
  CoTRA,
  lib.loc = libpath
)

cat("\nCoTRA installed successfully\n")
cat("Version: ")
print(packageVersion("CoTRA"))

cat("\nInstallation path:\n")
print(find.package("CoTRA"))
```

## Output folder

After launching CoTRA, open the Home page and select an output folder. CoTRA saves generated reports, figures, CSV tables, ZIP files, and session files into this user-selected folder. If no folder is selected, CoTRA uses `~/CoTRA_Results`.
