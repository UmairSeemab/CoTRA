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
