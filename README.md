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

## Install CoTRA without Admin rights
You can install CoTRA on your workstation (University, work-place) without having admin rights by following setup. Just copy it all and run all together.

```r
# ============================================================
# CoTRA - Windows Installation WITHOUT Administrator Rights
# ============================================================
#
# Run this entire script once in R or RStudio.
#
# It will:
#   - Create a personal R package library
#   - Save that library permanently using R_LIBS_USER
#   - Install required installation packages
#   - Install CoTRA and its dependencies
#
# Administrator rights are NOT required.
#
# CoTRA:
# https://github.com/UmairSeemab/CoTRA
# ============================================================


# ------------------------------------------------------------
# 1. Check Windows
# ------------------------------------------------------------

if (.Platform$OS.type != "windows") {
  stop(
    "This installation script is intended for Windows computers."
  )
}


# ------------------------------------------------------------
# 2. R version
# ------------------------------------------------------------

r_version <- paste0(
  R.version$major,
  ".",
  strsplit(R.version$minor, "\\.")[[1]][1]
)

cat("\nR version:\n")
cat(R.version.string, "\n")


# ------------------------------------------------------------
# 3. Create personal R library
# ------------------------------------------------------------

local_appdata <- Sys.getenv("LOCALAPPDATA")

if (!nzchar(local_appdata)) {
  local_appdata <- path.expand("~")
}

user_lib <- file.path(
  local_appdata,
  "R",
  "win-library",
  r_version
)

user_lib <- gsub(
  "\\\\",
  "/",
  user_lib
)

cat("\nPersonal R library:\n")
cat(user_lib, "\n")


if (!dir.exists(user_lib)) {

  cat("\nCreating personal R library...\n")

  dir.create(
    user_lib,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


if (!dir.exists(user_lib)) {

  stop(
    "\nERROR: Could not create personal R library:\n",
    user_lib,
    "\n\nCheck that your Windows user account can write ",
    "to its AppData directory."
  )
}


# ------------------------------------------------------------
# 4. Permanently set R_LIBS_USER
# ------------------------------------------------------------

user_home <- path.expand("~")

renviron_file <- file.path(
  user_home,
  ".Renviron"
)

renviron_file <- gsub(
  "\\\\",
  "/",
  renviron_file
)


cat("\nConfiguring permanent R library...\n")
cat(".Renviron location:\n")
cat(renviron_file, "\n")


if (file.exists(renviron_file)) {

  renviron_lines <- readLines(
    renviron_file,
    warn = FALSE
  )

} else {

  renviron_lines <- character(0)
}


# Remove any existing R_LIBS_USER entry
renviron_lines <- renviron_lines[
  !grepl(
    "^\\s*R_LIBS_USER\\s*=",
    renviron_lines
  )
]


# Add permanent personal library
renviron_lines <- c(
  renviron_lines,
  paste0(
    'R_LIBS_USER="',
    user_lib,
    '"'
  )
)


writeLines(
  renviron_lines,
  renviron_file
)


# ------------------------------------------------------------
# 5. Activate library NOW
# ------------------------------------------------------------

Sys.setenv(
  R_LIBS_USER = user_lib
)

.libPaths(
  unique(
    c(
      user_lib,
      .libPaths()
    )
  )
)


cat("\nCurrent R library paths:\n")

for (x in .libPaths()) {
  cat("  ", x, "\n")
}


# ------------------------------------------------------------
# 6. Confirm that library is writable
# ------------------------------------------------------------

test_file <- file.path(
  user_lib,
  "CoTRA_write_test.txt"
)


write_ok <- tryCatch({

  writeLines(
    "CoTRA write test",
    test_file
  )

  unlink(test_file)

  TRUE

}, error = function(e) {

  FALSE

})


if (!write_ok) {

  stop(
    "\nERROR: R cannot write to:\n",
    user_lib,
    "\n\nInstallation cannot continue."
  )
}


cat("\n[OK] Personal R library is writable.\n")


# ------------------------------------------------------------
# 7. General R installation settings
# ------------------------------------------------------------

options(
  repos = c(
    CRAN = "https://cloud.r-project.org"
  )
)

options(
  timeout = 2000
)

options(
  Ncpus = max(
    1,
    parallel::detectCores(logical = TRUE) - 1
  )
)


# ------------------------------------------------------------
# 8. Function for installing required CRAN packages
# ------------------------------------------------------------

install_required_cran <- function(pkg) {

  if (requireNamespace(
    pkg,
    quietly = TRUE
  )) {

    cat(
      "[OK]",
      pkg,
      "already installed:",
      as.character(packageVersion(pkg)),
      "\n"
    )

    return(invisible(TRUE))
  }


  cat(
    "\nInstalling required package:",
    pkg,
    "\n"
  )


  # Prefer Windows binary packages.
  # This avoids compilation when a binary is available.
  tryCatch({

    install.packages(
      pkg,
      lib = user_lib,
      repos = "https://cloud.r-project.org",
      type = "binary",
      dependencies = c(
        "Depends",
        "Imports",
        "LinkingTo"
      )
    )

  }, error = function(e) {

    cat(
      "\nInstallation error for",
      pkg,
      ":\n",
      conditionMessage(e),
      "\n"
    )
  })


  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {

    stop(
      "\nERROR: Required package '",
      pkg,
      "' could not be installed.\n\n",
      "CoTRA installation cannot continue."
    )
  }


  cat(
    "[OK]",
    pkg,
    "installed:",
    as.character(packageVersion(pkg)),
    "\n"
  )
}


# ------------------------------------------------------------
# 9. Install BiocManager
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)
cat(
  "Installing required installation packages\n"
)
cat(
  "============================================================\n\n"
)


install_required_cran(
  "BiocManager"
)


# ------------------------------------------------------------
# 10. Configure CRAN + Bioconductor repositories
# ------------------------------------------------------------

options(
  repos = BiocManager::repositories()
)


cat("\nBioconductor version detected:\n")
cat(
  as.character(
    BiocManager::version()
  ),
  "\n"
)


# ------------------------------------------------------------
# 11. Install remotes
# ------------------------------------------------------------

install_required_cran(
  "remotes"
)


# ------------------------------------------------------------
# 12. Install hdf5r explicitly
# ------------------------------------------------------------
#
# IMPORTANT:
# Seurat::Read10X_h5() requires hdf5r.
#
# CoTRA needs this package for:
#   10x Genomics HDF5 (.h5) input
#
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)
cat(
  "Installing HDF5 support for 10x Genomics .h5 files\n"
)
cat(
  "============================================================\n\n"
)


install_required_cran(
  "hdf5r"
)


# ------------------------------------------------------------
# 13. Verify hdf5r before installing CoTRA
# ------------------------------------------------------------

if (!requireNamespace(
  "hdf5r",
  quietly = TRUE
)) {

  stop(
    "\nERROR: hdf5r is unavailable.\n",
    "10x Genomics HDF5 (.h5) files cannot be imported.\n"
  )
}


cat(
  "\n[OK] hdf5r version: ",
  as.character(
    packageVersion("hdf5r")
  ),
  "\n",
  sep = ""
)


cat(
  "[OK] hdf5r location: ",
  find.package("hdf5r"),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 14. Install CoTRA
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "Installing CoTRA and dependencies\n"
)

cat(
  "============================================================\n\n"
)


tryCatch({

  remotes::install_github(
    "UmairSeemab/CoTRA",
    lib = user_lib,
    dependencies = TRUE,
    upgrade = "never",
    build_vignettes = FALSE,
    force = TRUE
  )

}, error = function(e) {

  stop(
    "\nCoTRA installation failed.\n\n",
    conditionMessage(e)
  )
})


# ------------------------------------------------------------
# 15. Verify CoTRA
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "Verifying CoTRA installation\n"
)

cat(
  "============================================================\n\n"
)


required_checks <- c(
  "CoTRA",
  "Seurat",
  "hdf5r"
)


check_results <- sapply(
  required_checks,
  function(pkg) {

    requireNamespace(
      pkg,
      quietly = TRUE
    )

  }
)


for (pkg in required_checks) {

  if (check_results[[pkg]]) {

    cat(
      "[OK] ",
      pkg,
      " ",
      as.character(
        packageVersion(pkg)
      ),
      "\n",
      sep = ""
    )

  } else {

    cat(
      "[MISSING] ",
      pkg,
      "\n",
      sep = ""
    )
  }
}


if (!all(check_results)) {

  missing_packages <- names(
    check_results
  )[!check_results]


  stop(
    "\nInstallation is incomplete.\n\n",
    "Missing package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\n"
  )
}


# ------------------------------------------------------------
# 16. Verify Seurat H5 functionality
# ------------------------------------------------------------

cat(
  "\nChecking 10x HDF5 support...\n"
)


if (!exists(
  "Read10X_h5",
  envir = asNamespace("Seurat"),
  inherits = FALSE
)) {

  stop(
    "\nERROR: Seurat::Read10X_h5() is unavailable."
  )
}


if (!requireNamespace(
  "hdf5r",
  quietly = TRUE
)) {

  stop(
    "\nERROR: hdf5r is unavailable."
  )
}


cat(
  "[OK] Seurat::Read10X_h5() is available.\n"
)

cat(
  "[OK] hdf5r is available.\n"
)

cat(
  "[OK] CoTRA is ready for 10x Genomics HDF5 (.h5) files.\n"
)


# ------------------------------------------------------------
# 17. Show installation locations
# ------------------------------------------------------------

cat(
  "\nInstalled package locations:\n\n"
)


cat(
  "CoTRA:\n",
  find.package("CoTRA"),
  "\n\n",
  sep = ""
)


cat(
  "Seurat:\n",
  find.package("Seurat"),
  "\n\n",
  sep = ""
)


cat(
  "hdf5r:\n",
  find.package("hdf5r"),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 18. Final message
# ------------------------------------------------------------

cat(
  "\n",
  "============================================================\n",
  "             CoTRA INSTALLATION COMPLETE\n",
  "============================================================\n\n",
  "Administrator rights were NOT required.\n\n",
  "Personal R library:\n",
  user_lib,
  "\n\n",
  "The library path has been permanently saved in:\n",
  renviron_file,
  "\n\n",
  "10x Genomics HDF5 (.h5) support:\n",
  "AVAILABLE\n\n",
  "Restart R or RStudio once.\n\n",
  "Then start CoTRA using:\n\n",
  "    library(CoTRA)\n",
  "    runCoTRA()\n\n",
  "You do NOT need to run this installation script again.\n\n",
  "============================================================\n",
  sep = ""
)
```

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
  version = "3.23",
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
install.packages("remotes")
remotes::install_github("UmairSeemab/CoTRA", dependencies = TRUE)
library(CoTRA)
CoTRA::install_cotra_dependencies()
cat("\nCoTRA installed successfully\n")
cat("Version: ")
print(packageVersion("CoTRA"))

cat("\nInstallation path:\n")
print(find.package("CoTRA"))
```

## Output folder

After launching CoTRA, open the Home page and select an output folder. CoTRA saves generated reports, figures, CSV tables, ZIP files, and session files into this user-selected folder. If no folder is selected, CoTRA uses `~/CoTRA_Results`.
