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

## Install CoTRA without user rights
You can install CoTRA on your workstation (University, work-place) without having admin rights by following setup. Just copy it all and run all together.

```r
# ============================================================
# CoTRA - Windows installation WITHOUT administrator rights
# ============================================================
#
# This script:
#   1. Creates a personal R library
#   2. Permanently sets R_LIBS_USER
#   3. Installs CoTRA dependencies into that library
#   4. Installs CoTRA from GitHub
#
# Administrator rights are NOT required.
#
# CoTRA:
# https://github.com/UmairSeemab/CoTRA
# ============================================================


# ------------------------------------------------------------
# 1. Check operating system
# ------------------------------------------------------------

if (.Platform$OS.type != "windows") {
  message(
    "This installer is designed primarily for Windows.\n",
    "Installation will continue using a user library."
  )
}


# ------------------------------------------------------------
# 2. Determine R version
# ------------------------------------------------------------

r_version <- paste0(
  R.version$major,
  ".",
  strsplit(R.version$minor, "\\.")[[1]][1]
)

message("R version: ", R.version.string)


# ------------------------------------------------------------
# 3. Create a user-owned R library
# ------------------------------------------------------------

local_appdata <- Sys.getenv("LOCALAPPDATA")

if (nzchar(local_appdata)) {

  user_lib <- file.path(
    local_appdata,
    "R",
    "win-library",
    r_version
  )

} else {

  # Fallback if LOCALAPPDATA is unavailable
  user_lib <- file.path(
    path.expand("~"),
    "R",
    "win-library",
    r_version
  )
}

# Convert Windows backslashes to forward slashes
user_lib <- gsub("\\\\", "/", user_lib)


message("\nCoTRA packages will be installed in:")
message(user_lib)


if (!dir.exists(user_lib)) {

  message("\nCreating personal R library...")

  dir.create(
    user_lib,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


if (!dir.exists(user_lib)) {
  stop(
    "\nCould not create the user library:\n",
    user_lib,
    "\n\nPlease check that you have write permission to your user profile."
  )
}


# ------------------------------------------------------------
# 4. Permanently configure R_LIBS_USER
# ------------------------------------------------------------

renviron_file <- file.path(
  path.expand("~"),
  ".Renviron"
)

renviron_file <- gsub("\\\\", "/", renviron_file)

message("\nConfiguring permanent user library...")
message(".Renviron file: ", renviron_file)


# Read existing .Renviron without deleting other settings
if (file.exists(renviron_file)) {

  renviron_lines <- readLines(
    renviron_file,
    warn = FALSE
  )

} else {

  renviron_lines <- character(0)
}


# Remove previous R_LIBS_USER setting, if present
renviron_lines <- renviron_lines[
  !grepl(
    "^\\s*R_LIBS_USER\\s*=",
    renviron_lines
  )
]


# Add permanent R_LIBS_USER
new_setting <- paste0(
  'R_LIBS_USER="',
  user_lib,
  '"'
)

renviron_lines <- c(
  renviron_lines,
  new_setting
)


writeLines(
  renviron_lines,
  renviron_file
)


message(
  "Permanent R_LIBS_USER configured successfully."
)


# ------------------------------------------------------------
# 5. Activate the library in CURRENT R session
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


message("\nCurrent R library paths:")

for (x in .libPaths()) {
  message("  ", x)
}


# ------------------------------------------------------------
# 6. Check write permission
# ------------------------------------------------------------

test_file <- file.path(
  user_lib,
  "CoTRA_write_test.txt"
)

write_test <- tryCatch(
  {

    writeLines(
      "CoTRA write test",
      test_file
    )

    unlink(test_file)

    TRUE

  },
  error = function(e) FALSE
)


if (!write_test) {

  stop(
    "\nR cannot write to:\n",
    user_lib,
    "\n\nInstallation cannot continue."
  )

} else {

  message(
    "\nUser library is writable."
  )
}


# ------------------------------------------------------------
# 7. Configure CRAN
# ------------------------------------------------------------

options(
  repos = c(
    CRAN = "https://cloud.r-project.org"
  )
)

options(
  timeout = 1000
)


# ------------------------------------------------------------
# 8. Install BiocManager
# ------------------------------------------------------------

if (!requireNamespace(
  "BiocManager",
  quietly = TRUE,
  lib.loc = user_lib
)) {

  message(
    "\nInstalling BiocManager..."
  )

  install.packages(
    "BiocManager",
    lib = user_lib,
    dependencies = TRUE
  )
}


# ------------------------------------------------------------
# 9. Configure Bioconductor repositories
# ------------------------------------------------------------

suppressPackageStartupMessages(
  library(
    BiocManager,
    lib.loc = user_lib
  )
)

options(
  repos = BiocManager::repositories()
)


# ------------------------------------------------------------
# 10. Install remotes
# ------------------------------------------------------------

if (!requireNamespace(
  "remotes",
  quietly = TRUE,
  lib.loc = user_lib
)) {

  message(
    "\nInstalling remotes..."
  )

  install.packages(
    "remotes",
    lib = user_lib,
    dependencies = TRUE
  )
}


# ------------------------------------------------------------
# 11. Install CoTRA
# ------------------------------------------------------------

message(
  "\n===================================================="
)

message(
  "Installing CoTRA and its dependencies..."
)

message(
  "====================================================\n"
)


remotes::install_github(
  "UmairSeemab/CoTRA",
  lib = user_lib,
  dependencies = TRUE,
  upgrade = "never",
  build_vignettes = FALSE
)


# ------------------------------------------------------------
# 12. Verify installation
# ------------------------------------------------------------

message(
  "\n===================================================="
)

message(
  "Checking CoTRA installation..."
)

message(
  "====================================================\n"
)


if (
  requireNamespace(
    "CoTRA",
    quietly = TRUE,
    lib.loc = user_lib
  )
) {

  cotra_location <- find.package(
    "CoTRA",
    lib.loc = user_lib
  )

  message(
    "SUCCESS: CoTRA has been installed."
  )

  message(
    "\nCoTRA location:"
  )

  message(
    cotra_location
  )

} else {

  stop(
    "\nCoTRA installation did not complete successfully.\n",
    "Please review the installation messages above."
  )
}


# ------------------------------------------------------------
# 13. Final instructions
# ------------------------------------------------------------

cat(
  "\n\n",
  "============================================================\n",
  "                 CoTRA INSTALLATION COMPLETE\n",
  "============================================================\n\n",
  "Administrator rights were not required.\n\n",
  "Personal R library:\n",
  user_lib,
  "\n\n",
  "The library path has been saved permanently in:\n",
  renviron_file,
  "\n\n",
  "You can now restart R/RStudio and run:\n\n",
  "    library(CoTRA)\n",
  "    runCoTRA()\n\n",
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
