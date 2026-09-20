# ============================================================
# CoTRA native R installation
# Windows, macOS, and Linux
# ============================================================

options(
  timeout = max(2000, getOption("timeout", 60))
)

cran_repo <- "https://cloud.r-project.org"
hard_dependencies <- c("Depends", "Imports", "LinkingTo")

message("R: ", R.version.string)
message("Platform: ", R.version$platform)
message("Architecture: ", R.version$arch)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages(
    "remotes",
    repos = cran_repo,
    dependencies = hard_dependencies
  )
}

# Install CoTRA itself first without its large Suggests dependency stack.
# CoTRA::install_cotra_dependencies() then installs the analysis stack in a
# controlled order and performs the macOS system-library checks when needed.
remotes::install_github(
  "UmairSeemab/CoTRA",
  dependencies = FALSE,
  upgrade = "never",
  force = TRUE
)

library(CoTRA)
CoTRA::install_cotra_dependencies(ask = FALSE, update = FALSE)

status <- CoTRA::check_cotra_dependencies(quiet = TRUE)
if (!isTRUE(status$ok)) {
  stop(
    "CoTRA installation is incomplete. Missing or unloadable packages: ",
    paste(status$missing, collapse = ", "),
    call. = FALSE
  )
}

message("")
message("CoTRA installation completed successfully.")
message("Start CoTRA with:")
message("  library(CoTRA)")
message("  CoTRA::runCoTRA()")
