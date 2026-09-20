options(timeout = 2000)
options(Ncpus = max(1L, parallel::detectCores(logical = TRUE) - 1L))

cran_repo <- Sys.getenv("CRAN", unset = "https://cloud.r-project.org")
hard_dependencies <- c("Depends", "Imports", "LinkingTo")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages(
    "remotes",
    repos = cran_repo,
    dependencies = hard_dependencies
  )
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = cran_repo,
    dependencies = hard_dependencies
  )
}

# Reuse the dependency definitions and installation logic shipped with CoTRA.
source("/opt/CoTRA/R/dependencies.R", local = TRUE)

install_cotra_dependencies(ask = FALSE, update = FALSE)

status <- check_cotra_dependencies(quiet = TRUE)
if (!isTRUE(status$ok)) {
  stop(
    "Container dependency installation is incomplete. Missing or unloadable: ",
    paste(status$missing, collapse = ", ")
  )
}

message("All declared CoTRA container dependencies are installed and loadable.")
