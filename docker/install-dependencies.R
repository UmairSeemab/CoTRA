options(timeout = 2000)
options(repos = c(CRAN = Sys.getenv("CRAN", unset = "https://cloud.r-project.org")))
options(Ncpus = max(1L, parallel::detectCores(logical = TRUE) - 1L))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Reuse CoTRA's own dependency declarations rather than maintaining a second list.
source("/opt/CoTRA/R/dependencies.R", local = TRUE)

cran <- unique(c(cotra_cran_packages(), "hdf5r"))
missing_cran <- cran[
  !vapply(cran, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, dependencies = TRUE)
}

bioc <- unique(cotra_bioc_packages())
missing_bioc <- bioc[
  !vapply(bioc, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

github <- cotra_github_packages()
for (pkg in names(github)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    remotes::install_github(
      github[[pkg]],
      dependencies = TRUE,
      upgrade = "never"
    )
  }
}

all_required <- unique(c(cran, bioc, names(github)))
missing <- all_required[
  !vapply(all_required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing) > 0) {
  stop(
    "Container dependency installation is incomplete. Missing: ",
    paste(missing, collapse = ", ")
  )
}

message("All declared CoTRA container dependencies are installed.")
