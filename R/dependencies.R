.cotra_loadable <- function(pkg) {
  isTRUE(
    suppressWarnings(
      suppressMessages(
        requireNamespace(pkg, quietly = TRUE)
      )
    )
  )
}

.cotra_installed_names <- function() {
  rownames(utils::installed.packages())
}

.cotra_prepend_path <- function(current, values, sep = .Platform$path.sep) {
  values <- unique(values[nzchar(values) & !is.na(values)])
  current_parts <- if (nzchar(current)) strsplit(current, sep, fixed = TRUE)[[1]] else character()
  paste(unique(c(values, current_parts)), collapse = sep)
}

.cotra_brew_path <- function() {
  candidates <- c(
    Sys.which("brew"),
    "/opt/homebrew/bin/brew",
    "/usr/local/bin/brew"
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  candidates[file.exists(candidates)][1] %||% ""
}

`%||%` <- function(x, y) {
  if (length(x) == 0L || is.na(x) || !nzchar(x)) y else x
}

.cotra_brew_prefix <- function(brew, formula = NULL) {
  if (!nzchar(brew)) return("")

  args <- if (is.null(formula)) "--prefix" else c("--prefix", formula)
  out <- suppressWarnings(
    system2(brew, args, stdout = TRUE, stderr = TRUE)
  )

  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) return("")
  if (length(out) == 0L) return("")

  trimws(out[[1]])
}

.cotra_prepare_macos_build <- function() {
  if (!identical(Sys.info()[["sysname"]], "Darwin")) {
    return(NULL)
  }

  xcode_select <- Sys.which("xcode-select")
  if (nzchar(xcode_select)) {
    xcode_check <- suppressWarnings(
      system2(xcode_select, "-p", stdout = TRUE, stderr = TRUE)
    )
    xcode_status <- attr(xcode_check, "status")
    if (!is.null(xcode_status) && xcode_status != 0L) {
      stop(
        paste0(
          "Apple Command Line Tools are required to compile some CoTRA dependencies.\n",
          "Open Terminal and run:\n\n",
          "  xcode-select --install\n\n",
          "After installation finishes, restart R/RStudio and run ",
          "CoTRA::install_cotra_dependencies() again."
        ),
        call. = FALSE
      )
    }
  }

  brew <- .cotra_brew_path()
  if (!nzchar(brew)) {
    stop(
      paste0(
        "Homebrew was not found. Some CoTRA dependencies on macOS require ",
        "HDF5 and OpenSSL when source compilation is needed.\n\n",
        "Install Homebrew from https://brew.sh, then run in Terminal:\n\n",
        "  brew install hdf5 pkgconf openssl@3\n\n",
        "Restart R/RStudio and run CoTRA::install_cotra_dependencies() again.\n",
        "Alternatively, use the recommended CoTRA Docker container."
      ),
      call. = FALSE
    )
  }

  brew_prefix <- .cotra_brew_prefix(brew)
  hdf5_prefix <- .cotra_brew_prefix(brew, "hdf5")
  openssl_prefix <- .cotra_brew_prefix(brew, "openssl@3")
  pkgconf_prefix <- .cotra_brew_prefix(brew, "pkgconf")

  missing_formulae <- character()
  if (!nzchar(hdf5_prefix)) missing_formulae <- c(missing_formulae, "hdf5")
  if (!nzchar(openssl_prefix)) missing_formulae <- c(missing_formulae, "openssl@3")
  if (!nzchar(pkgconf_prefix)) missing_formulae <- c(missing_formulae, "pkgconf")

  if (length(missing_formulae) > 0L) {
    stop(
      paste0(
        "Missing macOS system libraries required by CoTRA source packages: ",
        paste(missing_formulae, collapse = ", "),
        ".\n\nOpen Terminal and run:\n\n",
        "  brew install hdf5 pkgconf openssl@3\n\n",
        "Restart R/RStudio and run CoTRA::install_cotra_dependencies() again."
      ),
      call. = FALSE
    )
  }

  env_names <- c(
    "PATH", "PKG_CONFIG_PATH", "CPATH", "LIBRARY_PATH",
    "CPPFLAGS", "LDFLAGS", "R_MAKEVARS_USER"
  )
  old_env <- Sys.getenv(env_names, unset = NA_character_)
  names(old_env) <- env_names

  path_entries <- c(
    file.path(hdf5_prefix, "bin"),
    file.path(pkgconf_prefix, "bin"),
    file.path(brew_prefix, "bin")
  )

  pkgconfig_entries <- c(
    file.path(hdf5_prefix, "lib", "pkgconfig"),
    file.path(openssl_prefix, "lib", "pkgconfig"),
    file.path(pkgconf_prefix, "lib", "pkgconfig")
  )

  include_entries <- c(
    file.path(openssl_prefix, "include"),
    file.path(hdf5_prefix, "include")
  )

  library_entries <- c(
    file.path(openssl_prefix, "lib"),
    file.path(hdf5_prefix, "lib")
  )

  Sys.setenv(
    PATH = .cotra_prepend_path(Sys.getenv("PATH"), path_entries),
    PKG_CONFIG_PATH = .cotra_prepend_path(
      Sys.getenv("PKG_CONFIG_PATH"),
      pkgconfig_entries
    ),
    CPATH = .cotra_prepend_path(Sys.getenv("CPATH"), include_entries),
    LIBRARY_PATH = .cotra_prepend_path(
      Sys.getenv("LIBRARY_PATH"),
      library_entries
    )
  )

  old_cppflags <- Sys.getenv("CPPFLAGS")
  old_ldflags <- Sys.getenv("LDFLAGS")

  cppflags <- paste(
    paste0("-I", shQuote(include_entries)),
    collapse = " "
  )
  ldflags <- paste(
    paste0("-L", shQuote(library_entries)),
    collapse = " "
  )

  Sys.setenv(
    CPPFLAGS = trimws(paste(cppflags, old_cppflags)),
    LDFLAGS = trimws(paste(ldflags, old_ldflags))
  )

  # R package compilation reads Makevars rather than shell flags in some cases.
  # Create a temporary Makevars file that preserves the user's existing settings.
  existing_makevars <- Sys.getenv("R_MAKEVARS_USER")
  if (!nzchar(existing_makevars)) {
    existing_makevars <- file.path(path.expand("~"), ".R", "Makevars")
  }

  makevars_lines <- character()
  if (file.exists(existing_makevars)) {
    makevars_lines <- readLines(existing_makevars, warn = FALSE)
  }

  tmp_makevars <- tempfile("CoTRA-Makevars-")
  writeLines(
    c(
      makevars_lines,
      "",
      "# Temporary CoTRA macOS build flags",
      paste0("CPPFLAGS += -I", openssl_prefix, "/include -I", hdf5_prefix, "/include"),
      paste0("LDFLAGS += -L", openssl_prefix, "/lib -L", hdf5_prefix, "/lib")
    ),
    tmp_makevars
  )
  Sys.setenv(R_MAKEVARS_USER = tmp_makevars)

  pkg_config <- Sys.which("pkg-config")
  h5cc <- Sys.which("h5cc")

  if (!nzchar(pkg_config) || !nzchar(h5cc)) {
    stop(
      paste0(
        "CoTRA found Homebrew libraries but could not find pkg-config or h5cc.\n",
        "Run in Terminal:\n\n",
        "  brew reinstall hdf5 pkgconf\n\n",
        "Then restart R/RStudio and retry."
      ),
      call. = FALSE
    )
  }

  message(
    "Configured macOS build environment for ",
    R.version$arch,
    " using Homebrew HDF5 and OpenSSL."
  )

  list(old_env = old_env, makevars = tmp_makevars)
}

.cotra_restore_macos_build <- function(config) {
  if (is.null(config)) return(invisible(NULL))

  for (nm in names(config$old_env)) {
    value <- config$old_env[[nm]]
    if (is.na(value)) {
      Sys.unsetenv(nm)
    } else {
      do.call(Sys.setenv, setNames(list(value), nm))
    }
  }

  if (!is.null(config$makevars) && file.exists(config$makevars)) {
    unlink(config$makevars)
  }

  invisible(NULL)
}

cotra_cran_packages <- function() {
  c(
    "shiny", "shinyFiles", "shinydashboard", "bs4Dash", "shinyWidgets",
    "shinyjs", "plotly", "DT", "dplyr", "Seurat", "SeuratObject",
    "Matrix", "data.table", "ggplot2", "shinyalert", "knitr",
    "kableExtra", "scales", "stringr", "shinycssloaders", "svglite",
    "htmlwidgets", "msigdbr", "rmarkdown", "gridExtra", "cowplot",
    "sctransform", "corrplot", "igraph", "mgcv", "pagedown", "webshot2",
    "pheatmap", "RColorBrewer", "readr", "tidyr", "openxlsx", "hdf5r"
  )
}

cotra_bioc_packages <- function() {
  c(
    "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db",
    "enrichplot", "fgsea", "ReactomePA", "reactome.db", "pathview",
    "SingleCellExperiment", "scDblFinder", "MAST", "SummarizedExperiment",
    "slingshot", "UCell", "AUCell", "GSVA", "SingleR", "celldex",
    "alabaster.base", "ComplexHeatmap", "limma", "edgeR", "DESeq2",
    "AnnotationDbi", "biomaRt", "GenomeInfoDb", "IRanges", "S4Vectors",
    "BiocGenerics", "ensembldb", "EnsDb.Hsapiens.v86",
    "EnsDb.Mmusculus.v79", "EnsDb.Rnorvegicus.v79"
  )
}

cotra_github_packages <- function() {
  c(
    DoubletFinder = "chris-mcginnis-ucsf/DoubletFinder",
    BPCells = "bnprks/BPCells/r",
    monocle3 = "cole-trapnell-lab/monocle3",
    CellChat = "jinworks/CellChat"
  )
}

#' Check CoTRA dependencies
#'
#' @param quiet Logical. Suppress messages.
#' @export
check_cotra_dependencies <- function(quiet = FALSE) {
  pkgs <- unique(
    c(
      cotra_cran_packages(),
      cotra_bioc_packages(),
      names(cotra_github_packages())
    )
  )

  installed <- vapply(
    pkgs,
    .cotra_loadable,
    FUN.VALUE = logical(1)
  )
  missing <- names(installed)[!installed]

  result <- list(
    ok = length(missing) == 0L,
    installed = names(installed)[installed],
    missing = missing
  )

  if (!quiet) {
    if (result$ok) {
      message("All CoTRA dependencies are installed and loadable.")
    } else {
      message(
        "Missing or unloadable CoTRA dependencies: ",
        paste(missing, collapse = ", ")
      )
    }
  }

  result
}

#' Install CoTRA dependencies
#'
#' Installs CRAN, Bioconductor, and selected GitHub dependencies required by
#' the CoTRA app. On macOS, Homebrew HDF5/OpenSSL paths are configured when
#' source compilation is required.
#'
#' @param ask Logical. Ask before installing missing packages.
#' @param update Logical. Passed to BiocManager::install().
#' @export
install_cotra_dependencies <- function(ask = interactive(), update = FALSE) {
  options(timeout = max(2000, getOption("timeout", 60)))

  cran_repo <- "https://cloud.r-project.org"
  hard_dependencies <- c("Depends", "Imports", "LinkingTo")

  if (!requireNamespace("remotes", quietly = TRUE)) {
    utils::install.packages(
      "remotes",
      repos = cran_repo,
      dependencies = hard_dependencies
    )
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    utils::install.packages(
      "BiocManager",
      repos = cran_repo,
      dependencies = hard_dependencies
    )
  }

  cran <- cotra_cran_packages()
  bioc <- cotra_bioc_packages()
  github <- cotra_github_packages()

  missing_cran <- cran[!vapply(cran, .cotra_loadable, logical(1))]
  missing_bioc <- bioc[!vapply(bioc, .cotra_loadable, logical(1))]
  missing_github <- names(github)[
    !vapply(names(github), .cotra_loadable, logical(1))
  ]

  all_missing <- unique(c(missing_cran, missing_bioc, missing_github))

  mac_source_packages <- c(
    "alabaster.base", "celldex", "BPCells", "monocle3"
  )

  mac_config <- NULL
  if (
    identical(Sys.info()[["sysname"]], "Darwin") &&
      length(intersect(all_missing, mac_source_packages)) > 0L
  ) {
    mac_config <- .cotra_prepare_macos_build()
    on.exit(.cotra_restore_macos_build(mac_config), add = TRUE)
  }

  if (length(missing_cran) > 0L) {
    if (
      ask &&
        !utils::askYesNo(
          paste(
            "Install missing CRAN packages?",
            paste(missing_cran, collapse = ", ")
          )
        )
    ) {
      stop("CRAN dependency installation cancelled.", call. = FALSE)
    }

    utils::install.packages(
      missing_cran,
      repos = cran_repo,
      dependencies = hard_dependencies
    )
  }

  # Ensure Bioconductor's repositories are available instead of replacing them
  # with a CRAN-only repository option.
  old_repos <- getOption("repos")
  on.exit(options(repos = old_repos), add = TRUE)

  bioc_repos <- suppressWarnings(BiocManager::repositories())
  bioc_repos["CRAN"] <- cran_repo
  options(repos = bioc_repos)

  missing_bioc <- bioc[!vapply(bioc, .cotra_loadable, logical(1))]

  if (length(missing_bioc) > 0L) {
    if (
      ask &&
        !utils::askYesNo(
          paste(
            "Install missing Bioconductor packages?",
            paste(missing_bioc, collapse = ", ")
          )
        )
    ) {
      stop("Bioconductor dependency installation cancelled.", call. = FALSE)
    }

    installed_names <- .cotra_installed_names()
    not_installed <- setdiff(missing_bioc, installed_names)

    if (length(not_installed) > 0L) {
      BiocManager::install(
        not_installed,
        ask = FALSE,
        update = update
      )
    }

    # A package can exist in the library but fail to load because one of its
    # dependencies is missing. Recheck after installing missing dependencies,
    # then force-reinstall only packages that remain broken.
    remaining_bioc <- bioc[!vapply(bioc, .cotra_loadable, logical(1))]
    installed_names <- .cotra_installed_names()
    broken_bioc <- intersect(remaining_bioc, installed_names)

    if (length(broken_bioc) > 0L) {
      message(
        "Reinstalling installed but unloadable Bioconductor packages: ",
        paste(broken_bioc, collapse = ", ")
      )
      BiocManager::install(
        broken_bioc,
        ask = FALSE,
        update = update,
        force = TRUE
      )
    }
  }

  missing_github <- names(github)[
    !vapply(names(github), .cotra_loadable, logical(1))
  ]

  if (length(missing_github) > 0L) {
    if (
      ask &&
        !utils::askYesNo(
          paste(
            "Install missing GitHub packages?",
            paste(missing_github, collapse = ", ")
          )
        )
    ) {
      stop("GitHub dependency installation cancelled.", call. = FALSE)
    }

    # Order matters: BPCells must be available before monocle3.
    for (pkg in names(github)) {
      if (!pkg %in% missing_github) next

      message("Installing ", pkg, " from GitHub...")
      remotes::install_github(
        github[[pkg]],
        dependencies = NA,
        upgrade = "never"
      )
    }
  }

  result <- check_cotra_dependencies(quiet = TRUE)

  if (!result$ok) {
    extra <- ""
    if (identical(Sys.info()[["sysname"]], "Darwin")) {
      extra <- paste0(
        "\n\nFor native macOS installation, ensure Terminal has:\n",
        "  xcode-select --install\n",
        "  brew install hdf5 pkgconf openssl@3\n",
        "Then restart R/RStudio and retry."
      )
    }

    stop(
      paste0(
        "CoTRA dependency installation is incomplete. Missing or unloadable: ",
        paste(result$missing, collapse = ", "),
        extra
      ),
      call. = FALSE
    )
  }

  message("All CoTRA dependencies are installed and loadable.")
  invisible(result)
}
