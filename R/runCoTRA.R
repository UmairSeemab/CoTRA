#' Launch CoTRA
#'
#' Starts the CoTRA Shiny application from an installed package.
#'
#' @param app_dir Optional path to a writable app copy. If NULL, CoTRA creates a temporary writable copy.
#' @param install_missing Logical. Install missing R package dependencies before launch.
#' @param ask Logical. Ask before installing missing packages in interactive sessions.
#' @param launch.browser Passed to shiny::runApp().
#' @export
runCoTRA <- function(app_dir = NULL,
                     install_missing = TRUE,
                     ask = interactive(),
                     launch.browser = TRUE) {
  if (install_missing) {
    tryCatch(
      install_cotra_dependencies(ask = ask),
      error = function(e) {
        stop(
          paste0(
            "CoTRA dependency setup failed before the app could start.\n\n",
            conditionMessage(e),
            "\n\nYou can retry directly with:\n",
            "  CoTRA::install_cotra_dependencies(ask = FALSE)\n\n",
            "Container installation is recommended when native system-library ",
            "compilation is not desired."
          ),
          call. = FALSE
        )
      }
    )
  } else {
    deps <- check_cotra_dependencies(quiet = TRUE)
    if (length(deps$missing) > 0L) {
      stop(
        "Missing or unloadable CoTRA dependencies: ",
        paste(deps$missing, collapse = ", "),
        "\nRun CoTRA::install_cotra_dependencies(ask = FALSE) first.",
        call. = FALSE
      )
    }
  }

  pkg_app <- system.file("app", package = "CoTRA", mustWork = TRUE)

  if (is.null(app_dir)) {
    app_dir <- file.path(
      tempdir(),
      paste0("CoTRA_app_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    )
  }

  if (dir.exists(app_dir)) {
    unlink(app_dir, recursive = TRUE, force = TRUE)
  }

  ok <- file.copy(pkg_app, dirname(app_dir), recursive = TRUE)
  copied <- file.path(dirname(app_dir), basename(pkg_app))
  if (!ok || !dir.exists(copied)) {
    stop("Could not create a writable CoTRA app copy.", call. = FALSE)
  }

  if (!file.rename(copied, app_dir)) {
    stop("Could not prepare the writable CoTRA app copy.", call. = FALSE)
  }

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(app_dir)

  shiny::runApp(app_dir, launch.browser = launch.browser)
}
