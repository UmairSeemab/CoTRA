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
    install_cotra_dependencies(ask = ask)
  } else {
    deps <- check_cotra_dependencies(quiet = TRUE)
    if (length(deps$missing) > 0) {
      stop(
        "Missing CoTRA dependencies: ", paste(deps$missing, collapse = ", "),
        "\nRun CoTRA::install_cotra_dependencies() first.",
        call. = FALSE
      )
    }
  }

  pkg_app <- system.file("app", package = "CoTRA", mustWork = TRUE)

  if (is.null(app_dir)) {
    app_dir <- file.path(tempdir(), paste0("CoTRA_app_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  }

  if (dir.exists(app_dir)) {
    unlink(app_dir, recursive = TRUE, force = TRUE)
  }

  ok <- file.copy(pkg_app, dirname(app_dir), recursive = TRUE)
  copied <- file.path(dirname(app_dir), basename(pkg_app))
  if (!ok || !dir.exists(copied)) {
    stop("Could not create a writable CoTRA app copy.", call. = FALSE)
  }
  file.rename(copied, app_dir)

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(app_dir)

  shiny::runApp(app_dir, launch.browser = launch.browser)
}
