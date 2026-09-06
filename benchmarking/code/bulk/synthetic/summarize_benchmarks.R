#!/usr/bin/env Rscript

# Summarize CoTRA benchmark outputs into publication-ready CSV tables.
# Primary runtime = workflow_elapsed_seconds from the R benchmark scripts.
# Primary peak RAM = GNU /usr/bin/time Maximum resident set size when available,
# otherwise the Linux /proc VmHWM value recorded inside R.
#
# Example:
# Rscript --vanilla summarize_benchmarks.R --indir benchmark_results

options(stringsAsFactors = FALSE)

parse_cli <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token, call. = FALSE)
    key <- sub("^--", "", token)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

arg_value <- function(x, key, default = NULL) if (is.null(x[[key]])) default else x[[key]]

read_csv_set <- function(paths) {
  if (length(paths) == 0L) return(NULL)
  x <- lapply(paths, function(p) {
    d <- utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
    d$source_file <- p
    d
  })
  do.call(rbind, x)
}

parse_gnu_time_maxrss_mb <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  x <- readLines(path, warn = FALSE)
  hit <- grep("Maximum resident set size \\(kbytes\\)", x, value = TRUE)
  if (length(hit) == 0L) return(NA_real_)
  kb <- suppressWarnings(as.numeric(trimws(sub(".*:", "", hit[[1L]]))))
  if (is.na(kb)) NA_real_ else kb / 1024
}

attach_gnu_time <- function(df, directory) {
  if (is.null(df) || nrow(df) == 0L) return(df)
  df$gnu_time_peak_rss_mb <- vapply(df$benchmark_id, function(id) {
    parse_gnu_time_maxrss_mb(file.path(directory, paste0("resources_", id, ".txt")))
  }, numeric(1))
  df$peak_rss_mb_for_summary <- ifelse(
    is.finite(df$gnu_time_peak_rss_mb),
    df$gnu_time_peak_rss_mb,
    df$peak_rss_mb_process
  )
  df
}

summarize_groups <- function(df, group_cols) {
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  key <- interaction(df[group_cols], drop = TRUE, lex.order = TRUE)
  groups <- split(df, key)

  rows <- lapply(groups, function(g) {
    ok <- g[g$status == "success", , drop = FALSE]
    out <- g[1L, group_cols, drop = FALSE]
    out$measured_runs <- nrow(g)
    out$successful_runs <- nrow(ok)
    out$failed_runs <- nrow(g) - nrow(ok)

    vals <- ok$workflow_elapsed_seconds
    out$median_runtime_seconds <- if (length(vals)) stats::median(vals, na.rm = TRUE) else NA_real_
    out$runtime_Q1_seconds <- if (length(vals)) as.numeric(stats::quantile(vals, 0.25, na.rm = TRUE, names = FALSE)) else NA_real_
    out$runtime_Q3_seconds <- if (length(vals)) as.numeric(stats::quantile(vals, 0.75, na.rm = TRUE, names = FALSE)) else NA_real_
    out$runtime_IQR_seconds <- if (length(vals)) stats::IQR(vals, na.rm = TRUE) else NA_real_

    mem <- ok$peak_rss_mb_for_summary
    out$median_peak_rss_mb <- if (length(mem)) stats::median(mem, na.rm = TRUE) else NA_real_
    out$peak_rss_Q1_mb <- if (length(mem)) as.numeric(stats::quantile(mem, 0.25, na.rm = TRUE, names = FALSE)) else NA_real_
    out$peak_rss_Q3_mb <- if (length(mem)) as.numeric(stats::quantile(mem, 0.75, na.rm = TRUE, names = FALSE)) else NA_real_
    out$median_peak_rss_gb <- out$median_peak_rss_mb / 1024
    out
  })

  do.call(rbind, rows)
}

summarize_steps <- function(overall, step_dir, workflow, group_cols) {
  if (is.null(overall) || nrow(overall) == 0L) return(NULL)

  step_files <- list.files(step_dir, pattern = "^steps_.*\\.csv$", full.names = TRUE)
  if (length(step_files) == 0L) return(NULL)

  parts <- lapply(step_files, function(p) {
    id <- sub("^steps_(.*)\\.csv$", "\\1", basename(p))
    meta <- overall[overall$benchmark_id == id, , drop = FALSE]
    if (nrow(meta) != 1L || meta$status != "success") return(NULL)
    d <- utils::read.csv(p, stringsAsFactors = FALSE)
    if (nrow(d) == 0L) return(NULL)
    for (col in group_cols) d[[col]] <- meta[[col]][1L]
    d$replicate <- meta$replicate[1L]
    d
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) return(NULL)
  steps <- do.call(rbind, parts)

  key_cols <- c(group_cols, "step")
  key <- interaction(steps[key_cols], drop = TRUE, lex.order = TRUE)
  groups <- split(steps, key)
  do.call(rbind, lapply(groups, function(g) {
    out <- g[1L, key_cols, drop = FALSE]
    out$runs <- nrow(g)
    out$median_seconds <- stats::median(g$elapsed_seconds, na.rm = TRUE)
    out$Q1_seconds <- as.numeric(stats::quantile(g$elapsed_seconds, 0.25, na.rm = TRUE, names = FALSE))
    out$Q3_seconds <- as.numeric(stats::quantile(g$elapsed_seconds, 0.75, na.rm = TRUE, names = FALSE))
    out$IQR_seconds <- stats::IQR(g$elapsed_seconds, na.rm = TRUE)
    out
  }))
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
indir <- arg_value(args, "indir", "benchmark_results")
outdir <- arg_value(args, "outdir", file.path(indir, "summary"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bulk_dir <- file.path(indir, "bulk")
sc_dir <- file.path(indir, "scrna")

bulk_paths <- if (dir.exists(bulk_dir)) list.files(bulk_dir, pattern = "^overall_.*\\.csv$", full.names = TRUE) else character()
sc_paths <- if (dir.exists(sc_dir)) list.files(sc_dir, pattern = "^overall_.*\\.csv$", full.names = TRUE) else character()

bulk <- attach_gnu_time(read_csv_set(bulk_paths), bulk_dir)
sc <- attach_gnu_time(read_csv_set(sc_paths), sc_dir)

if (!is.null(bulk)) {
  utils::write.csv(bulk, file.path(outdir, "all_bulk_runs.csv"), row.names = FALSE, na = "")
  bulk_summary <- summarize_groups(bulk, c("method", "n_genes_input", "n_samples"))
  utils::write.csv(bulk_summary, file.path(outdir, "summary_bulk.csv"), row.names = FALSE, na = "")

  bulk_steps <- summarize_steps(bulk, bulk_dir, "bulk", c("method", "n_samples"))
  if (!is.null(bulk_steps)) {
    utils::write.csv(bulk_steps, file.path(outdir, "summary_bulk_steps.csv"), row.names = FALSE, na = "")
  }
}

if (!is.null(sc)) {
  utils::write.csv(sc, file.path(outdir, "all_scrna_runs.csv"), row.names = FALSE, na = "")
  sc_summary <- summarize_groups(sc, c("n_genes_input", "n_cells"))
  utils::write.csv(sc_summary, file.path(outdir, "summary_scrna.csv"), row.names = FALSE, na = "")

  sc_steps <- summarize_steps(sc, sc_dir, "scrna", c("n_cells"))
  if (!is.null(sc_steps)) {
    utils::write.csv(sc_steps, file.path(outdir, "summary_scrna_steps.csv"), row.names = FALSE, na = "")
  }
}

if (is.null(bulk) && is.null(sc)) {
  stop("No overall benchmark CSV files were found under: ", indir, call. = FALSE)
}

message("Summary written to: ", normalizePath(outdir, mustWork = FALSE))
