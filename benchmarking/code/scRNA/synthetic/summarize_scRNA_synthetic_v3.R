#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

indir <- if (length(commandArgs(trailingOnly = TRUE)) >= 1L) {
  commandArgs(trailingOnly = TRUE)[1L]
} else {
  "benchmark_results/scrna_synthetic_v3"
}

if (!dir.exists(indir)) stop("Directory not found: ", indir, call. = FALSE)

overall_files <- list.files(indir, "^overall_.*\\.csv$", full.names = TRUE)
if (!length(overall_files)) stop("No overall CSV files found.", call. = FALSE)

read_rss_mb <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  x <- readLines(path, warn = FALSE)
  hit <- grep("Maximum resident set size \\(kbytes\\)", x, value = TRUE)
  if (!length(hit)) return(NA_real_)
  kb <- suppressWarnings(as.numeric(trimws(sub(".*:", "", hit[[1L]]))))
  if (is.na(kb)) NA_real_ else kb / 1024
}

q <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))

runs <- do.call(rbind, lapply(overall_files, function(p) {
  d <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  d$gnu_peak_rss_mb <- vapply(
    d$benchmark_id,
    function(id) read_rss_mb(file.path(indir, paste0("resources_", id, ".txt"))),
    numeric(1)
  )
  d
}))

runs <- runs[order(runs$n_cells, runs$replicate), , drop = FALSE]
utils::write.csv(runs, file.path(indir, "all_scrna_synthetic_v3_runs.csv"),
                 row.names = FALSE, na = "")

groups <- split(runs, runs$n_cells)

summary_rows <- lapply(groups, function(g) {
  ok <- g[g$status == "success", , drop = FALSE]
  data.frame(
    n_cells = ok$n_cells[1L],
    runs = nrow(g),
    successful_runs = nrow(ok),
    core_median_s = median(ok$core_workflow_elapsed_seconds, na.rm = TRUE),
    core_q1_s = q(ok$core_workflow_elapsed_seconds, .25),
    core_q3_s = q(ok$core_workflow_elapsed_seconds, .75),
    marker_median_s = median(ok$marker_truth_elapsed_seconds, na.rm = TRUE),
    marker_q1_s = q(ok$marker_truth_elapsed_seconds, .25),
    marker_q3_s = q(ok$marker_truth_elapsed_seconds, .75),
    total_median_s = median(ok$total_core_plus_marker_seconds, na.rm = TRUE),
    total_q1_s = q(ok$total_core_plus_marker_seconds, .25),
    total_q3_s = q(ok$total_core_plus_marker_seconds, .75),
    peak_ram_median_gb = median(ok$gnu_peak_rss_mb, na.rm = TRUE) / 1024,
    peak_ram_q1_gb = q(ok$gnu_peak_rss_mb, .25) / 1024,
    peak_ram_q3_gb = q(ok$gnu_peak_rss_mb, .75) / 1024,
    inferred_louvain_clusters = median(ok$inferred_louvain_clusters, na.rm = TRUE),
    truth_marker_groups = median(ok$marker_truth_groups, na.rm = TRUE),
    marker_rows = median(ok$marker_rows, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})

summary_df <- do.call(rbind, summary_rows)
summary_df <- summary_df[order(summary_df$n_cells), , drop = FALSE]
utils::write.csv(summary_df,
                 file.path(indir, "summary_scrna_synthetic_v3.csv"),
                 row.names = FALSE, na = "")

step_files <- list.files(indir, "^steps_.*\\.csv$", full.names = TRUE)
parts <- list()

for (p in step_files) {
  id <- sub("^steps_(.*)\\.csv$", "\\1", basename(p))
  meta <- runs[runs$benchmark_id == id & runs$status == "success", , drop = FALSE]
  if (nrow(meta) != 1L) next
  d <- utils::read.csv(p, stringsAsFactors = FALSE)
  d$n_cells <- meta$n_cells[1L]
  d$replicate <- meta$replicate[1L]
  parts[[length(parts) + 1L]] <- d
}

if (length(parts)) {
  st <- do.call(rbind, parts)
  key <- interaction(st$n_cells, st$category, st$step, drop = TRUE, lex.order = TRUE)
  gs <- split(st, key)
  ss <- do.call(rbind, lapply(gs, function(g) {
    data.frame(
      n_cells = g$n_cells[1L],
      category = g$category[1L],
      step = g$step[1L],
      runs = nrow(g),
      median_seconds = median(g$elapsed_seconds, na.rm = TRUE),
      q1_seconds = q(g$elapsed_seconds, .25),
      q3_seconds = q(g$elapsed_seconds, .75),
      stringsAsFactors = FALSE
    )
  }))
  ss <- ss[order(ss$n_cells, ss$category, ss$step), , drop = FALSE]
  utils::write.csv(ss,
                   file.path(indir, "summary_scrna_synthetic_v3_steps.csv"),
                   row.names = FALSE, na = "")
}

print(summary_df)
