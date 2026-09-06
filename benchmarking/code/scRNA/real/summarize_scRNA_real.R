#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

indir <- if (length(commandArgs(trailingOnly = TRUE)) >= 1L) {
  commandArgs(trailingOnly = TRUE)[1L]
} else {
  "benchmark_results/scrna_real"
}

if (!dir.exists(indir)) stop("Directory not found: ", indir, call. = FALSE)

files <- list.files(indir, "^overall_.*\\.csv$", full.names = TRUE)
if (!length(files)) stop("No overall CSV files found.", call. = FALSE)

read_rss_mb <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  x <- readLines(path, warn = FALSE)
  hit <- grep("Maximum resident set size \\(kbytes\\)", x, value = TRUE)
  if (!length(hit)) return(NA_real_)
  kb <- suppressWarnings(as.numeric(trimws(sub(".*:", "", hit[[1L]]))))
  if (is.na(kb)) NA_real_ else kb / 1024
}

runs <- do.call(rbind, lapply(files, function(path) {
  d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  d$gnu_time_peak_rss_mb <- vapply(
    d$benchmark_id,
    function(id) read_rss_mb(file.path(indir, paste0("resources_", id, ".txt"))),
    numeric(1)
  )
  d$peak_rss_mb_for_summary <- ifelse(
    is.finite(d$gnu_time_peak_rss_mb),
    d$gnu_time_peak_rss_mb,
    d$peak_rss_mb_process
  )
  d
}))

runs <- runs[order(runs$n_cells, runs$replicate), , drop = FALSE]

utils::write.csv(
  runs,
  file.path(indir, "all_real_scrna_runs.csv"),
  row.names = FALSE,
  na = ""
)

groups <- split(runs, runs$n_cells)

summary_rows <- lapply(groups, function(g) {

  ok <- g[g$status == "success", , drop = FALSE]

  data.frame(
    n_genes = ok$n_genes_input[1L],
    n_cells = ok$n_cells[1L],
    n_WT = ok$n_WT[1L],
    n_rd10 = ok$n_rd10[1L],
    pcs_used = ok$pcs_used_downstream[1L],
    runs = nrow(g),
    successful_runs = nrow(ok),
    median_core_runtime_s = median(ok$core_workflow_elapsed_seconds, na.rm = TRUE),
    q1_core_runtime_s = as.numeric(quantile(ok$core_workflow_elapsed_seconds, 0.25, na.rm = TRUE, names = FALSE)),
    q3_core_runtime_s = as.numeric(quantile(ok$core_workflow_elapsed_seconds, 0.75, na.rm = TRUE, names = FALSE)),
    median_marker_runtime_s = median(ok$marker_elapsed_seconds, na.rm = TRUE),
    q1_marker_runtime_s = as.numeric(quantile(ok$marker_elapsed_seconds, 0.25, na.rm = TRUE, names = FALSE)),
    q3_marker_runtime_s = as.numeric(quantile(ok$marker_elapsed_seconds, 0.75, na.rm = TRUE, names = FALSE)),
    median_total_runtime_s = median(ok$total_core_plus_marker_seconds, na.rm = TRUE),
    q1_total_runtime_s = as.numeric(quantile(ok$total_core_plus_marker_seconds, 0.25, na.rm = TRUE, names = FALSE)),
    q3_total_runtime_s = as.numeric(quantile(ok$total_core_plus_marker_seconds, 0.75, na.rm = TRUE, names = FALSE)),
    median_peak_ram_gb = median(ok$peak_rss_mb_for_summary, na.rm = TRUE) / 1024,
    q1_peak_ram_gb = as.numeric(quantile(ok$peak_rss_mb_for_summary, 0.25, na.rm = TRUE, names = FALSE)) / 1024,
    q3_peak_ram_gb = as.numeric(quantile(ok$peak_rss_mb_for_summary, 0.75, na.rm = TRUE, names = FALSE)) / 1024,
    median_clusters = median(ok$inferred_louvain_clusters, na.rm = TRUE),
    median_marker_rows = median(ok$marker_rows, na.rm = TRUE),
    full_dataset_cluster_match = if (any(!is.na(ok$full_dataset_cluster_match))) {
      all(ok$full_dataset_cluster_match[!is.na(ok$full_dataset_cluster_match)])
    } else {
      NA
    },
    stringsAsFactors = FALSE
  )
})

summary_df <- do.call(rbind, summary_rows)
summary_df <- summary_df[order(summary_df$n_cells), , drop = FALSE]

utils::write.csv(
  summary_df,
  file.path(indir, "summary_real_scrna.csv"),
  row.names = FALSE,
  na = ""
)

# Step-level summary
step_files <- list.files(indir, "^steps_.*\\.csv$", full.names = TRUE)
parts <- list()

for (path in step_files) {
  id <- sub("^steps_(.*)\\.csv$", "\\1", basename(path))
  meta <- runs[runs$benchmark_id == id, , drop = FALSE]
  if (nrow(meta) != 1L || meta$status != "success") next

  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!nrow(d)) next

  d$n_cells <- meta$n_cells[1L]
  d$replicate <- meta$replicate[1L]
  parts[[length(parts) + 1L]] <- d
}

if (length(parts)) {
  steps <- do.call(rbind, parts)
  keys <- interaction(steps$n_cells, steps$category, steps$step, drop = TRUE, lex.order = TRUE)
  gs <- split(steps, keys)

  step_summary <- do.call(rbind, lapply(gs, function(g) {
    data.frame(
      n_cells = g$n_cells[1L],
      category = g$category[1L],
      step = g$step[1L],
      runs = nrow(g),
      median_seconds = median(g$elapsed_seconds, na.rm = TRUE),
      Q1_seconds = as.numeric(quantile(g$elapsed_seconds, 0.25, na.rm = TRUE, names = FALSE)),
      Q3_seconds = as.numeric(quantile(g$elapsed_seconds, 0.75, na.rm = TRUE, names = FALSE)),
      stringsAsFactors = FALSE
    )
  }))

  step_summary <- step_summary[order(step_summary$n_cells, step_summary$category, step_summary$step), , drop = FALSE]

  utils::write.csv(
    step_summary,
    file.path(indir, "summary_real_scrna_steps.csv"),
    row.names = FALSE,
    na = ""
  )
}

print(summary_df)
message("Summary written to: ", normalizePath(indir, mustWork = FALSE))
