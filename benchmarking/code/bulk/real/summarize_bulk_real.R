#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
indir <- if (length(args) >= 1L) args[[1L]] else "benchmark_results/bulk_real"

overall_files <- list.files(indir, pattern = "^overall_.*\\.csv$", full.names = TRUE)
if (!length(overall_files)) stop("No overall CSV files found in: ", indir)

overall <- do.call(rbind, lapply(overall_files, utils::read.csv, check.names = FALSE))

read_maxrss <- function(method, rep) {
  f <- file.path(indir, sprintf("resources_bulk_real_%s_rep%02d.txt", method, rep))
  if (!file.exists(f)) return(NA_real_)
  x <- readLines(f, warn = FALSE)
  hit <- grep("Maximum resident set size \\(kbytes\\):", x, value = TRUE)
  if (!length(hit)) return(NA_real_)
  kb <- suppressWarnings(as.numeric(sub(".*:\\s*", "", hit[[1L]])))
  kb / 1024 / 1024
}

overall$peak_ram_gb_gnutime <- mapply(
  read_maxrss,
  overall$method,
  overall$replicate
)

q1 <- function(x) stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE)
q3 <- function(x) stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE)

groups <- split(overall, overall$method)

summary <- do.call(rbind, lapply(groups, function(d) {
  data.frame(
    dataset = d$dataset[[1L]],
    contrast = d$contrast[[1L]],
    method = d$method[[1L]],
    n_genes_input = d$n_genes_input[[1L]],
    n_samples = d$n_samples[[1L]],
    reference_n = d$reference_n[[1L]],
    comparison_n = d$comparison_n[[1L]],
    runs = nrow(d),
    successful_runs = sum(d$status == "success"),
    median_runtime_s = stats::median(d$workflow_elapsed_seconds, na.rm = TRUE),
    q1_runtime_s = q1(d$workflow_elapsed_seconds),
    q3_runtime_s = q3(d$workflow_elapsed_seconds),
    median_peak_ram_gb = stats::median(d$peak_ram_gb_gnutime, na.rm = TRUE),
    min_peak_ram_gb = min(d$peak_ram_gb_gnutime, na.rm = TRUE),
    max_peak_ram_gb = max(d$peak_ram_gb_gnutime, na.rm = TRUE),
    median_genes_tested = stats::median(d$n_genes_tested, na.rm = TRUE),
    median_significant = stats::median(d$n_significant, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

summary <- summary[order(summary$method), ]
rownames(summary) <- NULL

utils::write.csv(
  overall,
  file.path(indir, "all_real_bulk_runs.csv"),
  row.names = FALSE
)
utils::write.csv(
  summary,
  file.path(indir, "summary_real_bulk.csv"),
  row.names = FALSE
)

print(summary, row.names = FALSE)
