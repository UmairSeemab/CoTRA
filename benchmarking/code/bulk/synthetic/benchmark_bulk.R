#!/usr/bin/env Rscript

# CoTRA bulk RNA-seq computational performance benchmark
#
# This script benchmarks the computational backend used by the current CoTRA
# bulk modules without measuring browser/Shiny rendering time.
#
# Current CoTRA logic reproduced here:
#   inst/app/modules/bulk/qc_visualization.R
#     - log2(count + 1)
#     - top 2,000 variable genes
#     - prcomp(center = TRUE, scale. = TRUE)
#     - Pearson sample correlation
#     - Euclidean hierarchical clustering
#   inst/app/modules/bulk/de_analysis.R
#     - DESeq2: genes with counts >=10 in >=2 samples, DESeq(), results()
#     - edgeR: filterByExpr(), calcNormFactors(), estimateDisp(), glmQLFit(),
#       glmQLFTest()
#
# Scalability data are generated before the timed workflow so simulation time
# does not contribute to the reported CoTRA-equivalent workflow runtime.
# Plot rendering and file export are deliberately excluded.
#
# Example:
# Rscript --vanilla benchmark_bulk.R \
#   --genes 20000 --samples 24 --method DESeq2 --rep 1 \
#   --seed 1234 --outdir benchmark_results/bulk

options(stringsAsFactors = FALSE)

parse_cli <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) {
      stop("Unexpected argument: ", token, call. = FALSE)
    }
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

arg_value <- function(x, key, default = NULL) {
  if (is.null(x[[key]])) default else x[[key]]
}

as_int <- function(x, name) {
  value <- suppressWarnings(as.integer(x))
  if (length(value) != 1L || is.na(value)) stop(name, " must be an integer.", call. = FALSE)
  value
}

as_num <- function(x, name) {
  value <- suppressWarnings(as.numeric(x))
  if (length(value) != 1L || is.na(value)) stop(name, " must be numeric.", call. = FALSE)
  value
}

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg, call. = FALSE)
  }
}

proc_status_kb <- function(field) {
  path <- "/proc/self/status"
  if (!file.exists(path)) return(NA_real_)
  txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  hit <- grep(paste0("^", field, ":"), txt, value = TRUE)
  if (length(hit) == 0L) return(NA_real_)
  suppressWarnings(as.numeric(sub(".*?:\\s*([0-9]+).*", "\\1", hit[[1L]])))
}

safe_pkg_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}

safe_git_commit <- function() {
  explicit <- Sys.getenv("COTRA_COMMIT", unset = "")
  if (nzchar(explicit)) return(explicit)

  repo <- Sys.getenv("COTRA_REPO", unset = "")
  if (!nzchar(repo) || !dir.exists(repo)) return(NA_character_)
  ans <- tryCatch(
    system2("git", c("-C", normalizePath(repo), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  if (length(ans) == 0L) NA_character_ else trimws(ans[[1L]])
}

write_session_info <- function(path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(capture.output(sessionInfo()), con)
}

measure_step <- function(step_name, expr, step_log_env) {
  invisible(gc())
  start <- proc.time()[["elapsed"]]
  result <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - start
  step_log_env$rows[[length(step_log_env$rows) + 1L]] <- data.frame(
    step = step_name,
    elapsed_seconds = as.numeric(elapsed),
    stringsAsFactors = FALSE
  )
  result
}

simulate_counts <- function(n_genes, n_samples, seed, de_fraction = 0.10) {
  if (n_samples < 4L || n_samples %% 2L != 0L) {
    stop("--samples must be an even integer >= 4 for a balanced two-group benchmark.", call. = FALSE)
  }
  if (n_genes < 100L) stop("--genes must be >= 100.", call. = FALSE)

  set.seed(seed)
  group <- factor(
    rep(c("control", "treatment"), each = n_samples / 2L),
    levels = c("control", "treatment")
  )

  # Heterogeneous expression means and dispersions produce a realistic range
  # of bulk RNA-seq count magnitudes while keeping generation reproducible.
  base_mean <- exp(stats::rnorm(n_genes, mean = log(80), sd = 1.15))
  dispersion <- pmax(0.02, pmin(1.5, exp(stats::rnorm(n_genes, mean = log(0.15), sd = 0.60))))
  library_factor <- exp(stats::rnorm(n_samples, mean = 0, sd = 0.18))

  n_de <- max(1L, round(n_genes * de_fraction))
  de_idx <- seq_len(n_de)
  lfc <- numeric(n_genes)
  lfc[de_idx] <- sample(c(-1.5, 1.5), n_de, replace = TRUE)

  mu <- outer(base_mean, library_factor)
  treatment_cols <- which(group == "treatment")
  mu[, treatment_cols] <- mu[, treatment_cols, drop = FALSE] * (2 ^ lfc)

  count_vector <- stats::rnbinom(
    n = n_genes * n_samples,
    mu = as.vector(mu),
    size = rep(1 / dispersion, times = n_samples)
  )

  counts <- matrix(count_vector, nrow = n_genes, ncol = n_samples)
  rownames(counts) <- sprintf("GENE_%05d", seq_len(n_genes))
  colnames(counts) <- sprintf("sample_%03d", seq_len(n_samples))
  storage.mode(counts) <- "integer"

  list(counts = counts, group = group)
}

run_qc <- function(counts, step_log) {
  top <- measure_step("qc_top_variable_genes", {
    log_mat <- log2(counts + 1)
    gene_variance <- apply(log_mat, 1L, stats::var, na.rm = TRUE)
    valid <- which(is.finite(gene_variance) & gene_variance > 0)
    if (length(valid) < 2L) stop("Too few genes with non-zero variance for QC.")
    ordered <- valid[order(gene_variance[valid], decreasing = TRUE)]
    keep <- ordered[seq_len(min(2000L, length(ordered)))]
    log_mat[keep, , drop = FALSE]
  }, step_log)

  pca <- measure_step("qc_pca", {
    stats::prcomp(t(top), center = TRUE, scale. = TRUE)
  }, step_log)

  cor_mat <- measure_step("qc_sample_correlation", {
    stats::cor(top, method = "pearson", use = "pairwise.complete.obs")
  }, step_log)

  hc <- measure_step("qc_hierarchical_clustering", {
    stats::hclust(stats::dist(t(top), method = "euclidean"))
  }, step_log)

  invisible(list(pca = pca, cor = cor_mat, hclust = hc))
}

run_deseq2 <- function(counts, group, step_log) {
  require_pkg("DESeq2")

  measure_step("differential_expression_DESeq2", {
    coldata <- data.frame(group = group, row.names = colnames(counts))
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(as.matrix(counts)),
      colData = coldata,
      design = ~ group
    )

    keep <- rowSums(DESeq2::counts(dds) >= 10) >= 2
    if (sum(keep) < 2L) stop("Too few genes passed the current CoTRA DESeq2 count filter.")

    dds <- dds[keep, ]
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    res <- DESeq2::results(dds, contrast = c("group", "treatment", "control"))
    df <- as.data.frame(res)
    df$gene <- rownames(df)
    df$lfc <- df$log2FoldChange
    df$padj <- df$padj
    df
  }, step_log)
}

run_edger <- function(counts, group, step_log) {
  require_pkg("edgeR")

  measure_step("differential_expression_edgeR", {
    y <- edgeR::DGEList(counts = round(as.matrix(counts)), group = group)
    keep <- edgeR::filterByExpr(y, group = group)
    if (sum(keep) < 2L) stop("Too few genes passed the current CoTRA edgeR filter.")

    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)

    design <- stats::model.matrix(~ group, data = data.frame(group = group))
    y <- edgeR::estimateDisp(y, design)
    fit <- edgeR::glmQLFit(y, design)
    test <- edgeR::glmQLFTest(fit, coef = 2)
    tab <- edgeR::topTags(test, n = Inf)$table

    df <- as.data.frame(tab)
    df$gene <- rownames(df)
    df$lfc <- df$logFC
    df$padj <- df$FDR
    df
  }, step_log)
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))

n_genes <- as_int(arg_value(args, "genes", "20000"), "--genes")
n_samples <- as_int(arg_value(args, "samples", "24"), "--samples")
method <- arg_value(args, "method", "DESeq2")
replicate_id <- as_int(arg_value(args, "rep", "1"), "--rep")
seed <- as_int(arg_value(args, "seed", "1234"), "--seed")
outdir <- arg_value(args, "outdir", file.path("benchmark_results", "bulk"))
padj_cutoff <- as_num(arg_value(args, "padj", "0.05"), "--padj")
lfc_cutoff <- as_num(arg_value(args, "lfc", "1"), "--lfc")

if (!method %in% c("DESeq2", "edgeR")) {
  stop("--method must be DESeq2 or edgeR.", call. = FALSE)
}
if (replicate_id < 0L) stop("--rep must be >= 0.", call. = FALSE)
if (padj_cutoff < 0 || padj_cutoff > 1) stop("--padj must be between 0 and 1.", call. = FALSE)
if (lfc_cutoff < 0) stop("--lfc must be >= 0.", call. = FALSE)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Keep the simulated dataset identical across methods and measured repetitions.
# The replicate ID is only a run label; it must not change the benchmark input.
effective_seed <- seed
sim <- simulate_counts(n_genes, n_samples, effective_seed)
counts <- sim$counts
group <- sim$group
rm(sim)
invisible(gc())

benchmark_id <- sprintf(
  "bulk_%s_%05dg_%03ds_rep%02d",
  method, n_genes, n_samples, replicate_id
)

message("Benchmark: ", benchmark_id)
message("Dimensions: ", nrow(counts), " genes x ", ncol(counts), " samples")
message("Method: ", method)

step_log <- new.env(parent = emptyenv())
step_log$rows <- list()

workflow_start <- proc.time()[["elapsed"]]
status <- "success"
error_message <- NA_character_
de_result <- NULL

tryCatch({
  run_qc(counts, step_log)
  de_result <- if (identical(method, "DESeq2")) {
    run_deseq2(counts, group, step_log)
  } else {
    run_edger(counts, group, step_log)
  }
}, error = function(e) {
  status <<- "failed"
  error_message <<- conditionMessage(e)
})

workflow_elapsed <- proc.time()[["elapsed"]] - workflow_start

steps <- if (length(step_log$rows)) {
  do.call(rbind, step_log$rows)
} else {
  data.frame(step = character(), elapsed_seconds = numeric())
}

n_tested <- if (!is.null(de_result)) nrow(de_result) else NA_integer_
n_significant <- if (!is.null(de_result)) {
  sum(
    is.finite(de_result$padj) &
      de_result$padj <= padj_cutoff &
      is.finite(de_result$lfc) &
      abs(de_result$lfc) >= lfc_cutoff,
    na.rm = TRUE
  )
} else {
  NA_integer_
}

peak_rss_kb <- proc_status_kb("VmHWM")
current_rss_kb <- proc_status_kb("VmRSS")

summary_df <- data.frame(
  benchmark_id = benchmark_id,
  workflow = "bulk",
  backend_scope = "CoTRA-equivalent computational backend; Shiny/browser rendering excluded",
  method = method,
  replicate = replicate_id,
  seed = effective_seed,
  n_genes_input = n_genes,
  n_samples = n_samples,
  group_1_n = sum(group == "control"),
  group_2_n = sum(group == "treatment"),
  qc_top_variable_genes = min(2000L, n_genes),
  padj_cutoff = padj_cutoff,
  abs_log2fc_cutoff = lfc_cutoff,
  workflow_elapsed_seconds = as.numeric(workflow_elapsed),
  n_genes_tested = n_tested,
  n_significant = n_significant,
  peak_rss_mb_process = peak_rss_kb / 1024,
  current_rss_mb_end = current_rss_kb / 1024,
  status = status,
  error_message = error_message,
  R_version = paste(R.version$major, R.version$minor, sep = "."),
  CoTRA_version = safe_pkg_version("CoTRA"),
  DESeq2_version = safe_pkg_version("DESeq2"),
  edgeR_version = safe_pkg_version("edgeR"),
  cotra_git_commit = safe_git_commit(),
  hostname = Sys.info()[["nodename"]],
  timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)

utils::write.csv(
  summary_df,
  file.path(outdir, paste0("overall_", benchmark_id, ".csv")),
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  steps,
  file.path(outdir, paste0("steps_", benchmark_id, ".csv")),
  row.names = FALSE,
  na = ""
)

write_session_info(file.path(outdir, paste0("session_", benchmark_id, ".txt")))

if (identical(status, "failed")) {
  message("FAILED: ", error_message)
  quit(save = "no", status = 1L)
}

message(sprintf("Completed in %.3f s", workflow_elapsed))
message(sprintf("Process peak RSS reported by /proc: %.1f MB", peak_rss_kb / 1024))
