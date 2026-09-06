#!/usr/bin/env Rscript

# CoTRA real-data bulk RNA-seq benchmark
# Benchmarks WT vs rd10 using the same computational logic as the current
# CoTRA bulk QC and DE modules. Shiny/browser rendering is excluded.

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

read_count_matrix <- function(path) {
  if (!file.exists(path)) stop("Input file does not exist: ", path, call. = FALSE)

  ext <- tolower(tools::file_ext(path))
  sep <- if (ext %in% c("tsv", "txt", "tab")) "\t" else ","

  dat <- utils::read.table(
    path,
    header = TRUE,
    sep = sep,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (ncol(dat) < 3L) stop("Expected a gene column plus at least two sample columns.", call. = FALSE)

  gene_ids <- dat[[1L]]
  if (anyNA(gene_ids) || any(!nzchar(as.character(gene_ids)))) {
    stop("Gene identifier column contains missing/empty values.", call. = FALSE)
  }
  if (anyDuplicated(gene_ids)) {
    stop("Gene identifiers are duplicated. Resolve duplicates before benchmarking.", call. = FALSE)
  }

  count_df <- dat[-1L]
  for (j in seq_along(count_df)) {
    count_df[[j]] <- suppressWarnings(as.numeric(count_df[[j]]))
  }

  counts <- as.matrix(count_df)
  rownames(counts) <- as.character(gene_ids)

  if (anyNA(counts)) stop("Count matrix contains missing or non-numeric values.", call. = FALSE)
  if (any(counts < 0)) stop("Count matrix contains negative values.", call. = FALSE)
  if (any(abs(counts - round(counts)) > .Machine$double.eps^0.5)) {
    stop("DESeq2/edgeR benchmark requires raw integer counts.", call. = FALSE)
  }

  storage.mode(counts) <- "integer"
  counts
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

run_deseq2 <- function(counts, group, reference, comparison, step_log) {
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

    res <- DESeq2::results(
      dds,
      contrast = c("group", comparison, reference)
    )

    df <- as.data.frame(res)
    df$gene <- rownames(df)
    df$lfc <- df$log2FoldChange
    df$padj <- df$padj
    df
  }, step_log)
}

run_edger <- function(counts, group, reference, comparison, step_log) {
  require_pkg("edgeR")

  measure_step("differential_expression_edgeR", {
    y <- edgeR::DGEList(counts = round(as.matrix(counts)), group = group)

    keep <- edgeR::filterByExpr(y, group = group)
    if (sum(keep) < 2L) stop("Too few genes passed the current CoTRA edgeR filter.")

    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)

    design <- stats::model.matrix(~ group)
    y <- edgeR::estimateDisp(y, design)
    fit <- edgeR::glmQLFit(y, design)

    coef_name <- paste0("group", comparison)
    coef_index <- match(coef_name, colnames(design))
    if (is.na(coef_index)) {
      stop("Could not identify edgeR coefficient for comparison: ", comparison)
    }

    test <- edgeR::glmQLFTest(fit, coef = coef_index)
    tab <- edgeR::topTags(test, n = Inf)$table

    df <- as.data.frame(tab)
    df$gene <- rownames(df)
    df$lfc <- df$logFC
    df$padj <- df$FDR
    df
  }, step_log)
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))

input <- arg_value(args, "input", NULL)
if (is.null(input)) stop("--input is required.", call. = FALSE)

method <- arg_value(args, "method", "DESeq2")
replicate_id <- as_int(arg_value(args, "rep", "1"), "--rep")
seed <- as_int(arg_value(args, "seed", "1234"), "--seed")
outdir <- arg_value(args, "outdir", file.path("benchmark_results", "bulk_real"))
reference <- arg_value(args, "reference", "WT")
comparison <- arg_value(args, "comparison", "rd10")
reference_regex <- arg_value(args, "reference-regex", "^WT")
comparison_regex <- arg_value(args, "comparison-regex", "^rd10")
padj_cutoff <- as_num(arg_value(args, "padj", "0.05"), "--padj")
lfc_cutoff <- as_num(arg_value(args, "lfc", "1"), "--lfc")

if (!method %in% c("DESeq2", "edgeR")) {
  stop("--method must be DESeq2 or edgeR.", call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Input loading is deliberately outside the timed workflow.
counts <- read_count_matrix(input)

ref_cols <- grep(reference_regex, colnames(counts), ignore.case = TRUE, value = TRUE)
cmp_cols <- grep(comparison_regex, colnames(counts), ignore.case = TRUE, value = TRUE)

if (length(ref_cols) == 0L) stop("No reference samples matched: ", reference_regex, call. = FALSE)
if (length(cmp_cols) == 0L) stop("No comparison samples matched: ", comparison_regex, call. = FALSE)
if (length(intersect(ref_cols, cmp_cols)) > 0L) stop("Reference and comparison regexes overlap.", call. = FALSE)

selected <- c(ref_cols, cmp_cols)
counts <- counts[, selected, drop = FALSE]

group <- factor(
  c(rep(reference, length(ref_cols)), rep(comparison, length(cmp_cols))),
  levels = c(reference, comparison)
)
names(group) <- selected

set.seed(seed)

benchmark_id <- sprintf(
  "bulk_real_%s_%s_vs_%s_%05dg_%03ds_rep%02d",
  method, comparison, reference, nrow(counts), ncol(counts), replicate_id
)

message("Benchmark: ", benchmark_id)
message("Input: ", normalizePath(input))
message("Dimensions: ", nrow(counts), " genes x ", ncol(counts), " samples")
message("Reference: ", reference, " (", length(ref_cols), ")")
message("Comparison: ", comparison, " (", length(cmp_cols), ")")
message("Contrast: ", comparison, " vs ", reference)
message("Positive log2FC = higher in ", comparison)
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
    run_deseq2(counts, group, reference, comparison, step_log)
  } else {
    run_edger(counts, group, reference, comparison, step_log)
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

summary_df <- data.frame(
  benchmark_id = benchmark_id,
  workflow = "bulk_real",
  dataset = basename(input),
  contrast = paste(comparison, "vs", reference),
  method = method,
  replicate = replicate_id,
  seed = seed,
  n_genes_input = nrow(counts),
  n_samples = ncol(counts),
  reference = reference,
  reference_n = length(ref_cols),
  comparison = comparison,
  comparison_n = length(cmp_cols),
  padj_cutoff = padj_cutoff,
  abs_log2fc_cutoff = lfc_cutoff,
  workflow_elapsed_seconds = as.numeric(workflow_elapsed),
  n_genes_tested = n_tested,
  n_significant = n_significant,
  peak_rss_mb_process = peak_rss_kb / 1024,
  status = status,
  error_message = error_message,
  R_version = paste(R.version$major, R.version$minor, sep = "."),
  DESeq2_version = safe_pkg_version("DESeq2"),
  edgeR_version = safe_pkg_version("edgeR"),
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

if (!is.null(de_result)) {
  de_out <- de_result
  de_out$significant <- with(
    de_out,
    is.finite(padj) & padj <= padj_cutoff &
      is.finite(lfc) & abs(lfc) >= lfc_cutoff
  )
  utils::write.csv(
    de_out,
    file.path(outdir, paste0("DE_", benchmark_id, ".csv")),
    row.names = FALSE,
    na = ""
  )
}

write_session_info(file.path(outdir, paste0("session_", benchmark_id, ".txt")))

if (identical(status, "failed")) {
  message("FAILED: ", error_message)
  quit(save = "no", status = 1L)
}

message(sprintf("Completed in %.3f s", workflow_elapsed))
message("Genes tested: ", n_tested)
message("Significant genes (padj <= ", padj_cutoff,
        ", |log2FC| >= ", lfc_cutoff, "): ", n_significant)
message(sprintf("Process peak RSS: %.1f MB", peak_rss_kb / 1024))
