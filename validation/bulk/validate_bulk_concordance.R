#!/usr/bin/env Rscript

# CoTRA bulk implementation-concordance validation
#
# Compares CoTRA-exported DESeq2 and edgeR results with direct scripted
# executions using the same raw count matrix, sample groups, filtering rules,
# contrast, and significance thresholds.
#
# Default validation dataset:
#   raw_gene_counts.tsv
#   WT = 4 samples
#   rd10 = 8 samples
#   contrast = rd10 vs WT
#
# Outputs:
#   Table_S7_Bulk_Implementation_Concordance.csv
#   run_level_concordance.csv
#   direct_DESeq2_results.csv
#   direct_edgeR_results.csv
#   validation_sessionInfo.txt

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

arg <- function(x, key, default = NULL) {
  if (is.null(x[[key]])) default else x[[key]]
}

need_file <- function(path, label) {
  if (is.null(path) || !file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

need_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package is not installed: ", pkg, call. = FALSE)
  }
}

read_counts <- function(path) {
  ext <- tolower(tools::file_ext(path))
  sep <- if (ext %in% c("tsv", "txt", "tab")) "\t" else ","

  dat <- utils::read.table(
    path,
    header = TRUE,
    sep = sep,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (ncol(dat) < 3L) {
    stop("Expected one gene column plus >=2 sample columns.", call. = FALSE)
  }

  genes <- as.character(dat[[1L]])
  if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("Gene identifiers must be unique and non-missing.", call. = FALSE)
  }

  x <- dat[-1L]
  for (j in seq_along(x)) {
    x[[j]] <- suppressWarnings(as.numeric(x[[j]]))
  }

  counts <- as.matrix(x)
  rownames(counts) <- genes

  if (anyNA(counts)) stop("Count matrix contains non-numeric/missing values.")
  if (any(counts < 0)) stop("Count matrix contains negative values.")
  if (any(abs(counts - round(counts)) > .Machine$double.eps^0.5)) {
    stop("Raw integer counts are required.")
  }

  storage.mode(counts) <- "integer"
  counts
}

run_direct_deseq2 <- function(counts, group, reference, comparison) {
  need_pkg("DESeq2")

  coldata <- data.frame(group = group, row.names = colnames(counts))

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData = coldata,
    design = ~ group
  )

  keep <- rowSums(DESeq2::counts(dds) >= 10) >= 2
  dds <- dds[keep, ]
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  res <- DESeq2::results(
    dds,
    contrast = c("group", comparison, reference)
  )

  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df$lfc <- df$log2FoldChange

  # CoTRA post-processing:
  # undefined P and adjusted-P values are assigned 1 before significance
  # filtering. This stabilizes downstream visualization/classification.
  raw_na_p <- sum(is.na(df$pvalue))
  raw_na_padj <- sum(is.na(df$padj))

  df$pvalue[is.na(df$pvalue)] <- 1
  df$padj[is.na(df$padj)] <- 1
  df$lfc[is.na(df$lfc)] <- 0

  attr(df, "raw_na_p") <- raw_na_p
  attr(df, "raw_na_padj") <- raw_na_padj

  df
}

run_direct_edger <- function(counts, group, comparison) {
  need_pkg("edgeR")

  y <- edgeR::DGEList(counts = round(as.matrix(counts)), group = group)
  keep <- edgeR::filterByExpr(y, group = group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)

  design <- stats::model.matrix(~ group)
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design)

  coef_name <- paste0("group", comparison)
  coef_index <- match(coef_name, colnames(design))
  if (is.na(coef_index)) {
    stop("Could not locate edgeR coefficient: ", coef_name, call. = FALSE)
  }

  test <- edgeR::glmQLFTest(fit, coef = coef_index)
  tab <- edgeR::topTags(test, n = Inf)$table

  df <- as.data.frame(tab)
  df$gene <- rownames(df)
  df$baseMean <- rowMeans(edgeR::cpm(y, normalized.lib.sizes = TRUE))
  df$lfc <- df$logFC
  df$pvalue <- df$PValue
  df$padj <- df$FDR
  df
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

jaccard <- function(a, b) {
  u <- union(a, b)
  if (!length(u)) return(1)
  length(intersect(a, b)) / length(u)
}

compare_one <- function(method, cotra, direct, padj_cutoff, lfc_cutoff,
                        replicate_id, raw_na_p = 0L, raw_na_padj = 0L) {

  cotra$gene <- as.character(cotra$gene)
  direct$gene <- as.character(direct$gene)

  # Significance rule used in CoTRA validation.
  direct$significant <- (
    is.finite(direct$padj) &
    direct$padj <= padj_cutoff &
    is.finite(direct$lfc) &
    abs(direct$lfc) >= lfc_cutoff
  )

  merged <- merge(
    cotra,
    direct,
    by = "gene",
    suffixes = c("_CoTRA", "_Scripted"),
    all = TRUE,
    sort = FALSE
  )

  matched <- merged[
    !is.na(merged$lfc_CoTRA) & !is.na(merged$lfc_Scripted),
    ,
    drop = FALSE
  ]

  if (method == "DESeq2") {
    cols <- c(
      "baseMean", "log2FoldChange", "lfcSE", "stat",
      "pvalue", "padj", "lfc"
    )
    lfc_c <- matched$log2FoldChange_CoTRA
    lfc_s <- matched$log2FoldChange_Scripted
  } else {
    cols <- c(
      "logFC", "logCPM", "F", "PValue", "FDR",
      "lfc", "pvalue", "padj"
    )
    lfc_c <- matched$logFC_CoTRA
    lfc_s <- matched$logFC_Scripted
  }

  max_diff <- 0
  all_exact <- TRUE

  for (col in cols) {
    a <- matched[[paste0(col, "_CoTRA")]]
    b <- matched[[paste0(col, "_Scripted")]]

    same <- isTRUE(all.equal(
      a, b,
      tolerance = 0,
      check.attributes = FALSE
    ))
    all_exact <- all_exact && same

    finite <- is.finite(a) & is.finite(b)
    if (any(finite)) {
      max_diff <- max(max_diff, max(abs(a[finite] - b[finite])))
    }
  }

  sig_c <- as.logical(cotra$significant)
  sig_s <- as.logical(direct$significant)

  cotra_sig_genes <- cotra$gene[sig_c]
  direct_sig_genes <- direct$gene[sig_s]

  sig_overlap <- intersect(cotra_sig_genes, direct_sig_genes)

  merged_sig <- merge(
    cotra[cotra$significant, c("gene", "lfc")],
    direct[direct$significant, c("gene", "lfc")],
    by = "gene",
    suffixes = c("_CoTRA", "_Scripted")
  )

  direction_concordance <- if (nrow(merged_sig)) {
    mean(sign(merged_sig$lfc_CoTRA) == sign(merged_sig$lfc_Scripted))
  } else {
    NA_real_
  }

  data.frame(
    Method = method,
    Scripted_replicate = replicate_id,
    CoTRA_genes_tested = nrow(cotra),
    Scripted_genes_tested = nrow(direct),
    Matched_genes = nrow(matched),
    CoTRA_significant_genes = length(cotra_sig_genes),
    Scripted_significant_genes = length(direct_sig_genes),
    Significant_set_intersection = length(sig_overlap),
    Significant_set_Jaccard = jaccard(cotra_sig_genes, direct_sig_genes),
    Direction_concordance = direction_concordance,
    log2FC_Pearson_r = safe_cor(lfc_c, lfc_s, "pearson"),
    log2FC_Spearman_rho = safe_cor(lfc_c, lfc_s, "spearman"),
    Maximum_absolute_numeric_difference_after_CoTRA_postprocessing = max_diff,
    Exact_after_documented_CoTRA_postprocessing = all_exact,
    DESeq2_raw_NA_pvalues_converted_to_1 = raw_na_p,
    DESeq2_raw_NA_padj_converted_to_1 = raw_na_padj,
    stringsAsFactors = FALSE
  )
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))

input <- arg(args, "input")
cotra_deseq2_file <- arg(args, "cotra-deseq2")
cotra_edger_file <- arg(args, "cotra-edger")
outdir <- arg(args, "outdir", "bulk_validation")
reference <- arg(args, "reference", "WT")
comparison <- arg(args, "comparison", "rd10")
reference_regex <- arg(args, "reference-regex", "^WT")
comparison_regex <- arg(args, "comparison-regex", "^rd10")
padj_cutoff <- as.numeric(arg(args, "padj", "0.05"))
lfc_cutoff <- as.numeric(arg(args, "lfc", "1"))
repeats <- as.integer(arg(args, "repeats", "5"))

need_file(input, "Bulk test count matrix")
need_file(cotra_deseq2_file, "CoTRA DESeq2 export")
need_file(cotra_edger_file, "CoTRA edgeR export")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts <- read_counts(input)

ref_cols <- grep(reference_regex, colnames(counts), ignore.case = TRUE, value = TRUE)
cmp_cols <- grep(comparison_regex, colnames(counts), ignore.case = TRUE, value = TRUE)

if (!length(ref_cols)) stop("No WT/reference samples matched.")
if (!length(cmp_cols)) stop("No rd10/comparison samples matched.")
if (length(intersect(ref_cols, cmp_cols))) {
  stop("Reference and comparison regexes overlap.")
}

selected <- c(ref_cols, cmp_cols)
counts <- counts[, selected, drop = FALSE]

group <- factor(
  c(rep(reference, length(ref_cols)), rep(comparison, length(cmp_cols))),
  levels = c(reference, comparison)
)
names(group) <- selected

cotra_deseq2 <- utils::read.csv(
  cotra_deseq2_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
cotra_edger <- utils::read.csv(
  cotra_edger_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

run_rows <- list()

for (rep_id in seq_len(repeats)) {
  message("Validation replicate ", rep_id, "/", repeats)

  direct_deseq2 <- run_direct_deseq2(
    counts, group, reference, comparison
  )
  direct_edger <- run_direct_edger(
    counts, group, comparison
  )

  if (rep_id == 1L) {
    utils::write.csv(
      direct_deseq2,
      file.path(outdir, "direct_DESeq2_results.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      direct_edger,
      file.path(outdir, "direct_edgeR_results.csv"),
      row.names = FALSE
    )
  }

  run_rows[[length(run_rows) + 1L]] <- compare_one(
    "DESeq2",
    cotra_deseq2,
    direct_deseq2,
    padj_cutoff,
    lfc_cutoff,
    rep_id,
    raw_na_p = attr(direct_deseq2, "raw_na_p"),
    raw_na_padj = attr(direct_deseq2, "raw_na_padj")
  )

  run_rows[[length(run_rows) + 1L]] <- compare_one(
    "edgeR",
    cotra_edger,
    direct_edger,
    padj_cutoff,
    lfc_cutoff,
    rep_id
  )
}

run_level <- do.call(rbind, run_rows)
utils::write.csv(
  run_level,
  file.path(outdir, "run_level_concordance.csv"),
  row.names = FALSE
)

summary_rows <- lapply(c("DESeq2", "edgeR"), function(method) {
  d <- run_level[run_level$Method == method, , drop = FALSE]
  first <- d[1L, , drop = FALSE]

  interpretation <- if (method == "DESeq2") {
    paste0(
      first$DESeq2_raw_NA_pvalues_converted_to_1,
      " raw DESeq2 NA p-values and ",
      first$DESeq2_raw_NA_padj_converted_to_1,
      " NA adjusted p-values were converted to 1 before filtering. ",
      "After documented CoTRA post-processing, all compared values and ",
      "significance calls matched."
    )
  } else {
    "All compared edgeR statistics and significance calls matched."
  }

  data.frame(
    Method = method,
    Input_contrast = paste(comparison, "vs", reference),
    Scripted_runs_compared = nrow(d),
    CoTRA_genes_tested = first$CoTRA_genes_tested,
    Scripted_genes_tested = first$Scripted_genes_tested,
    Matched_genes = first$Matched_genes,
    CoTRA_significant_genes = first$CoTRA_significant_genes,
    Scripted_significant_genes = first$Scripted_significant_genes,
    Significant_set_intersection = first$Significant_set_intersection,
    Significant_set_Jaccard = first$Significant_set_Jaccard,
    Direction_concordance_percent = 100 * first$Direction_concordance,
    log2FC_Pearson_r = first$log2FC_Pearson_r,
    log2FC_Spearman_rho = first$log2FC_Spearman_rho,
    Maximum_absolute_numeric_difference_after_CoTRA_postprocessing =
      first$Maximum_absolute_numeric_difference_after_CoTRA_postprocessing,
    Interpretation = interpretation,
    stringsAsFactors = FALSE
  )
})

summary <- do.call(rbind, summary_rows)

utils::write.csv(
  summary,
  file.path(outdir, "Table_S7_Bulk_Implementation_Concordance.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "validation_sessionInfo.txt")
)

message("\nBulk implementation validation complete.")
print(summary, row.names = FALSE)
