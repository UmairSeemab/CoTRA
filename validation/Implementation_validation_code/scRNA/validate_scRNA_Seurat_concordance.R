#!/usr/bin/env Rscript

# CoTRA scRNA-seq implementation validation
# Direct CoTRA versus standalone Seurat using the same two H5 inputs
# and the same fixed analysis settings.

options(stringsAsFactors = FALSE)

parse_cli <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
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
need_file <- function(x, label) {
  if (is.null(x) || !file.exists(x)) stop(label, " not found: ", x, call. = FALSE)
}
need_pkg <- function(x) {
  if (!requireNamespace(x, quietly = TRUE)) stop("Required package missing: ", x, call. = FALSE)
}

ari_base <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  choose2 <- function(x) x * (x - 1) / 2
  sum_nij <- sum(choose2(tab))
  sum_ai <- sum(choose2(rowSums(tab)))
  sum_bj <- sum(choose2(colSums(tab)))
  total <- choose2(n)
  expected <- (sum_ai * sum_bj) / total
  max_index <- 0.5 * (sum_ai + sum_bj)
  denom <- max_index - expected
  if (denom == 0) return(1)
  (sum_nij - expected) / denom
}

nmi_base <- function(a, b) {
  tab <- table(a, b)
  pxy <- tab / sum(tab)
  px <- rowSums(pxy)
  py <- colSums(pxy)

  nz <- pxy > 0
  mi <- sum(pxy[nz] * log(pxy[nz] / outer(px, py)[nz]))
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))

  if (hx == 0 && hy == 0) return(1)
  if (hx == 0 || hy == 0) return(0)
  mi / sqrt(hx * hy)
}

jaccard <- function(a, b) {
  a <- unique(as.character(a))
  b <- unique(as.character(b))
  u <- union(a, b)
  if (!length(u)) return(1)
  length(intersect(a, b)) / length(u)
}

extract_seurat <- function(x, depth = 0L) {
  if (inherits(x, "Seurat")) return(x)
  if (depth > 6L) return(NULL)
  if (is.list(x)) {
    for (i in seq_along(x)) {
      ans <- extract_seurat(x[[i]], depth + 1L)
      if (!is.null(ans)) return(ans)
    }
  }
  NULL
}

read_gex_h5 <- function(path) {
  x <- Seurat::Read10X_h5(path, use.names = TRUE, unique.features = TRUE)
  if (is.list(x)) {
    if ("Gene Expression" %in% names(x)) {
      x <- x[["Gene Expression"]]
    } else {
      x <- x[[1L]]
    }
  }
  x
}

prefix_cells <- function(mat, path) {
  sid <- sub("\\.h5$", "", basename(path), ignore.case = TRUE)
  colnames(mat) <- paste0(sid, "_", colnames(mat))
  mat
}

safe_cor <- function(a, b, method = "pearson") {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(stats::cor(a[ok], b[ok], method = method))
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))

wt_h5 <- arg(args, "wt-h5")
rd10_h5 <- arg(args, "rd10-h5")
cotra_clusters_file <- arg(args, "cotra-clusters")
cotra_markers_file <- arg(args, "cotra-markers")
cotra_session_file <- arg(args, "cotra-session")
outdir <- arg(args, "outdir", "validation_scRNA_Seurat")

need_file(wt_h5, "WT H5")
need_file(rd10_h5, "rd10 H5")
need_file(cotra_clusters_file, "CoTRA cluster assignments")
need_file(cotra_markers_file, "CoTRA marker table")
need_file(cotra_session_file, "CoTRA clustering session")
need_pkg("Seurat")
need_pkg("Matrix")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("Reading CoTRA exports...\n")
cotra_clusters <- utils::read.csv(
  cotra_clusters_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
cotra_markers <- utils::read.csv(
  cotra_markers_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("cell", "cluster_label") %in% names(cotra_clusters))) {
  stop("CoTRA cluster CSV must contain 'cell' and 'cluster_label'.", call. = FALSE)
}
if (!all(c("gene", "cluster", "avg_log2FC", "p_val_adj") %in% names(cotra_markers))) {
  stop("CoTRA marker CSV does not contain expected marker columns.", call. = FALSE)
}
if (anyDuplicated(cotra_clusters$cell)) {
  stop("Duplicate cell names in CoTRA cluster assignments.", call. = FALSE)
}

cat("Reading CoTRA clustering session...\n")
session_obj <- readRDS(cotra_session_file)
cotra_obj <- extract_seurat(session_obj)
if (is.null(cotra_obj)) {
  stop(
    "No Seurat object could be located recursively inside the CoTRA clustering session RDS.",
    call. = FALSE
  )
}
cat("Found CoTRA Seurat object: ",
    nrow(cotra_obj), " genes x ", ncol(cotra_obj), " cells\n", sep = "")

cat("Reading H5 files...\n")
wt <- prefix_cells(read_gex_h5(wt_h5), wt_h5)
rd10 <- prefix_cells(read_gex_h5(rd10_h5), rd10_h5)

common_genes <- intersect(rownames(wt), rownames(rd10))
wt <- wt[common_genes, , drop = FALSE]
rd10 <- rd10[common_genes, , drop = FALSE]

counts <- cbind(wt, rd10)

cat("Combined raw H5 input: ",
    nrow(counts), " genes x ", ncol(counts), " cells\n", sep = "")

missing_cells <- setdiff(cotra_clusters$cell, colnames(counts))
extra_cells <- setdiff(colnames(counts), cotra_clusters$cell)

if (length(missing_cells)) {
  stop(
    length(missing_cells),
    " CoTRA cells were not found in reconstructed H5 counts. First examples: ",
    paste(head(missing_cells), collapse = ", "),
    call. = FALSE
  )
}
if (length(extra_cells)) {
  stop(
    length(extra_cells),
    " reconstructed H5 cells were absent from CoTRA cluster assignments.",
    call. = FALSE
  )
}

# Match exact CoTRA cell order to remove ordering as a source of stochastic differences.
counts <- counts[, cotra_clusters$cell, drop = FALSE]

cat("Building direct standalone Seurat workflow...\n")
direct <- Seurat::CreateSeuratObject(
  counts = counts,
  project = "CoTRA_direct_validation",
  min.cells = 0,
  min.features = 0
)

direct$validation_sample <- ifelse(
  grepl("^GSM7474906_Wild_Type_non_treated_feature_bc_matrix_", colnames(direct)),
  "WT",
  "rd10"
)

direct <- Seurat::NormalizeData(
  direct,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

direct <- Seurat::FindVariableFeatures(
  direct,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

direct <- Seurat::ScaleData(
  direct,
  features = Seurat::VariableFeatures(direct),
  verbose = FALSE
)

direct <- Seurat::RunPCA(
  direct,
  features = Seurat::VariableFeatures(direct),
  npcs = 50,
  verbose = FALSE
)

direct <- Seurat::FindNeighbors(
  direct,
  dims = 1:7,
  k.param = 20,
  verbose = FALSE
)

direct <- Seurat::FindClusters(
  direct,
  resolution = 0.5,
  algorithm = 1,
  random.seed = 1234,
  verbose = FALSE
)

direct_clusters <- data.frame(
  cell = colnames(direct),
  direct_cluster = as.character(Seurat::Idents(direct)),
  stringsAsFactors = FALSE
)

utils::write.csv(
  direct_clusters,
  file.path(outdir, "direct_Seurat_cluster_assignments.csv"),
  row.names = FALSE
)

cat("Running direct Seurat marker analysis...\n")
direct_markers <- Seurat::FindAllMarkers(
  direct,
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  max.cells.per.ident = 500,
  verbose = FALSE
)

utils::write.csv(
  direct_markers,
  file.path(outdir, "direct_Seurat_all_cluster_markers.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------
# Core concordance metrics
# ----------------------------------------------------------
comparison <- merge(
  cotra_clusters[, c("cell", "cluster_label")],
  direct_clusters,
  by = "cell",
  all = FALSE,
  sort = FALSE
)

comparison$cluster_label <- as.character(comparison$cluster_label)
comparison$direct_cluster <- as.character(comparison$direct_cluster)

ari <- ari_base(comparison$cluster_label, comparison$direct_cluster)
nmi <- nmi_base(comparison$cluster_label, comparison$direct_cluster)

cont <- table(
  CoTRA = comparison$cluster_label,
  Direct_Seurat = comparison$direct_cluster
)
utils::write.csv(
  as.data.frame.matrix(cont),
  file.path(outdir, "cluster_contingency_matrix.csv")
)

# Map direct labels to CoTRA labels by maximum cell overlap.
mapping <- data.frame(
  direct_cluster = colnames(cont),
  mapped_CoTRA_cluster = vapply(
    seq_len(ncol(cont)),
    function(j) rownames(cont)[which.max(cont[, j])],
    character(1)
  ),
  overlap_cells = vapply(
    seq_len(ncol(cont)),
    function(j) max(cont[, j]),
    numeric(1)
  ),
  stringsAsFactors = FALSE
)
mapping$mapping_is_unique <- !duplicated(mapping$mapped_CoTRA_cluster) &
  !duplicated(mapping$mapped_CoTRA_cluster, fromLast = TRUE)

utils::write.csv(
  mapping,
  file.path(outdir, "cluster_label_mapping.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------
# HVG concordance
# ----------------------------------------------------------
cotra_hvg <- Seurat::VariableFeatures(cotra_obj)
direct_hvg <- Seurat::VariableFeatures(direct)

utils::write.csv(
  data.frame(gene = cotra_hvg),
  file.path(outdir, "CoTRA_variable_features.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(gene = direct_hvg),
  file.path(outdir, "direct_Seurat_variable_features.csv"),
  row.names = FALSE
)

hvg_intersection <- length(intersect(cotra_hvg, direct_hvg))
hvg_jaccard <- jaccard(cotra_hvg, direct_hvg)

# ----------------------------------------------------------
# PCA concordance
# Component signs are arbitrary; use absolute correlation.
# ----------------------------------------------------------
cotra_pca <- Seurat::Embeddings(cotra_obj, reduction = "pca")
direct_pca <- Seurat::Embeddings(direct, reduction = "pca")

common_pca_cells <- intersect(rownames(cotra_pca), rownames(direct_pca))
cotra_pca <- cotra_pca[common_pca_cells, , drop = FALSE]
direct_pca <- direct_pca[common_pca_cells, , drop = FALSE]

npc <- min(7L, ncol(cotra_pca), ncol(direct_pca))
pca_rows <- lapply(seq_len(npc), function(i) {
  r <- safe_cor(cotra_pca[, i], direct_pca[, i], method = "pearson")
  data.frame(
    PC = i,
    Pearson_r = r,
    absolute_Pearson_r = abs(r),
    stringsAsFactors = FALSE
  )
})
pca_summary <- do.call(rbind, pca_rows)

utils::write.csv(
  pca_summary,
  file.path(outdir, "PCA_concordance_PC1-PC7.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------
# Marker concordance after label mapping
# ----------------------------------------------------------
direct_markers_cmp <- direct_markers
direct_markers_cmp$direct_cluster <- as.character(direct_markers_cmp$cluster)

map_vec <- setNames(
  mapping$mapped_CoTRA_cluster,
  mapping$direct_cluster
)
direct_markers_cmp$mapped_CoTRA_cluster <- unname(
  map_vec[direct_markers_cmp$direct_cluster]
)

cotra_markers$cluster_chr <- as.character(cotra_markers$cluster)

cotra_all_keys <- paste(cotra_markers$cluster_chr, cotra_markers$gene, sep = "::")
direct_all_keys <- paste(
  direct_markers_cmp$mapped_CoTRA_cluster,
  direct_markers_cmp$gene,
  sep = "::"
)

marker_all_jaccard <- jaccard(cotra_all_keys, direct_all_keys)

cotra_sig <- cotra_markers[
  is.finite(cotra_markers$p_val_adj) &
    cotra_markers$p_val_adj <= 0.05,
  ,
  drop = FALSE
]
direct_sig <- direct_markers_cmp[
  is.finite(direct_markers_cmp$p_val_adj) &
    direct_markers_cmp$p_val_adj <= 0.05,
  ,
  drop = FALSE
]

cotra_sig_keys <- paste(cotra_sig$cluster_chr, cotra_sig$gene, sep = "::")
direct_sig_keys <- paste(
  direct_sig$mapped_CoTRA_cluster,
  direct_sig$gene,
  sep = "::"
)
marker_sig_jaccard <- jaccard(cotra_sig_keys, direct_sig_keys)

# log2FC agreement for overlapping marker gene-cluster pairs
cm <- cotra_markers[, c("cluster_chr", "gene", "avg_log2FC")]
names(cm)[3] <- "CoTRA_avg_log2FC"

dm <- direct_markers_cmp[, c("mapped_CoTRA_cluster", "gene", "avg_log2FC")]
names(dm) <- c("cluster_chr", "gene", "Direct_avg_log2FC")

marker_overlap <- merge(cm, dm, by = c("cluster_chr", "gene"))
marker_lfc_pearson <- safe_cor(
  marker_overlap$CoTRA_avg_log2FC,
  marker_overlap$Direct_avg_log2FC,
  method = "pearson"
)
marker_lfc_spearman <- safe_cor(
  marker_overlap$CoTRA_avg_log2FC,
  marker_overlap$Direct_avg_log2FC,
  method = "spearman"
)

utils::write.csv(
  marker_overlap,
  file.path(outdir, "overlapping_marker_log2FC_values.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------
# Summary
# ----------------------------------------------------------
summary <- data.frame(
  metric = c(
    "Input genes",
    "Input cells",
    "WT cells",
    "rd10 cells",
    "CoTRA clusters",
    "Direct Seurat clusters",
    "Adjusted Rand Index",
    "Normalized Mutual Information",
    "CoTRA HVGs",
    "Direct Seurat HVGs",
    "HVG intersection",
    "HVG Jaccard",
    "Median absolute PCA correlation PC1-PC7",
    "Minimum absolute PCA correlation PC1-PC7",
    "CoTRA marker rows",
    "Direct Seurat marker rows",
    "All marker gene-cluster pair Jaccard",
    "Significant marker gene-cluster pair Jaccard",
    "Overlapping marker pairs for log2FC comparison",
    "Marker avg_log2FC Pearson r",
    "Marker avg_log2FC Spearman rho"
  ),
  value = c(
    nrow(counts),
    ncol(counts),
    sum(direct$validation_sample == "WT"),
    sum(direct$validation_sample == "rd10"),
    length(unique(comparison$cluster_label)),
    length(unique(comparison$direct_cluster)),
    ari,
    nmi,
    length(cotra_hvg),
    length(direct_hvg),
    hvg_intersection,
    hvg_jaccard,
    stats::median(pca_summary$absolute_Pearson_r, na.rm = TRUE),
    min(pca_summary$absolute_Pearson_r, na.rm = TRUE),
    nrow(cotra_markers),
    nrow(direct_markers),
    marker_all_jaccard,
    marker_sig_jaccard,
    nrow(marker_overlap),
    marker_lfc_pearson,
    marker_lfc_spearman
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  summary,
  file.path(outdir, "S8_scRNA_implementation_validation_summary.csv"),
  row.names = FALSE
)

saveRDS(
  direct,
  file.path(outdir, "direct_Seurat_validation_object.rds"),
  compress = "xz"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "validation_sessionInfo.txt")
)

cat("\nValidation complete.\n")
cat("Output: ", normalizePath(outdir), "\n\n", sep = "")
print(summary, row.names = FALSE)
