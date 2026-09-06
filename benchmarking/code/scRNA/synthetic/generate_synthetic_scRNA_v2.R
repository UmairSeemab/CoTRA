#!/usr/bin/env Rscript

# ============================================================
# Strongly separated synthetic scRNA-seq count matrices
# for CoTRA computational scalability benchmarking.
#
# One 20,000-gene x 50,000-cell sparse matrix is generated,
# then nested 2,500 / 5,000 / 10,000 / 25,000 / 50,000-cell
# subsets are saved as independent RDS files.
#
# The signal is deliberately strong so that CoTRA's unchanged
# default Louvain resolution (0.5) yields multiple clusters.
# These data are for computational benchmarking only.
# ============================================================

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

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Matrix is required.", call. = FALSE)
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))

outdir <- arg_value(
  args,
  "outdir",
  "/home/umair/Downloads/CoTRA/benchmark/synthetic_scRNA"
)

n_genes <- as_int(arg_value(args, "genes", "20000"), "--genes")
n_cells <- as_int(arg_value(args, "cells", "50000"), "--cells")
n_clusters <- as_int(arg_value(args, "clusters", "12"), "--clusters")
markers_per_cluster <- as_int(
  arg_value(args, "markers-per-cluster", "150"),
  "--markers-per-cluster"
)
background_density <- as_num(
  arg_value(args, "background-density", "0.005"),
  "--background-density"
)
background_lambda <- as_num(
  arg_value(args, "background-lambda", "1.0"),
  "--background-lambda"
)
marker_lambda <- as_num(
  arg_value(args, "marker-lambda", "8.0"),
  "--marker-lambda"
)
housekeeping_genes <- as_int(
  arg_value(args, "housekeeping-genes", "100"),
  "--housekeeping-genes"
)
housekeeping_lambda <- as_num(
  arg_value(args, "housekeeping-lambda", "2.0"),
  "--housekeeping-lambda"
)
seed <- as_int(arg_value(args, "seed", "1234"), "--seed")

sizes <- c(2500L, 5000L, 10000L, 25000L, 50000L)
sizes <- sizes[sizes <= n_cells]

if (length(sizes) == 0L) stop("--cells must be at least 2500.", call. = FALSE)
if (n_clusters < 2L) stop("--clusters must be at least 2.", call. = FALSE)
if (markers_per_cluster < 20L) stop("--markers-per-cluster must be at least 20.", call. = FALSE)
if (background_density <= 0 || background_density >= 0.25) {
  stop("--background-density must be >0 and <0.25.", call. = FALSE)
}
if (marker_lambda <= background_lambda) {
  stop("--marker-lambda must exceed --background-lambda.", call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

set.seed(seed)

message("Generating strongly separated synthetic scRNA-seq benchmark matrix")
message("Genes: ", n_genes)
message("Cells: ", n_cells)
message("Clusters: ", n_clusters)
message("Markers per cluster: ", markers_per_cluster)
message("Background density: ", background_density)
message("Background lambda: ", background_lambda)
message("Marker lambda: ", marker_lambda)
message("Housekeeping genes: ", housekeeping_genes)
message("Housekeeping lambda: ", housekeeping_lambda)
message("Seed: ", seed)

# ------------------------------------------------------------
# Feature names
# Avoid underscores because Seurat replaces underscores in
# feature names.
# ------------------------------------------------------------
gene_names <- sprintf("Gene%05d", seq_len(n_genes))

n_mito <- min(20L, n_genes)
gene_names[seq_len(n_mito)] <- sprintf("mt-Syn%03d", seq_len(n_mito))

n_ribo <- min(80L, max(0L, n_genes - n_mito))
ribo_start <- n_mito + 1L
if (n_ribo > 0L) {
  ribo_idx <- seq.int(ribo_start, length.out = n_ribo)
  n_rpl <- ceiling(n_ribo / 2L)
  gene_names[ribo_idx[seq_len(n_rpl)]] <- sprintf("RplSyn%03d", seq_len(n_rpl))
  if (n_ribo > n_rpl) {
    gene_names[ribo_idx[(n_rpl + 1L):n_ribo]] <- sprintf(
      "RpsSyn%03d",
      seq_len(n_ribo - n_rpl)
    )
  }
}

house_start <- n_mito + n_ribo + 1L
house_end <- house_start + housekeeping_genes - 1L
if (house_end > n_genes) stop("Too many housekeeping genes.", call. = FALSE)

house_rows <- if (housekeeping_genes > 0L) seq.int(house_start, house_end) else integer()
if (length(house_rows)) {
  gene_names[house_rows] <- sprintf("House%03d", seq_len(length(house_rows)))
}

marker_start <- house_end + 1L
n_marker_genes <- n_clusters * markers_per_cluster
marker_end <- marker_start + n_marker_genes - 1L

if (marker_end > n_genes) {
  stop("Not enough genes for the requested cluster marker design.", call. = FALSE)
}

marker_rows_by_cluster <- vector("list", n_clusters)
for (k in seq_len(n_clusters)) {
  idx <- seq.int(
    marker_start + (k - 1L) * markers_per_cluster,
    length.out = markers_per_cluster
  )
  marker_rows_by_cluster[[k]] <- idx
  gene_names[idx] <- sprintf(
    "MarkerC%02d-%03d",
    k,
    seq_len(markers_per_cluster)
  )
}

cell_names <- sprintf("SyntheticCell-%06d", seq_len(n_cells))

# Balanced cyclic assignment ensures nested subsets remain balanced.
truth_cluster <- ((seq_len(n_cells) - 1L) %% n_clusters) + 1L

# ------------------------------------------------------------
# Sparse random background.
# ------------------------------------------------------------
message("Generating sparse background expression...")

background <- Matrix::rsparsematrix(
  nrow = n_genes,
  ncol = n_cells,
  density = background_density,
  rand.x = function(n) {
    as.numeric(stats::rpois(n, lambda = background_lambda) + 1L)
  }
)

# ------------------------------------------------------------
# Shared housekeeping expression.
# ------------------------------------------------------------
if (length(house_rows)) {
  message("Adding shared housekeeping expression...")

  n_house_entries <- as.integer(length(house_rows) * n_cells)

  house_matrix <- Matrix::sparseMatrix(
    i = rep(house_rows, times = n_cells),
    j = rep(seq_len(n_cells), each = length(house_rows)),
    x = as.numeric(
      stats::rpois(n_house_entries, lambda = housekeeping_lambda) + 1L
    ),
    dims = c(n_genes, n_cells),
    giveCsparse = TRUE
  )

  background <- background + house_matrix
  rm(house_matrix)
  invisible(gc())
}

# ------------------------------------------------------------
# Strong cluster-specific marker blocks.
# Each cluster gets a non-overlapping marker panel.
# ------------------------------------------------------------
message("Adding strong cluster-specific marker expression...")

total_marker_entries <- as.integer(n_cells * markers_per_cluster)

marker_i <- integer(total_marker_entries)
marker_j <- integer(total_marker_entries)
marker_x <- numeric(total_marker_entries)

cursor <- 1L

for (k in seq_len(n_clusters)) {

  cells_k <- which(truth_cluster == k)
  rows_k <- marker_rows_by_cluster[[k]]
  n_entries <- length(cells_k) * length(rows_k)

  idx <- seq.int(cursor, length.out = n_entries)

  marker_i[idx] <- rep(rows_k, times = length(cells_k))
  marker_j[idx] <- rep(cells_k, each = length(rows_k))
  marker_x[idx] <- as.numeric(
    stats::rpois(n_entries, lambda = marker_lambda) + 1L
  )

  cursor <- cursor + n_entries
}

marker_matrix <- Matrix::sparseMatrix(
  i = marker_i,
  j = marker_j,
  x = marker_x,
  dims = c(n_genes, n_cells),
  giveCsparse = TRUE
)

rm(marker_i, marker_j, marker_x)
invisible(gc())

counts <- Matrix::drop0(background + marker_matrix)

rm(background, marker_matrix)
invisible(gc())

rownames(counts) <- gene_names
colnames(counts) <- cell_names

message(
  "Full matrix nonzero entries: ",
  format(length(counts@x), big.mark = ",")
)

# ------------------------------------------------------------
# Save nested independent inputs.
# ------------------------------------------------------------
manifest <- list()

for (n in sizes) {

  message("Saving ", n, "-cell nested subset...")

  subset_counts <- counts[, seq_len(n), drop = FALSE]

  subset_path <- file.path(
    outdir,
    sprintf("counts_synthetic_%06dc.rds", n)
  )

  saveRDS(
    subset_counts,
    subset_path,
    compress = "xz"
  )

  tab <- table(truth_cluster[seq_len(n)])

  manifest[[length(manifest) + 1L]] <- data.frame(
    n_genes = nrow(subset_counts),
    n_cells = ncol(subset_counts),
    nonzero_entries = length(subset_counts@x),
    expected_clusters = n_clusters,
    min_cells_per_truth_cluster = min(tab),
    max_cells_per_truth_cluster = max(tab),
    file = normalizePath(subset_path, mustWork = FALSE),
    stringsAsFactors = FALSE
  )

  rm(subset_counts)
  invisible(gc())
}

manifest_df <- do.call(rbind, manifest)

utils::write.csv(
  manifest_df,
  file.path(outdir, "synthetic_scRNA_manifest.csv"),
  row.names = FALSE
)

truth_df <- data.frame(
  cell = cell_names,
  synthetic_cluster = sprintf("C%02d", truth_cluster),
  stringsAsFactors = FALSE
)

utils::write.csv(
  truth_df,
  file.path(outdir, "synthetic_scRNA_truth_clusters.csv"),
  row.names = FALSE
)

metadata <- data.frame(
  parameter = c(
    "purpose",
    "seed",
    "n_genes",
    "max_n_cells",
    "n_clusters",
    "markers_per_cluster",
    "background_density",
    "background_lambda",
    "marker_lambda",
    "housekeeping_genes",
    "housekeeping_lambda",
    "nested_subsets"
  ),
  value = c(
    "computational benchmarking only; not biological validation",
    seed,
    n_genes,
    n_cells,
    n_clusters,
    markers_per_cluster,
    background_density,
    background_lambda,
    marker_lambda,
    housekeeping_genes,
    housekeeping_lambda,
    paste(sizes, collapse = ",")
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  metadata,
  file.path(outdir, "synthetic_scRNA_generation_metadata.csv"),
  row.names = FALSE
)

message("Synthetic datasets created successfully.")
message("Output directory: ", normalizePath(outdir, mustWork = FALSE))
