#!/usr/bin/env Rscript

# ============================================================
# CoTRA scRNA-seq synthetic computational benchmark v3
#
# Core workflow is benchmarked exactly through graph clustering.
# Marker analysis is benchmarked separately using the known
# simulated truth labels, so marker runtime does not depend on
# whether Louvain reproduces the synthetic labels.
#
# Disk input and browser/Shiny rendering are excluded from the
# workflow timer. Use GNU /usr/bin/time -v externally for RAM.
# ============================================================

options(stringsAsFactors = FALSE)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  MKL_DYNAMIC = "FALSE",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

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

as_int <- function(x, name) {
  v <- suppressWarnings(as.integer(x))
  if (length(v) != 1L || is.na(v)) stop(name, " must be an integer.", call. = FALSE)
  v
}

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package is not installed: ", pkg, call. = FALSE)
}

safe_pkg_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}

proc_status_kb <- function(field) {
  path <- "/proc/self/status"
  if (!file.exists(path)) return(NA_real_)
  txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  hit <- grep(paste0("^", field, ":"), txt, value = TRUE)
  if (!length(hit)) return(NA_real_)
  suppressWarnings(as.numeric(sub(".*?:\\s*([0-9]+).*", "\\1", hit[[1L]])))
}

measure_step <- function(step_name, expr, logenv, category = "core") {
  invisible(gc())
  start <- proc.time()[["elapsed"]]
  ans <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - start
  logenv$rows[[length(logenv$rows) + 1L]] <- data.frame(
    category = category,
    step = step_name,
    elapsed_seconds = as.numeric(elapsed),
    stringsAsFactors = FALSE
  )
  ans
}

extract_counts <- function(path) {
  x <- readRDS(path)
  if (inherits(x, "Matrix") || is.matrix(x)) return(x)
  stop("Synthetic benchmark input must be a matrix/Matrix RDS.", call. = FALSE)
}

add_qc <- function(seu, species) {
  if (species == "mouse") {
    seu@meta.data$percent.mito <- Seurat::PercentageFeatureSet(seu, assay = "RNA", pattern = "^mt-")
    seu@meta.data$percent.ribo <- Seurat::PercentageFeatureSet(seu, assay = "RNA", pattern = "^Rp[ls]")
  } else if (species == "human") {
    seu@meta.data$percent.mito <- Seurat::PercentageFeatureSet(seu, assay = "RNA", pattern = "^MT-")
    seu@meta.data$percent.ribo <- Seurat::PercentageFeatureSet(seu, assay = "RNA", pattern = "^RP[LS]")
  } else {
    seu@meta.data$percent.mito <- NA_real_
    seu@meta.data$percent.ribo <- NA_real_
  }
  seu
}

write_session <- function(path) {
  writeLines(capture.output(sessionInfo()), path)
}

require_pkg("Seurat")
require_pkg("SeuratObject")
require_pkg("Matrix")

if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
  try(RhpcBLASctl::omp_set_num_threads(1L), silent = TRUE)
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))

input <- arg_value(args, "input", NULL)
if (is.null(input) || !file.exists(input)) stop("--input must point to an existing RDS.", call. = FALSE)

rep_id <- as_int(arg_value(args, "rep", "1"), "--rep")
seed <- as_int(arg_value(args, "seed", "1234"), "--seed")
species <- tolower(arg_value(args, "species", "mouse"))
truth_clusters <- as_int(arg_value(args, "truth-clusters", "12"), "--truth-clusters")
outdir <- arg_value(args, "outdir", "benchmark_results/scrna_synthetic")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Loading pre-sized raw count matrix: ", input)
counts <- extract_counts(input)
if (!inherits(counts, "Matrix")) counts <- Matrix::Matrix(counts, sparse = TRUE)
if (any(counts@x < 0)) stop("Counts contain negative values.", call. = FALSE)

n_genes <- nrow(counts)
n_cells <- ncol(counts)

id <- sprintf("scrna_synthetic_%06dc_rep%02d", n_cells, rep_id)

message("Benchmark: ", id)
message("Dimensions: ", n_genes, " genes x ", n_cells, " cells")
message("Species: ", species)
message("Synthetic truth groups for marker benchmark: ", truth_clusters)
message("Thread mode: single-thread constrained")

logenv <- new.env(parent = emptyenv())
logenv$rows <- list()

status <- "success"
error_message <- NA_character_
seu <- NULL
markers <- NULL
inferred_clusters <- NA_integer_
marker_truth_groups <- NA_integer_

core_start <- proc.time()[["elapsed"]]

tryCatch({
  seu <- measure_step("create_seurat_and_qc_metrics", {
    obj <- Seurat::CreateSeuratObject(counts = counts, assay = "RNA", project = "CoTRA_scRNA_benchmark")
    add_qc(obj, species)
  }, logenv)

  seu <- measure_step("normalize_lognormalize", {
    Seurat::NormalizeData(
      seu, assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  }, logenv)

  seu <- measure_step("find_variable_features", {
    Seurat::FindVariableFeatures(
      seu, assay = "RNA",
      selection.method = "vst",
      nfeatures = min(2000L, nrow(seu)),
      verbose = FALSE
    )
  }, logenv)

  features <- Seurat::VariableFeatures(seu, assay = "RNA")
  features <- unique(features[features %in% rownames(seu[["RNA"]])])
  if (length(features) < 50L) stop("Fewer than 50 variable features available.")

  seu <- measure_step("scale_data", {
    Seurat::ScaleData(seu, assay = "RNA", features = features, verbose = FALSE)
  }, logenv)

  npcs <- min(50L, ncol(seu) - 1L, length(features) - 1L)
  if (npcs < 30L) stop("Fewer than 30 PCs can be calculated.")

  seu <- measure_step("pca", {
    Seurat::RunPCA(seu, assay = "RNA", features = features, npcs = npcs, verbose = FALSE)
  }, logenv)

  dims <- seq_len(min(30L, npcs))

  set.seed(seed)
  seu <- measure_step("umap", {
    Seurat::RunUMAP(
      seu,
      reduction = "pca",
      dims = dims,
      n.neighbors = min(30L, ncol(seu) - 1L),
      min.dist = 0.3,
      metric = "cosine",
      reduction.name = "umap",
      reduction.key = "UMAP_",
      seed.use = seed,
      verbose = FALSE
    )
  }, logenv)

  seu <- measure_step("find_neighbors", {
    Seurat::FindNeighbors(
      seu,
      reduction = "pca",
      dims = dims,
      k.param = min(20L, ncol(seu) - 1L),
      verbose = FALSE
    )
  }, logenv)

  set.seed(seed)
  seu <- measure_step("find_clusters", {
    Seurat::FindClusters(
      seu,
      resolution = 0.5,
      algorithm = 1L,
      random.seed = seed,
      verbose = FALSE
    )
  }, logenv)

  inferred_clusters <- length(unique(as.character(Seurat::Idents(seu))))

}, error = function(e) {
  status <<- "failed"
  error_message <<- conditionMessage(e)
})

core_elapsed <- proc.time()[["elapsed"]] - core_start

marker_elapsed <- NA_real_

if (identical(status, "success")) {
  tryCatch({
    # Benchmark marker computation independently using known synthetic groups.
    truth <- factor(
      ((seq_len(ncol(seu)) - 1L) %% truth_clusters) + 1L,
      levels = seq_len(truth_clusters)
    )
    names(truth) <- colnames(seu)
    seu$synthetic_truth_cluster <- paste0("C", sprintf("%02d", as.integer(truth)))
    Seurat::Idents(seu) <- "synthetic_truth_cluster"
    marker_truth_groups <- length(unique(as.character(Seurat::Idents(seu))))

    set.seed(seed)
    marker_start <- proc.time()[["elapsed"]]

    markers <- measure_step("find_all_markers_truth_labels", {
      Seurat::FindAllMarkers(
        seu,
        assay = "RNA",
        only.pos = TRUE,
        test.use = "wilcox",
        logfc.threshold = 0.25,
        min.pct = 0.10,
        max.cells.per.ident = 500L,
        latent.vars = NULL,
        verbose = FALSE
      )
    }, logenv, category = "marker_truth_labels")

    marker_elapsed <- proc.time()[["elapsed"]] - marker_start

  }, error = function(e) {
    status <<- "failed"
    error_message <<- paste0("Marker benchmark failed: ", conditionMessage(e))
  })
}

total_elapsed <- core_elapsed + ifelse(is.finite(marker_elapsed), marker_elapsed, 0)

steps <- if (length(logenv$rows)) do.call(rbind, logenv$rows) else
  data.frame(category = character(), step = character(), elapsed_seconds = numeric())

peak_kb <- proc_status_kb("VmHWM")

summary <- data.frame(
  benchmark_id = id,
  workflow = "scrna_synthetic",
  replicate = rep_id,
  seed = seed,
  n_genes_input = n_genes,
  n_cells = n_cells,
  species = species,
  core_workflow_elapsed_seconds = as.numeric(core_elapsed),
  marker_truth_elapsed_seconds = as.numeric(marker_elapsed),
  total_core_plus_marker_seconds = as.numeric(total_elapsed),
  inferred_louvain_clusters = inferred_clusters,
  marker_truth_groups = marker_truth_groups,
  marker_rows = if (!is.null(markers)) nrow(markers) else NA_integer_,
  significant_marker_rows_padj_0_05 = if (!is.null(markers) && "p_val_adj" %in% names(markers))
    sum(is.finite(markers$p_val_adj) & markers$p_val_adj <= 0.05, na.rm = TRUE) else NA_integer_,
  hvg_used = if (!is.null(seu)) length(Seurat::VariableFeatures(seu, assay = "RNA")) else NA_integer_,
  pcs_used = if (exists("npcs", inherits = FALSE)) npcs else NA_integer_,
  peak_rss_mb_process = peak_kb / 1024,
  status = status,
  error_message = error_message,
  R_version = paste(R.version$major, R.version$minor, sep = "."),
  Seurat_version = safe_pkg_version("Seurat"),
  SeuratObject_version = safe_pkg_version("SeuratObject"),
  Matrix_version = safe_pkg_version("Matrix"),
  hostname = Sys.info()[["nodename"]],
  timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)

utils::write.csv(summary, file.path(outdir, paste0("overall_", id, ".csv")), row.names = FALSE, na = "")
utils::write.csv(steps, file.path(outdir, paste0("steps_", id, ".csv")), row.names = FALSE, na = "")

if (!is.null(markers)) {
  utils::write.csv(
    markers,
    file.path(outdir, paste0("markers_truth_", id, ".csv")),
    row.names = FALSE,
    na = ""
  )
}

write_session(file.path(outdir, paste0("session_", id, ".txt")))

if (status == "failed") {
  message("FAILED: ", error_message)
  quit(save = "no", status = 1L)
}

message(sprintf("Core workflow completed in %.3f s", core_elapsed))
message("Inferred Louvain clusters: ", inferred_clusters)
message(sprintf("Truth-label marker benchmark completed in %.3f s", marker_elapsed))
message("Truth groups used for markers: ", marker_truth_groups)
message("Marker rows: ", ifelse(is.null(markers), 0L, nrow(markers)))
message(sprintf("Process peak RSS from /proc: %.1f MB", peak_kb / 1024))
