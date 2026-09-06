#!/usr/bin/env Rscript

# ============================================================
# CoTRA real retinal scRNA-seq computational benchmark
#
# Uses 50 PCs calculated, but PCs 1:7 downstream for:
#   UMAP
#   FindNeighbors
#   Louvain clustering
#
# This matches the user's CoTRA case-study PC choice.
#
# Marker analysis uses the ACTUAL inferred Louvain clusters,
# not synthetic truth labels.
#
# Disk RDS loading is outside the workflow timer.
# GNU /usr/bin/time -v should be used externally for peak RAM.
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
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package not installed: ", pkg, call. = FALSE)
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
if (is.null(input) || !file.exists(input)) {
  stop("--input must point to an existing prepared real-scRNA RDS.", call. = FALSE)
}

rep_id <- as_int(arg_value(args, "rep", "1"), "--rep")
seed <- as_int(arg_value(args, "seed", "1234"), "--seed")
pcs_downstream <- as_int(arg_value(args, "pcs", "7"), "--pcs")
outdir <- arg_value(args, "outdir", "benchmark_results/scrna_real")

if (pcs_downstream < 2L || pcs_downstream > 50L) {
  stop("--pcs must be between 2 and 50.", call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Loading prepared real dataset: ", input)
payload <- readRDS(input)

if (!is.list(payload) || is.null(payload$counts) || is.null(payload$metadata)) {
  stop("Prepared RDS must contain counts and metadata.", call. = FALSE)
}

counts <- payload$counts
metadata <- payload$metadata

if (!inherits(counts, "Matrix")) counts <- Matrix::Matrix(counts, sparse = TRUE)

if (!identical(colnames(counts), rownames(metadata))) {
  metadata <- metadata[colnames(counts), , drop = FALSE]
}

n_genes <- nrow(counts)
n_cells <- ncol(counts)
n_wt <- sum(metadata$condition == "WT")
n_rd10 <- sum(metadata$condition == "rd10")

id <- sprintf("scrna_real_%06dc_rep%02d", n_cells, rep_id)

message("Benchmark: ", id)
message("Dimensions: ", n_genes, " genes x ", n_cells, " cells")
message("WT cells: ", n_wt)
message("rd10 cells: ", n_rd10)
message("PCs calculated: 50")
message("PCs used downstream: 1:", pcs_downstream)
message("Target full-dataset validation: 16 clusters (0-15)")
message("Thread mode: single-thread constrained")

logenv <- new.env(parent = emptyenv())
logenv$rows <- list()

status <- "success"
error_message <- NA_character_
seu <- NULL
markers <- NULL
n_clusters <- NA_integer_

core_start <- proc.time()[["elapsed"]]

tryCatch({

  seu <- measure_step("create_seurat_and_qc_metrics", {
    obj <- Seurat::CreateSeuratObject(
      counts = counts,
      assay = "RNA",
      project = "WT_rd10_retina",
      meta.data = metadata
    )

    obj@meta.data$percent.mito <- Seurat::PercentageFeatureSet(
      obj,
      assay = "RNA",
      pattern = "^mt-"
    )

    obj@meta.data$percent.ribo <- Seurat::PercentageFeatureSet(
      obj,
      assay = "RNA",
      pattern = "^Rp[ls]"
    )

    obj
  }, logenv)

  seu <- measure_step("normalize_lognormalize", {
    Seurat::NormalizeData(
      seu,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  }, logenv)

  seu <- measure_step("find_variable_features", {
    Seurat::FindVariableFeatures(
      seu,
      assay = "RNA",
      selection.method = "vst",
      nfeatures = min(2000L, nrow(seu)),
      verbose = FALSE
    )
  }, logenv)

  features <- Seurat::VariableFeatures(seu, assay = "RNA")
  features <- unique(features[features %in% rownames(seu[["RNA"]])])

  if (length(features) < 100L) stop("Too few variable features for PCA.", call. = FALSE)

  seu <- measure_step("scale_data", {
    Seurat::ScaleData(
      seu,
      assay = "RNA",
      features = features,
      verbose = FALSE
    )
  }, logenv)

  npcs <- min(50L, ncol(seu) - 1L, length(features) - 1L)

  if (npcs < pcs_downstream) {
    stop("Not enough PCs available for requested downstream PC count.", call. = FALSE)
  }

  seu <- measure_step("pca", {
    Seurat::RunPCA(
      seu,
      assay = "RNA",
      features = features,
      npcs = npcs,
      verbose = FALSE
    )
  }, logenv)

  dims <- seq_len(pcs_downstream)

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

  n_clusters <- length(unique(as.character(Seurat::Idents(seu))))

}, error = function(e) {
  status <<- "failed"
  error_message <<- conditionMessage(e)
})

core_elapsed <- proc.time()[["elapsed"]] - core_start

marker_elapsed <- NA_real_

if (identical(status, "success")) {

  if (n_clusters < 2L) {
    status <- "failed"
    error_message <- "Clustering produced fewer than two clusters."
  } else {
    tryCatch({

      Seurat::DefaultAssay(seu) <- "RNA"
      set.seed(seed)

      marker_start <- proc.time()[["elapsed"]]

      markers <- measure_step("find_all_markers_inferred_clusters", {
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
      }, logenv, category = "markers")

      marker_elapsed <- proc.time()[["elapsed"]] - marker_start

    }, error = function(e) {
      status <<- "failed"
      error_message <<- paste0("Marker analysis failed: ", conditionMessage(e))
    })
  }
}

total_elapsed <- core_elapsed + ifelse(is.finite(marker_elapsed), marker_elapsed, 0)

steps <- if (length(logenv$rows)) {
  do.call(rbind, logenv$rows)
} else {
  data.frame(category = character(), step = character(), elapsed_seconds = numeric())
}

peak_kb <- proc_status_kb("VmHWM")

cluster_table <- if (!is.null(seu) && n_clusters >= 1L) {
  tmp <- as.data.frame(table(
    cluster = as.character(Seurat::Idents(seu)),
    condition = seu$condition
  ))
  tmp[tmp$Freq > 0, , drop = FALSE]
} else {
  data.frame()
}

summary <- data.frame(
  benchmark_id = id,
  workflow = "scrna_real_retina",
  replicate = rep_id,
  seed = seed,
  n_genes_input = n_genes,
  n_cells = n_cells,
  n_WT = n_wt,
  n_rd10 = n_rd10,
  pcs_calculated = npcs,
  pcs_used_downstream = pcs_downstream,
  core_workflow_elapsed_seconds = as.numeric(core_elapsed),
  marker_elapsed_seconds = as.numeric(marker_elapsed),
  total_core_plus_marker_seconds = as.numeric(total_elapsed),
  inferred_louvain_clusters = n_clusters,
  full_dataset_expected_clusters = if (n_cells == (4478L + 5952L)) 16L else NA_integer_,
  full_dataset_cluster_match = if (n_cells == (4478L + 5952L)) n_clusters == 16L else NA,
  marker_rows = if (!is.null(markers)) nrow(markers) else NA_integer_,
  significant_marker_rows_padj_0_05 =
    if (!is.null(markers) && "p_val_adj" %in% names(markers)) {
      sum(is.finite(markers$p_val_adj) & markers$p_val_adj <= 0.05, na.rm = TRUE)
    } else {
      NA_integer_
    },
  hvg_used = if (!is.null(seu)) length(Seurat::VariableFeatures(seu, assay = "RNA")) else NA_integer_,
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

utils::write.csv(
  summary,
  file.path(outdir, paste0("overall_", id, ".csv")),
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  steps,
  file.path(outdir, paste0("steps_", id, ".csv")),
  row.names = FALSE,
  na = ""
)

if (nrow(cluster_table) > 0L) {
  utils::write.csv(
    cluster_table,
    file.path(outdir, paste0("cluster_condition_", id, ".csv")),
    row.names = FALSE,
    na = ""
  )
}

if (!is.null(markers)) {
  utils::write.csv(
    markers,
    file.path(outdir, paste0("markers_", id, ".csv")),
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
message("Inferred Louvain clusters: ", n_clusters)
message(sprintf("Marker analysis completed in %.3f s", marker_elapsed))
message("Marker rows: ", ifelse(is.null(markers), 0L, nrow(markers)))

if (n_cells == 10430L) {
  message(
    "Full-dataset 16-cluster validation: ",
    ifelse(n_clusters == 16L, "MATCH", "NO MATCH")
  )
}

message(sprintf("Process peak RSS from /proc: %.1f MB", peak_kb / 1024))
