#!/usr/bin/env Rscript

# ============================================================
# Reconstruct the exact GSE234797 WT + vehicle-rd10 benchmark
# using RAW integer counts from the two H5 files, while using
# GEO processed CellIDs/MetaData/Genes only to identify:
#   - retained processed cells
#   - retained processed genes
#
# TMB cells are excluded.
#
# Required final counts:
#   WT   = 4,478
#   rd10 = 5,952
#   total = 10,430
#
# This preparation is NOT timed.
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
  v <- suppressWarnings(as.integer(x))
  if (length(v) != 1L || is.na(v)) stop(name, " must be an integer.", call. = FALSE)
  v
}

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg, call. = FALSE)
  }
}

read_lines_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  x <- readLines(con, warn = FALSE)
  x <- trimws(x)
  x[nzchar(x)]
}

read_genes <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)

  dat <- utils::read.table(
    con,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fill = TRUE
  )

  # Use the last non-empty column as the gene symbol/name.
  for (j in rev(seq_len(ncol(dat)))) {
    x <- trimws(as.character(dat[[j]]))
    if (all(nzchar(x))) return(unique(x))
  }

  stop("Unable to derive processed gene names.", call. = FALSE)
}

extract_gene_expression <- function(x, label) {
  if (inherits(x, "Matrix") || is.matrix(x)) return(x)

  if (is.list(x)) {
    if ("Gene Expression" %in% names(x)) {
      return(x[["Gene Expression"]])
    }

    candidates <- which(vapply(
      x,
      function(z) inherits(z, "Matrix") || is.matrix(z),
      logical(1)
    ))

    if (length(candidates) == 1L) return(x[[candidates]])

    stop(
      label,
      ": unable to identify the Gene Expression matrix in H5.",
      call. = FALSE
    )
  }

  stop(label, ": unsupported H5 result.", call. = FALSE)
}

# Try to extract a standard 10x barcode such as
# AAACCCAAGACATACA-1 from a processed identifier.
extract_10x_barcode <- function(x) {

  x <- as.character(x)

  # Processed CellIDs in GSE234797 look like:
  # Wild_Type_non_treated-AAACCCAAGGTTGACG
  # Rd10_Female_vehicle-AAACCCAAGGTTGACG
  #
  # Raw Cell Ranger H5 barcodes normally look like:
  # AAACCCAAGGTTGACG-1
  #
  # Normalize both forms to the core 16-base barcode.

  out <- rep(NA_character_, length(x))

  # 1. Capture a 16-base A/C/G/T barcode at the end, optionally followed by -1, -2, etc.
  m <- regexec(
    "([ACGT]{16})(?:-[0-9]+)?$",
    x,
    perl = TRUE
  )
  hits <- regmatches(x, m)

  ok <- lengths(hits) >= 2L
  out[ok] <- vapply(hits[ok], function(z) z[[2L]], character(1))

  # 2. Fallback: strip sample prefix and Cell Ranger numeric suffix.
  need <- which(is.na(out))
  if (length(need)) {
    y <- x[need]
    y <- sub("^.*-", "", y)
    y <- sub("-[0-9]+$", "", y)
    y <- trimws(y)

    looks <- grepl("^[ACGT]{16}$", y)

    tmp <- rep(NA_character_, length(y))
    tmp[looks] <- y[looks]
    out[need] <- tmp
  }

  out
}
classify_rows <- function(meta, cell_ids) {

  # Search all metadata fields plus the processed cell ID.
  row_text <- apply(
    cbind(cell_id = cell_ids, meta),
    1L,
    function(z) paste(z, collapse = " | ")
  )

  wt <- grepl(
    "GSM7474906|Wild[_ -]?Type|Wild Type|non[_ -]?treated|(^|[^A-Za-z])WT([^A-Za-z]|$)",
    row_text,
    ignore.case = TRUE,
    perl = TRUE
  )

  vehicle <- grepl(
    "GSM7474907|vehicle",
    row_text,
    ignore.case = TRUE,
    perl = TRUE
  )

  tmb <- grepl(
    "GSM7474908|TMB",
    row_text,
    ignore.case = TRUE,
    perl = TRUE
  )

  out <- rep(NA_character_, length(row_text))

  out[wt & !vehicle & !tmb] <- "WT"
  out[vehicle & !tmb] <- "rd10"
  out[tmb] <- "TMB"

  out
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))

raw_dir <- arg_value(
  args,
  "raw-dir",
  "/home/umair/Downloads/CoTRA/data/scRNA"
)

processed_dir <- arg_value(
  args,
  "processed-dir",
  "/home/umair/Downloads/CoTRA/data/scRNA/GSE234797_processed"
)

outdir <- arg_value(
  args,
  "outdir",
  "/home/umair/Downloads/CoTRA/benchmark/real_scRNA_exact"
)

seed <- as_int(arg_value(args, "seed", "1234"), "--seed")

wt_h5 <- file.path(
  raw_dir,
  "GSM7474906_Wild_Type_non_treated_feature_bc_matrix.h5"
)

rd_h5 <- file.path(
  raw_dir,
  "GSM7474907_Rd10_Female_vehicle_feature_bc_matrix.h5"
)

cells_file <- file.path(
  processed_dir,
  "GSE234797_CellIDs.txt.gz"
)

genes_file <- file.path(
  processed_dir,
  "GSE234797_Genes.txt.gz"
)

meta_file <- file.path(
  processed_dir,
  "GSE234797_MetaData.csv.gz"
)

needed <- c(wt_h5, rd_h5, cells_file, genes_file, meta_file)
missing <- needed[!file.exists(needed)]

if (length(missing)) {
  stop(
    "Missing required files:\n",
    paste(missing, collapse = "\n"),
    call. = FALSE
  )
}

require_pkg("Seurat")
require_pkg("Matrix")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load processed CellIDs and metadata.
# ------------------------------------------------------------
message("Reading processed CellIDs...")

processed_ids <- read_lines_gz(cells_file)

message("Processed CellIDs: ", length(processed_ids))

message("Reading processed metadata...")

meta_con <- gzfile(meta_file, open = "rt")
meta <- utils::read.csv(
  meta_con,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
close(meta_con)

message(
  "Metadata: ",
  nrow(meta),
  " rows x ",
  ncol(meta),
  " columns"
)

message(
  "Metadata columns: ",
  paste(names(meta), collapse = " | ")
)

if (nrow(meta) != length(processed_ids)) {
  stop(
    "Metadata rows (", nrow(meta),
    ") do not equal CellIDs (", length(processed_ids), "). ",
    "The Series files cannot be safely aligned by row order.",
    call. = FALSE
  )
}

group <- classify_rows(meta, processed_ids)

message("Processed cell classification:")
message("  WT: ", sum(group == "WT", na.rm = TRUE))
message("  vehicle rd10: ", sum(group == "rd10", na.rm = TRUE))
message("  TMB: ", sum(group == "TMB", na.rm = TRUE))
message("  unclassified: ", sum(is.na(group)))

expected_wt <- 4478L
expected_rd <- 5952L

if (
  sum(group == "WT", na.rm = TRUE) != expected_wt ||
  sum(group == "rd10", na.rm = TRUE) != expected_rd
) {
  message("")
  message("Unique values from likely sample/group columns:")

  for (nm in names(meta)) {
    u <- unique(as.character(meta[[nm]]))
    u <- u[!is.na(u)]

    if (
      any(grepl(
        "WT|wild|vehicle|TMB|rd10|GSM747490",
        u,
        ignore.case = TRUE
      ))
    ) {
      message("  ", nm, ": ", paste(head(u, 30L), collapse = " | "))
    }
  }

  stop(
    "Could not recover the expected 4,478 WT and 5,952 vehicle-rd10 ",
    "cells from metadata. No benchmark files were written.",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# 2. Parse processed 10x barcodes.
# ------------------------------------------------------------
processed_barcode <- extract_10x_barcode(processed_ids)

target_wt_ids <- processed_ids[group == "WT"]
target_rd_ids <- processed_ids[group == "rd10"]

target_wt_bc <- processed_barcode[group == "WT"]
target_rd_bc <- processed_barcode[group == "rd10"]

message("Parsed standard 10x barcode from:")
message(
  "  WT: ",
  sum(!is.na(target_wt_bc)),
  "/",
  length(target_wt_bc)
)
message(
  "  rd10: ",
  sum(!is.na(target_rd_bc)),
  "/",
  length(target_rd_bc)
)

# ------------------------------------------------------------
# 3. Load the two raw H5 count matrices.
# ------------------------------------------------------------
message("Reading raw WT H5...")

wt <- extract_gene_expression(
  Seurat::Read10X_h5(
    wt_h5,
    use.names = TRUE,
    unique.features = TRUE
  ),
  "WT"
)

message("Reading raw vehicle-rd10 H5...")

rd <- extract_gene_expression(
  Seurat::Read10X_h5(
    rd_h5,
    use.names = TRUE,
    unique.features = TRUE
  ),
  "rd10"
)

if (!inherits(wt, "Matrix")) wt <- Matrix::Matrix(wt, sparse = TRUE)
if (!inherits(rd, "Matrix")) rd <- Matrix::Matrix(rd, sparse = TRUE)

message(
  "Raw WT: ",
  nrow(wt),
  " genes x ",
  ncol(wt),
  " cells"
)

message(
  "Raw rd10: ",
  nrow(rd),
  " genes x ",
  ncol(rd),
  " cells"
)

# ------------------------------------------------------------
# 4. Match processed IDs back to raw H5 barcodes.
# ------------------------------------------------------------
raw_wt_bc <- colnames(wt)
raw_rd_bc <- colnames(rd)

raw_wt_core <- extract_10x_barcode(raw_wt_bc)
raw_rd_core <- extract_10x_barcode(raw_rd_bc)

if (anyNA(raw_wt_core) || anyNA(raw_rd_core)) {
  stop(
    "Could not normalize every raw H5 barcode to a 16-base barcode.",
    call. = FALSE
  )
}

if (anyDuplicated(raw_wt_core)) {
  stop("Normalized WT raw barcodes are not unique.", call. = FALSE)
}

if (anyDuplicated(raw_rd_core)) {
  stop("Normalized rd10 raw barcodes are not unique.", call. = FALSE)
}

match_one_group <- function(
  processed_ids,
  processed_bc,
  raw_bc,
  raw_core,
  label
) {

  idx <- rep(NA_integer_, length(processed_ids))

  # Main mapping: compare normalized 16-base barcodes.
  good <- !is.na(processed_bc)
  idx[good] <- match(processed_bc[good], raw_core)

  # Exact raw barcode fallback, if ever needed.
  need <- is.na(idx)
  if (any(need)) {
    idx[need] <- match(processed_ids[need], raw_bc)
  }

  matched <- sum(!is.na(idx))

  message(
    label,
    " processed-to-raw matches: ",
    matched,
    "/",
    length(idx)
  )

  if (matched != length(idx)) {

    bad <- which(is.na(idx))

    message(
      label,
      " first unmatched processed IDs:"
    )

    message(
      paste(
        head(processed_ids[bad], 20L),
        collapse = "\n"
      )
    )

    stop(
      label,
      ": not every processed cell could be matched to its raw H5 barcode. ",
      "No benchmark files were written.",
      call. = FALSE
    )
  }

  if (anyDuplicated(idx)) {
    stop(
      label,
      ": duplicate raw-cell mappings detected.",
      call. = FALSE
    )
  }

  idx
}

wt_idx <- match_one_group(
  target_wt_ids,
  target_wt_bc,
  raw_wt_bc,
  raw_wt_core,
  "WT"
)

rd_idx <- match_one_group(
  target_rd_ids,
  target_rd_bc,
  raw_rd_bc,
  raw_rd_core,
  "rd10"
)

# ------------------------------------------------------------
# 5. Restrict raw counts to the processed gene set.
# ------------------------------------------------------------
message("Reading processed gene list...")

processed_genes <- read_genes(genes_file)

message(
  "Processed gene entries: ",
  length(processed_genes)
)

common_genes <- Reduce(
  intersect,
  list(
    processed_genes,
    rownames(wt),
    rownames(rd)
  )
)

message(
  "Processed genes found in both raw H5 matrices: ",
  length(common_genes)
)

if (length(common_genes) < 10000L) {
  stop(
    "Fewer than 10,000 processed genes matched the raw matrices. ",
    "Gene identifier format should be checked.",
    call. = FALSE
  )
}

wt_exact <- wt[
  common_genes,
  wt_idx,
  drop = FALSE
]

rd_exact <- rd[
  common_genes,
  rd_idx,
  drop = FALSE
]

# Prefix the raw barcodes to avoid duplicated 10x barcode names.
colnames(wt_exact) <- paste0(
  "WT_",
  colnames(wt_exact)
)

colnames(rd_exact) <- paste0(
  "rd10_",
  colnames(rd_exact)
)

counts <- cbind(
  wt_exact,
  rd_exact
)

metadata <- data.frame(
  cell = colnames(counts),
  condition = c(
    rep("WT", ncol(wt_exact)),
    rep("rd10", ncol(rd_exact))
  ),
  processed_cell_id = c(
    target_wt_ids,
    target_rd_ids
  ),
  stringsAsFactors = FALSE
)

rownames(metadata) <- metadata$cell

# ------------------------------------------------------------
# 6. Final hard validation.
# ------------------------------------------------------------
if (ncol(wt_exact) != expected_wt) {
  stop("Final WT count is not 4,478.", call. = FALSE)
}

if (ncol(rd_exact) != expected_rd) {
  stop("Final rd10 count is not 5,952.", call. = FALSE)
}

if (ncol(counts) != 10430L) {
  stop("Final combined cell count is not 10,430.", call. = FALSE)
}

if (length(counts@x)) {
  if (any(counts@x < 0)) {
    stop("Raw subset contains negative values.", call. = FALSE)
  }

  if (!all(abs(counts@x - round(counts@x)) < 1e-8)) {
    stop("Raw subset is not integer-valued.", call. = FALSE)
  }
}

message("")
message("Exact raw-count benchmark dataset reconstructed:")
message("  genes: ", nrow(counts))
message("  WT: ", ncol(wt_exact))
message("  rd10 vehicle: ", ncol(rd_exact))
message("  total: ", ncol(counts))
message("  values: non-negative integer counts")

# ------------------------------------------------------------
# 7. Create fixed nested stratified subsets.
# ------------------------------------------------------------
set.seed(seed)

wt_cells <- colnames(counts)[metadata$condition == "WT"]
rd_cells <- colnames(counts)[metadata$condition == "rd10"]

wt_order <- sample(
  wt_cells,
  length(wt_cells),
  replace = FALSE
)

rd_order <- sample(
  rd_cells,
  length(rd_cells),
  replace = FALSE
)

wt_fraction <- expected_wt / 10430

targets <- c(
  2500L,
  5000L,
  7500L,
  10000L,
  10430L
)

manifest <- list()

for (n in targets) {

  if (n == 10430L) {

    selected_wt <- wt_order
    selected_rd <- rd_order

  } else {

    nw <- as.integer(round(n * wt_fraction))
    nr <- n - nw

    nw <- min(nw, length(wt_order))
    nr <- min(nr, length(rd_order))

    selected_wt <- wt_order[seq_len(nw)]
    selected_rd <- rd_order[seq_len(nr)]
  }

  selected <- c(
    selected_wt,
    selected_rd
  )

  sub_counts <- counts[
    ,
    selected,
    drop = FALSE
  ]

  sub_meta <- metadata[
    selected,
    ,
    drop = FALSE
  ]

  payload <- list(
    counts = sub_counts,
    metadata = sub_meta,
    source = paste(
      "raw H5 counts subset using",
      "GSE234797 processed CellIDs/MetaData/Genes"
    ),
    seed = seed,
    n_cells = ncol(sub_counts),
    n_wt = sum(sub_meta$condition == "WT"),
    n_rd10 = sum(sub_meta$condition == "rd10")
  )

  path <- file.path(
    outdir,
    sprintf(
      "counts_real_%06dc.rds",
      ncol(sub_counts)
    )
  )

  saveRDS(
    payload,
    path,
    compress = "xz"
  )

  manifest[[length(manifest) + 1L]] <- data.frame(
    n_genes = nrow(sub_counts),
    n_cells = ncol(sub_counts),
    n_WT = sum(sub_meta$condition == "WT"),
    n_rd10 = sum(sub_meta$condition == "rd10"),
    WT_fraction = mean(sub_meta$condition == "WT"),
    rd10_fraction = mean(sub_meta$condition == "rd10"),
    nonzero_entries = length(sub_counts@x),
    file = normalizePath(
      path,
      mustWork = FALSE
    ),
    stringsAsFactors = FALSE
  )

  message(
    "Saved ",
    ncol(sub_counts),
    " cells: WT=",
    sum(sub_meta$condition == "WT"),
    ", rd10=",
    sum(sub_meta$condition == "rd10")
  )

  rm(
    sub_counts,
    sub_meta,
    payload
  )

  invisible(gc())
}

manifest_df <- do.call(
  rbind,
  manifest
)

utils::write.csv(
  manifest_df,
  file.path(
    outdir,
    "real_scRNA_exact_manifest.csv"
  ),
  row.names = FALSE
)

prep <- data.frame(
  item = c(
    "GEO_accession",
    "WT_raw_H5",
    "rd10_vehicle_raw_H5",
    "TMB_included",
    "processed_CellIDs_used_for_selection",
    "processed_MetaData_used_for_selection",
    "processed_Genes_used_for_selection",
    "final_genes",
    "final_WT_cells",
    "final_rd10_cells",
    "final_total_cells",
    "data_values",
    "seed"
  ),
  value = c(
    "GSE234797",
    normalizePath(wt_h5),
    normalizePath(rd_h5),
    "NO",
    normalizePath(cells_file),
    normalizePath(meta_file),
    normalizePath(genes_file),
    nrow(counts),
    expected_wt,
    expected_rd,
    10430L,
    "raw non-negative integer H5 counts",
    seed
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  prep,
  file.path(
    outdir,
    "real_scRNA_exact_preparation_metadata.csv"
  ),
  row.names = FALSE
)

message("")
message("Preparation completed successfully.")
message(
  "Output directory: ",
  normalizePath(outdir)
)
