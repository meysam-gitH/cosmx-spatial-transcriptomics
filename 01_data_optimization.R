# Phase 1, Step 1: Initialize Environment and Manage Memory
#
# The original methodology starts from a pre-built ~14GB Seurat object and
# strips its ~8GB transcript-coordinates slot. This mini project builds the
# Seurat object itself from CosMx instrument flat files (exprMat, metadata,
# per-molecule tx file), then applies the identical extract -> fwrite ->
# NULL-out -> gc(full=TRUE) memory-management pattern to the per-molecule
# transcript table, which plays the same role as the 8GB slot in the source
# methodology (largest single object, one row per detected transcript).

library(Seurat)
library(data.table)
library(Matrix)

source("R/00_config.R")

## 1. Load CosMx flat files ---------------------------------------------
message("Reading expression matrix, metadata, and FOV positions...")
expr_dt <- fread(paths$exprMat)
meta_dt <- fread(paths$metadata)
fov_pos <- fread(paths$fov_positions)

## 2. Build cell identifiers and align expression <-> metadata ----------
expr_dt[, cell_key := paste(fov, cell_ID, sep = "_")]
meta_dt[, cell_key := paste(fov, cell_ID, sep = "_")]

# cell_ID == 0 denotes transcripts not assigned to a segmented cell
# (background); drop those rows, matching a real QC-passed cell set.
expr_dt <- expr_dt[cell_ID != 0]
meta_dt <- meta_dt[cell_ID != 0]

common_cells <- intersect(expr_dt$cell_key, meta_dt$cell_key)
expr_dt <- expr_dt[match(common_cells, cell_key)]
meta_dt <- meta_dt[match(common_cells, cell_key)]
stopifnot(identical(expr_dt$cell_key, meta_dt$cell_key))

meta_dt[, unique_fov := paste(slide_name, sprintf("F%02d", fov), sep = "_")]
setDF(meta_dt)
rownames(meta_dt) <- meta_dt$cell_key

## 3. Separate true gene targets from technical control probes ----------
gene_cols <- setdiff(names(expr_dt), c("fov", "cell_ID", "cell_key"))
is_control <- Reduce(`|`, lapply(control_probe_patterns, grepl, x = gene_cols))
true_genes <- gene_cols[!is_control]
control_probes <- gene_cols[is_control]
message(sprintf(
  "Expression matrix: %d true gene targets, %d control probes (%s)",
  length(true_genes), length(control_probes),
  paste(unique(sub("[0-9]+$", "", control_probes)), collapse = ", ")
))

counts_genes <- t(as.matrix(expr_dt[, ..true_genes]))
colnames(counts_genes) <- expr_dt$cell_key
counts_genes <- Matrix(counts_genes, sparse = TRUE)

counts_controls <- t(as.matrix(expr_dt[, ..control_probes]))
colnames(counts_controls) <- expr_dt$cell_key

## 4. Construct the Seurat object (analog of READ_RDS of the big object) -
seurat_obj <- CreateSeuratObject(
  counts = counts_genes,
  assay = "RNA",
  meta.data = meta_dt,
  project = slide_name
)

# Keep the negative-probe counts alongside the object (needed later for
# InSituType's negative-probe background term), not inside the RNA assay.
seurat_obj@misc$negprobe_counts <- counts_controls

message(sprintf(
  "Seurat object built: %d genes x %d cells across %d FOVs",
  nrow(seurat_obj), ncol(seurat_obj), length(unique(meta_dt$fov))
))
saveRDS(seurat_obj, paths$seurat_raw)

## 4b. Basic per-cell QC filter (min counts / min genes) -----------------
# Mirrors the "188,000 cells that passed quality control" framing of the
# source methodology: drop empty/near-empty cells before normalization so
# SCTransform's per-cell log(UMI) covariate is always finite.
n_before <- ncol(seurat_obj)
keep_cells <- colnames(seurat_obj)[
  seurat_obj$nCount_RNA >= 20 & seurat_obj$nFeature_RNA >= 5
]
seurat_obj <- subset(seurat_obj, cells = keep_cells)
seurat_obj@misc$negprobe_counts <- seurat_obj@misc$negprobe_counts[, keep_cells, drop = FALSE]
message(sprintf(
  "QC filter (nCount_RNA >= 20, nFeature_RNA >= 5): kept %d / %d cells",
  ncol(seurat_obj), n_before
))

## 5. Attach the per-molecule transcript table (the "8GB slot" analog) --
message("Reading per-molecule transcript table...")
tx_dt <- fread(paths$tx_file)
tx_dt <- tx_dt[cell_ID != 0]
tx_dt[, cell_key := paste(fov, cell_ID, sep = "_")]
tx_dt <- tx_dt[cell_key %in% colnames(seurat_obj)]

seurat_obj@misc$transcript.coordinates <- tx_dt

size_with_transcripts <- format(object.size(seurat_obj), units = "MB")
message(sprintf("Object size WITH transcript coordinates: %s", size_with_transcripts))

## 6. Extract + save transcript coordinates efficiently, then drop them --
message("Writing transcript coordinates with data.table::fwrite (multi-threaded)...")
transcripts <- seurat_obj@misc$transcript.coordinates
fwrite(transcripts, file = paths$transcript_coords_out, compress = "gzip")

seurat_obj@misc$transcript.coordinates <- list(NULL)  # overwrite, do not merely NULL-assign
rm(transcripts, tx_dt)
invisible(gc(full = TRUE))

size_without_transcripts <- format(object.size(seurat_obj), units = "MB")
message(sprintf("Object size AFTER dropping transcript coordinates: %s", size_without_transcripts))

## 7. Remove any pre-existing reductions / normalized assays ------------
# None exist yet on a freshly built object, but this mirrors the source
# methodology so the same script works unchanged against a pre-computed
# object loaded via readRDS().
seurat_obj@reductions <- list()
if ("RNA_normalized" %in% Assays(seurat_obj)) {
  seurat_obj[["RNA_normalized"]] <- NULL
}
invisible(gc(full = TRUE))

saveRDS(seurat_obj, paths$seurat_optimized)
message(sprintf("Optimized object saved to %s", paths$seurat_optimized))
