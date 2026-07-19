# Phase 2, Step 3: Data Normalization (SCTransform)
#
# Normalizes raw RNA counts to remove technical (sequencing/imaging depth)
# variation while preserving biological variance, replacing basic
# log-normalization.

library(Seurat)

source("R/00_config.R")

seurat_obj <- readRDS(paths$seurat_optimized)

message("Running SCTransform...")
seurat_obj <- SCTransform(
  seurat_obj,
  assay = params$sct_assay,
  new.assay.name = params$sct_new_assay_name,
  vst.flavor = "v2",
  verbose = TRUE
)

DefaultAssay(seurat_obj) <- params$sct_new_assay_name
saveRDS(seurat_obj, paths$seurat_optimized)
message(sprintf(
  "SCTransform complete. New assay '%s' set as default; object re-saved to %s",
  params$sct_new_assay_name, paths$seurat_optimized
))
