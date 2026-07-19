# Phase 2, Step 4: Dimension Reduction (PCA and UMAP)

library(Seurat)

source("R/00_config.R")

seurat_obj <- readRDS(paths$seurat_optimized)
DefaultAssay(seurat_obj) <- params$sct_new_assay_name

message(sprintf("Running PCA (npcs = %d) on assay '%s'...", params$npcs, params$sct_new_assay_name))
seurat_obj <- RunPCA(seurat_obj, npcs = params$npcs, verbose = FALSE)

message(sprintf(
  "Running UMAP on dims %d:%d with repulsion.strength = %g...",
  min(params$umap_dims), max(params$umap_dims), params$umap_repulsion_strength
))
seurat_obj <- RunUMAP(
  seurat_obj,
  dims = params$umap_dims,
  repulsion.strength = params$umap_repulsion_strength,
  verbose = FALSE
)

saveRDS(seurat_obj, paths$seurat_optimized)
message(sprintf("PCA + UMAP complete. Object re-saved to %s", paths$seurat_optimized))
