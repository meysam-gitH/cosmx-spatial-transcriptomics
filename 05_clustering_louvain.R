# Phase 2, Step 5: Initial Clustering (Louvain Method)
#
# Preliminary, unsupervised cluster IDs from the PCA graph. These are later
# superseded by InSituType's expression + immunofluorescence-based cell
# typing (step 6), but are useful as a sanity check and for the UMAP plot.

library(Seurat)

source("R/00_config.R")

seurat_obj <- readRDS(paths$seurat_optimized)
DefaultAssay(seurat_obj) <- params$sct_new_assay_name

message(sprintf("Running FindNeighbors on dims 1:%d...", max(params$umap_dims)))
seurat_obj <- FindNeighbors(seurat_obj, dims = params$umap_dims, verbose = FALSE)

message(sprintf(
  "Running FindClusters (algorithm = %d [Louvain], resolution = %g)...",
  params$louvain_algorithm, params$louvain_resolution
))
seurat_obj <- FindClusters(
  seurat_obj,
  algorithm = params$louvain_algorithm,
  resolution = params$louvain_resolution,
  verbose = FALSE
)

n_clusters <- length(unique(Idents(seurat_obj)))
message(sprintf("Louvain clustering complete: %d clusters at resolution %g.", n_clusters, params$louvain_resolution))

saveRDS(seurat_obj, paths$seurat_optimized)

p <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, raster = FALSE) +
  ggplot2::ggtitle(sprintf("Louvain clusters (resolution = %g)", params$louvain_resolution))
ggplot2::ggsave(file.path(paths$results_dir, "01_umap_louvain_clusters.png"), p, width = 7, height = 6, dpi = 150)
message("Saved results/01_umap_louvain_clusters.png")
