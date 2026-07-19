# Phase 2, Step 6: Unsupervised Cell Annotation (InSituType)
#
# Runs NanoString's InSituType package for unsupervised cell typing from
# raw counts + a negative-probe background term + immunofluorescence-derived
# cohorting, producing cluster labels (A, B, C, ...) and a posterior
# probability ("confidence") per cell.

library(Seurat)
library(InSituType)

source("R/00_config.R")

seurat_obj <- readRDS(paths$seurat_optimized)

## 1. Raw (untransformed) counts, cells in rows -------------------------
raw_counts <- t(as.matrix(GetAssayData(seurat_obj, assay = "RNA", layer = "counts")))

## 2. Mean negative-probe count per cell (background / false-positive term)
negprobe_counts <- seurat_obj@misc$negprobe_counts  # probes x cells
negmean <- Matrix::colMeans(negprobe_counts)
negmean <- negmean[rownames(raw_counts)]

## 3. Immunofluorescence-derived cohorting (DAPI excluded) --------------
if_mat <- as.matrix(seurat_obj@meta.data[rownames(raw_counts), params$if_columns])
cohort <- InSituType::fastCohorting(if_mat, gaussian_transform = TRUE)

## 4. Unsupervised InSituType clustering ---------------------------------
# insitutype() is the single most expensive step in the pipeline (searches
# n_clusts several times over, ~7 minutes on the mini dataset) -- cache the
# raw result on disk so a bug in the code below never forces a re-run of the
# clustering itself.
insitutype_cache <- file.path(paths$work_dir, "insitutype_result.rds")
if (file.exists(insitutype_cache)) {
  message(sprintf("Loading cached InSituType result from %s", insitutype_cache))
  sup <- readRDS(insitutype_cache)
} else {
  message(sprintf(
    "Running insitutype() searching %d-%d clusters over %d cells x %d genes...",
    min(params$insitutype_n_clusters_range), max(params$insitutype_n_clusters_range),
    nrow(raw_counts), ncol(raw_counts)
  ))
  sup <- InSituType::insitutype(
    x = raw_counts,
    neg = negmean,
    cohort = cohort,
    n_clusts = params$insitutype_n_clusters_range
  )
  saveRDS(sup, insitutype_cache)
}

## 5. Attach cluster assignment + confidence to the Seurat object -------
# sup$clust: named character vector of cell -> assigned cluster label.
# sup$prob: named numeric vector, each cell's posterior probability of its
# *assigned* cluster (already a per-cell confidence score, not a matrix).
seurat_obj$insitu_type_class <- factor(sup$clust[colnames(seurat_obj)])
seurat_obj$insitu_type_confidence <- sup$prob[colnames(seurat_obj)]

message("InSituType cluster sizes:")
print(table(seurat_obj$insitu_type_class))

saveRDS(seurat_obj, paths$seurat_annotated)
message(sprintf("Annotated object saved to %s", paths$seurat_annotated))

## 6. Flight Path Plot: visual confidence check --------------------------
# InSituType ships its own flightpath_plot(), which lays cells out by their
# posterior-probability profile (cells near a cluster center = high
# confidence) and colors by cluster.
fp_plot <- InSituType::flightpath_plot(insitutype_result = sup)
ggplot2::ggsave(file.path(paths$results_dir, "02_flightpath_plot.png"), fp_plot, width = 7, height = 6, dpi = 150)
unlink("NBClust-Plots", recursive = TRUE)  # InSituType's own auto-saved copy; we keep only results/
message("Saved results/02_flightpath_plot.png")
