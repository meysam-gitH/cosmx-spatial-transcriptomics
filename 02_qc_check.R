# Phase 1, Step 2: Quality Control (QC) Check on RNA Assay
#
# Confirms the RNA count matrix that will be normalized contains only true
# gene targets -- no negative probes or other technical control probes.

library(Seurat)

source("R/00_config.R")

seurat_obj <- readRDS(paths$seurat_optimized)

rna_counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
gene_names <- rownames(rna_counts)

hits <- lapply(control_probe_patterns, function(pat) grep(pat, gene_names, value = TRUE))
names(hits) <- control_probe_patterns
n_hits <- sum(lengths(hits))

if (n_hits > 0) {
  stop(sprintf(
    "QC FAILED: %d control probe(s) found in RNA assay row names: %s",
    n_hits, paste(unlist(hits), collapse = ", ")
  ))
}

message(sprintf(
  "QC PASSED: RNA assay contains %d row names, none matching control probe patterns (%s).",
  length(gene_names), paste(control_probe_patterns, collapse = ", ")
))
message(sprintf(
  "Control probes were separated out during object construction and are retained in seurat_obj@misc$negprobe_counts (%d probes) for later use as the InSituType background term.",
  nrow(seurat_obj@misc$negprobe_counts)
))
