# Phase 3, Step 7: Spatial Data Structuring
#
# Converts Seurat metadata (cell locations + InSituType classes) into the
# three linked components spatialTIME needs: a Spatial List (one wide
# points data frame per FOV, one indicator column per cell class), a
# Summary DF (per-FOV cell-type counts/proportions), and a Clinical DF
# (per-FOV/sample-level covariates).

library(Seurat)
library(tidyverse)

source("R/00_config.R")

seurat_obj <- readRDS(paths$seurat_annotated)

class_levels <- sort(levels(droplevels(seurat_obj$insitu_type_class)))
message(sprintf("Cell classes for spatial analysis: %s", paste(class_levels, collapse = ", ")))

## 1. Spatial data frame (points), one row per cell ----------------------
spatial_df <- seurat_obj@meta.data %>%
  as_tibble(rownames = "cell_key") %>%
  select(unique_fov, CenterX_local_px, CenterY_local_px, insitu_type_class) %>%
  rename(x = CenterX_local_px, y = CenterY_local_px) %>%
  mutate(positive = 1)

## 2. Spread classes into one 0/1 indicator column per class -------------
spatial_wide <- spatial_df %>%
  pivot_wider(
    names_from = insitu_type_class,
    values_from = positive,
    values_fill = 0
  )
# Ensure every class has a column even if absent from a given FOV subset.
for (cl in class_levels) {
  if (!cl %in% names(spatial_wide)) spatial_wide[[cl]] <- 0
}

## 3. Spatial List: split by unique_fov -----------------------------------
spatial_list <- split(spatial_wide, f = spatial_wide$unique_fov)

## 4. Summary DF: per-FOV cell-type counts and proportions ---------------
summary_df <- spatial_list %>%
  lapply(function(df) {
    df %>%
      summarise(
        unique_fov = unique(unique_fov),
        total_cells = n(),
        across(all_of(class_levels), ~ sum(.x), .names = "counts_{.col}"),
      ) %>%
      mutate(across(starts_with("counts_"), ~ .x / total_cells, .names = "prop_{.col}"))
  }) %>%
  bind_rows()

## 5. Clinical DF: one row per FOV with sample-level QC metrics ----------
clinical_df <- seurat_obj@meta.data %>%
  as_tibble() %>%
  group_by(unique_fov) %>%
  summarise(
    fov = unique(fov),
    slide = slide_name,
    n_cells = n(),
    mean_transcripts_per_cell = mean(nCount_RNA),
    mean_genes_per_cell = mean(nFeature_RNA),
    .groups = "drop"
  )

saveRDS(
  list(spatial_list = spatial_list, summary_df = summary_df, clinical_df = clinical_df, class_levels = class_levels),
  file.path(paths$work_dir, "spatial_inputs.rds")
)

message(sprintf(
  "Spatial data prep complete: %d FOVs, %d cell classes.",
  length(spatial_list), length(class_levels)
))
print(summary_df)
