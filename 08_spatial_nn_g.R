# Phase 3, Step 8: Spatial Clustering Calculation
#
# Builds a spatialTIME `mif` object from the Spatial List / Summary DF /
# Clinical DF produced in step 7, then computes the univariate Nearest
# Neighbor G statistic for each InSituType cell class across a radius range
# well inside the smallest FOV dimension.

library(tidyverse)
library(spatialTIME)

source("R/00_config.R")

inputs <- readRDS(file.path(paths$work_dir, "spatial_inputs.rds"))
spatial_list <- inputs$spatial_list
summary_df   <- inputs$summary_df
clinical_df  <- inputs$clinical_df
class_levels <- inputs$class_levels

## 1. Create the spatialTIME mif object -----------------------------------
# unique_fov plays the role of BOTH patient.ID and sample.ID: in this mini
# dataset each field of view is treated as its own sample/patient unit, as
# specified by the source methodology.
mif <- create_mif(
  clinical_data = clinical_df,
  sample_data   = summary_df,
  spatial_list  = spatial_list,
  patient_id    = "unique_fov",
  sample_id     = "unique_fov"
)

## 2. Univariate Nearest Neighbor G(r) -------------------------------------
message(sprintf(
  "Running NN_G over r = %d..%d (n=%d radii), %d permutations, %d classes, %d workers...",
  min(params$nn_g_r_range), max(params$nn_g_r_range), length(params$nn_g_r_range),
  params$nn_g_n_permutations, length(class_levels), params$nn_g_n_workers
))
mif <- NN_G(
  mif = mif,
  mnames = class_levels,
  r_range = params$nn_g_r_range,
  num_permutations = params$nn_g_n_permutations,
  edge_correction = "rs",
  workers = params$nn_g_n_workers,
  overwrite = TRUE,
  xloc = "x",
  yloc = "y"
)

nn_g_results <- mif$derived$univariate_NN
write_csv(nn_g_results, file.path(paths$results_dir, "nn_g_results.csv"))
saveRDS(mif, file.path(paths$work_dir, "mif_object.rds"))

message(sprintf("Nearest Neighbor G complete: %d rows written to results/nn_g_results.csv", nrow(nn_g_results)))
