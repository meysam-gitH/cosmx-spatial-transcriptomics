# Shared paths and tunable parameters for the CosMx SMI pipeline.
#
# Data note: NanoString does not publicly host the Human Frontal Cortex CosMx
# WTx dataset referenced in the original methodology. This pipeline instead
# runs against a genuinely public CosMx SMI sample (NanoString "Lung9_Rep1",
# 2021 public FFPE lung release) so every step below operates on real
# instrument output, subset to a 2x2 block of fields of view (FOVs 1, 2, 5, 6;
# ~18,500 cells) to keep the mini-project fast to run end to end. The
# methodology -- memory optimization, SCTransform, PCA/UMAP/Louvain,
# InSituType annotation, spatialTIME Nearest Neighbor G -- is identical
# regardless of tissue of origin.

slide_name <- "Lung9_Rep1"

paths <- list(
  raw_dir       = "data/Lung9_Rep1/Lung9_Rep1-Flat_files_and_images",
  mini_dir      = "data/mini",
  exprMat       = "data/mini/mini_exprMat.csv",
  metadata      = "data/mini/mini_metadata.csv",
  fov_positions = "data/mini/mini_fov_positions.csv",
  tx_file       = "data/mini/mini_tx_file.csv",
  work_dir      = "data/work",
  results_dir   = "results",
  seurat_raw    = "data/work/seurat_raw.rds",
  seurat_optimized = "data/work/seurat_optimized.rds",
  seurat_annotated = "data/work/seurat_annotated.rds",
  transcript_coords_out = "data/work/transcript_locations.csv.gz"
)

dir.create(paths$work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)

params <- list(
  # Phase 2: normalization / dimension reduction / clustering
  sct_assay          = "RNA",
  sct_new_assay_name = "SCT",
  npcs               = 50,
  umap_dims          = 1:30,
  umap_repulsion_strength = 10,
  louvain_algorithm  = 1,      # 1 = original Louvain
  louvain_resolution = 0.3,

  # Phase 2: InSituType
  # Full CosMx WTx panel searches ~10-20 clusters; this mini 960-plex, ~18.5k
  # cell subset searches a smaller range so cluster search remains stable.
  insitutype_n_clusters_range = 4:8,
  # Immunofluorescence channels to hand to InSituType, DAPI excluded
  # (DAPI is a nuclear counterstain used for segmentation, not phenotyping).
  if_columns = c("Mean.MembraneStain", "Mean.PanCK", "Mean.CD45", "Mean.CD3"),

  # Phase 3: spatial statistics
  # Local FOV extent here is ~5430 x 3630 px, so R.range is kept well inside
  # the smaller (3630 px) dimension, matching the methodology's edge-
  # correction guidance.
  nn_g_r_range      = seq(0, 200, by = 10),
  nn_g_n_permutations = 50,
  nn_g_n_workers     = max(1, parallel::detectCores() - 1)
)

control_probe_patterns <- c("^NegPrb", "^SystemControl", "^FalseCode")

# Fixed-order, CVD-validated categorical palette (blue, green, magenta,
# yellow, aqua, orange, violet, red) used consistently for cell-class colors
# across all plots. Never cycled or reassigned by a subset/filter.
categorical_palette <- c(
  "#2a78d6", "#008300", "#e87ba4", "#eda100",
  "#1baf7a", "#eb6834", "#4a3aa7", "#e34948"
)
