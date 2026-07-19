# Runs the full CosMx SMI pipeline end to end, in methodology order.
# Usage: Rscript R/run_pipeline.R   (run from the cosmx-spatial-pipeline/ directory)

steps <- c(
  "R/01_data_optimization.R",
  "R/02_qc_check.R",
  "R/03_normalization_sctransform.R",
  "R/04_dimension_reduction.R",
  "R/05_clustering_louvain.R",
  "R/06_insitutype_annotation.R",
  "R/07_spatial_data_prep.R",
  "R/08_spatial_nn_g.R",
  "R/09_visualization.R"
)

for (step in steps) {
  cat(sprintf("\n==================== %s ====================\n", step))
  source(step)
}

cat("\nPipeline complete. See results/ for figures and tables.\n")
