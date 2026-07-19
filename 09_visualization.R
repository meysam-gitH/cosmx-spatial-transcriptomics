# Phase 3, Step 9: Interpretation and Visualization
#
# Plots the Degree of Clustering Permutation (Observed G - Permuted CSR)
# against radius r for every InSituType cell class. A positive curve means a
# cell type is found closer together than random chance (clustered /
# attracted); negative means farther apart than random chance (dispersed /
# repulsed); flat means consistent with complete spatial randomness (CSR).

library(tidyverse)

source("R/00_config.R")

nn_g <- read_csv(file.path(paths$results_dir, "nn_g_results.csv"), show_col_types = FALSE)
inputs <- readRDS(file.path(paths$work_dir, "spatial_inputs.rds"))
class_levels <- inputs$class_levels

palette_named <- setNames(categorical_palette[seq_along(class_levels)], class_levels)

## 1. Per-class mean Degree of Clustering Permutation across FOVs --------
nn_g_summary <- nn_g %>%
  group_by(Marker, r) %>%
  summarise(
    mean_doc = mean(`Degree of Clustering Permutation`, na.rm = TRUE),
    se_doc = sd(`Degree of Clustering Permutation`, na.rm = TRUE) / sqrt(sum(!is.na(`Degree of Clustering Permutation`))),
    .groups = "drop"
  )

p_overview <- ggplot(nn_g_summary, aes(r, mean_doc, color = Marker, fill = Marker)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#898781") +
  geom_ribbon(aes(ymin = mean_doc - se_doc, ymax = mean_doc + se_doc), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = palette_named) +
  scale_fill_manual(values = palette_named) +
  labs(
    title = "Degree of Clustering: Observed G - Permuted CSR",
    subtitle = "Positive = clustered (attracted); negative = dispersed (repulsed); mean +/- SE across FOVs",
    x = "Radius r (pixels)", y = "Degree of Clustering Permutation",
    color = "Cell class", fill = "Cell class"
  ) +
  theme_minimal()
ggsave(file.path(paths$results_dir, "03_degree_of_clustering_overview.png"), p_overview, width = 8, height = 5.5, dpi = 150)

## 2. Faceted per-class view (easier to read each curve in isolation) ----
p_facet <- ggplot(nn_g, aes(r, `Degree of Clustering Permutation`, group = unique_fov)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#898781") +
  geom_line(color = "#2a78d6", alpha = 0.6, linewidth = 0.6) +
  facet_wrap(~Marker) +
  labs(
    title = "Degree of Clustering Permutation by cell class and FOV",
    x = "Radius r (pixels)", y = "Observed G - Permuted G"
  ) +
  theme_minimal()
ggsave(file.path(paths$results_dir, "04_degree_of_clustering_by_fov.png"), p_facet, width = 9, height = 7, dpi = 150)

## 3. Summary table: clustering value at a fixed reference radius -------
ref_r <- 100
clustering_at_r <- nn_g_summary %>%
  filter(r == ref_r) %>%
  arrange(desc(mean_doc))
write_csv(clustering_at_r, file.path(paths$results_dir, "clustering_value_at_r100.csv"))

message("Saved results/03_degree_of_clustering_overview.png, results/04_degree_of_clustering_by_fov.png")
message(sprintf("Clustering ranking at r = %d:", ref_r))
print(clustering_at_r)
