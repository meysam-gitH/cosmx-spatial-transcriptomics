# CosMx SMI Spatial Transcriptomics Pipeline

An end-to-end R/Seurat pipeline for high-resolution, single-cell-resolution
spatial transcriptomics from NanoString's CosMx Spatial Molecular Imager
(SMI): memory-safe object handling, SCTransform normalization,
PCA/UMAP/Louvain clustering, unsupervised cell typing with **InSituType**,
and spatial point-pattern statistics with **spatialTIME** to test whether
each cell type is spatially clustered or dispersed relative to random
chance. Every step below ran against real CosMx instrument output and
produced the figures shown in [Results](#results).

## Objective

CosMx SMI captures subcellular-resolution gene expression while preserving
each cell's physical (x, y) location in tissue — something bulk and
droplet-based single-cell RNA-seq cannot do. This project builds a
reproducible pipeline that:

1. **Manages memory** for the very large objects CosMx experiments produce
   (per-molecule transcript tables can be gigabytes on their own).
2. **Phenotypes cells** — normalizes expression, reduces dimensionality, and
   assigns unsupervised cell-type labels from expression *and*
   immunofluorescence data.
3. **Quantifies spatial organization** — tests each cell type for spatial
   clustering (attraction) or dispersion (repulsion) using the univariate
   Nearest Neighbor G statistic, a formal point-pattern analysis method.

## A note on the dataset

The methodology this pipeline implements was originally described against a
CosMx **Human Frontal Cortex** dataset (6,078-gene WTx panel, ~188,000
cells). That specific dataset is not publicly downloadable. Rather than
simulate data, this project runs against NanoString's own **public** CosMx
SMI release instead: **Lung9_Rep1** (2021 non-small-cell lung cancer FFPE
dataset, 960-gene panel), subset to a 2x2 block of fields of view (FOVs 1,
2, 5, 6 — **17,622 cells** after QC) so the full pipeline runs in minutes on
a single machine rather than requiring 64GB RAM and hours of compute. The
methodology — memory optimization, SCTransform, PCA/UMAP/Louvain,
InSituType, spatialTIME NN-G — is identical regardless of tissue of origin
or panel size; only the biological narrative in the interpretation section
below is specific to lung tissue rather than brain.

## Pipeline

| Phase | Step | Script |
|---|---|---|
| 1. Data optimization | Build Seurat object from CosMx flat files; extract + `fwrite` the per-molecule transcript table; drop it from the object with `gc(full=TRUE)` | [`R/01_data_optimization.R`](R/01_data_optimization.R) |
| 1. QC | Confirm the RNA assay contains only true gene targets, no `NegPrb`/control probes | [`R/02_qc_check.R`](R/02_qc_check.R) |
| 2. Normalization | `SCTransform(assay="RNA", new.assay.name="SCT")` | [`R/03_normalization_sctransform.R`](R/03_normalization_sctransform.R) |
| 2. Dimension reduction | PCA (npcs=50) -> UMAP (dims=1:30, `repulsion.strength=10`) | [`R/04_dimension_reduction.R`](R/04_dimension_reduction.R) |
| 2. Louvain clustering | `FindNeighbors` + `FindClusters(algorithm=1, resolution=0.3)` | [`R/05_clustering_louvain.R`](R/05_clustering_louvain.R) |
| 2. Cell typing | NanoString **InSituType**: raw counts + negative-probe background + IF cohorting -> unsupervised clusters (searched 4-8) | [`R/06_insitutype_annotation.R`](R/06_insitutype_annotation.R) |
| 3. Spatial data prep | Metadata -> Spatial List (per FOV) + Summary DF + Clinical DF, linked by `unique_fov` | [`R/07_spatial_data_prep.R`](R/07_spatial_data_prep.R) |
| 3. Spatial statistics | **spatialTIME** `create_mif()` + `NN_G()`: univariate Nearest Neighbor G, r = 0-200px | [`R/08_spatial_nn_g.R`](R/08_spatial_nn_g.R) |
| 3. Interpretation | Degree of Clustering Permutation (Observed G - Permuted CSR) vs. radius | [`R/09_visualization.R`](R/09_visualization.R) |

Run the whole thing with `Rscript R/run_pipeline.R` from this directory.

## Setup

```bash
# 1. Install R + all dependencies (apt where possible, GitHub source builds
#    for InSituType, spatialTIME, and their few non-apt CRAN dependencies).
./scripts/install_dependencies.sh

# 2. Download and subset the public CosMx Lung9_Rep1 sample (~1.1GB download).
./scripts/download_and_subset_data.sh

# 3. Run the pipeline end to end.
Rscript R/run_pipeline.R
```

CRAN and Bioconductor were not reachable from the environment this was
built in, so `install_dependencies.sh` sources everything from Ubuntu's apt
archives (which mirror pre-built CRAN/Bioconductor binaries as
`r-cran-*`/`r-bioc-*` packages) plus `git clone` + `R CMD INSTALL` for the
handful of packages apt doesn't carry. One of those, `fastglm`, fails to
compile against modern Rcpp as released (a `List::create()` call exceeds
Rcpp's argument-count limit); [`scripts/patches/fastglm_fit_glm_hurdle.patch`](scripts/patches/fastglm_fit_glm_hurdle.patch)
fixes it by building the same result list element-by-element instead.

## Results

### Memory optimization (Phase 1)

Building the Seurat object and attaching the per-molecule transcript table
(the same role the ~8GB transcript-coordinates slot plays in the full-size
dataset) brought the object to **270.5 MB**. Writing that table out with
`data.table::fwrite` and dropping it from the object with
`gc(full=TRUE)` cut it to **32.7 MB** — an 8x reduction, the same
extract-then-drop pattern that takes the full dataset from ~21GB to ~14GB.

### Cellular phenotyping (Phase 2)

Louvain clustering at resolution 0.3 recovers 11 preliminary clusters:

![UMAP with Louvain clusters](results/01_umap_louvain_clusters.png)

InSituType's unsupervised search (4-8 clusters, raw counts + negative-probe
background + immunofluorescence cohorting on membrane stain, PanCK, CD45,
and CD3 — DAPI excluded) converged on **8 cell classes (a-h)**, ranging from
754 to 4,045 cells. The Flight Path Plot lays cells out by their posterior
cluster probability: cells near a cluster center are called with high
confidence, cells between centers represent ambiguous phenotypes. All 8
clusters here land at mean confidence 0.94-0.99 — tight, well-separated
calls:

![InSituType flight path plot](results/02_flightpath_plot.png)

### Spatial statistics (Phase 3)

The univariate Nearest Neighbor G statistic (Observed G minus permuted
Complete Spatial Randomness, `spatialTIME::NN_G()`, r = 0-200px in steps of
10, 50 permutations, edge-corrected) was computed for each of the 8 cell
classes across all 4 FOVs:

![Degree of clustering overview](results/03_degree_of_clustering_overview.png)
![Degree of clustering by FOV](results/04_degree_of_clustering_by_fov.png)

**Interpretation.** Every one of the 8 cell classes shows a *positive*
Degree of Clustering Permutation peaking around r = 70-100px before
decaying — i.e., in this tissue block, cells of the same InSituType class
consistently sit closer to each other than Complete Spatial Randomness
predicts, with no class showing net dispersion. Ranked by clustering
strength at r = 100px:

| Cell class | Degree of Clustering at r=100 |
|---|---|
| c | 0.554 |
| a | 0.454 |
| g | 0.374 |
| d | 0.351 |
| e | 0.314 |
| b | 0.293 |
| f | 0.229 |
| h | 0.201 |

(Full table: [`results/clustering_value_at_r100.csv`](results/clustering_value_at_r100.csv);
raw per-FOV/per-radius values: [`results/nn_g_results.csv`](results/nn_g_results.csv).)
This is consistent with tissue architecture in general — most cell types in
solid tissue occupy spatially coherent niches rather than being randomly
salted through the section — and with FFPE lung tissue specifically, where
tumor, immune, and stromal compartments form contiguous regions. A quieter
class like **h** clustering less strongly than **c** is exactly the kind of
per-class signal (a "clustering value at r" number, as flagged in the
Real-World Application section of the source methodology) that would be
compared across patients as a candidate prognostic feature in a clinical
cohort.

## Real-world application

This is the same analytical approach used to study the tumor immune
microenvironment: quantifying whether specific immune cell types cluster
around tumor cells (or each other) at a given radius turns "cell-cell
proximity" into a number that can be compared across patients — e.g., as a
candidate biomarker for immunotherapy response.

## Skills demonstrated

- Memory-safe handling of large Seurat S4 objects (extract -> persist ->
  drop -> `gc(full=TRUE)`) for datasets too large to hold comfortably in
  memory
- Single-cell normalization (SCTransform), dimensionality reduction
  (PCA/UMAP), and graph-based (Louvain) clustering with Seurat
- Unsupervised, multi-modal cell typing (expression + immunofluorescence)
  with NanoString's InSituType, including negative-probe background
  modeling
- Spatial point-pattern analysis (univariate Nearest Neighbor G against
  permuted Complete Spatial Randomness) with spatialTIME/spatstat
- Reshaping single-cell metadata into the three linked data structures
  (spatial list / summary / clinical) a spatial-statistics package expects
- Reproducible environment setup on a constrained network (apt-mirrored
  CRAN/Bioconductor binaries + GitHub source builds in place of direct
  CRAN/Bioconductor access), including diagnosing and patching a real C++
  compile failure in a third-party package
- Working directly with raw CosMx SMI instrument output (flat files: per-
  cell expression matrix, metadata, FOV positions, per-molecule transcript
  table) rather than a pre-packaged object

## Repository layout

```
cosmx-spatial-pipeline/
  R/                    pipeline scripts, run in numeric order (00-09)
  scripts/               install_dependencies.sh, download_and_subset_data.sh, patches/
  results/               figures + tables produced by the pipeline (checked in)
  data/                  downloaded/derived data (gitignored; regenerate with the download script)
  vendor/                 GitHub-sourced R package builds (gitignored; regenerate with the install script)
```

## References

- He, S. et al. "High-plex imaging of RNA and proteins at subcellular
  resolution in fixed tissue by spatial molecular imaging." *Nat.
  Biotechnol.* (2022) — CosMx SMI.
- Danaher, P. et al. InSituType: <https://github.com/Nanostring-Biostats/InSituType>
- Creed, J.H. et al. "spatialTIME and iTIME.io: A R package and Shiny
  application for visualization and analysis of immunofluorescence data."
  *Bioinformatics* (2021). <https://doi.org/10.1093/bioinformatics/btab757>
- NanoString CosMx SMI public FFPE dataset (Lung9_Rep1):
  <https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/>
