#!/usr/bin/env Rscript

# Input: outputs/PBMC_scRNAv3_harmony/checkpoints/PBMC_v3_harmony_integrated.RData
# Output: the unintegrated and Harmony UMAP PDFs, without rerunning analysis.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})
source(file.path(getwd(), "R", "pbmc_helpers.R"))

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
checkpoint_file <- file.path(
  project_dir, "outputs", "PBMC_scRNAv3_harmony", "checkpoints",
  "PBMC_v3_harmony_integrated.RData"
)
figure_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv3_harmony", "figures")
unintegrated_reduction <- "umap.unintegrated"
harmony_reduction <- "umap.harmony"
unintegrated_cluster_column <- "clusters.unintegrated"
harmony_cluster_column <- "clusters.harmony"
pbmc_validate_parameters(list(
  checkpoint_file = checkpoint_file, unintegrated_reduction = unintegrated_reduction,
  harmony_reduction = harmony_reduction, unintegrated_cluster_column = unintegrated_cluster_column,
  harmony_cluster_column = harmony_cluster_column
))
pbmc_make_dirs(figure_dir)

pbmc <- pbmc_load_rdata_object(checkpoint_file, "pbmc")
if (!inherits(pbmc, "Seurat")) stop("Checkpoint object `pbmc` is not a Seurat object.")
pbmc_require_reductions(pbmc, c(unintegrated_reduction, harmony_reduction))
pbmc_require_metadata(
  pbmc,
  c("sample_id", "day_group", unintegrated_cluster_column, harmony_cluster_column)
)

pbmc_save_pdf(
  pbmc_make_umap_panels(pbmc, unintegrated_reduction, unintegrated_cluster_column, "sample_id", "day_group"),
  file.path(figure_dir, "PBMC_unintegrated_UMAP.pdf"), 12, 18
)
pbmc_save_pdf(
  pbmc_make_umap_panels(pbmc, harmony_reduction, harmony_cluster_column, "sample_id", "day_group"),
  file.path(figure_dir, "PBMC_harmony_UMAP.pdf"), 12, 18
)
message("v3 Harmony plots refreshed. Results: ", figure_dir)
