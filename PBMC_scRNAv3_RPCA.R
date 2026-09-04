#!/usr/bin/env Rscript

# Input: original Seurat RDS files and v2 intersection keep lists.
# Output: unintegrated and RPCA UMAPs plus three compact RData checkpoints.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})
source(file.path(getwd(), "R", "pbmc_helpers.R"))
source(file.path(getwd(), "R", "pbmc_checkpoint_helpers.R"))

# Analysis parameters.  An NA value is an error; no Seurat default is hidden.
project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(project_dir, "inputs")
keep_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv2_intersection_summary")
out_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv3_RPCA")
figure_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")
checkpoint_dir <- file.path(out_dir, "checkpoints")
days <- c(0L, 1L, 3L, 5L)
measurements <- 1:5
biological_mouse_id <- "M1"
random_seed <- 1234L
n_variable_features <- 3000L
n_pcs <- 50L
dims_to_use <- 1:30
cluster_resolution <- 0.5
integration_workers <- 1L
normalization_method <- "LogNormalize"
scale_factor <- 10000
variable_feature_method <- "vst"
pbmc_validate_parameters(list(
  days = days, measurements = measurements, biological_mouse_id = biological_mouse_id,
  random_seed = random_seed, n_variable_features = n_variable_features, n_pcs = n_pcs,
  dims_to_use = dims_to_use, cluster_resolution = cluster_resolution,
  integration_workers = integration_workers, normalization_method = normalization_method,
  scale_factor = scale_factor, variable_feature_method = variable_feature_method
))
if (integration_workers < 1L || n_pcs < 2L || n_variable_features < 2L ||
    length(dims_to_use) < 2L || min(dims_to_use) < 1L) {
  stop("RPCA parameters are invalid.")
}
pbmc_make_dirs(c(figure_dir, table_dir, checkpoint_dir))
pbmc_require_packages(c("Seurat", "SeuratObject", "ggplot2", "patchwork", "dplyr"))
set.seed(random_seed)

manifest <- pbmc_make_manifest(
  input_dir, keep_dir, days, measurements, biological_mouse_id,
  "_high_quality_cells.txt"
)
pbmc_check_manifest(manifest, TRUE)
write.csv(manifest, file.path(table_dir, "PBMC_sample_manifest.csv"), row.names = FALSE)

# Rebuild one Seurat object from its v2-retained barcodes.
read_sample <- function(index) {
  row <- manifest[index, , drop = FALSE]
  keep <- pbmc_read_cells(row$keep_file[[1L]], paste0("v2 keep list for ", row$sample_id[[1L]]))
  object <- pbmc_read_seurat_counts(row$file[[1L]], keep, row$sample_id[[1L]], NULL)
  pbmc_add_manifest_metadata(object, row, paste0("Day", days))
}
sample_objects <- lapply(seq_len(nrow(manifest)), read_sample)
names(sample_objects) <- manifest$sample_id
write.csv(
  data.frame(sample_id = names(sample_objects),
             singlets = vapply(sample_objects, function(object) as.integer(ncol(object)), integer(1))),
  file.path(table_dir, "PBMC_v2_doubletfinder_singlet_counts.csv"), row.names = FALSE
)

pbmc <- pbmc_merge_samples(sample_objects, "LPS_PBMC_v3", "sample_id")
rm(sample_objects)
gc(verbose = FALSE)
analysis_parameters <- list(
  random_seed = random_seed, n_variable_features = n_variable_features, n_pcs = n_pcs,
  dims_to_use = dims_to_use, cluster_resolution = cluster_resolution,
  integration_workers = integration_workers, normalization_method = normalization_method,
  scale_factor = scale_factor, variable_feature_method = variable_feature_method,
  source_v2_intersection = normalizePath(keep_dir, mustWork = FALSE)
)
save_pbmc_rdata_checkpoint(
  pbmc, "PBMC_v3_after_v2_doubletfinder", checkpoint_dir,
  analysis_parameters, layers = "counts", reductions = character(), extra_metadata = list()
)

pbmc <- pbmc_preprocess(
  pbmc, normalization_method, scale_factor, variable_feature_method,
  n_variable_features, n_pcs, random_seed
)
if (max(dims_to_use) > ncol(Embeddings(pbmc, "pca"))) {
  stop("dims_to_use contains a PC that was not computed.")
}
pbmc_save_pdf(ElbowPlot(pbmc, ndims = ncol(Embeddings(pbmc, "pca"))),
              file.path(figure_dir, "PBMC_ElbowPlot.pdf"), 7, 5)

pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = dims_to_use)
pbmc <- FindClusters(pbmc, resolution = cluster_resolution, cluster.name = "clusters.unintegrated")
pbmc <- RunUMAP(
  pbmc, reduction = "pca", dims = dims_to_use,
  reduction.name = "umap.unintegrated", reduction.key = "UMAPunint_", seed.use = random_seed
)
pbmc_save_pdf(
  pbmc_make_umap_panels(pbmc, "umap.unintegrated", "clusters.unintegrated", "sample_id", "day_group"),
  file.path(figure_dir, "PBMC_unintegrated_UMAP.pdf"), 12, 18
)
save_pbmc_rdata_checkpoint(
  pbmc, "PBMC_v3_unintegrated_clustering", checkpoint_dir, analysis_parameters,
  layers = c("counts", "data"), reductions = c("pca", "umap.unintegrated"),
  extra_metadata = list()
)

if (length(unique(pbmc$sample_id)) < 2L) stop("RPCA requires at least two samples.")
old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
if (integration_workers == 1L) {
  future::plan(future::sequential)
} else {
  future::plan(future::multisession, workers = integration_workers)
}
pbmc <- IntegrateLayers(
  object = pbmc, method = RPCAIntegration, orig.reduction = "pca",
  new.reduction = "integrated.rpca", verbose = TRUE
)
pbmc <- FindNeighbors(pbmc, reduction = "integrated.rpca", dims = dims_to_use)
pbmc <- FindClusters(pbmc, resolution = cluster_resolution, cluster.name = "seurat_clusters")
pbmc <- RunUMAP(
  pbmc, reduction = "integrated.rpca", dims = dims_to_use,
  reduction.name = "umap.rpca", reduction.key = "UMAPrpca_", seed.use = random_seed
)
pbmc_save_pdf(
  DimPlot(pbmc, reduction = "umap.rpca", group.by = "seurat_clusters", label = TRUE) /
    DimPlot(pbmc, reduction = "umap.rpca", group.by = "day_group") /
    DimPlot(pbmc, reduction = "umap.rpca", group.by = "sample_id"),
  file.path(figure_dir, "PBMC_integrated_UMAP.pdf"), 12, 18
)
save_pbmc_rdata_checkpoint(
  pbmc, "PBMC_v3_rpca_integrated", checkpoint_dir, analysis_parameters,
  layers = "counts", reductions = c("pca", "umap.unintegrated", "integrated.rpca", "umap.rpca"),
  extra_metadata = list()
)
capture.output(sessionInfo(), file = file.path(table_dir, "sessionInfo.txt"))
message("v3 RPCA finished. Results: ", out_dir)
