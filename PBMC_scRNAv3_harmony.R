#!/usr/bin/env Rscript

# Input: original Seurat RDS files and v2 intersection keep lists.
# Output: unintegrated/Harmony UMAPs and compact Harmony checkpoints.
# Harmony corrects sample_id (technical library); day_group remains biological.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})
source(file.path(getwd(), "R", "pbmc_helpers.R"))
source(file.path(getwd(), "R", "pbmc_checkpoint_helpers.R"))

# Analysis parameters.  Set a parameter to NA to get an immediate error.
project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(project_dir, "inputs")
keep_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv2_intersection_summary")
out_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv3_harmony")
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
harmony_batch_variable <- "sample_id"
normalization_method <- "LogNormalize"
scale_factor <- 10000
variable_feature_method <- "vst"
pbmc_validate_parameters(list(
  days = days, measurements = measurements, biological_mouse_id = biological_mouse_id,
  random_seed = random_seed, n_variable_features = n_variable_features, n_pcs = n_pcs,
  dims_to_use = dims_to_use, cluster_resolution = cluster_resolution,
  harmony_batch_variable = harmony_batch_variable, normalization_method = normalization_method,
  scale_factor = scale_factor, variable_feature_method = variable_feature_method
))
if (n_pcs < 2L || n_variable_features < 2L || length(dims_to_use) < 2L || min(dims_to_use) < 1L) {
  stop("Harmony parameters are invalid.")
}
pbmc_make_dirs(c(figure_dir, table_dir, checkpoint_dir))
pbmc_require_packages(c("Seurat", "SeuratObject", "harmony", "ggplot2", "patchwork"))
set.seed(random_seed)

manifest <- pbmc_make_manifest(
  input_dir, keep_dir, days, measurements, biological_mouse_id,
  "_high_quality_cells.txt"
)
pbmc_check_manifest(manifest, TRUE)
write.csv(manifest, file.path(table_dir, "PBMC_sample_manifest.csv"), row.names = FALSE)

# Rebuild one Seurat object from the v2-intersection barcode list.
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
  file.path(table_dir, "PBMC_v2_intersection_retained_counts.csv"), row.names = FALSE
)

pbmc <- pbmc_merge_samples(sample_objects, "LPS_PBMC_v3", "sample_id")
rm(sample_objects)
gc(verbose = FALSE)
analysis_parameters <- list(
  random_seed = random_seed, n_variable_features = n_variable_features, n_pcs = n_pcs,
  dims_to_use = dims_to_use, cluster_resolution = cluster_resolution,
  harmony_batch_variable = harmony_batch_variable, normalization_method = normalization_method,
  scale_factor = scale_factor, variable_feature_method = variable_feature_method,
  source_v2_intersection = normalizePath(keep_dir, mustWork = FALSE)
)
save_pbmc_rdata_checkpoint(
  pbmc, "PBMC_v3_after_v2_intersection", checkpoint_dir,
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

if (!harmony_batch_variable %in% colnames(pbmc[[]])) {
  stop("Harmony batch column is missing: ", harmony_batch_variable)
}
if (length(unique(as.character(pbmc[[harmony_batch_variable]][, 1L]))) < 2L) {
  stop("Harmony requires at least two technical samples.")
}
old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
future::plan(future::sequential)
pbmc <- RunHarmony(
  object = pbmc, group.by.vars = harmony_batch_variable,
  reduction.use = "pca", dims.use = dims_to_use,
  reduction.save = "harmony", project.dim = FALSE,
  seed.use = random_seed, verbose = TRUE
)
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = dims_to_use)
pbmc <- FindClusters(pbmc, resolution = cluster_resolution, cluster.name = "clusters.harmony")
pbmc <- RunUMAP(
  pbmc, reduction = "harmony", dims = dims_to_use,
  reduction.name = "umap.harmony", reduction.key = "UMAPharmony_", seed.use = random_seed
)
pbmc_save_pdf(
  DimPlot(pbmc, reduction = "umap.harmony", group.by = "clusters.harmony", label = TRUE) /
    DimPlot(pbmc, reduction = "umap.harmony", group.by = "day_group") /
    DimPlot(pbmc, reduction = "umap.harmony", group.by = "sample_id"),
  file.path(figure_dir, "PBMC_harmony_UMAP.pdf"), 12, 18
)
save_pbmc_rdata_checkpoint(
  pbmc, "PBMC_v3_harmony_integrated", checkpoint_dir, analysis_parameters,
  layers = "counts", reductions = c("pca", "umap.unintegrated", "harmony", "umap.harmony"),
  extra_metadata = list()
)
capture.output(sessionInfo(), file = file.path(table_dir, "sessionInfo.txt"))
writeLines(c(
  "# PBMC v3 Harmony integration", "",
  paste0("- Batch variable: `", harmony_batch_variable, "` (technical library)."),
  paste0("- PCA/integration dimensions: ", min(dims_to_use), "-", max(dims_to_use), "."),
  paste0("- HVGs: ", n_variable_features, "; clustering resolution: ", cluster_resolution, "."),
  "- day_group is retained for visualization and is not corrected.", ""
), file.path(out_dir, "README.md"))
message("v3 Harmony finished. Results: ", out_dir)
