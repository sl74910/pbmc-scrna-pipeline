#!/usr/bin/env Rscript

# Input: v1 QC lists plus the original Seurat RDS files.
# Output: one DoubletFinder singlet list and one ElbowPlot per sample,
#         followed by QC_DoubletFinder_summary.csv.

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
})
source(file.path(getwd(), "R", "pbmc_helpers.R"))
source(file.path(getwd(), "R", "sigCell_QC.R"))

# All analysis choices are here.  Set any value to NA to stop before reading data.
project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(project_dir, "inputs")
v1_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv1_QC_diagnostics")
output_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv2_DoubletFinder")
figure_dir <- file.path(output_dir, "ElbowPlots")
n_variable_features <- 3000L
normalization_method <- "LogNormalize"
scale_factor <- 10000
variable_feature_method <- "vst"
doublet_dims <- 20L
pca_plot_dims <- 50L
precluster_resolution <- 0.5
doublet_rate_per_1000 <- 0.008
doublet_rate_max <- 0.20
doublet_pN <- 0.25
n_workers <- 10L
random_seed <- 1234L
pbmc_validate_parameters(list(
  n_variable_features = n_variable_features,
  normalization_method = normalization_method,
  scale_factor = scale_factor,
  variable_feature_method = variable_feature_method,
  doublet_dims = doublet_dims,
  pca_plot_dims = pca_plot_dims,
  precluster_resolution = precluster_resolution,
  doublet_rate_per_1000 = doublet_rate_per_1000,
  doublet_rate_max = doublet_rate_max,
  doublet_pN = doublet_pN,
  n_workers = n_workers,
  random_seed = random_seed
))
if (n_workers < 1L || doublet_dims < 2L || pca_plot_dims < 2L ||
    doublet_rate_max <= 0 || doublet_pN <= 0 || doublet_pN >= 1) {
  stop("DoubletFinder parameters are invalid.")
}
pbmc_make_dirs(c(output_dir, figure_dir))
pbmc_require_packages(c("Seurat", "SeuratObject", "DoubletFinder"))

files <- pbmc_sample_files(input_dir, "^PBMC_.*_Seurat\\.rds$")

# Read v1 barcodes, run DoubletFinder, and write this sample's outputs.
process_sample <- function(index) {
  sample_id <- files$sample_id[[index]]
  message("Processing ", sample_id)
  keep_file <- file.path(v1_dir, paste0(sample_id, "_high_quality_cells.txt"))
  keep <- pbmc_read_cells(keep_file, paste0("v1 QC list for ", sample_id))
  object <- pbmc_read_seurat_counts(files$file[[index]], keep, sample_id, NULL)
  cells_before <- ncol(object)

  result <- run_sigcell_doubletfinder(
    object = object,
    sample_id = sample_id,
    normalization_method = normalization_method,
    scale_factor = scale_factor,
    variable_feature_method = variable_feature_method,
    n_variable_features = min(n_variable_features, nrow(object)),
    n_pcs = min(pca_plot_dims, ncol(object) - 1L),
    dims_to_use = doublet_dims,
    precluster_resolution = precluster_resolution,
    doublet_rate_per_1000 = doublet_rate_per_1000,
    doublet_rate_max = doublet_rate_max,
    doublet_pN = doublet_pN
  )
  singlets <- colnames(result$object)[result$object$doublet_call == "Singlet"]
  if (!length(singlets)) stop("DoubletFinder removed every cell from ", sample_id)
  writeLines(singlets, file.path(output_dir, paste0(sample_id, "_high_quality_cells.txt")))
  pbmc_save_elbow_plot(
    result$object, result$summary$selected_pcs,
    file.path(figure_dir, paste0(sample_id, "_ElbowPlot.pdf")), pca_plot_dims
  )

  summary <- result$summary
  summary$cells_before_doublet <- cells_before
  summary$doublets_removed <- cells_before - length(singlets)
  summary$doublets_removed_percent <- round(100 * summary$doublets_removed / cells_before, 2L)
  summary$cells_after_doublet <- length(singlets)
  summary
}

results <- pbmc_apply_samples(
  seq_len(nrow(files)), process_sample, n_workers,
  c("Seurat", "SeuratObject", "DoubletFinder"), random_seed
)
summary <- do.call(rbind, results)
summary <- summary[match(files$sample_id, summary$sample_id), , drop = FALSE]
summary$estimated_doublet_rate_percent <- round(100 * summary$estimated_doublet_rate, 2L)
summary <- summary[, c(
  "sample_id", "cells_before_doublet", "doublets_removed", "doublets_removed_percent",
  "cells_after_doublet", "estimated_doublet_rate_percent", "selected_pK", "selected_pcs"
), drop = FALSE]
write.csv(summary, file.path(output_dir, "QC_DoubletFinder_summary.csv"), row.names = FALSE)
message("v2 DoubletFinder finished. Results: ", output_dir)
