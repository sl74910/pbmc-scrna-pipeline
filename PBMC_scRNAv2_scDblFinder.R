#!/usr/bin/env Rscript

# Input: v1 QC lists plus the original Seurat RDS files.
# Output: one scDblFinder singlet list, score table and PCA plot per sample,
#         followed by QC_scDblFinder_summary.csv.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(BiocParallel)
})
source(file.path(getwd(), "R", "pbmc_helpers.R"))
source(file.path(getwd(), "R", "sigCell_QC.R"))

# All analysis choices are explicit.  An NA value stops the script immediately.
project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(project_dir, "inputs")
v1_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv1_QC_diagnostics")
output_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv2_scDblFinder")
figure_dir <- file.path(output_dir, "ElbowPlots")
score_dir <- file.path(output_dir, "Scores")
n_variable_features <- 3000L
normalization_method <- "LogNormalize"
scale_factor <- 10000
variable_feature_method <- "vst"
doublet_dims <- 20L
pca_plot_dims <- 50L
precluster_resolution <- 0.5
doublet_rate_per_1000 <- 0.008
doublet_rate_max <- 0.20
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
  n_workers = n_workers,
  random_seed = random_seed
))
if (n_workers < 1L || doublet_dims < 2L || pca_plot_dims < 2L ||
    doublet_rate_max <= 0 || doublet_rate_per_1000 <= 0) {
  stop("scDblFinder parameters are invalid.")
}
pbmc_make_dirs(c(output_dir, figure_dir, score_dir))
pbmc_require_packages(c(
  "Seurat", "SeuratObject", "scDblFinder", "SingleCellExperiment", "BiocParallel"
))

files <- pbmc_sample_files(input_dir, "^PBMC_.*_Seurat\\.rds$")

# Read v1 barcodes, run scDblFinder, and write this sample's outputs.
process_sample <- function(index) {
  sample_id <- files$sample_id[[index]]
  message("Processing ", sample_id)
  keep_file <- file.path(v1_dir, paste0(sample_id, "_high_quality_cells.txt"))
  keep <- pbmc_read_cells(keep_file, paste0("v1 QC list for ", sample_id))
  old <- readRDS(files$file[[index]])
  if (!inherits(old, "Seurat")) stop("Input is not a Seurat object: ", files$file[[index]])
  if (!"RNA" %in% SeuratObject::Assays(old) ||
      !"counts" %in% SeuratObject::Layers(old[["RNA"]])) {
    stop("Input requires an RNA counts layer: ", files$file[[index]])
  }
  counts <- SeuratObject::LayerData(old, assay = "RNA", layer = "counts")
  keep <- intersect(keep, colnames(counts))
  if (!length(keep)) stop("No v1 QC barcodes match ", sample_id)
  counts <- counts[, keep, drop = FALSE]
  cells_before <- ncol(counts)

  result <- run_sigcell_scdblfinder(
    counts = counts,
    sample_id = sample_id,
    normalization_method = normalization_method,
    scale_factor = scale_factor,
    variable_feature_method = variable_feature_method,
    n_variable_features = n_variable_features,
    doublet_dims = doublet_dims,
    pca_plot_dims = pca_plot_dims,
    precluster_resolution = precluster_resolution,
    doublet_rate_per_1000 = doublet_rate_per_1000,
    doublet_rate_max = doublet_rate_max,
    random_seed = random_seed
  )
  singlets <- colnames(counts)[result$classification == "singlet"]
  if (!length(singlets)) stop("scDblFinder removed every cell from ", sample_id)
  writeLines(singlets, file.path(output_dir, paste0(sample_id, "_high_quality_cells.txt")))

  score_table <- data.frame(
    barcode = colnames(counts),
    scDblFinder_class = result$classification,
    scDblFinder_score = result$score,
    scDblFinder_cluster = if (is.null(result$clusters)) {
      rep(NA_character_, ncol(counts))
    } else {
      as.character(result$clusters[colnames(counts)])
    },
    stringsAsFactors = FALSE
  )
  write.csv(score_table, file.path(score_dir, paste0(sample_id, "_scDblFinder_scores.csv")), row.names = FALSE)
  pbmc_save_scdblfinder_elbow(
    result$pca_stdev, max(result$clustering_dims), sample_id,
    file.path(figure_dir, paste0(sample_id, "_ElbowPlot.pdf")), pca_plot_dims
  )

  doublets <- sum(result$classification == "doublet")
  data.frame(
    sample_id = sample_id,
    cells_before_doublet = cells_before,
    doublets_called = doublets,
    doublets_called_percent = round(100 * doublets / cells_before, 2L),
    singlets_retained = length(singlets),
    estimated_doublet_rate_percent = round(100 * result$dbr, 2L),
    scDblFinder_dims = min(doublet_dims, ncol(counts) - 1L),
    precluster_mode = result$precluster_mode,
    preclusters = result$n_preclusters,
    stringsAsFactors = FALSE
  )
}

results <- pbmc_apply_samples(
  seq_len(nrow(files)), process_sample, n_workers,
  c("Seurat", "SeuratObject", "scDblFinder", "SingleCellExperiment", "BiocParallel"),
  random_seed
)
summary <- do.call(rbind, results)
summary <- summary[match(files$sample_id, summary$sample_id), , drop = FALSE]
write.csv(summary, file.path(output_dir, "QC_scDblFinder_summary.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
message("v2 scDblFinder finished. Results: ", output_dir)
