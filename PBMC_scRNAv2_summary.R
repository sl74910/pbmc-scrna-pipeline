#!/usr/bin/env Rscript

# Input: v1, DoubletFinder, scDblFinder and intersection summary outputs.
# Output: one CSV comparison and one short Markdown report.

suppressPackageStartupMessages(library(Seurat))
source(file.path(getwd(), "R", "pbmc_helpers.R"))

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
doubletfinder_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv2_DoubletFinder")
scdblfinder_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv2_scDblFinder")
v1_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv1_QC_diagnostics")
intersection_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv2_intersection_summary")
output_file <- file.path(project_dir, "outputs", "PMBC_scRNAv2_method_comparison.csv")
markdown_file <- file.path(project_dir, "outputs", "PMBC_scRNAv2_method_comparison.md")

# These columns document the choices used by the two v2 scripts.
DoubletFinder_n_variable_features <- 3000L
DoubletFinder_normalization <- "LogNormalize"
DoubletFinder_scale_factor <- 10000
DoubletFinder_selection_method <- "vst"
DoubletFinder_pca_dims_for_plot <- 50L
DoubletFinder_pca_dims_for_detection <- 20L
DoubletFinder_precluster_resolution <- 0.5
DoubletFinder_doublet_rate_per_1000 <- 0.008
DoubletFinder_doublet_rate_max <- 0.20
DoubletFinder_pN <- 0.25
scDblFinder_nfeatures <- 3000L
scDblFinder_normalization <- "LogNormalize (custom processing)"
scDblFinder_scale_factor <- 10000
scDblFinder_selection_method <- "vst (Seurat preclustering)"
scDblFinder_pca_dims_for_plot <- 50L
scDblFinder_pca_dims_for_detection <- 20L
scDblFinder_precluster_resolution <- 0.5
scDblFinder_doublet_rate_per_1000 <- 0.008
scDblFinder_doublet_rate_max <- 0.20
pbmc_validate_parameters(list(
  DoubletFinder_n_variable_features = DoubletFinder_n_variable_features,
  DoubletFinder_normalization = DoubletFinder_normalization,
  DoubletFinder_scale_factor = DoubletFinder_scale_factor,
  DoubletFinder_selection_method = DoubletFinder_selection_method,
  DoubletFinder_pca_dims_for_plot = DoubletFinder_pca_dims_for_plot,
  DoubletFinder_pca_dims_for_detection = DoubletFinder_pca_dims_for_detection,
  DoubletFinder_precluster_resolution = DoubletFinder_precluster_resolution,
  DoubletFinder_doublet_rate_per_1000 = DoubletFinder_doublet_rate_per_1000,
  DoubletFinder_doublet_rate_max = DoubletFinder_doublet_rate_max,
  DoubletFinder_pN = DoubletFinder_pN,
  scDblFinder_nfeatures = scDblFinder_nfeatures,
  scDblFinder_normalization = scDblFinder_normalization,
  scDblFinder_scale_factor = scDblFinder_scale_factor,
  scDblFinder_selection_method = scDblFinder_selection_method,
  scDblFinder_pca_dims_for_plot = scDblFinder_pca_dims_for_plot,
  scDblFinder_pca_dims_for_detection = scDblFinder_pca_dims_for_detection,
  scDblFinder_precluster_resolution = scDblFinder_precluster_resolution,
  scDblFinder_doublet_rate_per_1000 = scDblFinder_doublet_rate_per_1000,
  scDblFinder_doublet_rate_max = scDblFinder_doublet_rate_max
))
pbmc_make_dirs(dirname(output_file))

read_summary <- function(file, label) {
  if (!file.exists(file)) stop("Missing ", label, " summary: ", file)
  result <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"sample_id" %in% names(result)) stop(label, " summary has no sample_id column.")
  result
}
df_summary <- read_summary(file.path(doubletfinder_dir, "QC_DoubletFinder_summary.csv"), "DoubletFinder")
sc_summary <- read_summary(file.path(scdblfinder_dir, "QC_scDblFinder_summary.csv"), "scDblFinder")
intersection_summary <- read_summary(
  file.path(intersection_dir, "PBMC_scRNAv2_intersection_summary.csv"), "intersection"
)
sample_ids <- sort(unique(c(df_summary$sample_id, sc_summary$sample_id)))

read_singlets <- function(directory, sample_id, label) {
  pbmc_read_cells(
    file.path(directory, paste0(sample_id, "_high_quality_cells.txt")),
    paste0(label, " list for ", sample_id)
  )
}
get_value <- function(data, sample_id, column) {
  if (!column %in% names(data)) return(NA)
  rows <- data[data$sample_id == sample_id, , drop = FALSE]
  if (!nrow(rows)) NA else rows[[column]][[1L]]
}

comparison <- do.call(rbind, lapply(sample_ids, function(sample_id) {
  input_cells <- pbmc_read_cells(
    file.path(v1_dir, paste0(sample_id, "_high_quality_cells.txt")),
    paste0("v1 QC list for ", sample_id)
  )
  df_removed <- setdiff(input_cells, read_singlets(doubletfinder_dir, sample_id, "DoubletFinder"))
  sc_removed <- setdiff(input_cells, read_singlets(scdblfinder_dir, sample_id, "scDblFinder"))
  both <- intersect(df_removed, sc_removed)
  union_removed <- union(df_removed, sc_removed)
  overall <- intersection_summary[intersection_summary$sample_id == sample_id, , drop = FALSE]
  if (!nrow(overall)) stop("Missing intersection summary for ", sample_id)
  initial <- as.integer(overall$cells_initial[[1L]])
  after <- as.integer(overall$cells_after_intersection_removal[[1L]])
  percent <- function(n, denominator) round(100 * n / max(1L, denominator), 2L)
  data.frame(
    sample_id = sample_id,
    cells_before_qc = initial,
    cells_after_qc_doublet = after,
    total_removed = initial - after,
    total_removed_percent = as.numeric(overall$total_removed_percent_of_initial[[1L]]),
    input_cells = length(input_cells),
    DoubletFinder_removed = length(df_removed),
    scDblFinder_removed = length(sc_removed),
    removed_intersection = length(both),
    removed_union = length(union_removed),
    DoubletFinder_only_removed = length(setdiff(df_removed, sc_removed)),
    scDblFinder_only_removed = length(setdiff(sc_removed, df_removed)),
    removed_intersection_percent_of_union = percent(length(both), length(union_removed)),
    removed_intersection_percent_of_DF = percent(length(both), length(df_removed)),
    removed_intersection_percent_of_scDblFinder = percent(length(both), length(sc_removed)),
    DF_estimated_doublet_rate_percent = get_value(df_summary, sample_id, "estimated_doublet_rate_percent"),
    scDblFinder_estimated_doublet_rate_percent = get_value(sc_summary, sample_id, "estimated_doublet_rate_percent"),
    DF_selected_pK = get_value(df_summary, sample_id, "selected_pK"),
    DF_selected_pcs = get_value(df_summary, sample_id, "selected_pcs"),
    scDblFinder_precluster_mode = get_value(sc_summary, sample_id, "precluster_mode"),
    scDblFinder_preclusters = get_value(sc_summary, sample_id, "preclusters"),
    scDblFinder_dims = get_value(sc_summary, sample_id, "scDblFinder_dims"),
    DoubletFinder_n_variable_features = DoubletFinder_n_variable_features,
    DoubletFinder_normalization = DoubletFinder_normalization,
    DoubletFinder_scale_factor = DoubletFinder_scale_factor,
    DoubletFinder_selection_method = DoubletFinder_selection_method,
    DoubletFinder_pca_dims_for_plot = DoubletFinder_pca_dims_for_plot,
    DoubletFinder_pca_dims_for_detection = DoubletFinder_pca_dims_for_detection,
    DoubletFinder_precluster_resolution = DoubletFinder_precluster_resolution,
    DoubletFinder_doublet_rate_per_1000 = DoubletFinder_doublet_rate_per_1000,
    DoubletFinder_doublet_rate_max = DoubletFinder_doublet_rate_max,
    DoubletFinder_pN = DoubletFinder_pN,
    scDblFinder_nfeatures = scDblFinder_nfeatures,
    scDblFinder_normalization = scDblFinder_normalization,
    scDblFinder_scale_factor = scDblFinder_scale_factor,
    scDblFinder_selection_method = scDblFinder_selection_method,
    scDblFinder_pca_dims_for_plot = scDblFinder_pca_dims_for_plot,
    scDblFinder_pca_dims_for_detection = scDblFinder_pca_dims_for_detection,
    scDblFinder_precluster_resolution = scDblFinder_precluster_resolution,
    scDblFinder_doublet_rate_per_1000 = scDblFinder_doublet_rate_per_1000,
    scDblFinder_doublet_rate_max = scDblFinder_doublet_rate_max,
    stringsAsFactors = FALSE
  )
}))
write.csv(comparison, output_file, row.names = FALSE, na = "")

fmt_pct <- function(value) sprintf("%.2f%%", as.numeric(value))
report <- c(
  "# PBMC v2 doublet-method comparison", "",
  "A cell is counted as removed when it is absent from that method's singlet list.", "",
  "| Sample | Input | DoubletFinder removed | scDblFinder removed | Both | Union | Jaccard |",
  "|---|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(comparison)), function(i) {
    row <- comparison[i, ]
    paste0("| ", row$sample_id, " | ", row$input_cells, " | ", row$DoubletFinder_removed,
           " | ", row$scDblFinder_removed, " | ", row$removed_intersection, " | ",
           row$removed_union, " | ", fmt_pct(100 * row$removed_intersection / max(1, row$removed_union)), " |")
  }, character(1)),
  "", "Parameters are recorded in the CSV columns; Jaccard = Both / Union."
)
writeLines(report, markdown_file)
print(comparison)
message("v2 comparison finished. Wrote: ", output_file)
