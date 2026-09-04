#!/usr/bin/env Rscript

# Input: v1 QC lists and the singlet lists from DoubletFinder/scDblFinder.
# Output: retain cells rejected by neither method, plus a summary table.

suppressPackageStartupMessages(library(Seurat))
source(file.path(getwd(), "R", "pbmc_helpers.R"))

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(project_dir, "inputs")
v1_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv1_QC_diagnostics")
doubletfinder_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv2_DoubletFinder")
scdblfinder_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv2_scDblFinder")
output_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv2_intersection_summary")
doubletfinder_summary_file <- file.path(doubletfinder_dir, "QC_DoubletFinder_summary.csv")
pbmc_validate_parameters(list(
  project_dir = project_dir,
  input_dir = input_dir,
  v1_dir = v1_dir,
  doubletfinder_dir = doubletfinder_dir,
  scdblfinder_dir = scdblfinder_dir,
  output_dir = output_dir,
  doubletfinder_summary_file = doubletfinder_summary_file
))
pbmc_make_dirs(output_dir)

files <- pbmc_sample_files(input_dir, "^PBMC_.*_Seurat\\.rds$")
if (!file.exists(doubletfinder_summary_file)) {
  stop("Missing DoubletFinder summary: ", doubletfinder_summary_file)
}
doubletfinder_summary <- read.csv(doubletfinder_summary_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample_id", "selected_pK") %in% names(doubletfinder_summary))) {
  stop("DoubletFinder summary must contain sample_id and selected_pK.")
}

get_df_pK <- function(sample_id) {
  row <- doubletfinder_summary[doubletfinder_summary$sample_id == sample_id, , drop = FALSE]
  if (!nrow(row)) stop("No DoubletFinder summary row for ", sample_id)
  pK <- suppressWarnings(as.numeric(as.character(row$selected_pK[[1L]])))
  if (!is.finite(pK)) stop("Invalid DoubletFinder selected_pK for ", sample_id)
  pK
}

read_method_cells <- function(directory, sample_id, label) {
  pbmc_read_cells(
    file.path(directory, paste0(sample_id, "_high_quality_cells.txt")),
    paste0(label, " list for ", sample_id)
  )
}

rows <- lapply(seq_len(nrow(files)), function(i) {
  sample_id <- files$sample_id[[i]]
  input <- readRDS(files$file[[i]])
  if (!inherits(input, "Seurat")) stop("Input is not a Seurat object: ", files$file[[i]])
  input_cells <- colnames(input)
  v1_cells <- intersect(
    pbmc_read_cells(file.path(v1_dir, paste0(sample_id, "_high_quality_cells.txt")),
                    paste0("v1 QC list for ", sample_id)),
    input_cells
  )
  if (!length(v1_cells)) stop("No v1 QC cells match input for ", sample_id)
  df_singlets <- intersect(read_method_cells(doubletfinder_dir, sample_id, "DoubletFinder"), v1_cells)
  sc_singlets <- intersect(read_method_cells(scdblfinder_dir, sample_id, "scDblFinder"), v1_cells)
  df_removed <- setdiff(v1_cells, df_singlets)
  sc_removed <- setdiff(v1_cells, sc_singlets)
  removed_both <- intersect(df_removed, sc_removed)
  retained <- setdiff(v1_cells, removed_both)
  if (!length(retained)) stop("Intersection removal removed every cell from ", sample_id)
  writeLines(retained, file.path(output_dir, paste0(sample_id, "_high_quality_cells.txt")))
  percent <- function(n, denominator) round(100 * n / max(1L, denominator), 2L)
  data.frame(
    sample_id = sample_id,
    cells_initial = length(input_cells),
    cells_before_v1_qc = length(input_cells),
    cells_after_v1_qc = length(v1_cells),
    cells_removed_v1_qc = length(input_cells) - length(v1_cells),
    percent_removed_v1_qc = percent(length(input_cells) - length(v1_cells), length(input_cells)),
    DoubletFinder_removed = length(df_removed),
    DF_selected_pK = get_df_pK(sample_id),
    scDblFinder_removed = length(sc_removed),
    doublet_removed_intersection = length(removed_both),
    total_removed_percent_of_initial = percent(length(input_cells) - length(retained), length(input_cells)),
    cells_after_intersection_removal = length(retained),
    cells_retained_percent_of_initial = percent(length(retained), length(input_cells)),
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, rows)
summary_file <- file.path(output_dir, "PBMC_scRNAv2_intersection_summary.csv")
write.csv(summary, summary_file, row.names = FALSE, na = "")
print(summary)
message("v2 intersection filtering finished. Results: ", output_dir)
