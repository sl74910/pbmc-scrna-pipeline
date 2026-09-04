#!/usr/bin/env Rscript

# 基于 marker基因对单细胞PBMC 注释
# 输入 Harmony 的结果，输出 marker 表、可编辑的 cluster 标签、注释对象和所有 UMAP 图。
# 首次运行会创建 PBMC_v4_cluster_annotation.csv。请编辑其中的 celltype 列，
# 然后重新运行以使用审核后的标签替代候选标签。

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})
source(file.path(getwd(), "R", "pbmc_helpers.R"))

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
checkpoint_file <- file.path(
  project_dir, "outputs", "PBMC_scRNAv3_harmony", "checkpoints",
  "PBMC_v3_harmony_integrated.RData"
)
out_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv4_annotation")
figure_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")
cluster_column <- "clusters.harmony"
umap_reduction <- "umap.harmony"
max_cells_per_cluster <- 5000L
marker_min_pct <- 0.25
marker_logfc_threshold <- 0.25
marker_return_threshold <- 0.05
top_marker_count <- 20L
normalization_method <- "LogNormalize"
scale_factor <- 10000
day_levels <- c("Day0", "Day1", "Day3", "Day5")
random_seed <- 1234L
pbmc_validate_parameters(list(
  checkpoint_file = checkpoint_file, cluster_column = cluster_column,
  umap_reduction = umap_reduction, max_cells_per_cluster = max_cells_per_cluster,
  marker_min_pct = marker_min_pct, marker_logfc_threshold = marker_logfc_threshold,
  marker_return_threshold = marker_return_threshold, top_marker_count = top_marker_count,
  normalization_method = normalization_method, scale_factor = scale_factor,
  day_levels = day_levels, random_seed = random_seed
))
if (max_cells_per_cluster < 1L || marker_min_pct <= 0 || marker_min_pct > 1 ||
    marker_logfc_threshold < 0 || marker_return_threshold <= 0 || top_marker_count < 1L) {
  stop("Marker parameters are invalid.")
}
pbmc_make_dirs(c(figure_dir, table_dir))
set.seed(random_seed)

pbmc <- pbmc_load_rdata_object(checkpoint_file, "pbmc")
if (!inherits(pbmc, "Seurat")) stop("Checkpoint object `pbmc` is not a Seurat object.")
pbmc_require_metadata(pbmc, c(
  cluster_column, "day_group", "sample_id", "mouse_timepoint_id", "mouse_id"
))
pbmc_require_reductions(pbmc, umap_reduction)
if (!"RNA" %in% Assays(pbmc)) stop("The checkpoint must contain an RNA assay.")
DefaultAssay(pbmc) <- "RNA"
Idents(pbmc) <- cluster_column
cluster_levels <- unique(as.character(pbmc[[cluster_column]][, 1L]))
cluster_levels <- cluster_levels[order(suppressWarnings(as.integer(cluster_levels)), cluster_levels)]

# marker 模块特意使用普通向量，便于初学者编辑。
lineage_markers <- list(
  T_cell = c("Cd3d", "Cd3e", "Trbc1", "Trbc2", "Lck", "Il7r", "Lef1", "Tcf7", "Cd4", "Cd8a", "Cd8b1"),
  B_cell = c("Cd79a", "Cd79b", "Ms4a1", "Cd19", "Cd22", "Cd74", "H2-Aa", "H2-Ab1", "Ighd", "Mzb1"),
  NK_cell = c("Nkg7", "Klrd1", "Ncr1", "Prf1", "Gzmb", "Gzmk", "Tyrobp", "Fcerg1", "Ccl5"),
  Monocyte = c("Lyz2", "Lst1", "Ctss", "Csf1r", "Adgre1", "Ms4a7", "Fcgr3", "Ly6c2", "Ccr2", "Vcan", "S100a8", "S100a9"),
  Neutrophil = c("Ly6g", "Mpo", "Camp", "Ngp", "Retnlg", "Ltf", "S100a8", "S100a9", "Wfdc21"),
  DC = c("Flt3", "Itgax", "Clec10a", "Clec9a", "Fscn1", "Cd74", "H2-Ab1", "Ciita", "Zbtb46"),
  pDC = c("Siglech", "Gzmb", "Bst2", "Tcf4", "Irf8", "Jchain"),
  Basophil = c("Mcpt8", "Ms4a2", "Fcer1a", "Cpa3", "Hdc"),
  Erythroid = c("Hbb-bs", "Hbb-bt", "Alas2", "Gata1", "Klf1", "Epb42", "Ermap"),
  Platelet = c("Pf4", "Ppbp", "Nrg1", "Itga2b", "Gp9"),
  Plasma_cell = c("Jchain", "Mzb1", "Sdc1", "Derl3", "Xbp1", "Igkc")
)
state_markers <- list(
  Cycling = c("Mki67", "Top2a", "Pcna", "Tuba1b", "Stmn1", "Ube2c", "Birc5"),
  IFN_responsive = c("Ifit1", "Ifit2", "Ifit3", "Isg15", "Mx1", "Oas1a", "Rsad2", "Ifitm3"),
  Cytotoxic = c("Nkg7", "Ccl5", "Prf1", "Gzmb", "Gzmk", "Ctsw"),
  Immediate_early = c("Fos", "Jun", "Dusp1", "Nfkbia", "Ier2", "Egr1", "Tnf", "Cxcl10")
)
all_marker_genes <- unique(unlist(c(lineage_markers, state_markers), use.names = FALSE))
marker_genes_present <- intersect(all_marker_genes, rownames(pbmc))
if (length(marker_genes_present) < 10L) stop("Too few marker genes match the object.")
marker_panel <- bind_rows(lapply(names(c(lineage_markers, state_markers)), function(group) {
  category <- if (group %in% names(lineage_markers)) "lineage" else "state"
  genes <- intersect(c(lineage_markers, state_markers)[[group]], rownames(pbmc))
  data.frame(category = rep(category, length(genes)), group = rep(group, length(genes)),
             gene = genes,
             stringsAsFactors = FALSE)
}))
write.csv(marker_panel, file.path(table_dir, "PBMC_v4_marker_panel.csv"), row.names = FALSE)
definitions <- data.frame(
  label = c(names(lineage_markers), names(state_markers), "No_dominant_state", "Review_required"),
  meaning = c(
    rep("Candidate lineage or transcriptional state; review marker table and UMAP.",
        length(lineage_markers) + length(state_markers)),
    "No state module has a positive centered score.",
    "Marker evidence is mixed or insufficient; review manually."
  ),
  canonical_markers = c(
    vapply(c(lineage_markers, state_markers), function(x) paste(head(x, 6L), collapse = ";"), character(1)),
    "", ""
  ),
  stringsAsFactors = FALSE
)
write.csv(definitions, file.path(table_dir, "PBMC_v4_celltype_definitions.csv"), row.names = FALSE)

# 计算 RNA 平均表达量前先合并 count 层。
if (length(grep("^counts\\.", Layers(pbmc[["RNA"]]))) > 1L) {
  pbmc <- JoinLayers(pbmc, assay = "RNA")
}
if (!"data" %in% Layers(pbmc[["RNA"]])) {
  pbmc <- NormalizeData(pbmc, normalization.method = normalization_method,
                        scale.factor = scale_factor, verbose = FALSE)
}
average_expression <- AverageExpression(
  pbmc, assays = "RNA", group.by = cluster_column, layer = "data", verbose = FALSE
)$RNA
average_expression <- as.matrix(average_expression)
colnames(average_expression) <- sub(paste0("^", cluster_column, "_"), "", colnames(average_expression))
colnames(average_expression) <- sub("^g(?=[0-9]+$)", "", colnames(average_expression), perl = TRUE)
if (!setequal(colnames(average_expression), cluster_levels)) stop("AverageExpression cluster names do not match.")
average_expression <- average_expression[, cluster_levels, drop = FALSE]
centered_expression <- sweep(average_expression, 1L, rowMeans(average_expression), "-")
lineage_scores <- pbmc_score_modules(centered_expression, lineage_markers)
state_scores <- pbmc_score_modules(centered_expression, state_markers)
lineage_best <- pbmc_best_module(lineage_scores)
state_best <- pbmc_best_module(state_scores)
state_best$candidate[!is.na(state_best$score) & state_best$score <= 0] <- "No_dominant_state"

candidate_table <- data.frame(
  cluster = cluster_levels,
  candidate_lineage = lineage_best[cluster_levels, "candidate"],
  lineage_score = lineage_best[cluster_levels, "score"],
  lineage_margin = lineage_best[cluster_levels, "margin"],
  candidate_state = state_best[cluster_levels, "candidate"],
  state_score = state_best[cluster_levels, "score"],
  state_margin = state_best[cluster_levels, "margin"],
  n_cells = as.integer(table(factor(as.character(pbmc[[cluster_column]][, 1L]), levels = cluster_levels))),
  stringsAsFactors = FALSE
)
write.csv(lineage_scores, file.path(table_dir, "PBMC_v4_cluster_lineage_scores.csv"))
write.csv(state_scores, file.path(table_dir, "PBMC_v4_cluster_state_scores.csv"))

message("Finding positive markers...")
cluster_markers <- FindAllMarkers(
  pbmc, assay = "RNA", only.pos = TRUE, min.pct = marker_min_pct,
  logfc.threshold = marker_logfc_threshold, max.cells.per.ident = max_cells_per_cluster,
  return.thresh = marker_return_threshold, densify = FALSE
)
write.csv(cluster_markers, file.path(table_dir, "PBMC_v4_cluster_markers.csv"), row.names = FALSE)
top_markers <- cluster_markers %>%
  mutate(cluster = as.character(cluster)) %>%
  group_by(cluster) %>% arrange(desc(avg_log2FC), p_val_adj, .by_group = TRUE) %>%
  slice_head(n = top_marker_count) %>% summarise(top_markers = paste(gene, collapse = ";"), .groups = "drop")
candidate_table <- candidate_table %>%
  left_join(top_markers, by = "cluster") %>%
  mutate(candidate_annotation = paste(candidate_state, candidate_lineage, sep = " + "),
         celltype = NA_character_, annotation_notes = "Review markers and UMAP; edit celltype.")

annotation_file <- file.path(table_dir, "PBMC_v4_cluster_annotation.csv")
if (file.exists(annotation_file)) {
  edited <- read.csv(annotation_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("cluster", "celltype") %in% names(edited))) {
    stop("Existing annotation file must contain cluster and celltype.")
  }
  candidate_table <- candidate_table %>%
    select(-celltype) %>%
    left_join(edited %>% select(cluster, celltype) %>% mutate(cluster = as.character(cluster)), by = "cluster")
}
write.csv(candidate_table, annotation_file, row.names = FALSE)

labels <- setNames(candidate_table$candidate_annotation, candidate_table$cluster)
manual <- setNames(candidate_table$celltype, candidate_table$cluster)
use_manual <- !is.na(manual) & nzchar(trimws(manual))
labels[use_manual] <- manual[use_manual]
cell_clusters <- as.character(pbmc[[cluster_column]][, 1L])
pbmc$celltype <- factor(unname(labels[cell_clusters]))
pbmc$candidate_lineage <- factor(unname(setNames(lineage_best$candidate, rownames(lineage_best))[cell_clusters]))
pbmc$candidate_state <- factor(unname(setNames(state_best$candidate, rownames(state_best))[cell_clusters]))

pbmc_save_pdf(
  DotPlot(pbmc, features = marker_genes_present, group.by = cluster_column) + RotatedAxis(),
  file.path(figure_dir, "PBMC_v4_cluster_marker_DotPlot.pdf"), 18, 10
)
pbmc_save_pdf(
  DimPlot(pbmc, reduction = umap_reduction, group.by = cluster_column, label = TRUE, repel = TRUE) /
    DimPlot(pbmc, reduction = umap_reduction, group.by = "celltype", label = TRUE, repel = TRUE),
  file.path(figure_dir, "PBMC_v4_annotation_UMAP.pdf"), 14, 16
)
# 同时保存按日期和样本拆分的 UMAP 图；不需要单独运行绘图脚本。
pbmc_make_annotation_plots(
  pbmc, umap_reduction, "celltype", figure_dir,
  "day_group", "sample_id", day_levels, "PBMC_v4_annotation"
)

counts <- pbmc[[]] %>% count(mouse_timepoint_id, mouse_id, day_group, celltype, name = "n")
proportions <- counts %>% group_by(mouse_timepoint_id) %>% mutate(proportion = n / sum(n)) %>% ungroup()
write.csv(counts, file.path(table_dir, "PBMC_v4_celltype_counts_by_mouse_timepoint.csv"), row.names = FALSE)
write.csv(proportions, file.path(table_dir, "PBMC_v4_celltype_proportion_by_mouse_timepoint.csv"), row.names = FALSE)
saveRDS(pbmc, file.path(out_dir, "PBMC_v4_annotated.rds"), compress = "gzip")
capture.output(sessionInfo(), file = file.path(table_dir, "sessionInfo.txt"))
writeLines(c(
  "# PBMC v4 注释", "",
  paste0("- Cluster column: `", cluster_column, "`; UMAP: `", umap_reduction, "`."),
  paste0("- Cells: ", ncol(pbmc), "; clusters: ", length(cluster_levels), "."),
  "- Edit `tables/PBMC_v4_cluster_annotation.csv` and rerun to apply reviewed labels.", ""
), file.path(out_dir, "README.md"))
message("v4 annotation finished. Results: ", out_dir)
