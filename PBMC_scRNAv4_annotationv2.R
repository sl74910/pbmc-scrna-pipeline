#!/usr/bin/env Rscript

# 按 DeepCellSeek 的 marker 规则生成可人工检查的 gene-by-cluster 表。
# 每个 cluster 一个 data.frame：先保留 FindAllMarkers 统计量，再列出所有 cluster 的平均表达。

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
})
source(file.path(getwd(), "R", "pbmc_helpers.R"))

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
checkpoint_file <- file.path(
  project_dir, "outputs", "PBMC_scRNAv3_harmony", "checkpoints",
  "PBMC_v3_harmony_integrated.RData"
)
out_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv4_annotationv2")
table_dir <- file.path(out_dir, "tables")
cluster_column <- "clusters.harmony"
marker_min_pct <- 0.25
marker_logfc_threshold <- 0.25
top_marker_count <- 50L
# DeepCellSeek::process_input_data() 使用的额外过滤阈值。
deepcellseek_p_adj_threshold <- 2.2e-16
normalization_method <- "LogNormalize"
scale_factor <- 10000

pbmc_make_dirs(c(out_dir, table_dir))
if (!file.exists(checkpoint_file)) {
  stop("Checkpoint file does not exist: ", checkpoint_file)
}

pbmc <- pbmc_load_rdata_object(checkpoint_file, "pbmc")
if (!inherits(pbmc, "Seurat")) stop("Checkpoint object `pbmc` is not a Seurat object.")
if (!cluster_column %in% colnames(pbmc[[]])) {
  stop("Cluster metadata column is missing: ", cluster_column)
}
if (!"RNA" %in% Assays(pbmc)) stop("The checkpoint must contain an RNA assay.")
DefaultAssay(pbmc) <- "RNA"
Idents(pbmc) <- cluster_column

# FindAllMarkers 需要单一 counts/data layer；检查点若只有 counts，则生成 normalized data。
if (length(grep("^counts\\.", Layers(pbmc[["RNA"]]))) > 1L) {
  pbmc <- JoinLayers(pbmc, assay = "RNA")
}
if (!"data" %in% Layers(pbmc[["RNA"]])) {
  pbmc <- NormalizeData(
    pbmc, assay = "RNA", normalization.method = normalization_method,
    scale.factor = scale_factor, verbose = FALSE
  )
}

message("Finding positive markers...")
markers_df <- FindAllMarkers(
  pbmc,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = marker_min_pct,
  logfc.threshold = marker_logfc_threshold,
  densify = FALSE
)
write.csv(
  markers_df,
  file.path(table_dir, "PBMC_v4_annotationv2_FindAllMarkers.csv"),
  row.names = FALSE
)

required_marker_columns <- c("gene", "cluster", "avg_log2FC", "p_val_adj")
missing_marker_columns <- setdiff(required_marker_columns, colnames(markers_df))
if (length(missing_marker_columns)) {
  stop("FindAllMarkers output is missing: ", paste(missing_marker_columns, collapse = ", "))
}

# 这一步严格复现 DeepCellSeek::process_input_data：显著性过滤后，
# 每个 cluster 按 avg_log2FC 降序取前 top_marker_count 个基因。
deepcellseek_markers <- markers_df %>%
  mutate(cluster = as.character(cluster)) %>%
  filter(p_val_adj < deepcellseek_p_adj_threshold, avg_log2FC > 0) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = top_marker_count) %>%
  ungroup()
write.csv(
  deepcellseek_markers,
  file.path(table_dir, "PBMC_v4_annotationv2_DeepCellSeek_top50_markers.csv"),
  row.names = FALSE
)

cluster_levels <- unique(as.character(pbmc[[cluster_column]][, 1L]))
cluster_levels <- cluster_levels[order(
  suppressWarnings(as.integer(cluster_levels)), cluster_levels, na.last = TRUE
)]
deepcellseek_markers$cluster <- factor(
  deepcellseek_markers$cluster, levels = cluster_levels
)

# 平均表达使用 RNA assay 的 normalized data；列名统一为实际 cluster ID。
average_expression <- AverageExpression(
  pbmc,
  assays = "RNA",
  group.by = cluster_column,
  layer = "data",
  verbose = FALSE
)$RNA
average_expression <- as.matrix(average_expression)
colnames(average_expression) <- sub(
  paste0("^", cluster_column, "_"), "", colnames(average_expression)
)
colnames(average_expression) <- sub("^g(?=[0-9]+$)", "", colnames(average_expression), perl = TRUE)
if (!setequal(colnames(average_expression), cluster_levels)) {
  stop("AverageExpression cluster names do not match metadata cluster names.")
}
average_expression <- average_expression[, cluster_levels, drop = FALSE]

marker_output_columns <- c(
  "cluster", "gene",
  setdiff(colnames(markers_df), c("cluster", "gene"))
)
marker_expression_list <- setNames(vector("list", length(cluster_levels)), cluster_levels)

for (target_cluster in cluster_levels) {
  target_markers <- deepcellseek_markers %>%
    filter(as.character(cluster) == target_cluster) %>%
    mutate(cluster = as.character(cluster))
  genes <- target_markers$gene
  genes <- genes[!duplicated(genes) & genes %in% rownames(average_expression)]
  target_markers <- target_markers[match(genes, target_markers$gene), , drop = FALSE]

  if (nrow(target_markers) < top_marker_count) {
    warning(
      "Cluster ", target_cluster, " has only ", nrow(target_markers),
      " markers after the exact DeepCellSeek filter (requested ", top_marker_count, ")."
    )
  }

  # 目标 cluster 的表达列放在最前，随后才是其他 cluster，便于人工逐群检查。
  expression_columns <- c(target_cluster, setdiff(cluster_levels, target_cluster))
  expression_table <- as.data.frame(
    average_expression[genes, expression_columns, drop = FALSE],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(expression_table) <- paste0("avg_expr_cluster_", expression_columns)

  marker_table <- target_markers[, marker_output_columns, drop = FALSE]
  marker_table <- marker_table[match(genes, marker_table$gene), , drop = FALSE]
  result <- cbind(marker_table, expression_table)
  rownames(result) <- NULL
  marker_expression_list[[target_cluster]] <- result

  write.csv(
    result,
    file.path(table_dir, paste0("PBMC_v4_cluster_", target_cluster, "_top50_gene_by_cluster.csv")),
    row.names = FALSE
  )
}

saveRDS(
  marker_expression_list,
  file.path(out_dir, "PBMC_v4_annotationv2_gene_by_cluster_top50.rds"),
  compress = "gzip"
)
message(
  "Saved ", length(marker_expression_list), " cluster tables to: ", table_dir,
  "\nSaved R list to: ", file.path(out_dir, "PBMC_v4_annotationv2_gene_by_cluster_top50.rds")
)
