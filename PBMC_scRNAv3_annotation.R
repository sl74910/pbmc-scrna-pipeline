#!/usr/bin/env Rscript

# PBMC v3：使用 Harmony 输出和 DeepCellSeek 自动进行 cluster 注释。

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(DeepCellSeek)
  library(ggplot2)
})
source(file.path(getwd(), "R", "pbmc_helpers.R"))

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
checkpoint_file <- file.path(
  project_dir, "outputs", "PBMC_scRNAv2_harmony",
  "PBMC_all_batches_harmony_Seurat.rds"
)
out_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv3_annotation")
figure_dir <- file.path(out_dir, "figures")
model <- Sys.getenv("DEEPCELLSEEK_MODEL", unset = "gpt-6-astra")
marker_min_pct <- 0.25
marker_logfc_threshold <- 0.25
day_levels <- c("Day0", "Day1", "Day3", "Day5")

# 与 PBMC_scRNAv4_DeepCellSeek.R 使用同一套外部接口配置；密钥只从环境变量读取。
Sys.setenv(DEEPCELLSEEK_REASONING_EFFORT = "high")
Sys.setenv(DEEPCELLSEEK_EXTERNAL_BASE_URL = "https://api.tryaigc.cn")
Sys.setenv(DEEPCELLSEEK_EXTERNAL_ENDPOINT_PATH = "/v1/responses")

pbmc_make_dirs(c(out_dir, figure_dir))
if (!file.exists(checkpoint_file)) {
  stop("找不到 Harmony 输出：", checkpoint_file)
}

pbmc <- readRDS(checkpoint_file)
if (!inherits(pbmc, "Seurat")) stop("Harmony 输出不是 Seurat 对象。")
if (!"RNA" %in% Assays(pbmc)) stop("Harmony 输出缺少 RNA assay。")
if (!"clusters.harmony" %in% colnames(pbmc[[]])) {
  stop("Harmony 输出缺少 clusters.harmony。")
}
pbmc_require_reductions(pbmc, "umap.harmony")
pbmc_require_metadata(pbmc, c("day_group", "sample_id"))
DefaultAssay(pbmc) <- "RNA"
Idents(pbmc) <- "clusters.harmony"

# 检查点通常保存分层 counts；FindAllMarkers 需要合并成一个 counts/data layer。
if (length(grep("^counts\\.", Layers(pbmc[["RNA"]]))) > 1L) {
  pbmc <- JoinLayers(pbmc, assay = "RNA")
}
if (!"data" %in% Layers(pbmc[["RNA"]])) {
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize",
                        scale.factor = 10000, verbose = FALSE)
}

message("正在寻找每个 cluster 的 marker...")
markers_df <- FindAllMarkers(
  pbmc,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = marker_min_pct,
  logfc.threshold = marker_logfc_threshold,
  densify = FALSE
)
saveRDS(
  markers_df,
  file.path(out_dir, "PBMC_v3_markers.rds"),
  compress = "gzip"
)
write.csv(
  markers_df,
  file.path(out_dir, "PBMC_v3_markers.csv"),
  row.names = FALSE
)

# DeepCellSeek 使用当前脚本所在环境中的 API 配置；不在这里重复写入密钥。
api_key_envs <- c("OPENAI_API_KEY", "DEEPCELLSEEK_EXTERNAL_API_KEY")
if (!any(nzchar(Sys.getenv(api_key_envs)))) {
  stop("请先设置 OPENAI_API_KEY 或 DEEPCELLSEEK_EXTERNAL_API_KEY。")
}

allowed_cell_types_file <- file.path(
  "/home/cylroot/proj_Immune/R_packLearn/DeepCellSeek/demo/inputs",
  "PeripheralBlood_celltype.rds"
)
if (!file.exists(allowed_cell_types_file)) {
  stop("找不到允许的细胞类型列表：", allowed_cell_types_file)
}
allowed_cell_types <- readRDS(allowed_cell_types_file)

celltype_results <- llm_celltype(
  input = markers_df,
  tissuename = "PBMC",
  species = "Mouse",
  model = model,
  topgenenumber = 30,
  wait_indefinitely = TRUE,
  allowed_cell_types = allowed_cell_types
)
if (is.null(names(celltype_results))) {
  stop("DeepCellSeek 返回结果没有 cluster 名称。")
}
cluster_ids <- unique(as.character(Idents(pbmc)))
missing_clusters <- setdiff(cluster_ids, names(celltype_results))
if (length(missing_clusters)) {
  stop("DeepCellSeek 缺少 cluster 标签：", paste(missing_clusters, collapse = ", "))
}

model_dir <- file.path(out_dir, model)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  celltype_results,
  file.path(model_dir, "PBMC_v3_celltype.rds"),
  compress = "xz"
)
write.csv(
  data.frame(
    cluster = names(celltype_results),
    celltype = unname(as.character(celltype_results)),
    stringsAsFactors = FALSE
  ),
  file.path(model_dir, "PBMC_v3_celltype.csv"),
  row.names = FALSE
)

# 将 DeepCellSeek 的 cluster 标签写回每个细胞。
pbmc$DeepCellSeek_celltype <- unname(
  celltype_results[as.character(Idents(pbmc))]
)
pbmc$celltype <- pbmc$DeepCellSeek_celltype

umap_plot <- DimPlot(
  pbmc,
  reduction = "umap.harmony",
  group.by = "DeepCellSeek_celltype",
  label = TRUE,
  repel = TRUE
) + NoLegend()
ggsave(
  file.path(model_dir, "PBMC_v3_DeepCellSeek_umap.png"),
  umap_plot, width = 12, height = 8, dpi = 150
)

pbmc$DeepCellSeek_plot_label <- factor(as.character(pbmc$DeepCellSeek_celltype))
pbmc_make_annotation_plots(
  pbmc,
  "umap.harmony",
  "DeepCellSeek_plot_label",
  figure_dir,
  "day_group",
  "sample_id",
  day_levels,
  "PBMC_v3_DeepCellSeek"
)

counts <- as.data.frame(
  table(
    day_group = pbmc$day_group,
    sample_id = pbmc$sample_id,
    celltype = pbmc$DeepCellSeek_celltype
  )
)
names(counts)[names(counts) == "Freq"] <- "n"
write.csv(counts, file.path(out_dir, "PBMC_v3_celltype_counts.csv"), row.names = FALSE)

final_file <- file.path(out_dir, "PBMC_v3_DeepCellSeek_annotated.rds")
saveRDS(pbmc, final_file, compress = "gzip")
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "# PBMC v3 DeepCellSeek annotation",
  "",
  "- Input: v2 Harmony Seurat object.",
  "- Cluster column: clusters.harmony.",
  "- UMAP: umap.harmony.",
  paste0("- Model: ", model, "."),
  paste0("- Output: `", basename(final_file), "`."),
  ""
), file.path(out_dir, "README.md"))
message("v3 DeepCellSeek annotation 完成：", final_file)
