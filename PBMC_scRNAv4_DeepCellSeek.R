#!/usr/bin/env Rscript

# 最简洁的 DeepCellSeek 注释示例。
# 请在项目根目录运行：source("PBMC_scRNAv4_DeepCellSeek.R")

library(Seurat)
library(DeepCellSeek)
library(ggplot2)
source(file.path(getwd(), "R", "pbmc_helpers.R"))



# 读取 Harmony 检查点；文件中包含对象 pbmc。
checkpoint_file <- "outputs/PBMC_scRNAv3_harmony/checkpoints/PBMC_v3_harmony_integrated.RData"
out_dir <- "outputs/PBMC_scRNAv4_DeepCellSeek"
figure_dir <- file.path(out_dir, "figures")
pbmc_make_dirs(c(out_dir, figure_dir))
pbmc <- pbmc_load_rdata_object(checkpoint_file, "pbmc")
DefaultAssay(pbmc) <- "RNA"
Idents(pbmc) <- "clusters.harmony"

# 检查点只保留 counts 时，先生成 data 层。
if (length(grep("^counts\\.", Layers(pbmc[["RNA"]]))) > 1L) {
  pbmc <- JoinLayers(pbmc, assay = "RNA")
}
if (!"data" %in% Layers(pbmc[["RNA"]])) {
  pbmc <- NormalizeData(pbmc, verbose = FALSE)
}

# 找出每个细胞群的标记基因。
markers_df <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

saveRDS(
  markers_df,
  file = "outputs/PBMC_scRNAv4_DeepCellSeek/PBMC_v4_markers.rds"
)
markers_df <- readRDS("outputs/PBMC_scRNAv4_DeepCellSeek/PBMC_v4_markers.rds")

# 默认追加使用外部 GPT-5.6-sol，并请求最大推理强度。
model <- "gpt-6-astra"
Sys.setenv(DEEPCELLSEEK_REASONING_EFFORT = "high")

# 如果不使用 OPENAI_API_KEY 环境变量，可改用：
Sys.setenv(OPENAI_API_KEY = "sk-c749ff1ffcf4be490006a58af3941c48382d61c516707c553c27e95e9b48993b")

Sys.setenv(DEEPCELLSEEK_EXTERNAL_BASE_URL = "https://api.tryaigc.cn")
# 若中转站要求 /v1/responses，可设置：
Sys.setenv(DEEPCELLSEEK_EXTERNAL_ENDPOINT_PATH = "/v1/responses")
api_key_envs <- switch(
  model,
  "kimi-k2.6" = "KIMI_API_KEY",
  "deepseek-v4-flash" = "DEEPSEEK_API_KEY",
  "gpt-5.6-sol" = c("OPENAI_API_KEY", "DEEPCELLSEEK_EXTERNAL_API_KEY"),
  "gpt-6-astra" = c("OPENAI_API_KEY", "DEEPCELLSEEK_EXTERNAL_API_KEY")
)
api_key_values <- Sys.getenv(api_key_envs, unset = "")
if (!any(nzchar(api_key_values))) {
  stop("请先设置以下任一变量：", paste(api_key_envs, collapse = " 或 "), "，再运行此 demo。")
}

allowed_cell_types_file <- file.path("/home/cylroot/proj_Immune/R_packLearn/DeepCellSeek/demo/inputs/PeripheralBlood_celltype.rds")
allowed_cell_types <- readRDS(allowed_cell_types_file)


# 这里使用自定义 OpenAI 兼容站点提供的 gpt-5.6-sol 模型。
celltype_results <- llm_celltype(
  input = markers_df,
  tissuename = "PBMC",
  species = "Mouse",
  model = model,
  topgenenumber = 30,
  wait_indefinitely = TRUE,
  allowed_cell_types = allowed_cell_types
)
dir.create(file.path(out_dir, model), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  celltype_results,
  file = file.path(out_dir, model,"PBMC_scRNAv4_celltype.rds"),
  compress = "xz"
)

# 将细胞群注释写回每个细胞并保存结果。
pbmc$DeepCellSeek_celltype <- unname(
  celltype_results[as.character(Idents(pbmc))]
)

# 保存主 UMAP 图。
plot <- DimPlot(
  pbmc,
  reduction = "umap.harmony",
  group.by = "DeepCellSeek_celltype",
  label = TRUE,
  repel = TRUE
) + NoLegend()
ggsave(file.path(out_dir, model, "PBMC_v4_DeepCellSeek_umap.png"),
       plot, width = 12, height = 8, dpi = 150)

# 保存按日期和样本拆分的 UMAP 图。
pbmc$DeepCellSeek_plot_label <- factor(as.character(pbmc$DeepCellSeek_celltype))
pbmc_make_annotation_plots(
  pbmc,
  "umap.harmony",
  "DeepCellSeek_plot_label",
  figure_dir,
  "day_group",
  "sample_id",
  c("Day0", "Day1", "Day3", "Day5"),
  "PBMC_v4_DeepCellSeek"
)

saveRDS(
  pbmc,
  file.path(out_dir, "PBMC_v4_DeepCellSeek_annotated.rds"),
  compress = "gzip"
)
