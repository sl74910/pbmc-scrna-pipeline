#!/usr/bin/env Rscript

# 最简洁的 DeepCellSeek 注释示例。
# 请在项目根目录运行：source("PBMC_scRNAv4_DeepCellSeek.R")

library(Seurat)
library(DeepCellSeek)
library(ggplot2)
source(file.path(getwd(), "R", "pbmc_helpers.R"))

# 运行前请在当前 shell 或 R 会话中设置 API 密钥；不要把真实密钥写进脚本。
# Shell：export OPENAI_API_KEY="你的密钥"
# R：   Sys.setenv(OPENAI_API_KEY = "你的密钥")
if (!nzchar(Sys.getenv("OPENAI_API_KEY"))) {
  stop("Please set OPENAI_API_KEY before running this script.")
}

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

# 这里使用自定义 OpenAI 兼容站点提供的 gpt-5.6-sol 模型。
celltype_results <- llm_celltype(
  input = markers_df,
  tissuename = "PBMC",
  species = "Mouse",
  model = "gpt-5.6-sol",
  topgenenumber = 10
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
ggsave(file.path(figure_dir, "PBMC_v4_DeepCellSeek_umap.png"),
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
