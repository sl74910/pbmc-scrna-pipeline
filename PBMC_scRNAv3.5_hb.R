#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
})

# 读取 v3 Harmony 检查点；其中已包含 umap.harmony。
checkpoint_file <- "outputs/PBMC_scRNAv3_harmony/checkpoints/PBMC_v3_harmony_integrated.RData"
out_dir <- "outputs/PBMC_scRNAv3.5_hb"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
load(checkpoint_file)
DefaultAssay(pbmc) <- "RNA"

# 检查点只保留 counts；合并分层后生成用于 FeaturePlot 的 data layer。
if (length(grep("^counts\\.", Layers(pbmc[["RNA"]]))) > 1L) {
  pbmc <- JoinLayers(pbmc, assay = "RNA")
}
hb_genes <- c("Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2")
pbmc[["percent_hb"]] <- PercentageFeatureSet(pbmc, features = hb_genes)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

erythrocyte_genes <- intersect(
  c("Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Alas2", "Gypa"),
  rownames(pbmc)
)
if (!length(erythrocyte_genes)) stop("No erythrocyte genes found in the object.")

erythrocyte_plot <- FeaturePlot(
  pbmc, features = erythrocyte_genes, reduction = "umap.harmony",
  cols = c("#268591", "#41498B", "#4A115E"),
  min.cutoff = 0, max.cutoff = NA, raster = FALSE
)
ggsave(
  file.path(out_dir, "erythrocyte_gene_featureplot.png"),
  erythrocyte_plot, width = 12, height = 8, dpi = 150
)

percent_hb_plot <- FeaturePlot(
  pbmc, features = "percent_hb", reduction = "umap.harmony",
  cols = c("#268591", "#41498B", "#4A115E"),
  min.cutoff = 0, max.cutoff = 5, raster = FALSE
)
ggsave(
  file.path(out_dir, "percent_hb_featureplot.png"),
  percent_hb_plot, width = 7, height = 6, dpi = 150
)
