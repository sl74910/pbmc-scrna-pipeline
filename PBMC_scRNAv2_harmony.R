#!/usr/bin/env Rscript

# PBMC v2：删除不需要的细胞，合并全部批次并进行 Harmony 整合。
# 输入是 v1 保存的 Seurat 对象列表；输出是一个包含全部批次的 Seurat 对象。

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_file <- file.path(
  project_dir, "outputs", "PBMC_scRNAv1_QC_diagnostics",
  "PBMC_all_samples_QC_Seurat.rds"
)
output_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv2_harmony")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

random_seed <- 1234L
n_variable_features <- 3000L
n_pcs <- 50L
dims_to_use <- 1:30
cluster_resolution <- 0.5
set.seed(random_seed)

if (!file.exists(input_file)) stop("缺少 v1 输入文件：", input_file)
objects <- readRDS(input_file)
if (!is.list(objects) || !length(objects) || is.null(names(objects))) {
  stop("v1 输入必须是带名称的 Seurat 对象列表。")
}
if (!all(vapply(objects, inherits, logical(1), what = "Seurat"))) {
  stop("v1 列表中的每个元素都必须是 Seurat 对象。")
}

# 每个批次重新建立对象，只保留三个筛选标记都为 TRUE 的细胞。
rebuild_sample <- function(object, sample_id) {
  needed <- c("qc_pass", "DoubletFinder_pass", "scDblFinder_pass")
  missing <- setdiff(needed, colnames(object[[]]))
  if (length(missing)) {
    stop(sample_id, " 缺少 metadata：", paste(missing, collapse = ", "))
  }

  metadata <- object[[]]
  keep <- with(metadata,
    !is.na(qc_pass) & qc_pass &
      !is.na(DoubletFinder_pass) & DoubletFinder_pass &
      !is.na(scDblFinder_pass) & scDblFinder_pass
  )
  keep_cells <- colnames(object)[keep]
  if (!length(keep_cells)) stop(sample_id, " 没有细胞通过筛选。")

  counts <- LayerData(object, assay = "RNA", layer = "counts")
  result <- CreateSeuratObject(
    counts = counts[, keep_cells, drop = FALSE],
    project = sample_id,
    meta.data = metadata[keep_cells, , drop = FALSE],
    min.cells = 0,
    min.features = 0
  )

  result$sample_id <- sample_id
  day <- sub("^PBMC_([0-9]+)day_.*$", "\\1", sample_id)
  measurement <- sub("^PBMC_[0-9]+day_([0-9]+)$", "\\1", sample_id)
  result$day <- as.integer(day)
  result$day_group <- factor(
    paste0("Day", day), levels = c("Day0", "Day1", "Day3", "Day5")
  )
  result$measurement <- as.integer(measurement)
  result$measurement_id <- paste0("R", measurement)
  result$mouse_id <- "M1"
  result$mouse_timepoint_id <- paste0("M1_Day", day)
  result$tissue <- "PBMC"
  result
}

sample_objects <- lapply(seq_along(objects), function(i) {
  rebuild_sample(objects[[i]], names(objects)[[i]])
})
names(sample_objects) <- names(objects)

write.csv(
  data.frame(
    sample_id = names(sample_objects),
    cells_before = vapply(objects, function(object) as.integer(ncol(object)), integer(1)),
    cells_after = vapply(sample_objects, function(object) as.integer(ncol(object)), integer(1)),
    stringsAsFactors = FALSE
  ),
  file.path(table_dir, "PBMC_cells_after_filtering.csv"),
  row.names = FALSE
)

# 从这里开始，pbmc 是一个 Seurat 对象，不再是 list。
if (length(sample_objects) == 1L) {
  pbmc <- sample_objects[[1L]]
} else {
  pbmc <- merge(
    x = sample_objects[[1L]],
    y = sample_objects[-1L],
    add.cell.ids = names(sample_objects),
    project = "LPS_PBMC_v2"
  )
}
rm(objects, sample_objects)
gc(verbose = FALSE)

# 这是删除不需要细胞后的起点文件，只保存原始 counts 和 metadata。
filtered_rds <- file.path(output_dir, "PBMC_all_batches_filtered_Seurat.rds")
saveRDS(pbmc, filtered_rds, compress = TRUE)

# 标准预处理：归一化、高变基因、缩放和 PCA。
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize",
                      scale.factor = 10000, verbose = FALSE)
pbmc <- FindVariableFeatures(
  pbmc, selection.method = "vst",
  nfeatures = min(n_variable_features, nrow(pbmc)), verbose = FALSE
)
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc), verbose = FALSE)
n_pcs <- min(n_pcs, ncol(pbmc) - 1L, length(VariableFeatures(pbmc)) - 1L)
if (n_pcs < 2L) stop("细胞或高变基因太少，无法进行 PCA。")
pbmc <- RunPCA(
  pbmc, features = VariableFeatures(pbmc), npcs = n_pcs,
  seed.use = random_seed, verbose = FALSE
)
dims_to_use <- dims_to_use[dims_to_use <= ncol(Embeddings(pbmc, "pca"))]
if (length(dims_to_use) < 2L) stop("PCA 维度太少，无法进行 Harmony。")

pdf(file.path(figure_dir, "PBMC_ElbowPlot.pdf"), width = 7, height = 5)
print(ElbowPlot(pbmc, ndims = ncol(Embeddings(pbmc, "pca"))))
dev.off()

# Harmony 按 sample_id 校正批次；counts 不会被 Harmony 修改。
if (length(unique(pbmc$sample_id)) < 2L) {
  stop("Harmony 至少需要两个批次。")
}
old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
future::plan(future::sequential)
pbmc <- RunHarmony(
  object = pbmc,
  group.by.vars = "sample_id",
  reduction.use = "pca",
  dims.use = dims_to_use,
  reduction.save = "harmony",
  project.dim = FALSE,
  seed.use = random_seed,
  verbose = TRUE
)

pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = dims_to_use, verbose = FALSE)
pbmc <- FindClusters(
  pbmc, resolution = cluster_resolution,
  cluster.name = "clusters.harmony", verbose = FALSE
)
pbmc <- RunUMAP(
  pbmc, reduction = "harmony", dims = dims_to_use,
  reduction.name = "umap.harmony", reduction.key = "UMAPharmony_",
  seed.use = random_seed, verbose = FALSE
)

pdf(file.path(figure_dir, "PBMC_harmony_UMAP.pdf"), width = 12, height = 18)
print(
  DimPlot(pbmc, reduction = "umap.harmony", group.by = "clusters.harmony", label = TRUE) /
    DimPlot(pbmc, reduction = "umap.harmony", group.by = "day_group") /
    DimPlot(pbmc, reduction = "umap.harmony", group.by = "sample_id")
)
dev.off()

# 这是后续注释和绘图使用的主文件：一个包含所有批次的 Seurat 对象。
final_rds <- file.path(output_dir, "PBMC_all_batches_harmony_Seurat.rds")
saveRDS(pbmc, final_rds, compress = TRUE)
capture.output(sessionInfo(), file = file.path(table_dir, "sessionInfo.txt"))
writeLines(c(
  "# PBMC v2 Harmony",
  "",
  "- Input: v1 QC Seurat list.",
  "- Cells kept: qc_pass & DoubletFinder_pass & scDblFinder_pass.",
  paste0("- Filtered starting Seurat RDS: `", basename(filtered_rds), "`."),
  "- Harmony batch variable: sample_id.",
  paste0("- PCA/Harmony dimensions: ", min(dims_to_use), "-", max(dims_to_use), "."),
  paste0("- Final Seurat RDS: `", basename(final_rds), "`."),
  ""
), file.path(output_dir, "README.md"))
message("v2 Harmony 完成：", final_rds)
