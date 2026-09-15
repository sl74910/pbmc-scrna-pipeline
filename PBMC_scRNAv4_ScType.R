#!/usr/bin/env Rscript

# ScType 注释单细胞RNA-seq
# 输入： Harmony 降维聚类后的分群和umap +
# ScType 数据库（基因list，对应各种类型的细胞应该表达或不应该表达）

# 输出： cluster 标签、注释对象和所有 UMAP 图。

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(HGNChelper)
  library(openxlsx)
})
source(file.path(getwd(), "R", "pbmc_helpers.R"))
# 读取工作目录
project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
# Harmony 整合后的 Seurat 对象文件
checkpoint_file <- file.path(
  project_dir, "outputs", "PBMC_scRNAv3_harmony", "checkpoints",
  "PBMC_v3_harmony_integrated.RData"
)
# 设置R包sc_type的输入文件目录，output输出目录
sc_type_dir <- file.path(project_dir, "down", "sc-type-master")
out_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv4_ScType")
figure_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")
# 设置读入seurat中存储cluster列的列名，指定使用harmony批次整合后的降维结果
cluster_column <- "clusters.harmony"
umap_reduction <- "umap.harmony"
# 设置ScType 使用“免疫系统”对应的 marker 基因list，他有很多list
sc_type_tissue <- "Immune system"
# 其中的一个图按照时间排序
day_levels <- c("Day0", "Day1", "Day3", "Day5")
# 达矩阵的归一化方式为 Seurat 的 LogNormalize，参数为10000，设置seed方便浮现
normalization_method <- "LogNormalize"
scale_factor <- 10000
random_seed <- 1234L
# 做参数合法性检查，有空参数，空文件报错
pbmc_validate_parameters(list(
  checkpoint_file = checkpoint_file, sc_type_dir = sc_type_dir,
  cluster_column = cluster_column, umap_reduction = umap_reduction,
  sc_type_tissue = sc_type_tissue, day_levels = day_levels,
  normalization_method = normalization_method, scale_factor = scale_factor,
  random_seed = random_seed
))
pbmc_make_dirs(c(out_dir, figure_dir, table_dir))
if (!file.exists(checkpoint_file)) stop("Missing Harmony checkpoint: ", checkpoint_file)
set.seed(random_seed)
# gene_sets_prepare.R：读取、清洗并准备正/负 marker 基因集。
# sctype_score_.R：根据 marker 基因集计算 ScType 评分。
for (file_name in c("gene_sets_prepare.R", "sctype_score_.R")) {
  file <- file.path(sc_type_dir, "R", file_name)
  if (!file.exists(file)) stop("Missing ScType file: ", file)
  source(file)
}
db_file <- file.path(sc_type_dir, "ScTypeDB_full_mousePB.xlsx")
if (!file.exists(db_file)) stop("Missing ScType database: ", db_file)
# 读入pbmc的Rdata，设置使用RNA_count
pbmc <- pbmc_load_rdata_object(checkpoint_file, "pbmc")
if (!inherits(pbmc, "Seurat")) stop("Checkpoint object `pbmc` is not a Seurat object.")
pbmc_require_metadata(pbmc, c(cluster_column, "day_group", "sample_id"))
pbmc_require_reductions(pbmc, umap_reduction)
if (!"RNA" %in% Assays(pbmc)) stop("The checkpoint must contain an RNA assay.")
DefaultAssay(pbmc) <- "RNA"

# 从 ScType Excel 数据库 db_file 中读取指定组织类型的 marker 基因。此时组织类型是 "Immune system"，有正向marker和反向marker
gs_list <- gene_sets_prepare(db_file, sc_type_tissue)
gs_list$gs_positive <- lapply(gs_list$gs_positive, function(x) unique(toupper(x)))
gs_list$gs_negative <- lapply(gs_list$gs_negative, function(x) unique(toupper(x)))
# 此时会得到261个marker基因
marker_keys <- unique(toupper(unlist(c(gs_list$gs_positive, gs_list$gs_negative))))
gene_names <- rownames(pbmc)
# 取交集，此时211个基因
marker_genes <- gene_names[match(marker_keys, toupper(gene_names))]
marker_genes <- unique(marker_genes[!is.na(marker_genes)])
if (length(marker_genes) < 10L) stop("Too few ScType markers match the object.")
# 如果现在的批次没有聚合，则聚合
if (length(grep("^counts\\.", Layers(pbmc[["RNA"]]))) > 1L) {
  pbmc <- JoinLayers(pbmc, assay = "RNA")
}
# 使用 Seurat 官方 NormalizeData 对完整 RNA counts 做归一化。
pbmc <- Seurat::NormalizeData(
  pbmc, assay = "RNA", normalization.method = normalization_method,
  scale.factor = scale_factor, verbose = FALSE
)
# 只从官方生成的 data layer 中提取 ScType marker。
expression <- SeuratObject::LayerData(pbmc, assay = "RNA", layer = "data")
expression <- expression[marker_genes, , drop = FALSE]
# 检查转换为大写后是否产生重复的基因名，如果有则对同名基因的表达值按行求和
rownames(expression) <- toupper(rownames(expression))
if (anyDuplicated(rownames(expression))) {
  expression <- rowsum(as.matrix(expression), rownames(expression), reorder = FALSE)
}
# 再次检查是否只保留了markergene
marker_genes <- intersect(toupper(marker_genes), rownames(expression))
# 基因表中取交集和我们的数据集
gs_list$gs_positive <- lapply(gs_list$gs_positive, intersect, y = marker_genes)
gs_list$gs_negative <- lapply(gs_list$gs_negative, intersect, y = marker_genes)
# 行mark，列细胞
expression <- expression[marker_genes, , drop = FALSE]
# 算出来每个基因在某个聚类里面所有细胞的平均表达
clusters <- as.character(pbmc[[cluster_column]][, 1L])
cluster_levels <- unique(clusters)
cluster_expression <- do.call(cbind, lapply(cluster_levels, function(cluster) {
  Matrix::rowMeans(expression[, clusters == cluster, drop = FALSE])
}))
colnames(cluster_expression) <- cluster_levels
# 保存 ScType 的输入表：marker 基因在行、cluster 在列，数值为各 cluster 的平均表达。
write.csv(
  as.data.frame(cluster_expression, check.names = FALSE),
  file.path(table_dir, "PBMC_v4_ScType_cluster_expression.csv"),
  row.names = TRUE
)
# 得到一个矩阵，关于每一个聚类对于每一种细胞类型的打分
scores <- sctype_score(
  scRNAseqData = as.matrix(cluster_expression), scaled = TRUE,
  gs = gs_list$gs_positive, gs2 = gs_list$gs_negative
)
if (!nrow(scores) || !ncol(scores)) stop("ScType produced no usable scores.")
# 保存 ScType 的完整评分表：细胞类型在行、cluster 在列。
write.csv(
  as.data.frame(scores, check.names = FALSE),
  file.path(table_dir, "PBMC_v4_ScType_scores.csv"),
  row.names = TRUE
)
# 只保存每个聚类的最高分的细胞类型
top <- do.call(rbind, lapply(seq_len(ncol(scores)), function(i) {
  sorted <- sort(scores[, i], decreasing = TRUE)
  data.frame(
    cluster = colnames(scores)[i], type = names(sorted)[1L], score = unname(sorted[1L]),
    stringsAsFactors = FALSE
  )
}))
top$ncells <- as.integer(table(factor(clusters, levels = top$cluster)))

pbmc_sctype <- pbmc
pbmc_sctype$sctype_label <- unname(setNames(top$type, top$cluster)[clusters])
write.csv(top, file.path(out_dir, "PBMC_v4_ScType_cluster_labels.csv"), row.names = FALSE)

# 保存主 UMAP 图，以及按日期和样本拆分的 UMAP 图。
umap_plot <- DimPlot(pbmc_sctype, reduction = umap_reduction, group.by = "sctype_label",
                     label = TRUE, repel = TRUE) + NoLegend()
ggsave(file.path(figure_dir, "PBMC_v4_ScType_umap.png"), umap_plot, width = 12, height = 8, dpi = 150)
pbmc_sctype$sc_type_plot_label <- factor(as.character(pbmc_sctype$sctype_label))
pbmc_make_annotation_plots(
  pbmc_sctype, umap_reduction, "sc_type_plot_label", figure_dir,
  "day_group", "sample_id", day_levels, "PBMC_v4_ScType"
)
save(pbmc_sctype, top, cluster_expression, scores,
     file = file.path(out_dir, "PBMC_v4_ScType_annotated.RData"), compress = TRUE)
message(
  "v4 ScType finished. Matched markers: ", length(marker_genes), "; clusters: ", nrow(top),
  "; tables: ", table_dir
)
