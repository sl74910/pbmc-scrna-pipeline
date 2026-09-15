#!/usr/bin/env Rscript

# ============================================================================
# PBMC 单细胞 RNA-seq v3：Harmony 批次校正、聚类和 UMAP 可视化
# ============================================================================
#
# 这个脚本承接 v2 的双细胞筛选和交集保留结果，主要完成以下工作：
#
# 1. 读取每个样本在 v2 中保留下来的细胞条形码；
# 2. 从原始 Seurat 文件中重新提取这些细胞，建立每个样本的 Seurat 对象；
# 3. 合并多个样本，并按照 sample_id 保留样本级 RNA layers；
# 4. 对表达矩阵进行归一化、筛选高变基因、缩放和 PCA；
# 5. 先使用 PCA 做一次未整合的聚类，作为批次校正前的参照；
# 6. 使用 Harmony 在 PCA 坐标上校正 sample_id 带来的技术差异；
# 7. 使用 Harmony 坐标重新计算邻居、聚类和 UMAP；
# 8. 保存图形、参数、日志和可以继续用于下游分析的检查点文件。
#
# 需要特别区分两个概念：
# - Harmony 处理的是细胞的低维坐标，结果保存在 harmony reduction 中；
# - RNA 的 counts/data layers 是表达矩阵的存储方式，JoinLayers() 只负责
#   合并这些 layers，不负责批次校正，也不会重新计算聚类。
#
# Input: original Seurat RDS files and v2 intersection keep lists.
# Output: unintegrated/Harmony UMAPs and compact Harmony checkpoints.
# Harmony corrects sample_id (technical library); day_group remains biological.

suppressPackageStartupMessages({
  # Seurat 和 SeuratObject：创建和操作单细胞对象、归一化、PCA、邻居和聚类。
  library(Seurat)
  library(SeuratObject)
  # harmony：根据指定的批次变量校正 PCA 低维坐标。
  library(harmony)
  # ggplot2 和 patchwork：绘图以及多个图形的组合。
  library(ggplot2)
  library(patchwork)
})
# 载入项目中统一编写的辅助函数：
# - pbmc_helpers.R：样本清单、读取对象、样本合并、预处理和绘图辅助函数；
# - pbmc_checkpoint_helpers.R：以较小体积保存 Seurat 检查点。
source(file.path(getwd(), "R", "pbmc_helpers.R"))
source(file.path(getwd(), "R", "pbmc_checkpoint_helpers.R"))

# ============================================================================
# 一、设置项目路径和分析参数
# ============================================================================
#
# 所有路径都相对于当前工作目录。运行脚本时，当前工作目录应当是 endotol
# 项目目录，这样 getwd() 才能正确找到 inputs、outputs 和 R 文件夹。
#
# 这些参数会被写入 checkpoint metadata 和 README，方便之后复查本次分析使用的
# 归一化方法、PCA 维度、Harmony 批次变量和聚类分辨率。若参数设置为 NA 或不
# 符合基本要求，脚本会尽早报错，而不是运行很久后才发现参数问题。
#
# Analysis parameters.  Set a parameter to NA to get an immediate error.
project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
# 原始 Seurat RDS 文件所在目录。
input_dir <- file.path(project_dir, "inputs")
# v2 生成的高质量细胞条形码列表所在目录。
keep_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv2_intersection_summary")
# v3 的输出目录，以及图形、表格和 checkpoint 子目录。
out_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv3_harmony")
figure_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")
checkpoint_dir <- file.path(out_dir, "checkpoints")
# days 和 measurements 用于构造样本清单。这里分析 Day0、Day1、Day3、Day5，
# 每个时间点包含 1 到 5 号测量样本；biological_mouse_id 指定生物学来源。
days <- c(0L, 1L, 3L, 5L)
measurements <- 1:5
biological_mouse_id <- "M1"
# 固定随机种子，使 PCA、UMAP、聚类和 Harmony 尽量可以重复。
random_seed <- 1234L
# 在 PCA 前选择的高变基因数量，以及需要计算的 PCA 维度数量。
n_variable_features <- 3000L
n_pcs <- 50L
# 后续邻居、聚类、UMAP 和 Harmony 实际使用前 30 个 PCA/Harmony 维度。
dims_to_use <- 1:30
# Louvain/Leiden 聚类的分辨率。分辨率越高，通常会得到更多较小的细胞群。
cluster_resolution <- 0.5
# Harmony 用这一列 metadata 判断技术批次。这里使用每个测序文库的 sample_id。
harmony_batch_variable <- "sample_id"
# RNA 归一化方法和每个细胞的目标总量。
normalization_method <- "LogNormalize"
scale_factor <- 10000
# FindVariableFeatures 使用的高变基因选择方法。
variable_feature_method <- "vst"
# 检查参数的类型、范围和必要字段是否完整。
pbmc_validate_parameters(list(
  days = days, measurements = measurements, biological_mouse_id = biological_mouse_id,
  random_seed = random_seed, n_variable_features = n_variable_features, n_pcs = n_pcs,
  dims_to_use = dims_to_use, cluster_resolution = cluster_resolution,
  harmony_batch_variable = harmony_batch_variable, normalization_method = normalization_method,
  scale_factor = scale_factor, variable_feature_method = variable_feature_method
))
# 再做一次维度相关检查，保证后面请求的 PCA 维度是合法的。
if (n_pcs < 2L || n_variable_features < 2L || length(dims_to_use) < 2L || min(dims_to_use) < 1L) {
  stop("Harmony parameters are invalid.")
}
# 创建输出目录；如果目录已经存在，不会删除其中已有文件。
pbmc_make_dirs(c(figure_dir, table_dir, checkpoint_dir))
# 确认所需 R 包已经安装，并固定随机数种子。
pbmc_require_packages(c("Seurat", "SeuratObject", "harmony", "ggplot2", "patchwork"))
set.seed(random_seed)

# ============================================================================
# 二、建立样本清单并检查输入文件
# ============================================================================
#
# pbmc_make_manifest() 会根据输入目录、v2 保留细胞列表和样本信息，建立一张
# manifest 表。表中包含原始 Seurat 文件、v2 保留细胞列表、sample_id 和日期等
# 信息。后面的读取步骤全部依赖这张表，因此先检查 manifest，可以尽早发现
# 文件缺失、样本名称不一致或样本元数据不完整的问题。
manifest <- pbmc_make_manifest(
  input_dir, keep_dir, days, measurements, biological_mouse_id,
  "_high_quality_cells.txt"
)
# 检查 manifest 中的文件是否存在、每个样本是否有对应的保留细胞列表。
pbmc_check_manifest(manifest, TRUE)
# 保存一份可直接查看的样本清单，便于核对本次究竟分析了哪些样本。
write.csv(manifest, file.path(table_dir, "PBMC_sample_manifest.csv"), row.names = FALSE)

# ============================================================================
# 三、根据 v2 保留细胞列表重建每个样本的 Seurat 对象
# ============================================================================
#
# v2 已经完成了质量控制、双细胞识别和交集筛选。这里不再把所有细胞读入后
# 重新筛选，而是直接读取每个样本对应的 v2 保留条形码，再从原始 Seurat 对象
# 中提取这些细胞。这样可以保证 v3 使用的细胞集合与 v2 最终结果一致。
#
# Rebuild one Seurat object from the v2-intersection barcode list.
read_sample <- function(index) {
  # 从 manifest 取出当前样本这一行；drop = FALSE 保证结果仍然是数据框。
  row <- manifest[index, , drop = FALSE]
  # 读取 v2 保留的细胞条形码。文件不存在或内容为空时，辅助函数会报错。
  keep <- pbmc_read_cells(row$keep_file[[1L]], paste0("v2 keep list for ", row$sample_id[[1L]]))
  # 从原始 Seurat 文件中只读取 keep 中列出的细胞和原始 counts。
  object <- pbmc_read_seurat_counts(row$file[[1L]], keep, row$sample_id[[1L]], NULL)
  # 将 manifest 中的样本、日期等信息写入每个细胞的 metadata。
  pbmc_add_manifest_metadata(object, row, paste0("Day", days))
}
# 对 manifest 中的每一行执行 read_sample()，生成一个样本一个 Seurat 对象。
sample_objects <- lapply(seq_len(nrow(manifest)), read_sample)
# 用 sample_id 作为列表名称，后面合并时会用这些名称给细胞添加前缀。
names(sample_objects) <- manifest$sample_id
# 统计每个样本最终保留下来的细胞数，作为 v2 到 v3 的输入核对表。
write.csv(
  data.frame(sample_id = names(sample_objects),
             singlets = vapply(sample_objects, function(object) as.integer(ncol(object)), integer(1))),
  file.path(table_dir, "PBMC_v2_intersection_retained_counts.csv"), row.names = FALSE
)

# ============================================================================
# 四、合并样本并保存 v2 筛选后的初始检查点
# ============================================================================
#
# pbmc_merge_samples() 会完成两件事：
# 1. 合并多个样本的 Seurat 对象，并通过 add.cell.ids 保证细胞名称不重复；
# 2. 当存在多个样本时，按照 sample_id 把 RNA assay 拆成多个 counts layers。
#
# 这里的 split layer 只是 Seurat v5 对多样本表达矩阵的存储方式，例如：
# counts.sample1、counts.sample2。它不是 Harmony 批次校正，也不代表细胞分群。
pbmc <- pbmc_merge_samples(sample_objects, "LPS_PBMC_v3", "sample_id")
# 合并完成后不再需要单独的样本对象，删除它们以释放内存。
rm(sample_objects)
# 主动进行垃圾回收，减少后续 PCA/Harmony 计算时的内存压力。
gc(verbose = FALSE)
# 把本次分析的重要参数集中记录下来，写入每个 checkpoint 的 metadata。
analysis_parameters <- list(
  random_seed = random_seed, n_variable_features = n_variable_features, n_pcs = n_pcs,
  dims_to_use = dims_to_use, cluster_resolution = cluster_resolution,
  harmony_batch_variable = harmony_batch_variable, normalization_method = normalization_method,
  scale_factor = scale_factor, variable_feature_method = variable_feature_method,
  source_v2_intersection = normalizePath(keep_dir, mustWork = FALSE)
)
# 保存一个“v2 筛选后、尚未预处理”的检查点。
# 这里只保存 counts，是为了保留最原始的表达数据并控制文件大小；后面的
# data、scale.data 和降维结果都可以在相应阶段重新计算。
save_pbmc_rdata_checkpoint(
  pbmc, "PBMC_v3_after_v2_intersection", checkpoint_dir,
  analysis_parameters, layers = "counts", reductions = character(), extra_metadata = list()
)

# ============================================================================
# 五、RNA 预处理：归一化、高变基因、缩放和 PCA
# ============================================================================
#
# pbmc_preprocess() 在辅助文件中封装了标准 Seurat 预处理流程：
# - NormalizeData：把每个细胞的 counts 转换为可比较的 log-normalized data；
# - FindVariableFeatures：找出最能反映细胞差异的高变基因；
# - ScaleData：对用于 PCA 的基因进行中心化和缩放；
# - RunPCA：把高维基因表达压缩成低维 PCA 坐标。
#
# 此时 PCA 仍然包含样本之间的技术差异，因此它是“未整合”的 PCA。
pbmc <- pbmc_preprocess(
  pbmc, normalization_method, scale_factor, variable_feature_method,
  n_variable_features, n_pcs, random_seed
)
# 如果请求的最大维度超过实际计算出的 PCA 维度，立即停止并提示参数错误。
if (max(dims_to_use) > ncol(Embeddings(pbmc, "pca"))) {
  stop("dims_to_use contains a PC that was not computed.")
}
# 保存 elbow plot，用于观察不同 PCA 维度的方差解释情况，辅助判断维度选择。
pbmc_save_pdf(ElbowPlot(pbmc, ndims = ncol(Embeddings(pbmc, "pca"))),
              file.path(figure_dir, "PBMC_ElbowPlot.pdf"), 7, 5)

# ============================================================================
# 六、基于未整合 PCA 的邻居、聚类和 UMAP
# ============================================================================
#
# 这一步故意不使用 Harmony，而是先用原始 PCA 坐标计算一次结果，作为批次
# 校正前的参照。之后可以比较 unintegrated UMAP 和 harmony UMAP，判断 Harmony
# 是否减少了 sample_id 造成的分离，以及是否保留了 day_group 等生物学变化。
#
# FindNeighbors：根据 PCA 坐标寻找每个细胞的近邻；
# FindClusters：根据近邻图识别细胞群；
# RunUMAP：把 PCA 的多维结构压缩到二维，便于绘图。
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = dims_to_use)
pbmc <- FindClusters(pbmc, resolution = cluster_resolution, cluster.name = "clusters.unintegrated")
pbmc <- RunUMAP(
  pbmc, reduction = "pca", dims = dims_to_use,
  reduction.name = "umap.unintegrated", reduction.key = "UMAPunint_", seed.use = random_seed
)
# 生成未整合 UMAP：分别查看聚类、样本和日期分组在图上的分布。
pbmc_save_pdf(
  pbmc_make_umap_panels(pbmc, "umap.unintegrated", "clusters.unintegrated", "sample_id", "day_group"),
  file.path(figure_dir, "PBMC_unintegrated_UMAP.pdf"), 12, 18
)
# 保存未整合结果。这里保留 counts 和 data，是为了让这个阶段的对象可以直接
# 用于检查表达和绘图；同时保存 PCA、UMAP 和未整合聚类结果。
save_pbmc_rdata_checkpoint(
  pbmc, "PBMC_v3_unintegrated_clustering", checkpoint_dir, analysis_parameters,
  layers = c("counts", "data"), reductions = c("pca", "umap.unintegrated"),
  extra_metadata = list()
)

# ============================================================================
# 七、检查 Harmony 的批次变量并执行低维批次校正
# ============================================================================
#
# Harmony 不是直接修改 RNA counts，也不是把不同样本的表达矩阵简单相加。
# 它读取 PCA 坐标，并根据 harmony_batch_variable 指定的 metadata 列调整细胞
# 在低维空间中的位置，使不同技术批次之间更容易比较，同时尽量保留真实的
# 细胞状态差异。这里校正的是 sample_id，day_group 保留为后续观察的生物学变量。
#
# 在运行 Harmony 前先检查：
# 1. sample_id 确实存在于细胞 metadata；
# 2. sample_id 至少包含两个不同样本，否则不存在可校正的批次差异。
if (!harmony_batch_variable %in% colnames(pbmc[[]])) {
  stop("Harmony batch column is missing: ", harmony_batch_variable)
}
if (length(unique(as.character(pbmc[[harmony_batch_variable]][, 1L]))) < 2L) {
  stop("Harmony requires at least two technical samples.")
}
# 记录调用 Harmony 前的 future 计划，并临时使用顺序执行，避免并行过程带来
# 不可控的内存占用或随机性。脚本退出时 on.exit() 会恢复原来的计划。
old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
future::plan(future::sequential)
# RunHarmony 的主要输入和输出：
# - reduction.use = "pca"：以 PCA 坐标作为待校正的输入；
# - group.by.vars = "sample_id"：按样本识别技术批次；
# - reduction.save = "harmony"：把校正后的坐标保存为 harmony reduction；
# - project.dim = FALSE：不把校正坐标重新投影到所有基因，而是直接使用 Harmony 坐标。
pbmc <- RunHarmony(
  object = pbmc, group.by.vars = harmony_batch_variable,
  reduction.use = "pca", dims.use = dims_to_use,
  reduction.save = "harmony", project.dim = FALSE,
  seed.use = random_seed, verbose = TRUE
)

# ============================================================================
# 八、基于 Harmony 坐标重新聚类和绘制整合后的 UMAP
# ============================================================================
#
# 下面三步与未整合分析类似，但输入从 reduction = "pca" 改为
# reduction = "harmony"。因此这里得到的 clusters.harmony 和 umap.harmony
# 才是本脚本用于后续注释的批次校正后结果。
#
# 注意：JoinLayers() 不需要在这里执行。JoinLayers 只改变 RNA layers 的存储
# 形式，而 Harmony、邻居、聚类和 UMAP 使用的是已经保存在对象中的低维坐标。
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = dims_to_use)
pbmc <- FindClusters(pbmc, resolution = cluster_resolution, cluster.name = "clusters.harmony")
pbmc <- RunUMAP(
  pbmc, reduction = "harmony", dims = dims_to_use,
  reduction.name = "umap.harmony", reduction.key = "UMAPharmony_", seed.use = random_seed
)
# 组合三张图：
# 1. 按 Harmony 聚类着色，检查细胞群结构；
# 2. 按 day_group 着色，观察时间/生物学变化；
# 3. 按 sample_id 着色，观察技术样本是否仍然强烈分离。
pbmc_save_pdf(
  DimPlot(pbmc, reduction = "umap.harmony", group.by = "clusters.harmony", label = TRUE) /
    DimPlot(pbmc, reduction = "umap.harmony", group.by = "day_group") /
    DimPlot(pbmc, reduction = "umap.harmony", group.by = "sample_id"),
  file.path(figure_dir, "PBMC_harmony_UMAP.pdf"), 12, 18
)
# ============================================================================
# 九、保存最终 Harmony 检查点和分析记录
# ============================================================================
#
# 最终 checkpoint 保存：
# - counts：原始 RNA counts；
# - pca：未整合 PCA；
# - umap.unintegrated：未整合 UMAP；
# - harmony：Harmony 校正后的低维坐标；
# - umap.harmony：基于 Harmony 的 UMAP。
#
# 这里 layers = "counts" 是有意的紧凑保存策略。由于输入样本在前面按 sample_id
# 分成了多个 counts layers，最终文件保留这些原始分样本 counts，但不保存可由
# counts 重新生成的 data/scale.data。v4 做 marker、ScType 或功能分析时，才会
# 根据需要 JoinLayers() 并重新 NormalizeData()。
save_pbmc_rdata_checkpoint(
  pbmc, "PBMC_v3_harmony_integrated", checkpoint_dir, analysis_parameters,
  layers = "counts", reductions = c("pca", "umap.unintegrated", "harmony", "umap.harmony"),
  extra_metadata = list()
)
# 保存 R 和所有已加载包的版本信息，便于日后复现或排查版本差异。
capture.output(sessionInfo(), file = file.path(table_dir, "sessionInfo.txt"))
# 写入一个简短的 README，记录批次变量、使用维度、高变基因数和聚类分辨率。
# 这里明确说明 day_group 只用于观察，并没有作为 Harmony 的校正变量。
writeLines(c(
  "# PBMC v3 Harmony integration", "",
  paste0("- Batch variable: `", harmony_batch_variable, "` (technical library)."),
  paste0("- PCA/integration dimensions: ", min(dims_to_use), "-", max(dims_to_use), "."),
  paste0("- HVGs: ", n_variable_features, "; clustering resolution: ", cluster_resolution, "."),
  "- day_group is retained for visualization and is not corrected.", ""
), file.path(out_dir, "README.md"))
# 所有步骤完成后，在终端打印输出目录，方便快速定位结果。
message("v3 Harmony finished. Results: ", out_dir)
