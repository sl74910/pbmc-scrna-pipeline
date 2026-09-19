#!/usr/bin/env Rscript

# PBMC QC 和双细胞检测：逐个样本计算指标，保存带有两种检测结果的 Seurat 对象列表。

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(DoubletFinder)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(BiocParallel)
})
source(file.path(getwd(), "R", "pbmc_plots.R"))

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- "/home/cylroot/proj_Immune/inputs/PBMC_count_seaurt"
output_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv1_QC_diagnostics")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 这里只设置 QC 阈值；nCount_RNA 和 nFeature_RNA 由 Seurat 自动生成。
qc <- list(
  min_features = 1000L,
  max_features = 6000L,
  max_counts = 30000L,
  max_percent_mt = 10,
  max_percent_hb = 5
)
hb.genes <- c("Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2")
# 为了可复现图中的散点抽样，固定随机种子；它不影响细胞筛选结果。
point_fraction <- 0.05
random_seed <- 1234L

files <- list.files(input_dir, pattern = "^PBMC_.*_Seurat\\.rds$", full.names = TRUE)
if (!length(files)) stop("No Seurat RDS files found in: ", input_dir)
# 正式分析：读取目录中的全部 20 个样本。
sample_id <- sub("_Seurat\\.rds$", "", basename(files))

# 每个输入文件对应一个独立的 Seurat 对象；QC 阶段不合并样本。
objects <- lapply(seq_along(files), function(i) {
  object <- readRDS(files[[i]])
  if (!inherits(object, "Seurat")) stop("Input is not a Seurat object: ", files[[i]])
  # 原始对象是 v3 Assay：LayerData(..., "data") <- NULL 不会删层，只会把 data 重置成 counts。
  # 用 counts 重建对象，才能真正丢掉 data、tsne、umap。
  object <- CreateSeuratObject(
    counts = GetAssayData(object, assay = "RNA", layer = "counts"),
    project = sample_id[[i]]
  )
  object$sample_id <- sample_id[[i]]
  object$day <- sub("^PBMC_([0-9]+)day_.*$", "Day\\1", sample_id[[i]])
  object$percent_mt <- PercentageFeatureSet(object, pattern = "^(mt-|MT-)")
  object$percent_hb <- PercentageFeatureSet(object, features = hb.genes)
  object
})
names(objects) <- sample_id

# 每个样本独立计算 QC 通过标记；结果保存在各自的 metadata 中。
objects <- lapply(objects, function(object) {
  meta <- object[[]]
  object[["qc_pass_pre_hb"]] <- with(meta,
    nFeature_RNA >= qc$min_features & nFeature_RNA <= qc$max_features &
      nCount_RNA <= qc$max_counts & percent_mt <= qc$max_percent_mt
  )
  object[["hb_pass"]] <- !is.na(meta$percent_hb) & meta$percent_hb <= qc$max_percent_hb
  object[["qc_pass"]] <- object$qc_pass_pre_hb & object$hb_pass
  object$qc_pass_pre_hb[is.na(object$qc_pass_pre_hb)] <- FALSE
  object$hb_pass[is.na(object$hb_pass)] <- FALSE
  object$qc_pass[is.na(object$qc_pass)] <- FALSE
  object
})
names(objects) <- sample_id

# 过滤前后 QC 图，以及每个样本过滤前后的细胞数。
save_qc(objects, output_dir, qc, point_fraction, random_seed)
# 保存完整 QC 对象列表：每个元素包含 RNA/counts、细胞 metadata 和基因 metadata。

# 最多同时处理 10 个样本；这个变量要在 Hb 和双细胞两个并行步骤之前定义。
n_workers <- 10L

# Hb 图：在 DoubletFinder 之前，使用“基本 QC 通过但尚未应用 Hb 阈值”的细胞计算 UMAP。
# 函数内部不需要 FindNeighbors/FindClusters；NormalizeData -> PCA -> RunUMAP
# 已经足够把 Hb 含量映射到同一套前后对照坐标上。
hb_plot_dir <- file.path(output_dir, "HB_UMAP")
dir.create(hb_plot_dir, recursive = TRUE, showWarnings = FALSE)
hb_umap_results <- parallel::mclapply(objects, function(object) {
  sample_name <- object$sample_id[[1L]]
  sample_code <- utf8ToInt(sample_name)
  sample_seed <- as.integer(
    (random_seed + sum(sample_code * seq_along(sample_code))) %% .Machine$integer.max
  )
  message("processing Hb UMAP：", sample_name)

  cells_before_hb <- colnames(object)[which(object$qc_pass_pre_hb)]
  if (length(cells_before_hb) < 100L) {
    stop(sample_name, "：基本 QC 通过的细胞少于 100 个，无法计算 Hb UMAP。")
  }
  hb_result <- plot_hb_umap(
    object[, cells_before_hb],
    hb_column = "percent_hb",
    keep_column = "hb_pass",
    sample_name = sample_name,
    hb_max_cutoff = qc$max_percent_hb,
    seed.use = sample_seed
  )
  # 函数返回的是 QC 子对象用于作图；这里把 summary 写回完整对象，
  # 这样后续 DoubletFinder/scDblFinder 仍然能看到所有原始细胞。
  object@misc$HB <- hb_result$object@misc$HB
  object@misc$HB$genes <- hb.genes
  object@misc$HB$threshold_percent <- qc$max_percent_hb
  list(object = object, before = hb_result$before, after = hb_result$after)
}, mc.cores = min(2L, n_workers), mc.preschedule = FALSE)
names(hb_umap_results) <- sample_id
hb_failed <- !vapply(
  hb_umap_results,
  function(result) is.list(result) && inherits(result$object, "Seurat"),
  logical(1)
)
if (any(hb_failed)) {
  stop("Hb UMAP 失败的样本：", paste(names(objects)[hb_failed], collapse = ", "))
}
hb_before_file <- file.path(hb_plot_dir, "HB_before_UMAP.pdf")
hb_after_file <- file.path(hb_plot_dir, "HB_after_UMAP.pdf")
hb_before_plots <- lapply(hb_umap_results, `[[`, "before")
hb_after_plots <- lapply(hb_umap_results, `[[`, "after")
hb_ncol <- min(5L, length(hb_before_plots))
hb_nrow <- as.integer(ceiling(length(hb_before_plots) / hb_ncol))
save_pdf(
  patchwork::wrap_plots(hb_before_plots, ncol = hb_ncol, nrow = hb_nrow, byrow = TRUE),
  hb_before_file, width = 4 * hb_ncol, height = 4 * hb_nrow
)
save_pdf(
  patchwork::wrap_plots(hb_after_plots, ncol = hb_ncol, nrow = hb_nrow, byrow = TRUE),
  hb_after_file, width = 4 * hb_ncol, height = 4 * hb_nrow
)
objects <- lapply(hb_umap_results, function(result) {
  object <- result$object
  object@misc$HB$plots <- list(before = hb_before_file, after = hb_after_file)
  object
})
names(objects) <- sample_id

# DoubletFinder
# 进程非常简单，首先列表里拿一个seurat数据集->取出qc后的细胞重建seurat->Normalize标准化->高变基因->Scale放缩(在采用PCA等降维技术前)->pca->聚类
df_params <- list(
  nfeatures = 3000L,               # 高变基因数量。
  normalization = "LogNormalize", # 按细胞总量归一化、缩放，再取 log1p。
  scale_factor = 10000,            # 对数变换前的总量缩放因子。
  selection_method = "vst",        # 高变基因筛选方法。
  pca_dims = 50L,                  # 最多计算 50 个主成分，供肘部图检查。
  detection_dims = 20L,            # 暂用前 20 个主成分，查看肘部图后可调整。
  resolution = 0.5,                # 预聚类分辨率，用于估计同型双细胞比例。
  rate_per_1000 = 0.008,           # 每千细胞的比例估算系数，需按建库信息核对。
  max_rate = 0.20,                 # 脚本设定的比例上限，并非包默认值。
  pN = 0.25                       # 人工双细胞占混合数据的比例，并非预期双细胞率。
)
df_plot_dir <- file.path(output_dir, "DoubletFinder")
dir.create(df_plot_dir, recursive = TRUE, showWarnings = FALSE)

# 每个进程处理一个 Seurat，返回带结果的对象，保留列表名称。

object <- parallel::mclapply(objects, function(object) {
  sample_name <- object$sample_id[[1L]]
  sample_code <- utf8ToInt(sample_name)
  sample_seed <- as.integer(
    (random_seed + sum(sample_code * seq_along(sample_code))) %% .Machine$integer.max
  )
  set.seed(sample_seed)
  message("processing DoubletFinder：", sample_name)

  # 提取通过 QC 的 counts 建立临时对象，原始对象的细胞和数据保持完整。
  cells <- colnames(object)[which(object$qc_pass)]
  if (length(cells) < 100L) stop(sample_name, "：通过 QC 的细胞少于 100 个。")
  counts <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  temp <- Seurat::CreateSeuratObject(counts[, cells, drop = FALSE], project = sample_name)
  temp <- Seurat::NormalizeData(temp, normalization.method = df_params$normalization,
                                scale.factor = df_params$scale_factor)
  temp <- Seurat::FindVariableFeatures(temp, selection.method = df_params$selection_method,
                                       nfeatures = min(df_params$nfeatures, nrow(temp)))
  # VariableFeatures(temp)          # 查看高变基因名单
  temp <- Seurat::ScaleData(temp, features = Seurat::VariableFeatures(temp))
  npcs <- min(df_params$pca_dims, ncol(temp) - 1L, length(Seurat::VariableFeatures(temp)) - 1L)
  if (npcs < 2L) stop(sample_name, "：可用于 PCA 的细胞或高变基因不足。")
  temp <- Seurat::RunPCA(temp, npcs = npcs, seed.use = sample_seed)
  available_pcs <- ncol(Seurat::Embeddings(temp, "pca"))
  PCs <- seq_len(min(df_params$detection_dims, available_pcs))

  # 保存肘部图，红线标记本次使用的 PCs；预聚类仅用于同型校正。
  elbow <- Seurat::ElbowPlot(temp, ndims = available_pcs) +
    ggplot2::geom_vline(xintercept = max(PCs), colour = "red", linetype = "dashed") +
    ggplot2::labs(title = sample_name, subtitle = paste0("DoubletFinder uses PCs 1:", max(PCs)))
  save_pdf(elbow, file.path(df_plot_dir, paste0(sample_name, "_ElbowPlot.pdf")), 7, 5)
  temp <- Seurat::FindNeighbors(temp, dims = PCs)
  temp <- Seurat::FindClusters(temp, resolution = df_params$resolution, random.seed = sample_seed)

  # 一、选 pK：尝试不同参数，得到每个 pK 的评分表。
  # sct = FALSE：使用普通标准化；GT = FALSE：没有已知的真实双细胞标签。
  sweep <- DoubletFinder::paramSweep(temp, PCs = PCs, sct = FALSE)
  sweep_stats <- DoubletFinder::summarizeSweep(sweep, GT = FALSE)
  # find.pK 会自动画图；这里不保存默认图，并保证结束后关闭绘图设备。
  grDevices::pdf(NULL)
  pK_scan <- tryCatch(DoubletFinder::find.pK(sweep_stats), finally = grDevices::dev.off())
  # 把表里的 pK 转成数字；评分无效的行不参与选择，完整表仍保留。
  pK_scan$pK <- as.numeric(as.character(pK_scan$pK))
  valid_rows <- is.finite(pK_scan$pK) & is.finite(pK_scan$BCmetric)
  valid_scan <- pK_scan[valid_rows, ]
  if (nrow(valid_scan) == 0L) stop(sample_name, "：未找到有效的 pK。")
  # which.max 返回最高分所在的行号，再从这一行取出 pK。
  best_row <- which.max(valid_scan$BCmetric)
  pK <- valid_scan$pK[best_row]

  # 保存你看到的 pK 曲线：横轴是候选 pK，纵轴是评分，红线是选中的值。
  pk_plot <- ggplot2::ggplot(valid_scan, ggplot2::aes(pK, BCmetric)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = pK, colour = "red", linetype = "dashed") +
    ggplot2::labs(title = sample_name, subtitle = paste0("Selected pK = ", pK)) +
    ggplot2::theme_bw()
  save_pdf(pk_plot, file.path(df_plot_dir, paste0(sample_name, "_pK_scan.pdf")), 7, 5)

  # 二、估计数量：例如 5000 个细胞，当前规则估计比例为 4%，数量为 200。
  # 这只是脚本的估算规则，比例系数和上限需按建库信息核对。
  n_cells <- ncol(temp)
  doublet_rate <- df_params$rate_per_1000 * n_cells / 1000
  doublet_rate <- min(df_params$max_rate, doublet_rate)
  expected_doublets <- round(n_cells * doublet_rate)
  # 按官方示例：估计同型比例，再校正预期双细胞数量。
  homotypic.prop <- DoubletFinder::modelHomotypic(temp$seurat_clusters)
  nExp <- round(expected_doublets * (1 - homotypic.prop))
  # 当前包用 1:nExp 选细胞，不能正确处理零个的情况。
  if (nExp == 0L) stop(sample_name, "：校正后 nExp 为 0，请检查预期比例和预聚类结果。")

  # 三、正式检测：用选好的 pK 打分，将排名靠前的 nExp 个细胞标为双细胞。
  temp <- DoubletFinder::doubletFinder(temp, PCs = PCs, pN = df_params$pN,
                                       pK = pK, nExp = nExp, sct = FALSE)
  # temp[[]] 取出细胞信息表；grep 找到列名以 DF.classifications 开头的分类列。
  cell_metadata <- temp[[]]
  call_column <- grep("^DF.classifications", colnames(cell_metadata), value = TRUE)
  calls <- as.character(cell_metadata[[call_column]])  # 每个细胞的 Singlet 或 Doublet 标签。

  # 四、写回原对象：先给全部细胞填 NA，再按 barcode 找到参与检测的细胞位置。
  pass <- rep(NA, ncol(object))  # rep 表示重复：生成与原细胞数相同长度的 NA。
  cell_positions <- match(colnames(temp), colnames(object))
  is_singlet <- calls == "Singlet"  # 标签是 Singlet 得到 TRUE，否则得到 FALSE。
  pass[cell_positions] <- is_singlet
  object$DoubletFinder_pass <- pass
  doublets_called <- sum(calls == "Doublet")
  # 每个细胞的标记放 metadata；整个样本的参数和简要汇总放 misc。
  object@misc$DoubletFinder <- list(pK = pK, nExp = nExp, PCs = PCs, pK_scan = pK_scan,
                                    parameters = df_params, seed = sample_seed,
                                    summary = list(
                                      cells_before_doublet = n_cells,
                                      doublets_called = doublets_called,
                                      doublets_called_percent = round(100 * doublets_called / n_cells, 2L),
                                      singlets_retained = sum(is_singlet)
                                    ))
  object  # 返回更新后的原 Seurat，供 mclapply 组成结果列表；全部 counts 保留。
}, mc.cores = n_workers, mc.preschedule = FALSE)
names(object) <- sample_id

# 子进程失败时停止，不用错误结果覆盖原始列表。
failed <- !vapply(object, inherits, logical(1), what = "Seurat")
if (any(failed)) stop("DoubletFinder 失败的样本：", paste(names(objects)[failed], collapse = ", "))
objects <- object


# scDblFinder：对每个样本的 QC 通过细胞进行第二种双细胞检测。
# 输入：上面 DoubletFinder 完成后的 `objects` 列表；每个元素仍是一个完整的 Seurat 对象。
# 输出：每个对象新增 scDblFinder 的分类、分数、是否保留的 metadata 和简要 summary。

# scDblFinder 官方默认 processing 使用 irlba。当前环境的 Matrix 1.6.5
# 与 irlba 存在 as_cholmod_sparse 兼容问题，因此只替换 PCA 后端；
# 人工双胞生成、kNN、分类和阈值仍由 scDblFinder 官方流程完成。
# 这个函数是 scDblFinder 的 processing 钩子：它只负责把输入计数转换成 PCA 坐标。
# 双细胞的生成、邻居搜索、分类和最终阈值判断仍由 scDblFinder::scDblFinder() 完成。
scdblfinder_processing <- function(e, dims) {
  # 按每个细胞的测序深度做 library-size normalization；这里仍然使用原始 counts。
  normalized <- scuttle::normalizeCounts(e)
  # 使用精确 PCA 后端，绕开当前 Matrix/irlba 的兼容问题。
  # seq_len(nrow(...)) 和 ntop=nrow(...) 表示这里不再额外筛选高变基因。
  pca <- scater::calculatePCA(
    normalized,
    ncomponents = dims,
    subset_row = seq_len(nrow(normalized)),
    ntop = nrow(normalized),
    BSPARAM = BiocSingular::ExactParam()
  )
  # 不同版本的 calculatePCA 可能返回矩阵或带有 `$x` 的结果列表；统一成矩阵。
  if (is.list(pca)) pca <- pca$x
  # PCA 每一行对应一个细胞；显式写回 barcode，后面才能安全映射回 Seurat。
  rownames(pca) <- colnames(e)
  pca
}

# mclapply 每次处理一个样本；函数内部的 `object` 是当前这个样本的完整 Seurat。
sc_results <- parallel::mclapply(objects, function(object) {
  sample_name <- object$sample_id[[1L]]
  # 基于总种子和样本名派生样本级种子：不同样本不同，但同一样本可复现。
  sample_code <- utf8ToInt(sample_name)
  sample_seed <- as.integer(
    (random_seed + sum(sample_code * seq_along(sample_code))) %% .Machine$integer.max
  )
  message("processing scDblFinder：", sample_name)

  # qc_pass 是前面 QC 阶段写入的逻辑列；只有 TRUE 的细胞进入双细胞检测。
  cells <- colnames(object)[which(object$qc_pass)]
  if (length(cells) < 100L) {
    stop(sample_name, "：通过 QC 的细胞少于 100 个。")
  }
  # 先从完整 Seurat 中取出 QC 子对象，再转换为 scDblFinder 官方支持的 SCE。
  # 转换会保留 RNA/counts 层；scDblFinder 仍然从原始 counts 开始计算。
  seurat_qc <- object[, cells]
  sce_input <- Seurat::as.SingleCellExperiment(seurat_qc)
  if (!"counts" %in% SummarizedExperiment::assayNames(sce_input)) {
    stop("Seurat 对象转换后没有 counts assay：", sample_name)
  }
  cells_before_doublet <- ncol(sce_input)
  # SerialParam 保证每个样本内部按固定种子执行随机步骤，结果可复现。
  bp <- BiocParallel::SerialParam(RNGseed = sample_seed)
  # 返回一个 SingleCellExperiment；分类和分数会写在 sce 的 colData 中。
  sce <- scDblFinder::scDblFinder(
    sce_input,
    processing = scdblfinder_processing,
    BPPARAM = bp,
    verbose = FALSE
  )
  # 先检查返回对象和细胞数量，避免错误结果继续覆盖当前对象。
  if (!inherits(sce, "SingleCellExperiment") || ncol(sce) != cells_before_doublet) {
    stop("scDblFinder 返回的对象无效：", sample_name)
  }
  # class 是 singlet/doublet；score 是每个细胞被判为 doublet 的数值分数。
  classification <- as.character(sce$scDblFinder.class)
  score <- as.numeric(sce$scDblFinder.score)
  if (length(classification) != cells_before_doublet ||
      length(score) != cells_before_doublet || anyNA(classification) || anyNA(score) ||
      !all(classification %in% c("doublet", "singlet"))) {
    stop("scDblFinder 返回的分类或分数长度不正确：", sample_name)
  }
  singlets <- colnames(sce)[classification == "singlet"]
  if (!length(singlets)) stop("scDblFinder removed every cell from ", sample_name)

  # 结果写回完整 Seurat：只给参与检测的 QC 通过细胞写值，其他细胞保持 NA。
  # NA 表示“没有运行 scDblFinder”，不是“被判定为 singlet”。
  call_values <- rep(NA_character_, ncol(object))
  score_values <- rep(NA_real_, ncol(object))
  pass_values <- rep(NA, ncol(object))
  # scDblFinder 可能改变内部对象，但 barcode 是细胞的唯一标识，必须按 barcode 对齐。
  positions <- match(colnames(sce), colnames(object))
  if (anyNA(positions) || anyDuplicated(positions)) {
    stop("scDblFinder 条形码无法唯一映射回 Seurat 对象：", sample_name)
  }
  call_values[positions] <- classification
  score_values[positions] <- score
  pass_values[positions] <- classification == "singlet"
  object$scDblFinder_call <- call_values
  object$scDblFinder_score <- score_values
  object$scDblFinder_pass <- pass_values
  doublets_called <- sum(classification == "doublet")
  # misc 保存本次运行的算法说明和参数，便于以后追溯。
  object@misc$scDblFinder <- list(
    input = "Seurat -> SingleCellExperiment",
    processing = "scuttle::normalizeCounts + scater::calculatePCA(ExactParam)",
    dbr = "auto",
    seed = sample_seed,
    summary = list(
      cells_before_doublet = cells_before_doublet,
      doublets_called = doublets_called,
      doublets_called_percent = round(100 * doublets_called / cells_before_doublet, 2L),
      singlets_retained = length(singlets)
    )
  )
  object
}, mc.cores = n_workers, mc.preschedule = FALSE)
names(sc_results) <- sample_id

# 只有每个样本都返回合法的 Seurat，才替换原来的 objects。
# 子进程失败时停止，不用错误结果覆盖 DoubletFinder 已完成的对象列表。
sc_failed <- !vapply(sc_results, inherits, logical(1), what = "Seurat")
if (any(sc_failed)) {
  stop("scDblFinder 失败的样本：", paste(names(objects)[sc_failed], collapse = ", "))
}
objects <- sc_results

# 一个 RDS 同时保留 qc_pass、DoubletFinder_pass 和 scDblFinder_pass。
qc_file <- file.path(output_dir, "PBMC_all_samples_QC_Seurat.rds")
saveRDS(objects, qc_file, compress = TRUE)

message("QC finished: ", qc_file)
