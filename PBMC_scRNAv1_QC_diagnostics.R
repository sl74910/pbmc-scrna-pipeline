#!/usr/bin/env Rscript

# PBMC 单细胞 RNA-seq 的基础质量控制（QC）。
#
# 输入：inputs/PBMC_*_Seurat.rds
# 输出：
#   - QC_all_samples_before.pdf：过滤前的 QC 指标分布；
#   - QC_all_samples_after.pdf ：过滤后的 QC 指标分布；
#   - *_high_quality_cells.txt：每个样本保留的细胞条形码；
#   - QC_cell_counts.csv         ：每个样本过滤前后的细胞数量。
#
# 图中的红色虚线就是下面定义的筛选阈值。箱线图使用全部细胞，散点仅
# 随机抽取一部分用于展示，以免细胞很多时 PDF 过大、难以阅读。

source(file.path(getwd(), "R", "pbmc_helpers.R"))
pbmc_require_packages(c("Seurat", "ggplot2", "patchwork"))

# ---- 1. 路径和 QC 规则 ---------------------------------------------------

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(project_dir, "inputs")
output_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv1_QC_diagnostics")

# 所有样本共用同一套阈值，便于横向比较。若要调整 QC，只修改这一段。
qc <- list(
  min_features = 1000L,   # 每个细胞至少检测到的基因数
  max_features = 6000L,   # 过高时可能是双细胞或异常细胞
  max_counts = 30000L,    # 过高的 UMI 总数通常提示异常捕获
  max_percent_mt = 10,    # 线粒体转录本比例上限（%）
  max_percent_hb = 5      # 血红蛋白基因转录本比例上限（%）
)

# 这些定义决定 percent_mt 和 percent_hb 两个 QC 指标如何计算。
mitochondrial_pattern <- "^(mt-|MT-)"
hemoglobin_genes <- c("Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2")

# 为了可复现图中的散点抽样，固定随机种子；它不影响细胞筛选结果。
point_fraction <- 0.05
random_seed <- 1234L

pbmc_validate_parameters(c(
  qc,
  list(
    mitochondrial_pattern = mitochondrial_pattern,
    hemoglobin_genes = hemoglobin_genes,
    point_fraction = point_fraction,
    random_seed = random_seed
  )
))
if (qc$min_features > qc$max_features ||
    qc$max_counts <= 0 ||
    qc$max_percent_mt < 0 ||
    qc$max_percent_hb < 0 ||
    point_fraction <= 0 ||
    point_fraction > 1) {
  stop("QC thresholds or point_fraction are invalid.")
}

pbmc_make_dirs(output_dir)
sample_files <- pbmc_sample_files(input_dir, "^PBMC_.*_Seurat\\.rds$")

# ---- 2. 读取样本并计算 QC 指标 ------------------------------------------

read_qc_sample <- function(file, sample_id) {
  object <- readRDS(file)
  if (!inherits(object, "Seurat")) {
    stop("Input is not a Seurat object: ", file)
  }

  # 保留样本和时间点信息，便于后续按样本或时间点作图、汇总。
  sample_label <- sub("^PBMC_", "", sample_id)
  object$orig.ident <- sample_label
  object$group <- sub("_.*$", "", sample_label)

  # 在 Seurat 元数据中新增 percent_mt 和 percent_hb；nFeature_RNA 与
  # nCount_RNA 已由原始 Seurat 对象提供。
  pbmc_add_qc_metadata(object, mitochondrial_pattern, hemoglobin_genes)
}

samples <- lapply(seq_len(nrow(sample_files)), function(index) {
  read_qc_sample(sample_files$file[[index]], sample_files$sample_id[[index]])
})
names(samples) <- sample_files$sample_id

# 每张小图对应一个样本，四个分面依次显示基因数、UMI 数、线粒体和
# 血红蛋白比例。两次调用共用同一函数，保证过滤前后图形完全可比。
make_qc_plots <- function(objects) {
  lapply(names(objects), function(sample_id) {
    pbmc_qc_plot(objects[[sample_id]], sample_id, qc, point_fraction)
  })
}

set.seed(random_seed)
before_plots <- make_qc_plots(samples)
pbmc_save_plot_grid(
  before_plots,
  file.path(output_dir, "QC_all_samples_before.pdf"),
  ncol = 5,
  nrow = 4,
  width = 32,
  height = 18
)

# ---- 3. 按图中阈值筛选细胞 ----------------------------------------------

# 筛选规则与 QC 图的红色虚线使用同一个 qc 对象，避免显示和实际筛选不一致。
filtered_samples <- lapply(samples, function(object) {
  pbmc_filter_qc(
    object,
    min_features = qc$min_features,
    max_features = qc$max_features,
    max_counts = qc$max_counts,
    max_percent_mt = qc$max_percent_mt,
    max_percent_hb = qc$max_percent_hb
  )
})

# ---- 4. 导出保留细胞及汇总结果 ------------------------------------------

write_retained_cells <- function(sample_id, object) {
  writeLines(
    colnames(object),
    file.path(output_dir, paste0(sample_id, "_high_quality_cells.txt"))
  )
}

# 每个文件一行一个细胞条形码，供后续双细胞识别脚本直接读取。
Map(write_retained_cells, names(filtered_samples), filtered_samples)

cells_before <- vapply(samples, ncol, integer(1))
cells_after <- vapply(filtered_samples, ncol, integer(1))
qc_summary <- data.frame(
  sample_id = names(samples),
  cells_before_qc = cells_before,
  cells_after_qc = cells_after,
  cells_removed_qc = cells_before - cells_after,
  percent_removed_qc = 100 * (cells_before - cells_after) / cells_before,
  stringsAsFactors = FALSE
)
write.csv(qc_summary, file.path(output_dir, "QC_cell_counts.csv"), row.names = FALSE)

# 过滤后仍使用相同阈值绘图，可直观看到被排除的低质量或异常细胞。
set.seed(random_seed)
after_plots <- make_qc_plots(filtered_samples)
pbmc_save_plot_grid(
  after_plots,
  file.path(output_dir, "QC_all_samples_after.pdf"),
  ncol = 5,
  nrow = 4,
  width = 32,
  height = 18
)

message("v1 QC finished. Results: ", output_dir)
