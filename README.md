# PBMC 单细胞 RNA-seq 分析流程

本仓库保存小鼠 PBMC 单细胞 RNA-seq 分析代码。当前主流程分为三个阶段：逐样本质量控制与双细胞检测、Harmony 批次整合、DeepCellSeek 细胞类型注释。输入为包含原始 RNA counts 的 Seurat RDS 文件，最终输出带细胞类型标签的合并 Seurat 对象、marker 表和可视化结果。

仓库仅保存代码和说明文档；原始数据、分析对象、图表、日志及 API 密钥保留在本地。

> 当前运行前提：v1 引用了 `R/pbmc_plots.R`，但当前仓库中缺少该文件。运行 v1 前需补齐其中的 `save_qc()`、`save_pdf()` 和 `plot_hb_umap()`。已有的 `R/pbmc_helpers.R` 不能直接替代这个文件。

## 流程概览

| 阶段 | 入口脚本 | 主要处理 | 输出对象 |
| --- | --- | --- | --- |
| v1 | [PBMC_scRNAv1_QC_diagnostics.R](PBMC_scRNAv1_QC_diagnostics.R) | QC 指标、Hb 过滤前后检查、DoubletFinder 和 scDblFinder 检测 | 按样本命名的 Seurat 对象列表，保留全部细胞及筛选标记 |
| v2 | [PBMC_scRNAv2_harmony.R](PBMC_scRNAv2_harmony.R) | 按标记筛选细胞，合并样本，标准化、PCA、Harmony、聚类和 UMAP | 包含全部保留细胞的单个 Seurat 对象 |
| v3 | [PBMC_scRNAv3_annotation.R](PBMC_scRNAv3_annotation.R) | 提取 cluster marker，调用 DeepCellSeek 注释，绘图并统计细胞数 | 写入细胞类型标签的 Seurat 对象 |

```text
各样本原始 Seurat RDS
  → v1：QC + Hb 检查 + 两种双细胞检测
  → 带筛选标记的 Seurat 列表
  → v2：筛选细胞 + 合并 + Harmony + 聚类 / UMAP
  → 合并后的 Harmony Seurat 对象
  → v3：cluster marker + DeepCellSeek 注释
  → 注释后的 Seurat 对象、细胞类型表和图形
```

原先独立的双细胞检测、Hb 检查和结果绘图已纳入上述阶段，当前运行顺序不再包含 RPCA、ScType 或人工 annotation CSV 回填步骤。

## 运行环境与依赖

当前代码使用 Seurat / SeuratObject 5.x 的 layer 接口。仓库没有锁定完整的 R 与依赖版本；v2、v3 会将实际运行环境写入各自输出目录的 `sessionInfo.txt`。

| 用途 | R 包 |
| --- | --- |
| Seurat 对象处理与绘图 | `Seurat`、`SeuratObject`、`ggplot2`、`patchwork` |
| v1 双细胞检测 | `DoubletFinder`、`scDblFinder`、`SingleCellExperiment`、`SummarizedExperiment`、`BiocParallel`、`scuttle`、`scater`、`BiocSingular` |
| v2 批次整合 | `harmony`、`future` |
| v3 自动注释 | `DeepCellSeek`，需支持当前脚本的外部接口配置与 `allowed_cell_types` 参数 |

v1 使用 `parallel::mclapply()`，当前按 Linux/macOS 的多进程方式编写：双细胞检测默认同时处理 10 个样本，Hb UMAP 最多同时处理 2 个样本。可在 v1 中修改 `n_workers` 控制并发；Windows 下需将其设为 `1L`。

## 输入与配置

所有脚本都以 `getwd()` 作为项目根目录，需在本仓库根目录运行。输入路径、QC 阈值和分析参数直接写在脚本中，未提供统一配置文件或命令行参数。

### 输入数据

每个文件对应一个独立样本，必须能由 `readRDS()` 读取为 Seurat 对象，并包含 `RNA` assay 的原始 `counts`。文件名使用：

```text
PBMC_<时间点>day_<测量编号>_Seurat.rds

例如：PBMC_0day_1_Seurat.rds
      PBMC_3day_5_Seurat.rds
```

时间点和测量编号应为数字。v1 实际读取所有匹配 `^PBMC_.*_Seurat\.rds$` 的文件；v2 则依赖上述完整命名格式解析样本信息。

v1 当前的 `input_dir` 为：

```text
/home/cylroot/proj_Immune/inputs/PBMC_count_seaurt
```

若将数据放在仓库内的 `inputs/`，需先将 v1 的 `input_dir` 改为 `file.path(project_dir, "inputs")`；脚本不会自动切换到该目录。v1 会从 counts 重建对象，重新计算 QC，原输入对象的已有降维结果和其他 metadata 不会自动沿用。

### 样本信息与物种

当前脚本按小鼠 PBMC 设置：

- v1 的线粒体基因匹配规则为 `^(mt-|MT-)`，Hb 基因为 `Hbb-bs`、`Hbb-bt`、`Hba-a1`、`Hba-a2`。
- v2 的 `day_group` 固定为 `Day0`、`Day1`、`Day3`、`Day5`，并写入 `mouse_id = "M1"`、`tissue = "PBMC"`；`measurement` 和 `measurement_id` 从文件名生成。
- v3 使用 `species = "Mouse"`、`tissuename = "PBMC"`，绘图时间点顺序同上。

更换物种、时间点或动物设计时，需要同步修改这些设置。文件名末尾的测量编号不会自动转换为不同的动物 ID。

### DeepCellSeek 配置

v3 在计算 marker 后调用外部 API。运行前在本地设置以下任意一个密钥环境变量：

```bash
export OPENAI_API_KEY="你的密钥"
# 或设置 DEEPCELLSEEK_EXTERNAL_API_KEY
```

| 配置 | 当前值 | 修改方式 |
| --- | --- | --- |
| 模型 | `gpt-6-astra` | 设置 `DEEPCELLSEEK_MODEL` 环境变量 |
| 推理强度 | `high` | 修改 v3 中的 `DEEPCELLSEEK_REASONING_EFFORT` |
| Base URL | `https://api.tryaigc.cn` | 修改 v3 中的 `DEEPCELLSEEK_EXTERNAL_BASE_URL` |
| Endpoint | `/v1/responses` | 修改 v3 中的 `DEEPCELLSEEK_EXTERNAL_ENDPOINT_PATH` |
| 允许的细胞类型 | 本地 RDS 文件 | 修改 v3 中的 `allowed_cell_types_file` |

推理强度、Base URL 和 Endpoint 由脚本内的 `Sys.setenv()` 直接赋值，运行前设置同名环境变量会被覆盖。允许的细胞类型文件当前位于：

```text
/home/cylroot/proj_Immune/R_packLearn/DeepCellSeek/demo/inputs/PeripheralBlood_celltype.rds
```

该文件和 DeepCellSeek 包不随本仓库提供，迁移环境时需准备并调整路径。密钥只从环境变量读取，不应写入代码或版本库。

## 运行顺序

补齐 v1 绘图依赖，安装所需 R 包，并确认输入路径、样本设置和 v3 API 配置后，在仓库根目录依次执行；每一步成功结束后再运行下一步：

```bash
Rscript PBMC_scRNAv1_QC_diagnostics.R
Rscript PBMC_scRNAv2_harmony.R
Rscript PBMC_scRNAv3_annotation.R
```

各阶段读取上一步固定路径下的 RDS，可从已有的有效阶段结果继续运行。重复执行会重新计算该阶段并覆盖同名输出，没有自动跳过已完成步骤的机制。

## 各阶段的处理规则

### v1：QC、Hb 检查与双细胞检测

每个样本独立计算 QC 指标，当前保留阈值为：

| 指标 | 条件 |
| --- | --- |
| 检测到的基因数 `nFeature_RNA` | `1000 ≤ nFeature_RNA ≤ 6000` |
| UMI 总数 `nCount_RNA` | `≤ 30000` |
| 线粒体比例 `percent_mt` | `≤ 10%` |
| 血红蛋白比例 `percent_hb` | `≤ 5%` |

`qc_pass_pre_hb` 记录前三项是否通过，`hb_pass` 记录 Hb 条件，`qc_pass` 为二者同时通过。Hb UMAP 使用基本 QC 通过、尚未应用 Hb 阈值的细胞计算，在同一套坐标上展示过滤前后结果。

DoubletFinder 和 scDblFinder 均在各样本的 `qc_pass = TRUE` 细胞上运行。scDblFinder 不会先剔除 DoubletFinder 判定的双细胞，因此两种方法作用于同一批 QC 通过细胞。

- **DoubletFinder**：`LogNormalize`，缩放因子 10000，最多 3000 个高变基因和 50 个 PC，检测默认使用前 20 个 PC；预聚类分辨率 0.5，`pN = 0.25`。按最高 `BCmetric` 选择 `pK`，预期双细胞率为 `min(0.20, 0.008 × QC通过细胞数 / 1000)`，再进行同型双细胞校正。该比例规则需结合实际建库信息核对。
- **scDblFinder**：自动估计双细胞率；当前 `processing` 钩子使用 `scuttle::normalizeCounts()` 和 `scater::calculatePCA(..., BSPARAM = BiocSingular::ExactParam())`，其余检测步骤交由 scDblFinder 执行。

v1 输出保留全部原始细胞的 counts，并写入 `DoubletFinder_pass`、`scDblFinder_call`、`scDblFinder_score`、`scDblFinder_pass` 等 metadata。未参与双细胞检测的细胞对应字段为 `NA`。Hb、DoubletFinder 和 scDblFinder 的参数或汇总分别保存在各样本对象的 `@misc$HB`、`@misc$DoubletFinder` 和 `@misc$scDblFinder` 中。

每个样本的基本 QC 通过细胞及完整 QC 通过细胞均需至少 100 个；DoubletFinder 同型校正后的 `nExp` 为 0 时也会停止。v1、v2 的基础随机种子为 `1234`，v1 按样本名派生检测种子。

### v2：筛选、合并与 Harmony 整合

v2 读取 v1 的 Seurat 列表，仅保留以下三个标记均为 `TRUE` 且非 `NA` 的细胞：

```r
qc_pass & DoubletFinder_pass & scDblFinder_pass
```

即保留两种方法共同判定的 singlet，任一方法判为 doublet 的细胞都会被移除。随后用筛选后的 counts 和 metadata 重建各样本对象、合并，并保存整合前的起点 RDS。

整合采用 `LogNormalize`（缩放因子 10000）、最多 3000 个高变基因、最多 50 个 PC。Harmony 按 `sample_id` 校正批次，默认使用前 30 个 PC；可用 PC 不足时自动缩减。该阶段至少需要两个样本批次。

邻居图、分辨率为 0.5 的聚类和 UMAP 均基于 Harmony 坐标。结果名称为 `harmony`、`clusters.harmony` 和 `umap.harmony`；Harmony 不修改原始 counts。

### v3：DeepCellSeek 注释与绘图

v3 以 `clusters.harmony` 为 cluster 标签、`umap.harmony` 为绘图坐标，在 RNA assay 上按需合并 layers、补充标准化数据，再调用 `FindAllMarkers()` 提取正向 marker：`only.pos = TRUE`、`min.pct = 0.25`、`logfc.threshold = 0.25`。

marker 表作为 `llm_celltype()` 的输入，设置 `topgenenumber = 30`，使用允许的细胞类型列表，并启用 `wait_indefinitely = TRUE`。脚本检查返回标签是否覆盖全部 cluster，然后将标签映射回每个细胞，写入 `DeepCellSeek_celltype` 和 `celltype`。

输出包含整体注释 UMAP、按时间点和样本拆分的 UMAP，以及按 `day_group × sample_id × celltype` 统计的细胞数表。细胞类型计数表由 `table()` 展开，包含计数为 0 的组合。

## 主要输出

以下为各阶段完整运行后写入的主要文件；所有路径均相对于仓库根目录。

```text
outputs/
├── PBMC_scRNAv1_QC_diagnostics/
│   ├── PBMC_all_samples_QC_Seurat.rds
│   ├── HB_UMAP/
│   │   ├── HB_before_UMAP.pdf
│   │   └── HB_after_UMAP.pdf
│   └── DoubletFinder/
│       ├── <sample_id>_ElbowPlot.pdf
│       └── <sample_id>_pK_scan.pdf
├── PBMC_scRNAv2_harmony/
│   ├── PBMC_all_batches_filtered_Seurat.rds
│   ├── PBMC_all_batches_harmony_Seurat.rds
│   ├── figures/
│   │   ├── PBMC_ElbowPlot.pdf
│   │   └── PBMC_harmony_UMAP.pdf
│   ├── tables/
│   │   ├── PBMC_cells_after_filtering.csv
│   │   └── sessionInfo.txt
│   └── README.md
└── PBMC_scRNAv3_annotation/
    ├── PBMC_v3_markers.rds
    ├── PBMC_v3_markers.csv
    ├── <model>/
    │   ├── PBMC_v3_celltype.rds
    │   ├── PBMC_v3_celltype.csv
    │   └── PBMC_v3_DeepCellSeek_umap.png
    ├── figures/
    │   ├── PBMC_v3_DeepCellSeek_UMAP_by_day.pdf
    │   └── PBMC_v3_DeepCellSeek_UMAP_by_batch_day.pdf
    ├── PBMC_v3_celltype_counts.csv
    ├── PBMC_v3_DeepCellSeek_annotated.rds
    ├── sessionInfo.txt
    └── README.md
```

v1 还调用 `save_qc()` 输出 QC 图和过滤前后细胞数，具体文件名取决于尚需补齐的绘图函数。v3 仅将 cluster 标签表和整体 UMAP 放入模型子目录；marker、分面图、计数表及最终 Seurat 对象使用共用路径，切换模型后重跑会覆盖这些文件。

后续分析可直接读取最终注释对象：

```r
pbmc <- readRDS(
  "outputs/PBMC_scRNAv3_annotation/PBMC_v3_DeepCellSeek_annotated.rds"
)
table(pbmc$celltype)
```

## 辅助代码与版本管理

[R/pbmc_helpers.R](R/pbmc_helpers.R) 提供通用校验、Seurat 处理和绘图函数，当前由 v3 加载。[R/pbmc_checkpoint_helpers.R](R/pbmc_checkpoint_helpers.R) 与 [R/sigCell_QC.R](R/sigCell_QC.R) 保留了检查点及 QC 辅助函数，但上述三个入口脚本目前没有加载它们。

`.gitignore` 采用白名单方式保留 R 脚本和 Markdown 说明；`inputs/`、`outputs/` 仅保留各自的目录说明。实际数据和生成结果不纳入版本管理。
