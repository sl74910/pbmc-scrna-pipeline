# PBMC 单细胞 RNA-seq 分析流程

本仓库整理了我在**河南省医学科学院感染与免疫研究所**开展 PBMC 单细胞 RNA 测序分析时使用的 R 脚本。公开版本仅包含分析代码和运行说明；原始数据、Seurat 对象、分析结果、日志及密钥均不包含在仓库中。

## 内容

| 阶段 | 脚本 | 内容 |
| --- | --- | --- |
| v1 | `PBMC_scRNAv1_QC_diagnostics.R` | 质量控制及高质量细胞筛选。 |
| v2 | `PBMC_scRNAv2_DoubletFinder.R` | DoubletFinder 双细胞识别。 |
| v2 | `PBMC_scRNAv2_scDblFinder.R` | scDblFinder 双细胞识别。 |
| v2 | `PBMC_scRNAv2_intersection_summary.R` | 合并两种双细胞识别结果。 |
| v2 | `PBMC_scRNAv2_summary.R` | 双细胞识别结果汇总。 |
| v3 | `PBMC_scRNAv3_harmony.R` | Harmony 批次整合、降维与聚类。 |
| v3 | `PBMC_scRNAv3_RPCA.R` | RPCA 整合路线，可替代 Harmony。 |
| v3 | `PBMC_scRNAv3_harmony_replot.R` | 读取 Harmony 检查点并重新绘图。 |
| v3.5 | `PBMC_scRNAv3.5_hb.R` | 血红蛋白相关 QC 检查。 |
| v4 | `PBMC_scRNAv4_annotation.R` | 基于 marker 基因的人工审核注释。 |
| v4 | `PBMC_scRNAv4_ScType.R` | ScType 自动注释。 |
| v4 | `PBMC_scRNAv4_DeepCellSeek.R` | DeepCellSeek 模型辅助注释。 |

`R/` 目录存放主流程共享的函数。

## 运行前准备

1. 使用 R 4.3 或更高版本，并按需要安装 `Seurat`、`SeuratObject`、`ggplot2`、`patchwork`、`dplyr`、`DoubletFinder`、`harmony`、`scDblFinder`、`SingleCellExperiment`、`BiocParallel`、`HGNChelper` 和 `openxlsx`。
2. 将本地 Seurat RDS 输入文件放入 `inputs/`，文件名需符合 `PBMC_<时间点>day_<重复>_Seurat.rds`。
3. ScType 和 DeepCellSeek 是外部依赖：请按其上游文档在本地安装。ScType 数据库默认放在本地 `down/sc-type-master/`，该目录不会提交到 GitHub。
4. DeepCellSeek 运行前在本地设置 API 密钥，绝不将密钥写入脚本或提交到 GitHub：

```bash
export OPENAI_API_KEY="你的密钥"
```

## 推荐运行顺序

在仓库根目录执行：

```bash
Rscript PBMC_scRNAv1_QC_diagnostics.R
Rscript PBMC_scRNAv2_DoubletFinder.R
Rscript PBMC_scRNAv2_scDblFinder.R
Rscript PBMC_scRNAv2_intersection_summary.R
Rscript PBMC_scRNAv2_summary.R
Rscript PBMC_scRNAv3_harmony.R
Rscript PBMC_scRNAv3_harmony_replot.R
Rscript PBMC_scRNAv3.5_hb.R
Rscript PBMC_scRNAv4_annotation.R
Rscript PBMC_scRNAv4_ScType.R
```

`PBMC_scRNAv3_RPCA.R` 是 Harmony 的替代整合方案，不需要和 Harmony 流程同时运行。人工注释脚本首次运行后需要审核并修改生成的 cluster annotation CSV，再次运行以使用最终标签。需要使用 DeepCellSeek 时，在设置好 API 密钥后单独运行：

```bash
Rscript PBMC_scRNAv4_DeepCellSeek.R
```

## 发布范围

`.gitignore` 已配置为仅允许 R 脚本、Markdown 说明文件和 `inputs/`、`outputs/` 的目录说明被 Git 跟踪。实际数据和结果始终保留在本地。

单位署名、数据公开范围、第三方软件引用和许可证请以河南省医学科学院感染与免疫研究所的审核结果为准。
