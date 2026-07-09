# HNSCC / OSCC SASH1 多组学复现流程

本仓库用于复现一套以 SASH1 为核心的头颈鳞状细胞癌 / 口腔鳞状细胞癌（HNSCC / OSCC）多组学分析流程。

项目包含两条主线：

1. GSE215403 单细胞 RNA-seq 分析
2. GSE252265 空间转录组分析

整体目标是从公开数据出发，完成细胞类型注释、恶性上皮细胞识别、CopyKAT 辅助验证、Monocle3 全局细胞状态轨迹构建，以及 SASH1 / MYH11 / EMP1 / COL1A1 等核心基因在单细胞和空间转录组中的表达展示。

## 1. 仓库结构

scripts/：正式复现流程代码
data/：原始数据与下载文件
results/objects/：分析过程中生成的 RDS 对象
results/tables/：CSV / TXT 结果表
results/figures/：PDF 结果图
config/：项目配置文件
README.md：本说明文件

大文件如 .h5、.rds、.tar.gz、.csv.gz、.mtx.gz 等通过 Git LFS 管理。

## 2. 数据集

### 2.1 GSE215403

GSE215403 用于单细胞 RNA-seq 分析。流程包括原始 10x 矩阵读取、质控、Seurat 标准分析、细胞类型注释、肿瘤相关上皮候选群定义、CopyKAT 分析、恶性细胞审查，以及 Monocle3 全局轨迹构建。

使用的样本包括：

OSCC
scB1
scB2
scB5
scB7
scB8
scB9
scB10
scB12
scB13
scB14
scB15

### 2.2 GSE252265

GSE252265 用于空间转录组分析。流程包括 Visium 表达矩阵读取、spot 坐标匹配、spot-level QC，以及 SASH1 / COL1A1 / EMP1 / MYH11 空间表达展示。

## 3. 运行环境

主要 R 包包括：

Seurat
SeuratObject
Matrix
dplyr
tidyr
ggplot2
patchwork
CopyKAT
monocle3
igraph
pheatmap
hdf5r

建议使用：

R >= 4.4
Seurat >= 5
monocle3 >= 1.3.1

## 4. 脚本说明

### 01_GSE215403_download_and_inventory.R

功能：

* 下载或整理 GSE215403 原始 10x 数据；
* 建立样本目录；
* 生成原始输入文件清单；
* 为后续读取矩阵做准备。

运行前需要检查：

* project_dir
* raw_dir
* sample_dir

### 02_GSE215403_read10x_and_merge.R

功能：

* 读取每个样本的 10x 表达矩阵；
* 建立 Seurat 对象；
* 合并多个样本；
* 添加样本 ID 和基础 metadata；
* 输出合并后的原始 Seurat 对象。

推荐运行位置：

* 服务器优先；
* 高内存本地电脑也可运行。

### 03_QC_reproduction_candidate.R

功能：

* 计算每个细胞的基因数、UMI 数和线粒体比例；
* 根据设定阈值过滤细胞；
* 输出 QC 前后统计表；
* 输出 QC 可视化图；
* 保存 QC 后 Seurat 对象。

推荐运行位置：

* 服务器优先；
* 本地高内存环境可运行。

### 04_standard_Seurat_multi_resolution.R

功能：

* LogNormalize；
* 高变基因选择；
* ScaleData；
* PCA；
* UMAP；
* 多分辨率聚类；
* 输出不同分辨率下的 cluster 结果；
* 保存标准 Seurat 对象。

后续主要使用：

* RNA_snn_res.0.2

推荐运行位置：

* 服务器优先。

### 05_major_cell_type_annotation.R

功能：

* 根据 marker 基因对主要 cluster 进行细胞类型注释；
* 输出主细胞群 UMAP；
* 输出 marker DotPlot；
* 展示核心基因 SASH1 / MYH11 / EMP1 / COL1A1 在主要细胞群中的表达。

主要细胞群包括：

Cytotoxic_T_NKT
Macrophage
Differentiated_Tumor
Cycling_Tumor
CT_Antigen_Tumor
Treg
Tumor_Epithelial
Fibroblast_CAF
B_Cell
Plasma_Cell
Blood_Endothelial
Mast_Cell
pDC
Lymphatic_Endothelial
Salivary_Epithelial

推荐运行位置：

* 本地可运行。

### 06_tumor_related_epithelial_candidates.R

功能：

* 基于 05c 注释定义肿瘤相关上皮候选 cluster；
* 标记 salivary epithelial reference；
* 输出候选群 UMAP；
* 输出核心基因 DotPlot；
* 为 CopyKAT 和后续恶性细胞判定做准备。

初始肿瘤相关上皮候选群为：

cluster 2：Differentiated_Tumor
cluster 3：Cycling_Tumor
cluster 4：CT_Antigen_Tumor
cluster 6：Tumor_Epithelial
cluster 11：Tumor_Epithelial

推荐运行位置：

* 本地可运行。

### 07_CopyKAT_per_sample.R

配套脚本：

07_CopyKAT_per_sample_batch.sh

功能：

* 对每个样本分别运行 CopyKAT；
* 判断 diploid / aneuploid 状态；
* 输出每个样本的 CopyKAT 结果；
* 为严格恶性细胞定义提供 CNV 支持。

推荐运行位置：

* 服务器优先。

### 08_final_malignant_call.R

配套脚本：

08_final_malignant_call_batch.sh

功能：

* 整合各样本 CopyKAT 输出；
* 将 CopyKAT aneuploid 结果映射回全细胞 Seurat 对象；
* 定义 CopyKAT 支持的严格恶性细胞；
* 输出最终恶性细胞 UMAP 和统计表。

严格恶性细胞定义：

tumor-related epithelial candidate + CopyKAT aneuploid

推荐运行位置：

* 本地可运行。

### 09_strict_malignant_composition.R

配套脚本：

09_strict_malignant_composition_batch.sh

功能：

* 统计严格恶性细胞在不同 cluster 中的分布；
* 统计严格恶性细胞在不同样本中的分布；
* 输出组成表和辅助图。

推荐运行位置：

* 本地可运行。

### 10_manual_review_epithelial_core.R

功能：

* 对严格恶性候选细胞进行二次审查；
* 根据 marker 表达排除明显免疫、CAF、髓系混杂 cluster；
* 保留更可信的上皮核心群；
* 输出人工审查表；
* 输出核心上皮状态相关图和表。

推荐运行位置：

* 本地可运行。

### 11a_malignant_focus_cluster_audit.R

功能：

* 审查 cluster 2 / 3 / 4 / 6 / 11；
* 整合 CopyKAT 支持比例；
* 检查样本组成；
* 检查 marker 表达；
* 检查在 Monocle3 vertex-bin 中的分布；
* 为最终选择 Core_Malignant_Focused 提供依据。

本流程采用：

Core_Malignant_Focused = cluster 6 + cluster 11

同时保留：

Extended_Malignant_Focused = cluster 4 + cluster 6 + cluster 11

推荐运行位置：

* 本地可运行。

### 11b_core_malignant_global_trajectory.R

功能：

* 基于全细胞对象构建 Monocle3 global cellular-state graph；
* 提取 principal graph edge 和 vertex；
* 定义 Core 与 Extended 两套 malignant-focused 区域；
* 在全局 graph 上展示 SASH1 / MYH11 / EMP1 / COL1A1；
* 输出主图和敏感性分析图；
* 保存用于后续重画的轻量化坐标和表格。

主图使用：

Core_Malignant_Focused = cluster 6 + cluster 11

敏感性分析使用：

Extended_Malignant_Focused = cluster 4 + cluster 6 + cluster 11

推荐运行位置：

* 服务器优先；
* 高内存本地电脑可运行。

### 12_spatial_download_QC_gene_maps.R

功能：

* 下载 GSE252265 空间转录组数据；
* 读取 Visium H5 表达矩阵；
* 读取 spot 坐标；
* 检查 barcode 与坐标匹配；
* 建立空间 Seurat 对象；
* 计算 spot-level QC；
* 输出 SASH1 / COL1A1 / EMP1 / MYH11 空间表达图；
* 输出用于后续空间关系分析的 metadata。

运行前需要检查：

* project_dir
* raw_dir
* object_dir
* table_dir
* figure_dir

推荐运行位置：

* 本地可运行；
* 数据下载较慢时可在服务器运行。

## 5. 主要结果文件

### 5.1 单细胞主要图

results/figures/05_major_cell_populations_UMAP_final.pdf
results/figures/06_tumor_candidate_and_salivary_reference_UMAP.pdf
results/figures/08_final_malignant_call_UMAP.pdf
results/figures/11b_core_vs_extended_malignant_focus_comparison.pdf
results/figures/11b_core_malignant_focused_gene_overlays.pdf
results/figures/11b_extended_malignant_focused_gene_overlays.pdf

### 5.2 单细胞主要表格

results/tables/08_CopyKAT_status_summary.csv
results/tables/09_malignant_cell_composition_summary.csv
results/tables/10_manual_review_cluster_decision.csv
results/tables/11a_candidate_cluster_integrated_review_summary.csv
results/tables/11b_malignant_focus_definition_summary.csv
results/tables/11b_core_extended_gene_overlay_summary.csv

### 5.3 空间转录组主要结果

results/figures/12_spatial_QC_UMI_distribution.pdf
results/figures/12_spatial_QC_detected_gene_distribution.pdf
results/figures/12_SASH1_COL1A1_EMP1_MYH11_spatial_expression.pdf
results/tables/12_spatial_spot_QC_summary.csv
results/tables/12_SASH1_COL1A1_spot_detection_summary.csv

## 6. 最终分析逻辑

单细胞部分核心逻辑：

全细胞 Seurat 聚类与注释
→ 肿瘤相关上皮候选 cluster 定义
→ CopyKAT 辅助识别 aneuploid 恶性细胞
→ 审查候选 malignant-focused cluster
→ 选定 cluster 6 + 11 作为 Core_Malignant_Focused
→ 在全细胞 Monocle3 global graph 上展示核心基因表达

主图定义：

Core_Malignant_Focused = Seurat cluster 6 + 11

补充敏感性定义：

Extended_Malignant_Focused = Seurat cluster 4 + 6 + 11

## 7. 使用说明

所有脚本默认项目路径为：

~/Desktop/HNSCC_SASH1_reproduction

在其他电脑或服务器运行时，请修改脚本中的 project_dir。

克隆仓库后请先确认 Git LFS 已安装：

git lfs install
git lfs pull

查看结果可直接进入：

results/figures/
results/tables/