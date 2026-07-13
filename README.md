# HNSCC/OSCC SASH1多组学论文复现

本仓库用于复现HNSCC/OSCC中SASH1相关研究的单细胞RNA-seq、空间转录组和基因调控网络虚拟扰动分析。

当前01–17流程已经覆盖单细胞数据处理、人工注释、CopyKAT辅助恶性判断、恶性focus群体定义、Fig.4-like核心基因表达、Fig.5-like轨迹、GSE252265空间分析、Fig.6-like空间结构域以及scTenifold虚拟KO/OE。

需要明确的是：本仓库已完成论文的单细胞、空间转录组和虚拟扰动相关模块，但尚未完整覆盖bulk RNA-seq、机器学习特征筛选、TCGA预后模型和外部验证队列。

## 当前完成状态

- GSE215403单细胞RNA-seq原始数据下载和10x矩阵整理
- 单细胞QC、Harmony整合、PCA、UMAP和多resolution聚类
- 基于marker和cluster差异基因的人工细胞类型注释
- 上皮及肿瘤相关细胞群体复核
- CopyKAT辅助恶性细胞判断
- 核心恶性focus clusters：6和11
- 扩展恶性focus cluster：4
- 候选恶性相关肿瘤clusters：2和3
- Fig.4-like单细胞核心基因表达图
- Fig.5-like恶性细胞轨迹和核心基因表达图
- GSE252265空间转录组下载、坐标匹配和spot QC
- Fig.6-like空间结构域近似注释和核心基因空间图
- SASH1-high与COL1A1-high共定位、最近邻距离和置换检验
- SASH1、COL1A1、EMP1和MYH11虚拟KO与虚拟OE
- 8个scTenifold分析成功，0失败，0跳过

## 数据集

### GSE215403

用于单细胞RNA-seq分析。原始补充文件和整理后的10x矩阵分别位于：

```text
data/raw/scRNA_GSE215403/
data/processed/scRNA_GSE215403/
```

### GSE252265

用于空间转录组分析。表达矩阵、barcode、feature和空间坐标位于：

```text
data/raw/GSE252265/
```

## 仓库目录结构

```text
config/
data/
  raw/
  processed/
results/
  figures/
  tables/
  objects/
scripts/
scripts_utils/
split_file_manifest/
README.md
HNSCC_SASH1_reproduction.Rproj
```

- `config/`：样本metadata和项目配置。
- `data/raw/`：GEO原始补充文件和空间原始输入。
- `data/processed/`：按样本整理后的10x表达矩阵。
- `scripts/`：01–17主分析脚本和对应batch脚本。
- `scripts_utils/`：大型RDS重建和README生成工具。
- `results/figures/`：PDF和PNG结果图。
- `results/tables/`：CSV统计结果、运行参数、输出检查和sessionInfo。
- `results/objects/`：Seurat对象、轨迹对象、网络对象和大型对象分片。
- `split_file_manifest/`：大型分片文件的重建记录。

## 01–17主流程

| 脚本 | 作用 |
|---|---|
| `scripts/01_download_and_prepare_scRNA.R` | 下载并整理GSE215403单细胞RNA-seq原始数据，核对GEO补充文件，建立按样本组织的10x表达矩阵目录。 |
| `scripts/02_read_and_QC_scRNA.R` | 读取各样本10x表达矩阵，创建Seurat对象，计算nFeature、nCount和线粒体比例，完成基础质控检查并保存过滤前对象。 |
| `scripts/03_QC_reproduction_candidate.R` | 按照论文复现目标重新评估质控阈值，对过滤前后细胞数和QC指标进行比较，保存候选QC对象和统计结果。 |
| `scripts/04_standard_Seurat_PCA_UMAP_resolution_scan.R` | 执行标准Seurat流程，包括归一化、高变基因选择、ScaleData、PCA、Harmony、邻居图、UMAP以及多resolution聚类扫描。 |
| `scripts/05_manual_annotation_and_target_gene_summary.R` | 基于经典marker、cluster差异基因和样本构成进行人工细胞类型注释，并总结SASH1、COL1A1、EMP1和MYH11表达。 |
| `scripts/06_malignant_candidate_diagnostic.R` | 针对上皮和肿瘤相关cluster进行进一步诊断，结合marker、状态评分、样本构成和轨迹信息筛选恶性候选群体。 |
| `scripts/07_CopyKAT_malignant_call.R` | 使用CopyKAT推断细胞拷贝数状态，将aneuploid预测作为恶性细胞判断的辅助证据。 |
| `scripts/07_CopyKAT_malignant_call_batch.sh` | 脚本07的批处理启动脚本，用于长时间运行CopyKAT并记录日志。 |
| `scripts/08_finalize_malignant_call.R` | 整合人工注释、cluster信息、上皮状态和CopyKAT结果，生成最终恶性细胞标签。 |
| `scripts/08_finalize_malignant_call_batch.sh` | 脚本08的批处理启动脚本。 |
| `scripts/09_malignant_cell_composition_check.R` | 检查恶性细胞在不同样本、cluster和细胞类型中的数量及比例，评估后续轨迹分析的样本可用性。 |
| `scripts/09_malignant_cell_composition_check_batch.sh` | 脚本09的批处理启动脚本。 |
| `scripts/10_manual_review_epithelial_core.R` | 对上皮细胞和核心肿瘤相关cluster进行人工复核，比较不同状态评分、样本组成和核心基因表达。 |
| `scripts/11a_malignant_focus_cluster_audit.R` | 审查核心、扩展和候选恶性focus clusters，综合CopyKAT支持率、marker、样本组成和轨迹位置形成最终cluster层面判断。 |
| `scripts/11b_core_extended_malignant_overlay_rebuild.R` | 重建核心恶性focus和扩展恶性focus覆盖标签，整理Monocle3轨迹坐标、主图边和顶点信息，供脚本16使用。 |
| `scripts/12_spatial_download_QC_gene_maps.R` | 下载并读取GSE252265空间转录组表达矩阵和坐标文件，完成spot质控、barcode匹配检查及核心基因空间表达图。 |
| `scripts/13_spatial_SASH1_COL1A1_neighborhood_analysis.R` | 定义SASH1-high和COL1A1-high spots，计算重叠比例、最近邻距离及空间置换检验。 |
| `scripts/14_spatial_domain_annotation_and_core_gene_maps.R` | 根据空间marker module score近似注释Fibrotic Stroma、Inflammatory Zone、不同肿瘤状态和CSC-like Niche，并生成Fig.6-like图。 |
| `scripts/15_scRNA_core_gene_expression_Figure4_like.R` | 在主要细胞类型层面生成SASH1、COL1A1、EMP1和MYH11的UMAP、DotPlot、FeaturePlot和VlnPlot，形成Fig.4-like组合图。 |
| `scripts/15b_scRNA_core_gene_expression_Figure4_cluster_like.R` | 在cluster与cell type组合层面生成更细粒度的Fig.4-like核心基因表达图。 |
| `scripts/16_scRNA_malignant_trajectory_Figure5_like.R` | 将恶性focus标签、轨迹坐标和核心基因表达结合，生成Fig.5-like轨迹、gene-high状态和拟时序或拟时序代理趋势图。 |
| `scripts/17_scTenifold_virtual_KO_OE.R` | 选取cluster4、6和11中的7353个恶性focus细胞，使用scTenifoldKnk和scTenifoldNet对SASH1、COL1A1、EMP1和MYH11执行虚拟KO与虚拟OE分析。 |
| `scripts/17_scTenifold_virtual_KO_OE_batch.sh` | 脚本17的批处理启动脚本，适合后台运行耗时较长的scTenifold网络分析。 |

## 推荐运行顺序

在仓库根目录依次运行：

```bash
Rscript scripts/01_download_and_prepare_scRNA.R
Rscript scripts/02_read_and_QC_scRNA.R
Rscript scripts/03_QC_reproduction_candidate.R
Rscript scripts/04_standard_Seurat_PCA_UMAP_resolution_scan.R
Rscript scripts/05_manual_annotation_and_target_gene_summary.R
Rscript scripts/06_malignant_candidate_diagnostic.R
Rscript scripts/07_CopyKAT_malignant_call.R
Rscript scripts/08_finalize_malignant_call.R
Rscript scripts/09_malignant_cell_composition_check.R
Rscript scripts/10_manual_review_epithelial_core.R
Rscript scripts/11a_malignant_focus_cluster_audit.R
Rscript scripts/11b_core_extended_malignant_overlay_rebuild.R
Rscript scripts/12_spatial_download_QC_gene_maps.R
Rscript scripts/13_spatial_SASH1_COL1A1_neighborhood_analysis.R
Rscript scripts/14_spatial_domain_annotation_and_core_gene_maps.R
Rscript scripts/15_scRNA_core_gene_expression_Figure4_like.R
Rscript scripts/15b_scRNA_core_gene_expression_Figure4_cluster_like.R
Rscript scripts/16_scRNA_malignant_trajectory_Figure5_like.R
Rscript scripts/17_scTenifold_virtual_KO_OE.R
```

耗时较长的步骤可使用：

```bash
bash scripts/07_CopyKAT_malignant_call_batch.sh
bash scripts/08_finalize_malignant_call_batch.sh
bash scripts/09_malignant_cell_composition_check_batch.sh
bash scripts/17_scTenifold_virtual_KO_OE_batch.sh
```

## 恶性细胞focus定义

```text
核心恶性focus clusters：6和11
扩展恶性focus cluster：4
候选恶性相关肿瘤clusters：2和3
```

该定义综合考虑人工注释、上皮/肿瘤marker、样本组成、轨迹位置和CopyKAT支持率。
CopyKAT仅作为辅助证据，不作为单独的恶性判定标准。

脚本17使用cluster4、6和11，共7353个细胞作为虚拟扰动输入。

## 空间转录组分析说明

当前空间分析使用GSE252265公开表达矩阵和坐标文件。
由于公开文件缺少标准Seurat空间对象所需的完整H&E图像和明确的按样本图像拆分信息，当前分析采用基于坐标的空间表达和marker-score结构域近似注释。

当前限制：

```text
空间样本ID统一记为All_spots
没有H&E病理图像叠加
空间结构域属于marker-score近似注释
不等同于原文作者基于病理图像完成的人工区域标注
```

当前空间分析显示SASH1-high和COL1A1-high spots重叠有限；在现有置换检验框架下，没有观察到明显空间共定位富集。

## scTenifold虚拟KO/OE分析

脚本17使用：

```text
scTenifoldKnk 1.1
scTenifoldNet 1.4
RSpectra 0.16.2
```

目标基因：

```text
SASH1
COL1A1
EMP1
MYH11
```

每个基因分别进行虚拟KO和虚拟OE，输出全部差异调控结果、显著结果、火山图、Top基因条形图及OE网络对象。

| 基因 | 扰动 | 状态 | 表达细胞数 | 表达比例 | 显著基因数 |
|---|---|---|---:|---:|---:|
| SASH1 | KO | Success | 922 | 12.54% | 0 |
| SASH1 | OE | Success | 922 | 12.54% | 0 |
| COL1A1 | KO | Success | 2799 | 38.07% | 0 |
| COL1A1 | OE | Success | 2799 | 38.07% | 0 |
| EMP1 | KO | Success | 2520 | 34.27% | 0 |
| EMP1 | OE | Success | 2520 | 34.27% | 0 |
| MYH11 | KO | Success | 386 | 5.25% | 0 |
| MYH11 | OE | Success | 386 | 5.25% | 2 |

当前结果中，仅MYH11虚拟OE检测到2个FDR显著差异调控基因。
其余分析在当前细胞选择、网络参数和FDR阈值下未达到显著性。

虚拟KO/OE属于基因调控网络计算模拟，不等同于CRISPR敲除、转染过表达、动物模型或其他湿实验验证。无FDR显著结果也不等于相应基因没有生物学作用。

## 大型RDS对象和Git LFS

仓库中的大型输入文件、Seurat对象、网络对象和RDS分片由Git LFS管理。

完整RDS对象包括：

```text
results/objects/12_spatial_tissue_spots_Seurat.rds
results/objects/13_spatial_SASH1_COL1A1_high_annotated_Seurat.rds
results/objects/14_spatial_domain_annotated_paper_style_Seurat.rds
results/objects/17_scTenifold_analysis_input_Seurat.rds
results/objects/17_scTenifold_OE_networks/
```

采用分片方式保存的对象包括：

```text
results/objects/11_global_trajectory_Seurat.rds.parts/
results/objects/15_Figure4_like_scRNA_core_gene_expression_Seurat.rds.parts/
results/objects/15b_Figure4_like_scRNA_cluster_core_gene_expression_Seurat.rds.parts/
results/objects/16_Figure5_like_malignant_trajectory_core_gene_Seurat.rds.parts/
```

重建大型RDS：

```bash
bash scripts_utils/reconstruct_split_rds_objects.sh
```

通用分片文件重建：

```bash
bash scripts_utils/reconstruct_split_files.sh
```

克隆后下载LFS对象：

```bash
git lfs install
git clone https://github.com/ssshuoshuo/HNSCC_SASH1_multiomics_reproduction.git
cd HNSCC_SASH1_multiomics_reproduction
git lfs pull
```

## 可重复性说明

各阶段尽量保存以下信息：

- 分析参数和实际输入对象
- PDF或PNG结果图
- CSV统计表
- Seurat或网络RDS对象
- `sessionInfo.txt`
- 输出文件存在性检查表

部分脚本包含本地绝对路径。换电脑或移动项目后，需要修改脚本顶部的`project_dir`。

脚本17保留了缺失依赖自动安装逻辑。Apple Silicon macOS首次运行前还需要安装：

```text
Apple Command Line Tools
GNU Fortran 14.2
```

已经安装依赖后，可以把脚本17中的：

```r
install_missing_packages <- TRUE
```

改为：

```r
install_missing_packages <- FALSE
```

以避免每次运行都重复检查依赖。

## 项目范围和限制

- 当前仓库重点复现单细胞RNA-seq、空间转录组和虚拟扰动模块。
- 尚未纳入bulk RNA-seq差异分析、机器学习筛选、TCGA预后模型和外部队列验证。
- 恶性细胞注释综合依赖marker、人工判断、cluster结构和CopyKAT。
- 空间结构域为marker-score近似注释。
- 公开空间数据缺少完整H&E图像和明确样本拆分。
- 拟时序或拟时序代理结果不能替代真实时间序列。
- 虚拟KO/OE不能替代真实湿实验。

## 01–17脚本与输出文件对应关系

以下清单根据当前Git已追踪文件自动生成。

### 阶段01：单细胞原始数据下载、文件核对和10x目录整理

对应代码：

- `scripts/01_download_and_prepare_scRNA.R`：下载并整理GSE215403单细胞RNA-seq原始数据，核对GEO补充文件，建立按样本组织的10x表达矩阵目录。

当前已追踪输出：

- `results/tables/01_10x_folder_check.csv`：由脚本01生成，保存文件、基因或分析条件检查结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/01_core_package_check.csv`：由脚本01生成，保存文件、基因或分析条件检查结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/01_raw_file_list.csv`：由脚本01生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/01_raw_sample_file_check.csv`：由脚本01生成，保存文件、基因或分析条件检查结果，文件格式为CSV，可用于后续统计或绘图。

### 阶段02：单细胞表达矩阵读取和基础QC

对应代码：

- `scripts/02_read_and_QC_scRNA.R`：读取各样本10x表达矩阵，创建Seurat对象，计算nFeature、nCount和线粒体比例，完成基础质控检查并保存过滤前对象。

当前已追踪输出：

- `results/figures/02_QC_scatter_nCount_vs_feature.pdf`：由脚本02生成的质量控制图。
- `results/figures/02_QC_scatter_nCount_vs_mt.pdf`：由脚本02生成的质量控制图。
- `results/figures/02_QC_violin_by_sample.pdf`：由脚本02生成的质量控制图。
- `results/figures/02_QC_violin_by_sample.png`：由脚本02生成的质量控制图。
- `results/objects/02_raw_before_QC_filtering.rds`：脚本02生成的Seurat或分析中间对象，供后续脚本继续读取。
- `results/tables/02_QC_threshold_suggestions.csv`：由脚本02生成，保存分析中使用的阈值，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/02_sample_QC_summary.csv`：由脚本02生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/02_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。
- `results/tables/02_target_gene_check.csv`：由脚本02生成，保存文件、基因或分析条件检查结果，文件格式为CSV，可用于后续统计或绘图。

### 阶段03：论文复现候选QC过滤

对应代码：

- `scripts/03_QC_reproduction_candidate.R`：按照论文复现目标重新评估质控阈值，对过滤前后细胞数和QC指标进行比较，保存候选QC对象和统计结果。

当前已追踪输出：

- `results/figures/03_QC_before_after_reproduction_candidate.pdf`：由脚本03生成的质量控制图。
- `results/figures/03_QC_before_after_reproduction_candidate.png`：由脚本03生成的质量控制图。
- `results/figures/03_QC_violin_after_filtering.pdf`：由脚本03生成的质量控制图。
- `results/figures/03_cell_number_before_after_QC.pdf`：由脚本03生成的质量控制图。
- `results/figures/03_cell_number_before_after_QC.png`：由脚本03生成的质量控制图。
- `results/figures/03_cell_number_old03_vs_reproduction_candidate.pdf`：由脚本03生成的分析结果图。
- `results/objects/03_QC_reproduction_candidate.rds`：脚本03生成的Seurat或分析中间对象，供后续脚本继续读取。
- `results/tables/03_QC_comparison_old03_vs_reproduction_candidate.csv`：由脚本03生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/03_QC_filter_summary_by_sample.csv`：由脚本03生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/03_QC_filter_thresholds_by_sample.csv`：由脚本03生成，保存分析中使用的阈值，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/03_QC_reproduction_candidate_parameters.csv`：由脚本03生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/03_QC_reproduction_candidate_summary.csv`：由脚本03生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/03_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。

### 阶段04：Seurat标准降维、Harmony整合和聚类

对应代码：

- `scripts/04_standard_Seurat_PCA_UMAP_resolution_scan.R`：执行标准Seurat流程，包括归一化、高变基因选择、ScaleData、PCA、Harmony、邻居图、UMAP以及多resolution聚类扫描。

当前已追踪输出：

- `results/figures/04_Harmony_UMAP_sample_and_cluster.pdf`：由脚本04生成的UMAP降维可视化。
- `results/figures/04_Harmony_UMAP_sample_and_cluster.png`：由脚本04生成的UMAP降维可视化。
- `results/figures/04_PCA_elbow_plot.pdf`：由脚本04生成的PCA相关图。
- `results/figures/04_PCA_elbow_plot.png`：由脚本04生成的PCA相关图。
- `results/figures/04_UMAP_cluster_resolution_0.2.pdf`：由脚本04生成的UMAP降维可视化。
- `results/figures/04_UMAP_cluster_resolution_0.3.pdf`：由脚本04生成的UMAP降维可视化。
- `results/figures/04_UMAP_cluster_resolution_0.5.pdf`：由脚本04生成的UMAP降维可视化。
- `results/figures/04_UMAP_sample_and_primary_cluster.pdf`：由脚本04生成的UMAP降维可视化。
- `results/figures/04_UMAP_sample_and_primary_cluster.png`：由脚本04生成的UMAP降维可视化。
- `results/objects/04_standard_Seurat_multi_resolution.rds`：脚本04生成的Seurat或分析中间对象，供后续脚本继续读取。
- `results/tables/04_cluster_by_sample_cell_numbers.csv`：由脚本04生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/04_cluster_summary_resolution_0.2.csv`：由脚本04生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/04_cluster_summary_resolution_0.3.csv`：由脚本04生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/04_cluster_summary_resolution_0.5.csv`：由脚本04生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/04_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。
- `results/tables/04_variable_features.csv`：由脚本04生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/04_variable_features_2000.csv`：由脚本04生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。

### 阶段05：人工细胞类型注释和目标基因表达总结

对应代码：

- `scripts/05_manual_annotation_and_target_gene_summary.R`：基于经典marker、cluster差异基因和样本构成进行人工细胞类型注释，并总结SASH1、COL1A1、EMP1和MYH11表达。

当前已追踪输出：

- `results/figures/05_UMAP_cluster_and_manual_celltype.pdf`：由脚本05生成的UMAP降维可视化。
- `results/figures/05_UMAP_cluster_and_manual_celltype.png`：由脚本05生成的UMAP降维可视化。
- `results/figures/05_UMAP_cluster_and_manual_celltype_polished.pdf`：由脚本05生成的UMAP降维可视化。
- `results/figures/05_UMAP_cluster_and_manual_celltype_polished.png`：由脚本05生成的UMAP降维可视化。
- `results/figures/05_UMAP_cluster_and_preliminary_annotation.pdf`：由脚本05生成的UMAP降维可视化。
- `results/figures/05_UMAP_cluster_and_preliminary_annotation.png`：由脚本05生成的UMAP降维可视化。
- `results/figures/05_canonical_marker_DotPlot.pdf`：由脚本05生成的marker或目标基因DotPlot。
- `results/figures/05_canonical_marker_DotPlot.png`：由脚本05生成的marker或目标基因DotPlot。
- `results/figures/05_diagnostic_FeaturePlot_Epithelial_markers.pdf`：由脚本05生成的基因表达FeaturePlot。
- `results/figures/05_diagnostic_FeaturePlot_Immune_markers.pdf`：由脚本05生成的基因表达FeaturePlot。
- `results/figures/05_diagnostic_FeaturePlot_Pericyte_mural_markers.pdf`：由脚本05生成的基因表达FeaturePlot。
- `results/figures/05_diagnostic_FeaturePlot_Stromal_vascular_markers.pdf`：由脚本05生成的基因表达FeaturePlot。
- `results/figures/05_diagnostic_cluster_by_sample_heatmap_resolution_0.2.pdf`：由脚本05生成的热图。
- `results/figures/05_diagnostic_major_lineage_DotPlot_resolution_0.2.pdf`：由脚本05生成的marker或目标基因DotPlot。
- `results/figures/05_diagnostic_major_lineage_DotPlot_resolution_0.2.png`：由脚本05生成的marker或目标基因DotPlot。
- `results/figures/05_diagnostic_target_genes_UMAP_resolution_0.2.pdf`：由脚本05生成的UMAP降维可视化。
- `results/figures/05_major_cell_populations_UMAP_final.pdf`：由脚本05生成的UMAP降维可视化。
- `results/figures/05_major_cell_populations_UMAP_final.png`：由脚本05生成的UMAP降维可视化。
- `results/figures/05_manual_celltype_UMAP_paper_style.pdf`：由脚本05生成的UMAP降维可视化。
- `results/figures/05_manual_celltype_UMAP_paper_style.png`：由脚本05生成的UMAP降维可视化。
- `results/figures/05_target_genes_DotPlot_by_manual_celltype.pdf`：由脚本05生成的marker或目标基因DotPlot。
- `results/figures/05_target_genes_UMAP.pdf`：由脚本05生成的UMAP降维可视化。
- `results/figures/05_target_genes_UMAP.png`：由脚本05生成的UMAP降维可视化。
- `results/figures/05_target_genes_UMAP_quantile_scaled.pdf`：由脚本05生成的UMAP降维可视化。
- `results/figures/05_target_genes_VlnPlot_by_manual_celltype.pdf`：由脚本05生成的基因表达小提琴图。
- `results/figures/05_target_genes_by_preliminary_celltype.pdf`：由脚本05生成的分析结果图。
- `results/figures/05_top5_markers_heatmap.pdf`：由脚本05生成的热图。
- `results/objects/05_diagnostic_manual_annotation_diagnostic.rds`：脚本05生成的Seurat或分析中间对象，供后续脚本继续读取。
- `results/objects/05_manual_annotated_before_malignant_call.rds`：脚本05生成的Seurat或分析中间对象，供后续脚本继续读取。
- `results/objects/05_manual_annotated_plot_ready.rds`：脚本05生成的Seurat或分析中间对象，供后续脚本继续读取。
- `results/tables/05_DotPlot_marker_genes_used.csv`：由脚本05生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/05_all_cluster_markers.csv`：由脚本05生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/05_cluster_preliminary_annotation_template.csv`：由脚本05生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/05_diagnostic_all_markers_resolution_0.2.csv`：由脚本05生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/05_diagnostic_cluster_QC_summary_resolution_0.2.csv`：由脚本05生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/05_diagnostic_cluster_by_sample_resolution_0.2.csv`：由脚本05生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/05_diagnostic_major_lineage_markers_used.csv`：由脚本05生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/05_diagnostic_manual_annotation_template_resolution_0.2.csv`：由脚本05生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/05_diagnostic_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。
- `results/tables/05_diagnostic_top30_markers_resolution_0.2.csv`：由脚本05生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/05_manual_annotation_table_resolution_0.2.csv`：由脚本05生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/05_marker_panel_gene_check.csv`：由脚本05生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/05_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。
- `results/tables/05_target_gene_expression_by_cluster.csv`：由脚本05生成，保存基因表达结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/05_target_gene_expression_by_manual_celltype.csv`：由脚本05生成，保存基因表达结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/05_top15_markers_by_cluster.csv`：由脚本05生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。

### 阶段06：上皮及恶性候选细胞诊断

对应代码：

- `scripts/06_malignant_candidate_diagnostic.R`：针对上皮和肿瘤相关cluster进行进一步诊断，结合marker、状态评分、样本构成和轨迹信息筛选恶性候选群体。

当前已追踪输出：

- `results/figures/06_epithelial_pseudotime_UMAP.pdf`：由脚本06生成的UMAP降维可视化。
- `results/figures/06_epithelial_pseudotime_UMAP.png`：由脚本06生成的UMAP降维可视化。
- `results/figures/06_epithelial_salivary_tumor_marker_DotPlot.pdf`：由脚本06生成的marker或目标基因DotPlot。
- `results/figures/06_epithelial_state_UMAP.pdf`：由脚本06生成的UMAP降维可视化。
- `results/figures/06_epithelial_state_scores_along_pseudotime.pdf`：由脚本06生成的分析结果图。
- `results/figures/06_epithelial_subset_UMAP_cluster_and_sample.pdf`：由脚本06生成的UMAP降维可视化。
- `results/figures/06_epithelial_subset_UMAP_cluster_and_sample.png`：由脚本06生成的UMAP降维可视化。
- `results/figures/06_gene_expression_trends_along_pseudotime.pdf`：由脚本06生成的基因表达图。
- `results/figures/06_target_genes_DotPlot_by_diagnostic_status.pdf`：由脚本06生成的marker或目标基因DotPlot。
- `results/figures/06_target_genes_VlnPlot_by_diagnostic_status.pdf`：由脚本06生成的基因表达小提琴图。
- `results/figures/06_target_genes_by_epithelial_cluster.pdf`：由脚本06生成的分析结果图。
- `results/figures/06_target_genes_in_epithelial_like_UMAP.pdf`：由脚本06生成的UMAP降维可视化。
- `results/figures/06_tumor_candidate_and_salivary_reference_UMAP.pdf`：由脚本06生成的UMAP降维可视化。
- `results/figures/06_tumor_candidate_and_salivary_reference_UMAP.png`：由脚本06生成的UMAP降维可视化。
- `results/figures/06c_chromosome_aware_CNV_expression_heatmap.pdf`：由脚本06生成的热图。
- `results/figures/06c_chromosome_aware_CNV_expression_scores.pdf`：由脚本06生成的基因表达图。
- `results/figures/06c_chromosome_level_deviation_heatmap.pdf`：由脚本06生成的热图。
- `results/objects/06_malignant_candidate_diagnostic.rds`：脚本06生成的Seurat或分析中间对象，供后续脚本继续读取。
- `results/tables/06_diagnostic_gene_check.csv`：由脚本06生成，保存文件、基因或分析条件检查结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/06_epithelial_cell_level_pseudotime.csv`：由脚本06生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/06_epithelial_cluster_cell_numbers.csv`：由脚本06生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/06_epithelial_cluster_state_scores.csv`：由脚本06生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/06_epithelial_state_marker_genes_used.csv`：由脚本06生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/06_malignant_candidate_by_sample.csv`：由脚本06生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/06_malignant_candidate_cluster_summary.csv`：由脚本06生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/06_pseudotime_root_cluster.csv`：由脚本06生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/06_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。
- `results/tables/06_target_gene_expression_by_diagnostic_status.csv`：由脚本06生成，保存基因表达结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/06_target_gene_pseudobulk_by_sample.csv`：由脚本06生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/06c_chromosome_aware_CNV_expression_scores.csv`：由脚本06生成，保存基因表达结果，文件格式为CSV，可用于后续统计或绘图。

### 阶段07：CopyKAT拷贝数推断和恶性支持分析

对应代码：

- `scripts/07_CopyKAT_malignant_call.R`：使用CopyKAT推断细胞拷贝数状态，将aneuploid预测作为恶性细胞判断的辅助证据。
- `scripts/07_CopyKAT_malignant_call_batch.sh`：脚本07的批处理启动脚本，用于长时间运行CopyKAT并记录日志。

当前已追踪输出：

- `results/figures/07_CopyKAT_aneuploid_tumor_candidate_UMAP.pdf`：由脚本07生成的UMAP降维可视化。
- `results/figures/07_CopyKAT_prediction_UMAP.pdf`：由脚本07生成的UMAP降维可视化。
- `results/figures/07a_cluster6_strict_malignant_UMAP_sample_and_subcluster.pdf`：由脚本07生成的UMAP降维可视化。
- `results/figures/07a_cluster6_strict_malignant_core_gene_UMAP.pdf`：由脚本07生成的UMAP降维可视化。
- `results/figures/07a_cluster6_strict_malignant_core_gene_by_sample.pdf`：由脚本07生成的分析结果图。
- `results/figures/07b_refined_epithelial_malignant_UMAP_sample_and_subcluster.pdf`：由脚本07生成的UMAP降维可视化。
- `results/figures/07b_refined_epithelial_malignant_core_gene_UMAP.pdf`：由脚本07生成的UMAP降维可视化。
- `results/figures/07b_strict_malignant_epithelial_refinement_UMAP.pdf`：由脚本07生成的UMAP降维可视化。
- `results/figures/07b_strict_malignant_lineage_module_scores_UMAP.pdf`：由脚本07生成的UMAP降维可视化。
- `results/figures/07c_core_gene_expression_by_malignant_state.pdf`：由脚本07生成的基因表达图。
- `results/figures/07c_malignant_state_cluster_and_annotation_UMAP.pdf`：由脚本07生成的UMAP降维可视化。
- `results/figures/07c_malignant_state_module_scores_UMAP.pdf`：由脚本07生成的UMAP降维可视化。
- `results/objects/07_CopyKAT_malignant_call.rds`：脚本07生成的Seurat或分析中间对象，供后续脚本继续读取。
- `results/tables/07_CopyKAT_run_summary.csv`：由脚本07生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/07_CopyKAT_tumor_candidate_prediction_by_sample_cluster.csv`：由脚本07生成，保存CopyKAT预测或支持率结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/07_CopyKAT_tumor_candidate_prediction_overall.csv`：由脚本07生成，保存CopyKAT预测或支持率结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/07a_cluster6_internal_subcluster_top20_markers.csv`：由脚本07生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/07a_cluster6_shared_state_selected_samples.csv`：由脚本07生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/07a_cluster6_top20_sample_diagnostic_markers.csv`：由脚本07生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/07b_refined_epithelial_malignant_by_sample_cluster.csv`：由脚本07生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/07b_refinement_summary.csv`：由脚本07生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/07c_core_gene_expression_by_sample_and_malignant_state.csv`：由脚本07生成，保存基因表达结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/07c_internal_malignant_state_markers_limited_features.csv`：由脚本07生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/07c_internal_malignant_state_top20_markers.csv`：由脚本07生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/07c_malignant_state_cluster_summary.csv`：由脚本07生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/07c_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。
- `results/tables/07c_state_gene_sets_available.csv`：由脚本07生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。

### 阶段08：最终恶性细胞标签整合

对应代码：

- `scripts/08_finalize_malignant_call.R`：整合人工注释、cluster信息、上皮状态和CopyKAT结果，生成最终恶性细胞标签。
- `scripts/08_finalize_malignant_call_batch.sh`：脚本08的批处理启动脚本。

当前已追踪输出：

- `results/figures/08_CopyKAT_aneuploid_fraction_by_cluster.pdf`：由脚本08生成的CopyKAT预测或支持率图。
- `results/figures/08_CopyKAT_aneuploid_fraction_by_sample_cluster.pdf`：由脚本08生成的CopyKAT预测或支持率图。
- `results/figures/08_core_gene_expression_by_final_status.pdf`：由脚本08生成的基因表达图。
- `results/figures/08_final_malignant_status_UMAP.pdf`：由脚本08生成的UMAP降维可视化。
- `results/figures/08a_epithelial_core_sample_cluster_heatmap.pdf`：由脚本08生成的热图。
- `results/figures/08a_within_sample_trajectory_candidate_UMAP_panel.pdf`：由脚本08生成的UMAP降维可视化。
- `results/objects/08_final_malignant_call.rds`：脚本08生成的Seurat或分析中间对象，供后续脚本继续读取。
- `results/tables/08_aneuploid_fraction_by_cluster.csv`：由脚本08生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/08_aneuploid_fraction_by_sample_cluster.csv`：由脚本08生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/08_core_gene_expression_by_final_status.csv`：由脚本08生成，保存基因表达结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/08a_epithelial_core_cell_metadata_for_trajectory.csv`：由脚本08生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/08a_epithelial_core_sample_cluster_composition.csv`：由脚本08生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/08a_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。
- `results/tables/08a_within_sample_trajectory_feasibility_summary.csv`：由脚本08生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/08c_output_file_check.csv`：检查该阶段预期输出文件或目录是否真实存在。
- `results/tables/08c_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。

### 阶段09：恶性细胞样本和cluster组成检查

对应代码：

- `scripts/09_malignant_cell_composition_check.R`：检查恶性细胞在不同样本、cluster和细胞类型中的数量及比例，评估后续轨迹分析的样本可用性。
- `scripts/09_malignant_cell_composition_check_batch.sh`：脚本09的批处理启动脚本。

当前已追踪输出：

- `results/figures/09_strict_malignant_cell_composition_by_cluster.pdf`：由脚本09生成的细胞组成图。
- `results/figures/09_strict_malignant_cell_composition_by_sample.pdf`：由脚本09生成的细胞组成图。
- `results/figures/09_strict_malignant_cells_by_sample_cluster.pdf`：由脚本09生成的分析结果图。
- `results/tables/09_pseudotime_sample_recommendation.csv`：由脚本09生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/09_strict_malignant_cells_per_sample_summary.csv`：由脚本09生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。

### 阶段10：上皮核心群体人工复核

对应代码：

- `scripts/10_manual_review_epithelial_core.R`：对上皮细胞和核心肿瘤相关cluster进行人工复核，比较不同状态评分、样本组成和核心基因表达。

当前已追踪输出：

- `results/figures/10_core_gene_expression_by_epithelial_core_state.pdf`：由脚本10生成的基因表达图。
- `results/figures/10_epithelial_core_UMAP_sample_cluster_program.pdf`：由脚本10生成的UMAP降维可视化。
- `results/figures/10_epithelial_core_relative_program_scores_UMAP.pdf`：由脚本10生成的UMAP降维可视化。
- `results/figures/10_manual_review_cluster_decision_UMAP.pdf`：由脚本10生成的UMAP降维可视化。
- `results/objects/10_malignant_epithelial_state_characterization.rds`：脚本10生成的Seurat或分析中间对象，供后续脚本继续读取。
- `results/objects/10_manual_review_epithelial_core.rds`：脚本10生成的Seurat或分析中间对象，供后续脚本继续读取。
- `results/tables/10_core_gene_expression_by_sample_cluster_and_program.csv`：由脚本10生成，保存基因表达结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/10_epithelial_core_by_sample_cluster_and_program.csv`：由脚本10生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/10_epithelial_core_relative_program_summary.csv`：由脚本10生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/10_manual_review_cell_summary.csv`：由脚本10生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/10_manual_review_cluster_decision.csv`：由脚本10生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/10_relative_program_gene_sets_available.csv`：由脚本10生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/10_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。

### 阶段11：恶性focus cluster审查和轨迹输入重建

对应代码：

- `scripts/11a_malignant_focus_cluster_audit.R`：审查核心、扩展和候选恶性focus clusters，综合CopyKAT支持率、marker、样本组成和轨迹位置形成最终cluster层面判断。
- `scripts/11b_core_extended_malignant_overlay_rebuild.R`：重建核心恶性focus和扩展恶性focus覆盖标签，整理Monocle3轨迹坐标、主图边和顶点信息，供脚本16使用。

当前已追踪输出：

- `results/figures/11_global_cellular_trajectory_vertex_bins_paper_style.pdf`：由脚本11生成的轨迹或拟时序相关图。
- `results/figures/11_global_cellular_trajectory_vertex_groups.pdf`：由脚本11生成的轨迹或拟时序相关图。
- `results/figures/11_global_trajectory_major_cell_type.pdf`：由脚本11生成的轨迹或拟时序相关图。
- `results/figures/11_global_trajectory_malignant_focused_cells.pdf`：由脚本11生成的轨迹或拟时序相关图。
- `results/figures/11_malignant_focused_gene_trajectory_overlays.pdf`：由脚本11生成的轨迹或拟时序相关图。
- `results/figures/11_malignant_focused_gene_trajectory_overlays_paper_style.pdf`：由脚本11生成的轨迹或拟时序相关图。
- `results/figures/11a_candidate_cluster_CopyKAT_support.pdf`：由脚本11生成的CopyKAT预测或支持率图。
- `results/figures/11a_candidate_cluster_marker_DotPlot.pdf`：由脚本11生成的marker或目标基因DotPlot。
- `results/figures/11a_candidate_cluster_sample_composition_heatmap.pdf`：由脚本11生成的热图。
- `results/figures/11a_candidate_cluster_vertex_bin_distribution.pdf`：由脚本11生成的分析结果图。
- `results/figures/11a_candidate_malignant_clusters_UMAP.pdf`：由脚本11生成的UMAP降维可视化。
- `results/figures/11b_core_malignant_focused_gene_overlays.pdf`：由脚本11生成的分析结果图。
- `results/figures/11b_core_vs_extended_malignant_focus_comparison.pdf`：由脚本11生成的分析结果图。
- `results/figures/11b_extended_malignant_focused_gene_overlays.pdf`：由脚本11生成的分析结果图。
- `results/objects/11_global_trajectory_Seurat.rds.parts/11_global_trajectory_Seurat.rds.part_aa`：大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。
- `results/objects/11_global_trajectory_Seurat.rds.parts/11_global_trajectory_Seurat.rds.part_ab`：大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。
- `results/tables/11_global_trajectory_cell_metadata.csv`：由脚本11生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11_global_trajectory_cell_metadata_with_vertex_bins.csv`：由脚本11生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11_malignant_focused_gene_overlay_summary.csv`：由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11_output_file_check.csv`：检查该阶段预期输出文件或目录是否真实存在。
- `results/tables/11_paper_faithful_cell_type_and_malignant_focus_summary.csv`：由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11_paper_style_malignant_focused_gene_overlay_summary.csv`：由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11_paper_style_vertex_bin_cell_type_summary.csv`：由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11_paper_style_vertex_bin_definition.csv`：由脚本11生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11_principal_graph_edge_coordinates.csv`：由脚本11生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11_principal_graph_vertex_coordinates.csv`：由脚本11生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。
- `results/tables/11_vertex_group_cell_type_summary.csv`：由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11a_candidate_cluster_CopyKAT_support_summary.csv`：由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11a_candidate_cluster_integrated_review_summary.csv`：由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11a_candidate_cluster_sample_composition.csv`：由脚本11生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11a_candidate_cluster_size_summary.csv`：由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11a_candidate_cluster_vertex_bin_distribution.csv`：由脚本11生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11a_candidate_malignant_cluster_cell_metadata.csv`：由脚本11生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11a_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。
- `results/tables/11b_core_extended_gene_overlay_summary.csv`：由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11b_malignant_focus_definition_summary.csv`：由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11b_monocle3_cell_umap_coordinates_and_focus_labels.csv`：由脚本11生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11b_output_file_check.csv`：检查该阶段预期输出文件或目录是否真实存在。
- `results/tables/11b_principal_graph_edge_coordinates.csv`：由脚本11生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11b_principal_graph_vertex_coordinates.csv`：由脚本11生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/11b_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。

### 阶段12：空间转录组下载、QC和核心基因空间图

对应代码：

- `scripts/12_spatial_download_QC_gene_maps.R`：下载并读取GSE252265空间转录组表达矩阵和坐标文件，完成spot质控、barcode匹配检查及核心基因空间表达图。

当前已追踪输出：

- `results/figures/12_SASH1_COL1A1_EMP1_MYH11_spatial_expression.pdf`：由脚本12生成的空间坐标或空间表达图。
- `results/figures/12_spatial_QC_UMI_distribution.pdf`：由脚本12生成的空间坐标或空间表达图。
- `results/figures/12_spatial_QC_detected_gene_distribution.pdf`：由脚本12生成的空间坐标或空间表达图。
- `results/figures/12_spatial_QC_percent_mt_distribution.pdf`：由脚本12生成的空间坐标或空间表达图。
- `results/objects/12_spatial_tissue_spots_Seurat.rds`：空间转录组分析生成的Seurat中间对象，保存表达矩阵、空间坐标和新增metadata。
- `results/tables/12_SASH1_COL1A1_spot_detection_summary.csv`：由脚本12生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/12_barcode_coordinate_match_summary.csv`：由脚本12生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/12_input_file_inventory.csv`：由脚本12生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/12_output_file_check.csv`：检查该阶段预期输出文件或目录是否真实存在。
- `results/tables/12_raw_tissue_positions_table.csv`：由脚本12生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/12_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。
- `results/tables/12_spatial_spot_QC_cell_metadata.csv`：由脚本12生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/12_spatial_spot_QC_summary.csv`：由脚本12生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/12_spatial_tissue_spot_expression_metadata.csv`：由脚本12生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。

### 阶段13：SASH1-high与COL1A1-high空间关系分析

对应代码：

- `scripts/13_spatial_SASH1_COL1A1_neighborhood_analysis.R`：定义SASH1-high和COL1A1-high spots，计算重叠比例、最近邻距离及空间置换检验。

当前已追踪输出：

- `results/figures/13_SASH1_COL1A1_high_colocalization_barplot.pdf`：由脚本13生成的空间共定位统计图。
- `results/figures/13_SASH1_COL1A1_high_neighborhood_permutation_test.pdf`：由脚本13生成的空间置换检验图。
- `results/figures/13_SASH1_high_COL1A1_high_spatial_overlay.pdf`：由脚本13生成的空间坐标或空间表达图。
- `results/figures/13_core_gene_spatial_expression_panel.pdf`：由脚本13生成的空间坐标或空间表达图。
- `results/objects/13_spatial_SASH1_COL1A1_high_annotated_Seurat.rds`：空间转录组分析生成的Seurat中间对象，保存表达矩阵、空间坐标和新增metadata。
- `results/tables/13_SASH1_COL1A1_high_colocalization_summary.csv`：由脚本13生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/13_SASH1_COL1A1_high_neighborhood_permutation_summary.csv`：由脚本13生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/13_SASH1_COL1A1_high_thresholds_by_sample.csv`：由脚本13生成，保存分析中使用的阈值，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/13_output_file_check.csv`：检查该阶段预期输出文件或目录是否真实存在。
- `results/tables/13_permutation_distribution_All_spots.csv`：由脚本13生成，保存空间置换检验结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/13_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。
- `results/tables/13_spatial_gene_expression_with_coordinates.csv`：由脚本13生成，保存基因表达结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/13_spatial_high_status_metadata.csv`：由脚本13生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。

### 阶段14：空间结构域近似注释和Fig.6-like图

对应代码：

- `scripts/14_spatial_domain_annotation_and_core_gene_maps.R`：根据空间marker module score近似注释Fibrotic Stroma、Inflammatory Zone、不同肿瘤状态和CSC-like Niche，并生成Fig.6-like图。

当前已追踪输出：

- `results/figures/14_spatial_domain_COL1A1_expression.pdf`：由脚本14生成的空间坐标或空间表达图。
- `results/figures/14_spatial_domain_EMP1_expression.pdf`：由脚本14生成的空间坐标或空间表达图。
- `results/figures/14_spatial_domain_Figure6_like_panel.pdf`：由脚本14生成的空间坐标或空间表达图。
- `results/figures/14_spatial_domain_MYH11_expression.pdf`：由脚本14生成的空间坐标或空间表达图。
- `results/figures/14_spatial_domain_SASH1_expression.pdf`：由脚本14生成的空间坐标或空间表达图。
- `results/figures/14_spatial_domain_annotation.pdf`：由脚本14生成的空间坐标或空间表达图。
- `results/figures/14_spatial_domain_tissue_layout.pdf`：由脚本14生成的空间坐标或空间表达图。
- `results/objects/14_spatial_domain_annotated_paper_style_Seurat.rds`：空间转录组分析生成的Seurat中间对象，保存表达矩阵、空间坐标和新增metadata。
- `results/tables/14_cluster_domain_marker_score_summary.csv`：由脚本14生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/14_cluster_paper_style_domain_annotation.csv`：由脚本14生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/14_domain_core_gene_expression_summary.csv`：由脚本14生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/14_domain_marker_genes_found.csv`：由脚本14生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/14_output_file_check.csv`：检查该阶段预期输出文件或目录是否真实存在。
- `results/tables/14_possible_spatial_image_files.csv`：由脚本14生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/14_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。
- `results/tables/14_spatial_domain_and_core_gene_metadata.csv`：由脚本14生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。

### 阶段15：Fig.4-like单细胞核心基因表达

对应代码：

- `scripts/15_scRNA_core_gene_expression_Figure4_like.R`：在主要细胞类型层面生成SASH1、COL1A1、EMP1和MYH11的UMAP、DotPlot、FeaturePlot和VlnPlot，形成Fig.4-like组合图。
- `scripts/15b_scRNA_core_gene_expression_Figure4_cluster_like.R`：在cluster与cell type组合层面生成更细粒度的Fig.4-like核心基因表达图。

当前已追踪输出：

- `results/figures/15_Figure4_core_gene_FeaturePlot_panel.pdf`：由脚本15生成的基因表达FeaturePlot。
- `results/figures/15_Figure4_core_gene_VlnPlot_panel.pdf`：由脚本15生成的基因表达小提琴图。
- `results/figures/15_Figure4_like_scRNA_core_gene_expression_panel.pdf`：由脚本15生成的基因表达图。
- `results/figures/15_Figure4a_scRNA_celltype_UMAP.pdf`：由脚本15生成的UMAP降维可视化。
- `results/figures/15_Figure4a_scRNA_celltype_UMAP.png`：由脚本15生成的UMAP降维可视化。
- `results/figures/15_Figure4b_core_gene_DotPlot_by_celltype.pdf`：由脚本15生成的marker或目标基因DotPlot。
- `results/figures/15b_Figure4_FeaturePlot_COL1A1.pdf`：由脚本15生成的基因表达FeaturePlot。
- `results/figures/15b_Figure4_FeaturePlot_EMP1.pdf`：由脚本15生成的基因表达FeaturePlot。
- `results/figures/15b_Figure4_FeaturePlot_MYH11.pdf`：由脚本15生成的基因表达FeaturePlot。
- `results/figures/15b_Figure4_FeaturePlot_SASH1.pdf`：由脚本15生成的基因表达FeaturePlot。
- `results/figures/15b_Figure4_VlnPlot_COL1A1_by_cluster_celltype.pdf`：由脚本15生成的基因表达小提琴图。
- `results/figures/15b_Figure4_VlnPlot_EMP1_by_cluster_celltype.pdf`：由脚本15生成的基因表达小提琴图。
- `results/figures/15b_Figure4_VlnPlot_MYH11_by_cluster_celltype.pdf`：由脚本15生成的基因表达小提琴图。
- `results/figures/15b_Figure4_VlnPlot_SASH1_by_cluster_celltype.pdf`：由脚本15生成的基因表达小提琴图。
- `results/figures/15b_Figure4_cluster_celltype_UMAP.pdf`：由脚本15生成的UMAP降维可视化。
- `results/figures/15b_Figure4_core_gene_FeaturePlot_panel.pdf`：由脚本15生成的基因表达FeaturePlot。
- `results/figures/15b_Figure4_core_gene_VlnPlot_by_cluster_celltype_panel.pdf`：由脚本15生成的基因表达小提琴图。
- `results/figures/15b_Figure4_like_scRNA_core_gene_expression_cluster_panel.pdf`：由脚本15生成的基因表达图。
- `results/figures/15b_Figure4a_scRNA_celltype_UMAP.pdf`：由脚本15生成的UMAP降维可视化。
- `results/figures/15b_Figure4a_scRNA_celltype_UMAP.png`：由脚本15生成的UMAP降维可视化。
- `results/figures/15b_Figure4b_core_gene_DotPlot_by_cluster_celltype.pdf`：由脚本15生成的marker或目标基因DotPlot。
- `results/objects/15_Figure4_like_scRNA_core_gene_expression_Seurat.rds.parts/15_Figure4_like_scRNA_core_gene_expression_Seurat.rds.part_aa`：大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。
- `results/objects/15_Figure4_like_scRNA_core_gene_expression_Seurat.rds.parts/15_Figure4_like_scRNA_core_gene_expression_Seurat.rds.part_ab`：大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。
- `results/objects/15b_Figure4_like_scRNA_cluster_core_gene_expression_Seurat.rds.parts/15b_Figure4_like_scRNA_cluster_core_gene_expression_Seurat.rds.part_aa`：大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。
- `results/objects/15b_Figure4_like_scRNA_cluster_core_gene_expression_Seurat.rds.parts/15b_Figure4_like_scRNA_cluster_core_gene_expression_Seurat.rds.part_ab`：大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。
- `results/tables/15_Figure4_cell_metadata_with_core_gene_expression.csv`：由脚本15生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/15_Figure4_celltype_core_gene_expression_summary.csv`：由脚本15生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/15_output_file_check.csv`：检查该阶段预期输出文件或目录是否真实存在。
- `results/tables/15_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。
- `results/tables/15b_Figure4_cell_metadata_with_core_gene_expression.csv`：由脚本15生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/15b_Figure4_celltype_core_gene_expression_summary.csv`：由脚本15生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/15b_Figure4_cluster_celltype_core_gene_expression_summary.csv`：由脚本15生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/15b_output_file_check.csv`：检查该阶段预期输出文件或目录是否真实存在。
- `results/tables/15b_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。

### 阶段16：Fig.5-like恶性细胞轨迹

对应代码：

- `scripts/16_scRNA_malignant_trajectory_Figure5_like.R`：将恶性focus标签、轨迹坐标和核心基因表达结合，生成Fig.5-like轨迹、gene-high状态和拟时序或拟时序代理趋势图。

当前已追踪输出：

- `results/figures/16_Figure5_core_gene_high_status_trajectory_panel.pdf`：由脚本16生成的轨迹或拟时序相关图。
- `results/figures/16_Figure5_core_gene_trajectory_expression_panel.pdf`：由脚本16生成的轨迹或拟时序相关图。
- `results/figures/16_Figure5_like_malignant_trajectory_core_gene_panel.pdf`：由脚本16生成的轨迹或拟时序相关图。
- `results/figures/16_Figure5_supplementary_core_gene_pseudotime_or_proxy_trend.pdf`：由脚本16生成的分析结果图。
- `results/figures/16_Figure5_trajectory_expression_COL1A1.pdf`：由脚本16生成的轨迹或拟时序相关图。
- `results/figures/16_Figure5_trajectory_expression_EMP1.pdf`：由脚本16生成的轨迹或拟时序相关图。
- `results/figures/16_Figure5_trajectory_expression_MYH11.pdf`：由脚本16生成的轨迹或拟时序相关图。
- `results/figures/16_Figure5_trajectory_expression_SASH1.pdf`：由脚本16生成的轨迹或拟时序相关图。
- `results/figures/16_Figure5_trajectory_high_status_COL1A1.pdf`：由脚本16生成的轨迹或拟时序相关图。
- `results/figures/16_Figure5_trajectory_high_status_EMP1.pdf`：由脚本16生成的轨迹或拟时序相关图。
- `results/figures/16_Figure5_trajectory_high_status_MYH11.pdf`：由脚本16生成的轨迹或拟时序相关图。
- `results/figures/16_Figure5_trajectory_high_status_SASH1.pdf`：由脚本16生成的轨迹或拟时序相关图。
- `results/figures/16_Figure5a_trajectory_status.pdf`：由脚本16生成的轨迹或拟时序相关图。
- `results/figures/16_Figure5b_trajectory_celltype.pdf`：由脚本16生成的轨迹或拟时序相关图。
- `results/figures/16_Figure5c_supplementary_pseudotime_or_proxy.pdf`：由脚本16生成的分析结果图。
- `results/objects/16_Figure5_like_malignant_trajectory_core_gene_Seurat.rds.parts/16_Figure5_like_malignant_trajectory_core_gene_Seurat.rds.part_aa`：大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。
- `results/objects/16_Figure5_like_malignant_trajectory_core_gene_Seurat.rds.parts/16_Figure5_like_malignant_trajectory_core_gene_Seurat.rds.part_ab`：大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。
- `results/tables/16_Figure5_core_gene_high_thresholds.csv`：由脚本16生成，保存分析中使用的阈值，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/16_Figure5_gene_high_status_distribution_summary.csv`：由脚本16生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/16_Figure5_status_source_summary.csv`：由脚本16生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/16_Figure5_trajectory_metadata_with_core_gene_expression.csv`：由脚本16生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/16_Figure5_trajectory_status_core_gene_summary.csv`：由脚本16生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/16_available_trajectory_annotation_columns.csv`：由脚本16生成，保存轨迹分析相关结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/16_output_file_check.csv`：检查该阶段预期输出文件或目录是否真实存在。
- `results/tables/16_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。

### 阶段17：scTenifold虚拟KO和虚拟OE

对应代码：

- `scripts/17_scTenifold_virtual_KO_OE.R`：选取cluster4、6和11中的7353个恶性focus细胞，使用scTenifoldKnk和scTenifoldNet对SASH1、COL1A1、EMP1和MYH11执行虚拟KO与虚拟OE分析。
- `scripts/17_scTenifold_virtual_KO_OE_batch.sh`：脚本17的批处理启动脚本，适合后台运行耗时较长的scTenifold网络分析。

当前已追踪输出：

- `results/objects/17_scTenifold_analysis_input_Seurat.rds`：脚本17实际用于虚拟KO/OE分析的7353个cluster4、6和11细胞Seurat对象。
- `results/tables/17_output_file_check.csv`：检查该阶段预期输出文件或目录是否真实存在。
- `results/tables/17_run_status_summary.csv`：汇总运行输入、细胞数量以及成功、失败和跳过的分析数量。
- `results/tables/17_runtime_config.txt`：记录脚本实际使用的输入对象、细胞筛选来源、目标基因、网络参数和软件版本。
- `results/tables/17_scTenifold_KO_OE_all_results.csv`：由脚本17生成，保存全部差异调控或统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/17_scTenifold_KO_OE_summary.csv`：由脚本17生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。
- `results/tables/17_sessionInfo.txt`：记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。

## Git已追踪文件详细清单

本节由`scripts_utils/generate_detailed_readme.py`自动生成。
只列出当前Git实际追踪的文件，因此不会把`.DS_Store`、临时日志或未追踪文件写入README。

### 根目录

| 文件 | 作用 |
|---|---|
| `.gitattributes` | 定义Git LFS追踪规则，用于管理RDS、压缩包和其他大型文件。 |
| `.gitignore` | 定义不进入Git版本控制的临时文件、本地缓存和可重建文件。 |
| `HNSCC_SASH1_reproduction.Rproj` | RStudio项目文件，用于从项目根目录打开和运行分析。 |
| `README.md` | 本项目的总体说明、运行顺序、分析结果、文件清单和限制。 |

### config/

| 文件 | 作用 |
|---|---|
| `config/GSE215403_sample_metadata.csv` | 样本metadata或项目配置文件。 |

### data/raw/

| 文件 | 作用 |
|---|---|
| `data/raw/scRNA_GSE215403/GSE215403_RAW.tar` | GEO提供的补充原始文件压缩包，由Git LFS管理。 |

### data/raw/GSE252265/

| 文件 | 作用 |
|---|---|
| `data/raw/GSE252265/GSE252265_RAW.tar` | GEO提供的补充原始文件压缩包，由Git LFS管理。 |
| `data/raw/GSE252265/GSE252265_aggr_tissue_positions.csv.gz` | 空间spot组织坐标表，用于将barcode映射到空间位置。 |
| `data/raw/GSE252265/GSE252265_aggregation.csv.gz` | GSE252265聚合信息文件。 |
| `data/raw/GSE252265/GSE252265_barcodes.tsv.gz` | 从GEO补充文件整理得到的样本barcode列表。 |
| `data/raw/GSE252265/GSE252265_features.tsv.gz` | 从GEO补充文件整理得到的样本feature注释表。 |
| `data/raw/GSE252265/GSE252265_filtered_feature_bc_matrix.h5` | 过滤后的feature-barcode表达矩阵，HDF5格式。 |
| `data/raw/GSE252265/GSE252265_matrix.mtx.gz` | 从GEO补充文件整理得到的样本稀疏表达矩阵。 |

### data/processed/

| 文件 | 作用 |
|---|---|
| `data/processed/scRNA_GSE215403/10x_by_sample/OSCC/barcodes.tsv.gz` | 10x格式barcode列表，对应表达矩阵中的细胞列。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/OSCC/features.tsv.gz` | 10x格式feature列表，记录基因ID、基因符号和feature类型。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/OSCC/matrix.mtx.gz` | 10x格式稀疏表达矩阵，行为基因，列为细胞。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB1/barcodes.tsv.gz` | 10x格式barcode列表，对应表达矩阵中的细胞列。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB1/features.tsv.gz` | 10x格式feature列表，记录基因ID、基因符号和feature类型。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB1/matrix.mtx.gz` | 10x格式稀疏表达矩阵，行为基因，列为细胞。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB10/barcodes.tsv.gz` | 10x格式barcode列表，对应表达矩阵中的细胞列。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB10/features.tsv.gz` | 10x格式feature列表，记录基因ID、基因符号和feature类型。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB10/matrix.mtx.gz` | 10x格式稀疏表达矩阵，行为基因，列为细胞。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB12/barcodes.tsv.gz` | 10x格式barcode列表，对应表达矩阵中的细胞列。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB12/features.tsv.gz` | 10x格式feature列表，记录基因ID、基因符号和feature类型。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB12/matrix.mtx.gz` | 10x格式稀疏表达矩阵，行为基因，列为细胞。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB13/barcodes.tsv.gz` | 10x格式barcode列表，对应表达矩阵中的细胞列。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB13/features.tsv.gz` | 10x格式feature列表，记录基因ID、基因符号和feature类型。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB13/matrix.mtx.gz` | 10x格式稀疏表达矩阵，行为基因，列为细胞。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB14/barcodes.tsv.gz` | 10x格式barcode列表，对应表达矩阵中的细胞列。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB14/features.tsv.gz` | 10x格式feature列表，记录基因ID、基因符号和feature类型。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB14/matrix.mtx.gz` | 10x格式稀疏表达矩阵，行为基因，列为细胞。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB15/barcodes.tsv.gz` | 10x格式barcode列表，对应表达矩阵中的细胞列。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB15/features.tsv.gz` | 10x格式feature列表，记录基因ID、基因符号和feature类型。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB15/matrix.mtx.gz` | 10x格式稀疏表达矩阵，行为基因，列为细胞。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB2/barcodes.tsv.gz` | 10x格式barcode列表，对应表达矩阵中的细胞列。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB2/features.tsv.gz` | 10x格式feature列表，记录基因ID、基因符号和feature类型。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB2/matrix.mtx.gz` | 10x格式稀疏表达矩阵，行为基因，列为细胞。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB5/barcodes.tsv.gz` | 10x格式barcode列表，对应表达矩阵中的细胞列。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB5/features.tsv.gz` | 10x格式feature列表，记录基因ID、基因符号和feature类型。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB5/matrix.mtx.gz` | 10x格式稀疏表达矩阵，行为基因，列为细胞。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB7/barcodes.tsv.gz` | 10x格式barcode列表，对应表达矩阵中的细胞列。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB7/features.tsv.gz` | 10x格式feature列表，记录基因ID、基因符号和feature类型。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB7/matrix.mtx.gz` | 10x格式稀疏表达矩阵，行为基因，列为细胞。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB8/barcodes.tsv.gz` | 10x格式barcode列表，对应表达矩阵中的细胞列。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB8/features.tsv.gz` | 10x格式feature列表，记录基因ID、基因符号和feature类型。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB8/matrix.mtx.gz` | 10x格式稀疏表达矩阵，行为基因，列为细胞。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB9/barcodes.tsv.gz` | 10x格式barcode列表，对应表达矩阵中的细胞列。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB9/features.tsv.gz` | 10x格式feature列表，记录基因ID、基因符号和feature类型。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB9/matrix.mtx.gz` | 10x格式稀疏表达矩阵，行为基因，列为细胞。 |
| `data/processed/scRNA_GSE215403/GSM6634869_OSCC_barcodes.tsv.gz` | 从GEO补充文件整理得到的样本barcode列表。 |
| `data/processed/scRNA_GSE215403/GSM6634869_OSCC_features.tsv.gz` | 从GEO补充文件整理得到的样本feature注释表。 |
| `data/processed/scRNA_GSE215403/GSM6634869_OSCC_matrix.mtx.gz` | 从GEO补充文件整理得到的样本稀疏表达矩阵。 |
| `data/processed/scRNA_GSE215403/GSM6634870_scB1_barcodes.tsv.gz` | 从GEO补充文件整理得到的样本barcode列表。 |
| `data/processed/scRNA_GSE215403/GSM6634870_scB1_features.tsv.gz` | 从GEO补充文件整理得到的样本feature注释表。 |
| `data/processed/scRNA_GSE215403/GSM6634870_scB1_matrix.mtx.gz` | 从GEO补充文件整理得到的样本稀疏表达矩阵。 |
| `data/processed/scRNA_GSE215403/GSM6634871_scB2_barcodes.tsv.gz` | 从GEO补充文件整理得到的样本barcode列表。 |
| `data/processed/scRNA_GSE215403/GSM6634871_scB2_features.tsv.gz` | 从GEO补充文件整理得到的样本feature注释表。 |
| `data/processed/scRNA_GSE215403/GSM6634871_scB2_matrix.mtx.gz` | 从GEO补充文件整理得到的样本稀疏表达矩阵。 |
| `data/processed/scRNA_GSE215403/GSM6634872_scB5_barcodes.tsv.gz` | 从GEO补充文件整理得到的样本barcode列表。 |
| `data/processed/scRNA_GSE215403/GSM6634872_scB5_features.tsv.gz` | 从GEO补充文件整理得到的样本feature注释表。 |
| `data/processed/scRNA_GSE215403/GSM6634872_scB5_matrix.mtx.gz` | 从GEO补充文件整理得到的样本稀疏表达矩阵。 |
| `data/processed/scRNA_GSE215403/GSM6634873_scB7_barcodes.tsv.gz` | 从GEO补充文件整理得到的样本barcode列表。 |
| `data/processed/scRNA_GSE215403/GSM6634873_scB7_features.tsv.gz` | 从GEO补充文件整理得到的样本feature注释表。 |
| `data/processed/scRNA_GSE215403/GSM6634873_scB7_matrix.mtx.gz` | 从GEO补充文件整理得到的样本稀疏表达矩阵。 |
| `data/processed/scRNA_GSE215403/GSM6634874_scB8_barcodes.tsv.gz` | 从GEO补充文件整理得到的样本barcode列表。 |
| `data/processed/scRNA_GSE215403/GSM6634874_scB8_features.tsv.gz` | 从GEO补充文件整理得到的样本feature注释表。 |
| `data/processed/scRNA_GSE215403/GSM6634874_scB8_matrix.mtx.gz` | 从GEO补充文件整理得到的样本稀疏表达矩阵。 |
| `data/processed/scRNA_GSE215403/GSM6634875_scB9_barcodes.tsv.gz` | 从GEO补充文件整理得到的样本barcode列表。 |
| `data/processed/scRNA_GSE215403/GSM6634875_scB9_features.tsv.gz` | 从GEO补充文件整理得到的样本feature注释表。 |
| `data/processed/scRNA_GSE215403/GSM6634875_scB9_matrix.mtx.gz` | 从GEO补充文件整理得到的样本稀疏表达矩阵。 |
| `data/processed/scRNA_GSE215403/GSM6634876_scB10_barcodes.tsv.gz` | 从GEO补充文件整理得到的样本barcode列表。 |
| `data/processed/scRNA_GSE215403/GSM6634876_scB10_features.tsv.gz` | 从GEO补充文件整理得到的样本feature注释表。 |
| `data/processed/scRNA_GSE215403/GSM6634876_scB10_matrix.mtx.gz` | 从GEO补充文件整理得到的样本稀疏表达矩阵。 |
| `data/processed/scRNA_GSE215403/GSM6634877_scB12_barcodes.tsv.gz` | 从GEO补充文件整理得到的样本barcode列表。 |
| `data/processed/scRNA_GSE215403/GSM6634877_scB12_features.tsv.gz` | 从GEO补充文件整理得到的样本feature注释表。 |
| `data/processed/scRNA_GSE215403/GSM6634877_scB12_matrix.mtx.gz` | 从GEO补充文件整理得到的样本稀疏表达矩阵。 |
| `data/processed/scRNA_GSE215403/GSM6634878_scB13_barcodes.tsv.gz` | 从GEO补充文件整理得到的样本barcode列表。 |
| `data/processed/scRNA_GSE215403/GSM6634878_scB13_features.tsv.gz` | 从GEO补充文件整理得到的样本feature注释表。 |
| `data/processed/scRNA_GSE215403/GSM6634878_scB13_matrix.mtx.gz` | 从GEO补充文件整理得到的样本稀疏表达矩阵。 |
| `data/processed/scRNA_GSE215403/GSM6634879_scB14_barcodes.tsv.gz` | 从GEO补充文件整理得到的样本barcode列表。 |
| `data/processed/scRNA_GSE215403/GSM6634879_scB14_features.tsv.gz` | 从GEO补充文件整理得到的样本feature注释表。 |
| `data/processed/scRNA_GSE215403/GSM6634879_scB14_matrix.mtx.gz` | 从GEO补充文件整理得到的样本稀疏表达矩阵。 |
| `data/processed/scRNA_GSE215403/GSM6634880_scB15_barcodes.tsv.gz` | 从GEO补充文件整理得到的样本barcode列表。 |
| `data/processed/scRNA_GSE215403/GSM6634880_scB15_features.tsv.gz` | 从GEO补充文件整理得到的样本feature注释表。 |
| `data/processed/scRNA_GSE215403/GSM6634880_scB15_matrix.mtx.gz` | 从GEO补充文件整理得到的样本稀疏表达矩阵。 |

### scripts/

| 文件 | 作用 |
|---|---|
| `scripts/01_download_and_prepare_scRNA.R` | 下载并整理GSE215403单细胞RNA-seq原始数据，核对GEO补充文件，建立按样本组织的10x表达矩阵目录。 |
| `scripts/02_read_and_QC_scRNA.R` | 读取各样本10x表达矩阵，创建Seurat对象，计算nFeature、nCount和线粒体比例，完成基础质控检查并保存过滤前对象。 |
| `scripts/03_QC_reproduction_candidate.R` | 按照论文复现目标重新评估质控阈值，对过滤前后细胞数和QC指标进行比较，保存候选QC对象和统计结果。 |
| `scripts/04_standard_Seurat_PCA_UMAP_resolution_scan.R` | 执行标准Seurat流程，包括归一化、高变基因选择、ScaleData、PCA、Harmony、邻居图、UMAP以及多resolution聚类扫描。 |
| `scripts/05_manual_annotation_and_target_gene_summary.R` | 基于经典marker、cluster差异基因和样本构成进行人工细胞类型注释，并总结SASH1、COL1A1、EMP1和MYH11表达。 |
| `scripts/06_malignant_candidate_diagnostic.R` | 针对上皮和肿瘤相关cluster进行进一步诊断，结合marker、状态评分、样本构成和轨迹信息筛选恶性候选群体。 |
| `scripts/07_CopyKAT_malignant_call.R` | 使用CopyKAT推断细胞拷贝数状态，将aneuploid预测作为恶性细胞判断的辅助证据。 |
| `scripts/07_CopyKAT_malignant_call_batch.sh` | 脚本07的批处理启动脚本，用于长时间运行CopyKAT并记录日志。 |
| `scripts/08_finalize_malignant_call.R` | 整合人工注释、cluster信息、上皮状态和CopyKAT结果，生成最终恶性细胞标签。 |
| `scripts/08_finalize_malignant_call_batch.sh` | 脚本08的批处理启动脚本。 |
| `scripts/09_malignant_cell_composition_check.R` | 检查恶性细胞在不同样本、cluster和细胞类型中的数量及比例，评估后续轨迹分析的样本可用性。 |
| `scripts/09_malignant_cell_composition_check_batch.sh` | 脚本09的批处理启动脚本。 |
| `scripts/10_manual_review_epithelial_core.R` | 对上皮细胞和核心肿瘤相关cluster进行人工复核，比较不同状态评分、样本组成和核心基因表达。 |
| `scripts/11a_malignant_focus_cluster_audit.R` | 审查核心、扩展和候选恶性focus clusters，综合CopyKAT支持率、marker、样本组成和轨迹位置形成最终cluster层面判断。 |
| `scripts/11b_core_extended_malignant_overlay_rebuild.R` | 重建核心恶性focus和扩展恶性focus覆盖标签，整理Monocle3轨迹坐标、主图边和顶点信息，供脚本16使用。 |
| `scripts/12_spatial_download_QC_gene_maps.R` | 下载并读取GSE252265空间转录组表达矩阵和坐标文件，完成spot质控、barcode匹配检查及核心基因空间表达图。 |
| `scripts/13_spatial_SASH1_COL1A1_neighborhood_analysis.R` | 定义SASH1-high和COL1A1-high spots，计算重叠比例、最近邻距离及空间置换检验。 |
| `scripts/14_spatial_domain_annotation_and_core_gene_maps.R` | 根据空间marker module score近似注释Fibrotic Stroma、Inflammatory Zone、不同肿瘤状态和CSC-like Niche，并生成Fig.6-like图。 |
| `scripts/15_scRNA_core_gene_expression_Figure4_like.R` | 在主要细胞类型层面生成SASH1、COL1A1、EMP1和MYH11的UMAP、DotPlot、FeaturePlot和VlnPlot，形成Fig.4-like组合图。 |
| `scripts/15b_scRNA_core_gene_expression_Figure4_cluster_like.R` | 在cluster与cell type组合层面生成更细粒度的Fig.4-like核心基因表达图。 |
| `scripts/16_scRNA_malignant_trajectory_Figure5_like.R` | 将恶性focus标签、轨迹坐标和核心基因表达结合，生成Fig.5-like轨迹、gene-high状态和拟时序或拟时序代理趋势图。 |
| `scripts/17_scTenifold_virtual_KO_OE.R` | 选取cluster4、6和11中的7353个恶性focus细胞，使用scTenifoldKnk和scTenifoldNet对SASH1、COL1A1、EMP1和MYH11执行虚拟KO与虚拟OE分析。 |
| `scripts/17_scTenifold_virtual_KO_OE_batch.sh` | 脚本17的批处理启动脚本，适合后台运行耗时较长的scTenifold网络分析。 |

### scripts_utils/

| 文件 | 作用 |
|---|---|
| `scripts_utils/reconstruct_split_files.sh` | 根据manifest重建通用大型分片文件。 |
| `scripts_utils/reconstruct_split_rds_objects.sh` | 重建脚本11、15、15b和16产生的大型分片RDS对象。 |

### split_file_manifest/

| 文件 | 作用 |
|---|---|
| `split_file_manifest/oversized_files.txt` | 记录大型文件分片路径、大小或重建信息。 |

### results/figures/

| 文件 | 作用 |
|---|---|
| `results/figures/02_QC_scatter_nCount_vs_feature.pdf` | 由脚本02生成的质量控制图。 |
| `results/figures/02_QC_scatter_nCount_vs_mt.pdf` | 由脚本02生成的质量控制图。 |
| `results/figures/02_QC_violin_by_sample.pdf` | 由脚本02生成的质量控制图。 |
| `results/figures/02_QC_violin_by_sample.png` | 由脚本02生成的质量控制图。 |
| `results/figures/03_QC_before_after_reproduction_candidate.pdf` | 由脚本03生成的质量控制图。 |
| `results/figures/03_QC_before_after_reproduction_candidate.png` | 由脚本03生成的质量控制图。 |
| `results/figures/03_QC_violin_after_filtering.pdf` | 由脚本03生成的质量控制图。 |
| `results/figures/03_cell_number_before_after_QC.pdf` | 由脚本03生成的质量控制图。 |
| `results/figures/03_cell_number_before_after_QC.png` | 由脚本03生成的质量控制图。 |
| `results/figures/03_cell_number_old03_vs_reproduction_candidate.pdf` | 由脚本03生成的分析结果图。 |
| `results/figures/04_Harmony_UMAP_sample_and_cluster.pdf` | 由脚本04生成的UMAP降维可视化。 |
| `results/figures/04_Harmony_UMAP_sample_and_cluster.png` | 由脚本04生成的UMAP降维可视化。 |
| `results/figures/04_PCA_elbow_plot.pdf` | 由脚本04生成的PCA相关图。 |
| `results/figures/04_PCA_elbow_plot.png` | 由脚本04生成的PCA相关图。 |
| `results/figures/04_UMAP_cluster_resolution_0.2.pdf` | 由脚本04生成的UMAP降维可视化。 |
| `results/figures/04_UMAP_cluster_resolution_0.3.pdf` | 由脚本04生成的UMAP降维可视化。 |
| `results/figures/04_UMAP_cluster_resolution_0.5.pdf` | 由脚本04生成的UMAP降维可视化。 |
| `results/figures/04_UMAP_sample_and_primary_cluster.pdf` | 由脚本04生成的UMAP降维可视化。 |
| `results/figures/04_UMAP_sample_and_primary_cluster.png` | 由脚本04生成的UMAP降维可视化。 |
| `results/figures/05_UMAP_cluster_and_manual_celltype.pdf` | 由脚本05生成的UMAP降维可视化。 |
| `results/figures/05_UMAP_cluster_and_manual_celltype.png` | 由脚本05生成的UMAP降维可视化。 |
| `results/figures/05_UMAP_cluster_and_manual_celltype_polished.pdf` | 由脚本05生成的UMAP降维可视化。 |
| `results/figures/05_UMAP_cluster_and_manual_celltype_polished.png` | 由脚本05生成的UMAP降维可视化。 |
| `results/figures/05_UMAP_cluster_and_preliminary_annotation.pdf` | 由脚本05生成的UMAP降维可视化。 |
| `results/figures/05_UMAP_cluster_and_preliminary_annotation.png` | 由脚本05生成的UMAP降维可视化。 |
| `results/figures/05_canonical_marker_DotPlot.pdf` | 由脚本05生成的marker或目标基因DotPlot。 |
| `results/figures/05_canonical_marker_DotPlot.png` | 由脚本05生成的marker或目标基因DotPlot。 |
| `results/figures/05_diagnostic_FeaturePlot_Epithelial_markers.pdf` | 由脚本05生成的基因表达FeaturePlot。 |
| `results/figures/05_diagnostic_FeaturePlot_Immune_markers.pdf` | 由脚本05生成的基因表达FeaturePlot。 |
| `results/figures/05_diagnostic_FeaturePlot_Pericyte_mural_markers.pdf` | 由脚本05生成的基因表达FeaturePlot。 |
| `results/figures/05_diagnostic_FeaturePlot_Stromal_vascular_markers.pdf` | 由脚本05生成的基因表达FeaturePlot。 |
| `results/figures/05_diagnostic_cluster_by_sample_heatmap_resolution_0.2.pdf` | 由脚本05生成的热图。 |
| `results/figures/05_diagnostic_major_lineage_DotPlot_resolution_0.2.pdf` | 由脚本05生成的marker或目标基因DotPlot。 |
| `results/figures/05_diagnostic_major_lineage_DotPlot_resolution_0.2.png` | 由脚本05生成的marker或目标基因DotPlot。 |
| `results/figures/05_diagnostic_target_genes_UMAP_resolution_0.2.pdf` | 由脚本05生成的UMAP降维可视化。 |
| `results/figures/05_major_cell_populations_UMAP_final.pdf` | 由脚本05生成的UMAP降维可视化。 |
| `results/figures/05_major_cell_populations_UMAP_final.png` | 由脚本05生成的UMAP降维可视化。 |
| `results/figures/05_manual_celltype_UMAP_paper_style.pdf` | 由脚本05生成的UMAP降维可视化。 |
| `results/figures/05_manual_celltype_UMAP_paper_style.png` | 由脚本05生成的UMAP降维可视化。 |
| `results/figures/05_target_genes_DotPlot_by_manual_celltype.pdf` | 由脚本05生成的marker或目标基因DotPlot。 |
| `results/figures/05_target_genes_UMAP.pdf` | 由脚本05生成的UMAP降维可视化。 |
| `results/figures/05_target_genes_UMAP.png` | 由脚本05生成的UMAP降维可视化。 |
| `results/figures/05_target_genes_UMAP_quantile_scaled.pdf` | 由脚本05生成的UMAP降维可视化。 |
| `results/figures/05_target_genes_VlnPlot_by_manual_celltype.pdf` | 由脚本05生成的基因表达小提琴图。 |
| `results/figures/05_target_genes_by_preliminary_celltype.pdf` | 由脚本05生成的分析结果图。 |
| `results/figures/05_top5_markers_heatmap.pdf` | 由脚本05生成的热图。 |
| `results/figures/06_epithelial_pseudotime_UMAP.pdf` | 由脚本06生成的UMAP降维可视化。 |
| `results/figures/06_epithelial_pseudotime_UMAP.png` | 由脚本06生成的UMAP降维可视化。 |
| `results/figures/06_epithelial_salivary_tumor_marker_DotPlot.pdf` | 由脚本06生成的marker或目标基因DotPlot。 |
| `results/figures/06_epithelial_state_UMAP.pdf` | 由脚本06生成的UMAP降维可视化。 |
| `results/figures/06_epithelial_state_scores_along_pseudotime.pdf` | 由脚本06生成的分析结果图。 |
| `results/figures/06_epithelial_subset_UMAP_cluster_and_sample.pdf` | 由脚本06生成的UMAP降维可视化。 |
| `results/figures/06_epithelial_subset_UMAP_cluster_and_sample.png` | 由脚本06生成的UMAP降维可视化。 |
| `results/figures/06_gene_expression_trends_along_pseudotime.pdf` | 由脚本06生成的基因表达图。 |
| `results/figures/06_target_genes_DotPlot_by_diagnostic_status.pdf` | 由脚本06生成的marker或目标基因DotPlot。 |
| `results/figures/06_target_genes_VlnPlot_by_diagnostic_status.pdf` | 由脚本06生成的基因表达小提琴图。 |
| `results/figures/06_target_genes_by_epithelial_cluster.pdf` | 由脚本06生成的分析结果图。 |
| `results/figures/06_target_genes_in_epithelial_like_UMAP.pdf` | 由脚本06生成的UMAP降维可视化。 |
| `results/figures/06_tumor_candidate_and_salivary_reference_UMAP.pdf` | 由脚本06生成的UMAP降维可视化。 |
| `results/figures/06_tumor_candidate_and_salivary_reference_UMAP.png` | 由脚本06生成的UMAP降维可视化。 |
| `results/figures/06c_chromosome_aware_CNV_expression_heatmap.pdf` | 由脚本06生成的热图。 |
| `results/figures/06c_chromosome_aware_CNV_expression_scores.pdf` | 由脚本06生成的基因表达图。 |
| `results/figures/06c_chromosome_level_deviation_heatmap.pdf` | 由脚本06生成的热图。 |
| `results/figures/07_CopyKAT_aneuploid_tumor_candidate_UMAP.pdf` | 由脚本07生成的UMAP降维可视化。 |
| `results/figures/07_CopyKAT_prediction_UMAP.pdf` | 由脚本07生成的UMAP降维可视化。 |
| `results/figures/07a_cluster6_strict_malignant_UMAP_sample_and_subcluster.pdf` | 由脚本07生成的UMAP降维可视化。 |
| `results/figures/07a_cluster6_strict_malignant_core_gene_UMAP.pdf` | 由脚本07生成的UMAP降维可视化。 |
| `results/figures/07a_cluster6_strict_malignant_core_gene_by_sample.pdf` | 由脚本07生成的分析结果图。 |
| `results/figures/07b_refined_epithelial_malignant_UMAP_sample_and_subcluster.pdf` | 由脚本07生成的UMAP降维可视化。 |
| `results/figures/07b_refined_epithelial_malignant_core_gene_UMAP.pdf` | 由脚本07生成的UMAP降维可视化。 |
| `results/figures/07b_strict_malignant_epithelial_refinement_UMAP.pdf` | 由脚本07生成的UMAP降维可视化。 |
| `results/figures/07b_strict_malignant_lineage_module_scores_UMAP.pdf` | 由脚本07生成的UMAP降维可视化。 |
| `results/figures/07c_core_gene_expression_by_malignant_state.pdf` | 由脚本07生成的基因表达图。 |
| `results/figures/07c_malignant_state_cluster_and_annotation_UMAP.pdf` | 由脚本07生成的UMAP降维可视化。 |
| `results/figures/07c_malignant_state_module_scores_UMAP.pdf` | 由脚本07生成的UMAP降维可视化。 |
| `results/figures/08_CopyKAT_aneuploid_fraction_by_cluster.pdf` | 由脚本08生成的CopyKAT预测或支持率图。 |
| `results/figures/08_CopyKAT_aneuploid_fraction_by_sample_cluster.pdf` | 由脚本08生成的CopyKAT预测或支持率图。 |
| `results/figures/08_core_gene_expression_by_final_status.pdf` | 由脚本08生成的基因表达图。 |
| `results/figures/08_final_malignant_status_UMAP.pdf` | 由脚本08生成的UMAP降维可视化。 |
| `results/figures/08a_epithelial_core_sample_cluster_heatmap.pdf` | 由脚本08生成的热图。 |
| `results/figures/08a_within_sample_trajectory_candidate_UMAP_panel.pdf` | 由脚本08生成的UMAP降维可视化。 |
| `results/figures/09_strict_malignant_cell_composition_by_cluster.pdf` | 由脚本09生成的细胞组成图。 |
| `results/figures/09_strict_malignant_cell_composition_by_sample.pdf` | 由脚本09生成的细胞组成图。 |
| `results/figures/09_strict_malignant_cells_by_sample_cluster.pdf` | 由脚本09生成的分析结果图。 |
| `results/figures/10_core_gene_expression_by_epithelial_core_state.pdf` | 由脚本10生成的基因表达图。 |
| `results/figures/10_epithelial_core_UMAP_sample_cluster_program.pdf` | 由脚本10生成的UMAP降维可视化。 |
| `results/figures/10_epithelial_core_relative_program_scores_UMAP.pdf` | 由脚本10生成的UMAP降维可视化。 |
| `results/figures/10_manual_review_cluster_decision_UMAP.pdf` | 由脚本10生成的UMAP降维可视化。 |
| `results/figures/11_global_cellular_trajectory_vertex_bins_paper_style.pdf` | 由脚本11生成的轨迹或拟时序相关图。 |
| `results/figures/11_global_cellular_trajectory_vertex_groups.pdf` | 由脚本11生成的轨迹或拟时序相关图。 |
| `results/figures/11_global_trajectory_major_cell_type.pdf` | 由脚本11生成的轨迹或拟时序相关图。 |
| `results/figures/11_global_trajectory_malignant_focused_cells.pdf` | 由脚本11生成的轨迹或拟时序相关图。 |
| `results/figures/11_malignant_focused_gene_trajectory_overlays.pdf` | 由脚本11生成的轨迹或拟时序相关图。 |
| `results/figures/11_malignant_focused_gene_trajectory_overlays_paper_style.pdf` | 由脚本11生成的轨迹或拟时序相关图。 |
| `results/figures/11a_candidate_cluster_CopyKAT_support.pdf` | 由脚本11生成的CopyKAT预测或支持率图。 |
| `results/figures/11a_candidate_cluster_marker_DotPlot.pdf` | 由脚本11生成的marker或目标基因DotPlot。 |
| `results/figures/11a_candidate_cluster_sample_composition_heatmap.pdf` | 由脚本11生成的热图。 |
| `results/figures/11a_candidate_cluster_vertex_bin_distribution.pdf` | 由脚本11生成的分析结果图。 |
| `results/figures/11a_candidate_malignant_clusters_UMAP.pdf` | 由脚本11生成的UMAP降维可视化。 |
| `results/figures/11b_core_malignant_focused_gene_overlays.pdf` | 由脚本11生成的分析结果图。 |
| `results/figures/11b_core_vs_extended_malignant_focus_comparison.pdf` | 由脚本11生成的分析结果图。 |
| `results/figures/11b_extended_malignant_focused_gene_overlays.pdf` | 由脚本11生成的分析结果图。 |
| `results/figures/12_SASH1_COL1A1_EMP1_MYH11_spatial_expression.pdf` | 由脚本12生成的空间坐标或空间表达图。 |
| `results/figures/12_spatial_QC_UMI_distribution.pdf` | 由脚本12生成的空间坐标或空间表达图。 |
| `results/figures/12_spatial_QC_detected_gene_distribution.pdf` | 由脚本12生成的空间坐标或空间表达图。 |
| `results/figures/12_spatial_QC_percent_mt_distribution.pdf` | 由脚本12生成的空间坐标或空间表达图。 |
| `results/figures/13_SASH1_COL1A1_high_colocalization_barplot.pdf` | 由脚本13生成的空间共定位统计图。 |
| `results/figures/13_SASH1_COL1A1_high_neighborhood_permutation_test.pdf` | 由脚本13生成的空间置换检验图。 |
| `results/figures/13_SASH1_high_COL1A1_high_spatial_overlay.pdf` | 由脚本13生成的空间坐标或空间表达图。 |
| `results/figures/13_core_gene_spatial_expression_panel.pdf` | 由脚本13生成的空间坐标或空间表达图。 |
| `results/figures/14_spatial_domain_COL1A1_expression.pdf` | 由脚本14生成的空间坐标或空间表达图。 |
| `results/figures/14_spatial_domain_EMP1_expression.pdf` | 由脚本14生成的空间坐标或空间表达图。 |
| `results/figures/14_spatial_domain_Figure6_like_panel.pdf` | 由脚本14生成的空间坐标或空间表达图。 |
| `results/figures/14_spatial_domain_MYH11_expression.pdf` | 由脚本14生成的空间坐标或空间表达图。 |
| `results/figures/14_spatial_domain_SASH1_expression.pdf` | 由脚本14生成的空间坐标或空间表达图。 |
| `results/figures/14_spatial_domain_annotation.pdf` | 由脚本14生成的空间坐标或空间表达图。 |
| `results/figures/14_spatial_domain_tissue_layout.pdf` | 由脚本14生成的空间坐标或空间表达图。 |
| `results/figures/15_Figure4_core_gene_FeaturePlot_panel.pdf` | 由脚本15生成的基因表达FeaturePlot。 |
| `results/figures/15_Figure4_core_gene_VlnPlot_panel.pdf` | 由脚本15生成的基因表达小提琴图。 |
| `results/figures/15_Figure4_like_scRNA_core_gene_expression_panel.pdf` | 由脚本15生成的基因表达图。 |
| `results/figures/15_Figure4a_scRNA_celltype_UMAP.pdf` | 由脚本15生成的UMAP降维可视化。 |
| `results/figures/15_Figure4a_scRNA_celltype_UMAP.png` | 由脚本15生成的UMAP降维可视化。 |
| `results/figures/15_Figure4b_core_gene_DotPlot_by_celltype.pdf` | 由脚本15生成的marker或目标基因DotPlot。 |
| `results/figures/15b_Figure4_FeaturePlot_COL1A1.pdf` | 由脚本15生成的基因表达FeaturePlot。 |
| `results/figures/15b_Figure4_FeaturePlot_EMP1.pdf` | 由脚本15生成的基因表达FeaturePlot。 |
| `results/figures/15b_Figure4_FeaturePlot_MYH11.pdf` | 由脚本15生成的基因表达FeaturePlot。 |
| `results/figures/15b_Figure4_FeaturePlot_SASH1.pdf` | 由脚本15生成的基因表达FeaturePlot。 |
| `results/figures/15b_Figure4_VlnPlot_COL1A1_by_cluster_celltype.pdf` | 由脚本15生成的基因表达小提琴图。 |
| `results/figures/15b_Figure4_VlnPlot_EMP1_by_cluster_celltype.pdf` | 由脚本15生成的基因表达小提琴图。 |
| `results/figures/15b_Figure4_VlnPlot_MYH11_by_cluster_celltype.pdf` | 由脚本15生成的基因表达小提琴图。 |
| `results/figures/15b_Figure4_VlnPlot_SASH1_by_cluster_celltype.pdf` | 由脚本15生成的基因表达小提琴图。 |
| `results/figures/15b_Figure4_cluster_celltype_UMAP.pdf` | 由脚本15生成的UMAP降维可视化。 |
| `results/figures/15b_Figure4_core_gene_FeaturePlot_panel.pdf` | 由脚本15生成的基因表达FeaturePlot。 |
| `results/figures/15b_Figure4_core_gene_VlnPlot_by_cluster_celltype_panel.pdf` | 由脚本15生成的基因表达小提琴图。 |
| `results/figures/15b_Figure4_like_scRNA_core_gene_expression_cluster_panel.pdf` | 由脚本15生成的基因表达图。 |
| `results/figures/15b_Figure4a_scRNA_celltype_UMAP.pdf` | 由脚本15生成的UMAP降维可视化。 |
| `results/figures/15b_Figure4a_scRNA_celltype_UMAP.png` | 由脚本15生成的UMAP降维可视化。 |
| `results/figures/15b_Figure4b_core_gene_DotPlot_by_cluster_celltype.pdf` | 由脚本15生成的marker或目标基因DotPlot。 |
| `results/figures/16_Figure5_core_gene_high_status_trajectory_panel.pdf` | 由脚本16生成的轨迹或拟时序相关图。 |
| `results/figures/16_Figure5_core_gene_trajectory_expression_panel.pdf` | 由脚本16生成的轨迹或拟时序相关图。 |
| `results/figures/16_Figure5_like_malignant_trajectory_core_gene_panel.pdf` | 由脚本16生成的轨迹或拟时序相关图。 |
| `results/figures/16_Figure5_supplementary_core_gene_pseudotime_or_proxy_trend.pdf` | 由脚本16生成的分析结果图。 |
| `results/figures/16_Figure5_trajectory_expression_COL1A1.pdf` | 由脚本16生成的轨迹或拟时序相关图。 |
| `results/figures/16_Figure5_trajectory_expression_EMP1.pdf` | 由脚本16生成的轨迹或拟时序相关图。 |
| `results/figures/16_Figure5_trajectory_expression_MYH11.pdf` | 由脚本16生成的轨迹或拟时序相关图。 |
| `results/figures/16_Figure5_trajectory_expression_SASH1.pdf` | 由脚本16生成的轨迹或拟时序相关图。 |
| `results/figures/16_Figure5_trajectory_high_status_COL1A1.pdf` | 由脚本16生成的轨迹或拟时序相关图。 |
| `results/figures/16_Figure5_trajectory_high_status_EMP1.pdf` | 由脚本16生成的轨迹或拟时序相关图。 |
| `results/figures/16_Figure5_trajectory_high_status_MYH11.pdf` | 由脚本16生成的轨迹或拟时序相关图。 |
| `results/figures/16_Figure5_trajectory_high_status_SASH1.pdf` | 由脚本16生成的轨迹或拟时序相关图。 |
| `results/figures/16_Figure5a_trajectory_status.pdf` | 由脚本16生成的轨迹或拟时序相关图。 |
| `results/figures/16_Figure5b_trajectory_celltype.pdf` | 由脚本16生成的轨迹或拟时序相关图。 |
| `results/figures/16_Figure5c_supplementary_pseudotime_or_proxy.pdf` | 由脚本16生成的分析结果图。 |

### results/figures/17_scTenifold_KO/

| 文件 | 作用 |
|---|---|
| `results/figures/17_scTenifold_KO/Gene_COL1A1/COL1A1_KO_top_genes_barplot.pdf` | 项目生成的Top差异调控基因条形图。 |
| `results/figures/17_scTenifold_KO/Gene_COL1A1/COL1A1_KO_volcano.pdf` | 项目生成的差异调控火山图。 |
| `results/figures/17_scTenifold_KO/Gene_EMP1/EMP1_KO_top_genes_barplot.pdf` | 项目生成的Top差异调控基因条形图。 |
| `results/figures/17_scTenifold_KO/Gene_EMP1/EMP1_KO_volcano.pdf` | 项目生成的差异调控火山图。 |
| `results/figures/17_scTenifold_KO/Gene_MYH11/MYH11_KO_top_genes_barplot.pdf` | 项目生成的Top差异调控基因条形图。 |
| `results/figures/17_scTenifold_KO/Gene_MYH11/MYH11_KO_volcano.pdf` | 项目生成的差异调控火山图。 |
| `results/figures/17_scTenifold_KO/Gene_SASH1/SASH1_KO_top_genes_barplot.pdf` | 项目生成的Top差异调控基因条形图。 |
| `results/figures/17_scTenifold_KO/Gene_SASH1/SASH1_KO_volcano.pdf` | 项目生成的差异调控火山图。 |

### results/figures/17_scTenifold_OE/

| 文件 | 作用 |
|---|---|
| `results/figures/17_scTenifold_OE/Gene_COL1A1/COL1A1_OE_top_genes_barplot.pdf` | 项目生成的Top差异调控基因条形图。 |
| `results/figures/17_scTenifold_OE/Gene_COL1A1/COL1A1_OE_volcano.pdf` | 项目生成的差异调控火山图。 |
| `results/figures/17_scTenifold_OE/Gene_EMP1/EMP1_OE_top_genes_barplot.pdf` | 项目生成的Top差异调控基因条形图。 |
| `results/figures/17_scTenifold_OE/Gene_EMP1/EMP1_OE_volcano.pdf` | 项目生成的差异调控火山图。 |
| `results/figures/17_scTenifold_OE/Gene_MYH11/MYH11_OE_top_genes_barplot.pdf` | 项目生成的Top差异调控基因条形图。 |
| `results/figures/17_scTenifold_OE/Gene_MYH11/MYH11_OE_volcano.pdf` | 项目生成的差异调控火山图。 |
| `results/figures/17_scTenifold_OE/Gene_SASH1/SASH1_OE_top_genes_barplot.pdf` | 项目生成的Top差异调控基因条形图。 |
| `results/figures/17_scTenifold_OE/Gene_SASH1/SASH1_OE_volcano.pdf` | 项目生成的差异调控火山图。 |

### results/tables/

| 文件 | 作用 |
|---|---|
| `results/tables/01_10x_folder_check.csv` | 由脚本01生成，保存文件、基因或分析条件检查结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/01_core_package_check.csv` | 由脚本01生成，保存文件、基因或分析条件检查结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/01_raw_file_list.csv` | 由脚本01生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/01_raw_sample_file_check.csv` | 由脚本01生成，保存文件、基因或分析条件检查结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/02_QC_threshold_suggestions.csv` | 由脚本02生成，保存分析中使用的阈值，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/02_sample_QC_summary.csv` | 由脚本02生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/02_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/02_target_gene_check.csv` | 由脚本02生成，保存文件、基因或分析条件检查结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/03_QC_comparison_old03_vs_reproduction_candidate.csv` | 由脚本03生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/03_QC_filter_summary_by_sample.csv` | 由脚本03生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/03_QC_filter_thresholds_by_sample.csv` | 由脚本03生成，保存分析中使用的阈值，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/03_QC_reproduction_candidate_parameters.csv` | 由脚本03生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/03_QC_reproduction_candidate_summary.csv` | 由脚本03生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/03_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/04_cluster_by_sample_cell_numbers.csv` | 由脚本04生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/04_cluster_summary_resolution_0.2.csv` | 由脚本04生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/04_cluster_summary_resolution_0.3.csv` | 由脚本04生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/04_cluster_summary_resolution_0.5.csv` | 由脚本04生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/04_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/04_variable_features.csv` | 由脚本04生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/04_variable_features_2000.csv` | 由脚本04生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/05_DotPlot_marker_genes_used.csv` | 由脚本05生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/05_all_cluster_markers.csv` | 由脚本05生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/05_cluster_preliminary_annotation_template.csv` | 由脚本05生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/05_diagnostic_all_markers_resolution_0.2.csv` | 由脚本05生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/05_diagnostic_cluster_QC_summary_resolution_0.2.csv` | 由脚本05生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/05_diagnostic_cluster_by_sample_resolution_0.2.csv` | 由脚本05生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/05_diagnostic_major_lineage_markers_used.csv` | 由脚本05生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/05_diagnostic_manual_annotation_template_resolution_0.2.csv` | 由脚本05生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/05_diagnostic_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/05_diagnostic_top30_markers_resolution_0.2.csv` | 由脚本05生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/05_manual_annotation_table_resolution_0.2.csv` | 由脚本05生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/05_marker_panel_gene_check.csv` | 由脚本05生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/05_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/05_target_gene_expression_by_cluster.csv` | 由脚本05生成，保存基因表达结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/05_target_gene_expression_by_manual_celltype.csv` | 由脚本05生成，保存基因表达结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/05_top15_markers_by_cluster.csv` | 由脚本05生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/06_diagnostic_gene_check.csv` | 由脚本06生成，保存文件、基因或分析条件检查结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/06_epithelial_cell_level_pseudotime.csv` | 由脚本06生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/06_epithelial_cluster_cell_numbers.csv` | 由脚本06生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/06_epithelial_cluster_state_scores.csv` | 由脚本06生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/06_epithelial_state_marker_genes_used.csv` | 由脚本06生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/06_malignant_candidate_by_sample.csv` | 由脚本06生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/06_malignant_candidate_cluster_summary.csv` | 由脚本06生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/06_pseudotime_root_cluster.csv` | 由脚本06生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/06_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/06_target_gene_expression_by_diagnostic_status.csv` | 由脚本06生成，保存基因表达结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/06_target_gene_pseudobulk_by_sample.csv` | 由脚本06生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/06c_chromosome_aware_CNV_expression_scores.csv` | 由脚本06生成，保存基因表达结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/07_CopyKAT_run_summary.csv` | 由脚本07生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/07_CopyKAT_tumor_candidate_prediction_by_sample_cluster.csv` | 由脚本07生成，保存CopyKAT预测或支持率结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/07_CopyKAT_tumor_candidate_prediction_overall.csv` | 由脚本07生成，保存CopyKAT预测或支持率结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/07a_cluster6_internal_subcluster_top20_markers.csv` | 由脚本07生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/07a_cluster6_shared_state_selected_samples.csv` | 由脚本07生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/07a_cluster6_top20_sample_diagnostic_markers.csv` | 由脚本07生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/07b_refined_epithelial_malignant_by_sample_cluster.csv` | 由脚本07生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/07b_refinement_summary.csv` | 由脚本07生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/07c_core_gene_expression_by_sample_and_malignant_state.csv` | 由脚本07生成，保存基因表达结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/07c_internal_malignant_state_markers_limited_features.csv` | 由脚本07生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/07c_internal_malignant_state_top20_markers.csv` | 由脚本07生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/07c_malignant_state_cluster_summary.csv` | 由脚本07生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/07c_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/07c_state_gene_sets_available.csv` | 由脚本07生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/08_aneuploid_fraction_by_cluster.csv` | 由脚本08生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/08_aneuploid_fraction_by_sample_cluster.csv` | 由脚本08生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/08_core_gene_expression_by_final_status.csv` | 由脚本08生成，保存基因表达结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/08a_epithelial_core_cell_metadata_for_trajectory.csv` | 由脚本08生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/08a_epithelial_core_sample_cluster_composition.csv` | 由脚本08生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/08a_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/08a_within_sample_trajectory_feasibility_summary.csv` | 由脚本08生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/08c_output_file_check.csv` | 检查该阶段预期输出文件或目录是否真实存在。 |
| `results/tables/08c_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/09_pseudotime_sample_recommendation.csv` | 由脚本09生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/09_strict_malignant_cells_per_sample_summary.csv` | 由脚本09生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/10_core_gene_expression_by_sample_cluster_and_program.csv` | 由脚本10生成，保存基因表达结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/10_epithelial_core_by_sample_cluster_and_program.csv` | 由脚本10生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/10_epithelial_core_relative_program_summary.csv` | 由脚本10生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/10_manual_review_cell_summary.csv` | 由脚本10生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/10_manual_review_cluster_decision.csv` | 由脚本10生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/10_relative_program_gene_sets_available.csv` | 由脚本10生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/10_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/11_global_trajectory_cell_metadata.csv` | 由脚本11生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11_global_trajectory_cell_metadata_with_vertex_bins.csv` | 由脚本11生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11_malignant_focused_gene_overlay_summary.csv` | 由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11_output_file_check.csv` | 检查该阶段预期输出文件或目录是否真实存在。 |
| `results/tables/11_paper_faithful_cell_type_and_malignant_focus_summary.csv` | 由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11_paper_style_malignant_focused_gene_overlay_summary.csv` | 由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11_paper_style_vertex_bin_cell_type_summary.csv` | 由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11_paper_style_vertex_bin_definition.csv` | 由脚本11生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11_principal_graph_edge_coordinates.csv` | 由脚本11生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11_principal_graph_vertex_coordinates.csv` | 由脚本11生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/11_vertex_group_cell_type_summary.csv` | 由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11a_candidate_cluster_CopyKAT_support_summary.csv` | 由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11a_candidate_cluster_integrated_review_summary.csv` | 由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11a_candidate_cluster_sample_composition.csv` | 由脚本11生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11a_candidate_cluster_size_summary.csv` | 由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11a_candidate_cluster_vertex_bin_distribution.csv` | 由脚本11生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11a_candidate_malignant_cluster_cell_metadata.csv` | 由脚本11生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11a_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/11b_core_extended_gene_overlay_summary.csv` | 由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11b_malignant_focus_definition_summary.csv` | 由脚本11生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11b_monocle3_cell_umap_coordinates_and_focus_labels.csv` | 由脚本11生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11b_output_file_check.csv` | 检查该阶段预期输出文件或目录是否真实存在。 |
| `results/tables/11b_principal_graph_edge_coordinates.csv` | 由脚本11生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11b_principal_graph_vertex_coordinates.csv` | 由脚本11生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/11b_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/12_SASH1_COL1A1_spot_detection_summary.csv` | 由脚本12生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/12_barcode_coordinate_match_summary.csv` | 由脚本12生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/12_input_file_inventory.csv` | 由脚本12生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/12_output_file_check.csv` | 检查该阶段预期输出文件或目录是否真实存在。 |
| `results/tables/12_raw_tissue_positions_table.csv` | 由脚本12生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/12_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/12_spatial_spot_QC_cell_metadata.csv` | 由脚本12生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/12_spatial_spot_QC_summary.csv` | 由脚本12生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/12_spatial_tissue_spot_expression_metadata.csv` | 由脚本12生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/13_SASH1_COL1A1_high_colocalization_summary.csv` | 由脚本13生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/13_SASH1_COL1A1_high_neighborhood_permutation_summary.csv` | 由脚本13生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/13_SASH1_COL1A1_high_thresholds_by_sample.csv` | 由脚本13生成，保存分析中使用的阈值，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/13_output_file_check.csv` | 检查该阶段预期输出文件或目录是否真实存在。 |
| `results/tables/13_permutation_distribution_All_spots.csv` | 由脚本13生成，保存空间置换检验结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/13_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/13_spatial_gene_expression_with_coordinates.csv` | 由脚本13生成，保存基因表达结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/13_spatial_high_status_metadata.csv` | 由脚本13生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/14_cluster_domain_marker_score_summary.csv` | 由脚本14生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/14_cluster_paper_style_domain_annotation.csv` | 由脚本14生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/14_domain_core_gene_expression_summary.csv` | 由脚本14生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/14_domain_marker_genes_found.csv` | 由脚本14生成，保存marker基因或marker表达统计，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/14_output_file_check.csv` | 检查该阶段预期输出文件或目录是否真实存在。 |
| `results/tables/14_possible_spatial_image_files.csv` | 由脚本14生成，保存该分析阶段生成的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/14_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/14_spatial_domain_and_core_gene_metadata.csv` | 由脚本14生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/15_Figure4_cell_metadata_with_core_gene_expression.csv` | 由脚本15生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/15_Figure4_celltype_core_gene_expression_summary.csv` | 由脚本15生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/15_output_file_check.csv` | 检查该阶段预期输出文件或目录是否真实存在。 |
| `results/tables/15_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/15b_Figure4_cell_metadata_with_core_gene_expression.csv` | 由脚本15生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/15b_Figure4_celltype_core_gene_expression_summary.csv` | 由脚本15生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/15b_Figure4_cluster_celltype_core_gene_expression_summary.csv` | 由脚本15生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/15b_output_file_check.csv` | 检查该阶段预期输出文件或目录是否真实存在。 |
| `results/tables/15b_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/16_Figure5_core_gene_high_thresholds.csv` | 由脚本16生成，保存分析中使用的阈值，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/16_Figure5_gene_high_status_distribution_summary.csv` | 由脚本16生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/16_Figure5_status_source_summary.csv` | 由脚本16生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/16_Figure5_trajectory_metadata_with_core_gene_expression.csv` | 由脚本16生成，保存细胞或spot层面的metadata，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/16_Figure5_trajectory_status_core_gene_summary.csv` | 由脚本16生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/16_available_trajectory_annotation_columns.csv` | 由脚本16生成，保存轨迹分析相关结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/16_output_file_check.csv` | 检查该阶段预期输出文件或目录是否真实存在。 |
| `results/tables/16_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |
| `results/tables/17_output_file_check.csv` | 检查该阶段预期输出文件或目录是否真实存在。 |
| `results/tables/17_run_status_summary.csv` | 汇总运行输入、细胞数量以及成功、失败和跳过的分析数量。 |
| `results/tables/17_runtime_config.txt` | 记录脚本实际使用的输入对象、细胞筛选来源、目标基因、网络参数和软件版本。 |
| `results/tables/17_scTenifold_KO_OE_all_results.csv` | 由脚本17生成，保存全部差异调控或统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/17_scTenifold_KO_OE_summary.csv` | 由脚本17生成，汇总该分析阶段的主要统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/17_sessionInfo.txt` | 记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。 |

### results/tables/17_scTenifold_KO/

| 文件 | 作用 |
|---|---|
| `results/tables/17_scTenifold_KO/Gene_COL1A1/COL1A1_KO_all_results.csv` | 保存全部差异调控或统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/17_scTenifold_KO/Gene_COL1A1/COL1A1_KO_significant_results.csv` | 保存达到显著性阈值的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/17_scTenifold_KO/Gene_EMP1/EMP1_KO_all_results.csv` | 保存全部差异调控或统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/17_scTenifold_KO/Gene_EMP1/EMP1_KO_significant_results.csv` | 保存达到显著性阈值的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/17_scTenifold_KO/Gene_MYH11/MYH11_KO_all_results.csv` | 保存全部差异调控或统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/17_scTenifold_KO/Gene_MYH11/MYH11_KO_significant_results.csv` | 保存达到显著性阈值的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/17_scTenifold_KO/Gene_SASH1/SASH1_KO_all_results.csv` | 保存全部差异调控或统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/17_scTenifold_KO/Gene_SASH1/SASH1_KO_significant_results.csv` | 保存达到显著性阈值的结果，文件格式为CSV，可用于后续统计或绘图。 |

### results/tables/17_scTenifold_OE/

| 文件 | 作用 |
|---|---|
| `results/tables/17_scTenifold_OE/Gene_COL1A1/COL1A1_OE_all_results.csv` | 保存全部差异调控或统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/17_scTenifold_OE/Gene_COL1A1/COL1A1_OE_significant_results.csv` | 保存达到显著性阈值的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/17_scTenifold_OE/Gene_EMP1/EMP1_OE_all_results.csv` | 保存全部差异调控或统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/17_scTenifold_OE/Gene_EMP1/EMP1_OE_significant_results.csv` | 保存达到显著性阈值的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/17_scTenifold_OE/Gene_MYH11/MYH11_OE_all_results.csv` | 保存全部差异调控或统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/17_scTenifold_OE/Gene_MYH11/MYH11_OE_significant_results.csv` | 保存达到显著性阈值的结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/17_scTenifold_OE/Gene_SASH1/SASH1_OE_all_results.csv` | 保存全部差异调控或统计结果，文件格式为CSV，可用于后续统计或绘图。 |
| `results/tables/17_scTenifold_OE/Gene_SASH1/SASH1_OE_significant_results.csv` | 保存达到显著性阈值的结果，文件格式为CSV，可用于后续统计或绘图。 |

### results/objects/

| 文件 | 作用 |
|---|---|
| `results/objects/02_raw_before_QC_filtering.rds` | 脚本02生成的Seurat或分析中间对象，供后续脚本继续读取。 |
| `results/objects/03_QC_reproduction_candidate.rds` | 脚本03生成的Seurat或分析中间对象，供后续脚本继续读取。 |
| `results/objects/04_standard_Seurat_multi_resolution.rds` | 脚本04生成的Seurat或分析中间对象，供后续脚本继续读取。 |
| `results/objects/05_diagnostic_manual_annotation_diagnostic.rds` | 脚本05生成的Seurat或分析中间对象，供后续脚本继续读取。 |
| `results/objects/05_manual_annotated_before_malignant_call.rds` | 脚本05生成的Seurat或分析中间对象，供后续脚本继续读取。 |
| `results/objects/05_manual_annotated_plot_ready.rds` | 脚本05生成的Seurat或分析中间对象，供后续脚本继续读取。 |
| `results/objects/06_malignant_candidate_diagnostic.rds` | 脚本06生成的Seurat或分析中间对象，供后续脚本继续读取。 |
| `results/objects/07_CopyKAT_malignant_call.rds` | 脚本07生成的Seurat或分析中间对象，供后续脚本继续读取。 |
| `results/objects/08_final_malignant_call.rds` | 脚本08生成的Seurat或分析中间对象，供后续脚本继续读取。 |
| `results/objects/10_malignant_epithelial_state_characterization.rds` | 脚本10生成的Seurat或分析中间对象，供后续脚本继续读取。 |
| `results/objects/10_manual_review_epithelial_core.rds` | 脚本10生成的Seurat或分析中间对象，供后续脚本继续读取。 |
| `results/objects/11_global_trajectory_Seurat.rds.parts/11_global_trajectory_Seurat.rds.part_aa` | 大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。 |
| `results/objects/11_global_trajectory_Seurat.rds.parts/11_global_trajectory_Seurat.rds.part_ab` | 大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。 |
| `results/objects/12_spatial_tissue_spots_Seurat.rds` | 空间转录组分析生成的Seurat中间对象，保存表达矩阵、空间坐标和新增metadata。 |
| `results/objects/13_spatial_SASH1_COL1A1_high_annotated_Seurat.rds` | 空间转录组分析生成的Seurat中间对象，保存表达矩阵、空间坐标和新增metadata。 |
| `results/objects/14_spatial_domain_annotated_paper_style_Seurat.rds` | 空间转录组分析生成的Seurat中间对象，保存表达矩阵、空间坐标和新增metadata。 |
| `results/objects/15_Figure4_like_scRNA_core_gene_expression_Seurat.rds.parts/15_Figure4_like_scRNA_core_gene_expression_Seurat.rds.part_aa` | 大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。 |
| `results/objects/15_Figure4_like_scRNA_core_gene_expression_Seurat.rds.parts/15_Figure4_like_scRNA_core_gene_expression_Seurat.rds.part_ab` | 大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。 |
| `results/objects/15b_Figure4_like_scRNA_cluster_core_gene_expression_Seurat.rds.parts/15b_Figure4_like_scRNA_cluster_core_gene_expression_Seurat.rds.part_aa` | 大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。 |
| `results/objects/15b_Figure4_like_scRNA_cluster_core_gene_expression_Seurat.rds.parts/15b_Figure4_like_scRNA_cluster_core_gene_expression_Seurat.rds.part_ab` | 大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。 |
| `results/objects/16_Figure5_like_malignant_trajectory_core_gene_Seurat.rds.parts/16_Figure5_like_malignant_trajectory_core_gene_Seurat.rds.part_aa` | 大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。 |
| `results/objects/16_Figure5_like_malignant_trajectory_core_gene_Seurat.rds.parts/16_Figure5_like_malignant_trajectory_core_gene_Seurat.rds.part_ab` | 大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。 |
| `results/objects/17_scTenifold_analysis_input_Seurat.rds` | 脚本17实际用于虚拟KO/OE分析的7353个cluster4、6和11细胞Seurat对象。 |

### results/objects/17_scTenifold_OE_networks/

| 文件 | 作用 |
|---|---|
| `results/objects/17_scTenifold_OE_networks/Gene_COL1A1/COL1A1_OE_networks.rds` | 脚本17保存的WT网络、虚拟过表达网络和manifold alignment结果。 |
| `results/objects/17_scTenifold_OE_networks/Gene_EMP1/EMP1_OE_networks.rds` | 脚本17保存的WT网络、虚拟过表达网络和manifold alignment结果。 |
| `results/objects/17_scTenifold_OE_networks/Gene_MYH11/MYH11_OE_networks.rds` | 脚本17保存的WT网络、虚拟过表达网络和manifold alignment结果。 |
| `results/objects/17_scTenifold_OE_networks/Gene_SASH1/SASH1_OE_networks.rds` | 脚本17保存的WT网络、虚拟过表达网络和manifold alignment结果。 |

## README自动更新

新增、删除或重命名Git追踪文件后，可重新运行：

```bash
python3 scripts_utils/generate_detailed_readme.py
```

该脚本会重新扫描`git ls-files`并更新本README中的脚本说明、脚本与输出对应关系以及完整文件清单。
