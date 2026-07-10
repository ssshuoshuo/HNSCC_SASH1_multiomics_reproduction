# HNSCC/OSCC SASH1多组学论文复现

本仓库用于复现HNSCC/OSCC中SASH1相关多组学研究的单细胞转录组和空间转录组分析模块。

当前版本聚焦单细胞RNA-seq和空间转录组相关复现，包括细胞注释、恶性细胞候选群体判断、核心基因表达图、恶性轨迹图、空间基因表达图、空间结构域近似注释，以及SASH1/COL1A1空间邻近关系分析。

## 当前已完成内容

- 单细胞RNA-seq数据下载、质控、降维、聚类和人工注释
- 上皮/肿瘤细胞群体复核
- CopyKAT辅助恶性细胞判断
- 核心恶性focus clusters和扩展恶性focus clusters定义
- Fig.4-like单细胞核心基因表达图
- Fig.5-like恶性细胞轨迹和核心基因表达图
- GSE252265空间转录组下载和质控
- Fig.6-like空间结构域和核心基因表达图
- SASH1-high和COL1A1-high空间共定位、最近邻距离和置换检验

## 仓库文件夹说明

- `config/`：配置文件目录，用于保存项目参数、路径或可复用配置。
- `data/`：输入数据目录，包含单细胞和空间转录组分析所需的原始或整理后数据。
- `data/raw/GSE252265/`：GSE252265空间转录组原始输入文件目录，用于脚本12–14。
- `scripts/`：主分析脚本目录，01–16按顺序构成当前复现流程。
- `scripts_utils/`：辅助脚本目录，目前用于split文件重建等非主流程工具。
- `split_file_manifest/`：大文件切分后的manifest目录，用于记录分片文件和重建信息。
- `results/`：分析结果目录，包含图、表和中间对象。
- `results/figures/`：生成的图像结果，包括Fig.4-like、Fig.5-like、Fig.6-like和空间邻近分析图。
- `results/tables/`：生成的表格结果，包括QC统计、表达汇总、空间邻近检验、轨迹metadata和sessionInfo。
- `results/objects/`：生成的中间对象目录，主要保存Seurat对象、Monocle轨迹对象或大文件分片。
- `tools/`：辅助工具目录，保存项目运行过程中需要的工具或本地辅助文件。

## 01–16主流程脚本说明

| 脚本 | 作用 |
|---|---|
| `scripts/01_download_and_prepare_scRNA.R` | 下载并整理GSE215403单细胞RNA-seq原始数据，建立后续分析所需的输入目录和基础文件。 |
| `scripts/02_read_and_QC_scRNA.R` | 读取单细胞表达矩阵，创建Seurat对象，进行基础质控、过滤和标准化前处理。 |
| `scripts/03_QC_reproduction_candidate.R` | 复现并检查论文单细胞质控指标，输出候选QC阈值、细胞数统计和相关图表。 |
| `scripts/04_standard_Seurat_PCA_UMAP_resolution_scan.R` | 执行标准Seurat流程，包括归一化、高变基因、PCA、UMAP和不同resolution下的聚类扫描。 |
| `scripts/05_manual_annotation_and_target_gene_summary.R` | 基于marker和聚类结果进行人工细胞类型注释，并总结SASH1等目标基因在各细胞类型中的表达。 |
| `scripts/06_malignant_candidate_diagnostic.R` | 对上皮/肿瘤相关cluster进行诊断，筛选可能的恶性细胞候选群体。 |
| `scripts/07_CopyKAT_malignant_call.R` | 运行CopyKAT推断拷贝数变化，用于辅助判断恶性细胞。 |
| `scripts/07_CopyKAT_malignant_call_batch.sh` | CopyKAT分析的批处理脚本，适合长时间运行或分批运行。 |
| `scripts/08_finalize_malignant_call.R` | 整合人工注释、cluster信息和CopyKAT结果，生成最终恶性细胞标签。 |
| `scripts/08_finalize_malignant_call_batch.sh` | 最终恶性细胞注释流程的批处理脚本。 |
| `scripts/09_malignant_cell_composition_check.R` | 检查恶性细胞在样本、cluster和细胞类型中的组成分布。 |
| `scripts/09_malignant_cell_composition_check_batch.sh` | 恶性细胞组成检查的批处理脚本。 |
| `scripts/10_manual_review_epithelial_core.R` | 对上皮细胞和核心肿瘤相关cluster进行人工复核。 |
| `scripts/11a_malignant_focus_cluster_audit.R` | 审查核心恶性focus cluster和候选恶性相关cluster，输出审查表和诊断图。 |
| `scripts/11b_core_extended_malignant_overlay_rebuild.R` | 重建核心/扩展恶性细胞状态覆盖图，并生成后续轨迹分析需要的输入文件。 |
| `scripts/12_spatial_download_QC_gene_maps.R` | 下载并读取GSE252265空间转录组数据，完成空间spot质控和SASH1、COL1A1、EMP1、MYH11空间表达图。 |
| `scripts/13_spatial_SASH1_COL1A1_neighborhood_analysis.R` | 分析SASH1-high和COL1A1-high spots的共定位、最近邻距离和空间置换检验。 |
| `scripts/14_spatial_domain_annotation_and_core_gene_maps.R` | 基于marker score近似注释空间结构域，并生成Fig.6-like空间结构域和核心基因表达图。 |
| `scripts/15_scRNA_core_gene_expression_Figure4_like.R` | 生成Fig.4-like单细胞核心基因表达图，重点展示主要细胞类型层面的SASH1、COL1A1、EMP1和MYH11表达。 |
| `scripts/15b_scRNA_core_gene_expression_Figure4_cluster_like.R` | 生成cluster/cell type级别的Fig.4-like核心基因表达图，更接近论文中按cluster展示表达的形式。 |
| `scripts/16_scRNA_malignant_trajectory_Figure5_like.R` | 生成Fig.5-like恶性细胞轨迹、核心基因表达和gene-high状态图。 |

## 运行顺序

建议按照脚本编号顺序运行。部分步骤耗时较长，CopyKAT相关脚本可以使用对应batch脚本。

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
```

## 恶性细胞注释策略

当前复现采用以下恶性细胞focus定义：

```text
核心恶性focus clusters：6和11
扩展恶性focus cluster：4
候选恶性相关肿瘤clusters：2和3
```

CopyKAT结果作为辅助证据之一，用于支持恶性细胞群体判断。

## 空间转录组分析说明

当前空间分析基于GSE252265提供的表达矩阵和空间坐标文件进行。由于公开文件中没有标准Seurat空间对象所需的完整H&E图像和按样本拆分的图像信息，因此目前采用基于坐标的空间表达分析和marker-score空间结构域近似注释。

当前限制：

```text
空间样本ID目前统一记为All_spots
没有H&E图像叠加
空间结构域为marker-score近似注释，不等同于原文作者基于病理图像的人工结构域注释
```

## 当前主要结果

- SASH1、COL1A1、EMP1和MYH11已经在主要细胞类型、cluster级别细胞类型和恶性轨迹空间中完成可视化。
- SASH1-high和COL1A1-high spots在当前空间坐标分析中重叠有限。
- 在当前置换检验框架下，SASH1-high和COL1A1-high没有表现出明显空间共定位富集。

## 大型中间RDS文件说明

脚本12–16生成的大型Seurat对象已经通过Git LFS保存。由于GitHub LFS单文件大小限制，部分超过2GB的RDS对象采用分片方式上传。

已作为完整RDS上传的对象：

```text
results/objects/12_spatial_tissue_spots_Seurat.rds
results/objects/13_spatial_SASH1_COL1A1_high_annotated_Seurat.rds
results/objects/14_spatial_domain_annotated_paper_style_Seurat.rds
```

已作为分片文件上传的对象：

```text
results/objects/15_Figure4_like_scRNA_core_gene_expression_Seurat.rds.parts/
results/objects/15b_Figure4_like_scRNA_cluster_core_gene_expression_Seurat.rds.parts/
results/objects/16_Figure5_like_malignant_trajectory_core_gene_Seurat.rds.parts/
```

分片RDS可以用以下脚本重建：

```bash
bash scripts_utils/reconstruct_split_rds_objects.sh
```

重建后会得到：

```text
results/objects/15_Figure4_like_scRNA_core_gene_expression_Seurat.rds
results/objects/15b_Figure4_like_scRNA_cluster_core_gene_expression_Seurat.rds
results/objects/16_Figure5_like_malignant_trajectory_core_gene_Seurat.rds
```


## GitHub文件清单与作用

本节自动列出当前Git已追踪文件，并说明每个文件或文件夹的作用。

### 根目录

| 文件 | 作用 |
|---|---|
| `.gitattributes` | Git LFS追踪规则，用于管理RDS、压缩包等大文件。 |
| `.gitignore` | Git忽略规则，避免上传临时文件、日志文件和暂不上传的大型中间RDS。 |
| `HNSCC_SASH1_reproduction.Rproj` | RStudio项目文件，便于在RStudio中打开整个复现项目。 |
| `README.md` | 项目说明文件，记录仓库结构、运行顺序、当前完成度、文件作用和注意事项。 |

### config/

配置文件目录，用于保存项目参数、路径或可复用配置。

| 文件 | 作用 |
|---|---|
| `config/GSE215403_sample_metadata.csv` | 项目配置文件。 |

### data/

输入数据目录，包含单细胞和空间转录组分析所需的原始或整理后数据。

| 文件 | 作用 |
|---|---|
| `data/processed/scRNA_GSE215403/10x_by_sample/OSCC/barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/OSCC/features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/OSCC/matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB1/barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB1/features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB1/matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB10/barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB10/features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB10/matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB12/barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB12/features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB12/matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB13/barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB13/features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB13/matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB14/barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB14/features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB14/matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB15/barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB15/features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB15/matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB2/barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB2/features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB2/matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB5/barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB5/features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB5/matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB7/barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB7/features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB7/matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB8/barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB8/features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB8/matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB9/barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB9/features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/10x_by_sample/scB9/matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634869_OSCC_barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634869_OSCC_features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634869_OSCC_matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634870_scB1_barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634870_scB1_features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634870_scB1_matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634871_scB2_barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634871_scB2_features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634871_scB2_matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634872_scB5_barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634872_scB5_features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634872_scB5_matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634873_scB7_barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634873_scB7_features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634873_scB7_matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634874_scB8_barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634874_scB8_features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634874_scB8_matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634875_scB9_barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634875_scB9_features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634875_scB9_matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634876_scB10_barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634876_scB10_features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634876_scB10_matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634877_scB12_barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634877_scB12_features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634877_scB12_matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634878_scB13_barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634878_scB13_features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634878_scB13_matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634879_scB14_barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634879_scB14_features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634879_scB14_matrix.mtx.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634880_scB15_barcodes.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634880_scB15_features.tsv.gz` | 项目文件。 |
| `data/processed/scRNA_GSE215403/GSM6634880_scB15_matrix.mtx.gz` | 项目文件。 |
| `data/raw/scRNA_GSE215403/GSE215403_RAW.tar` | 项目文件。 |

### data/raw/GSE252265/

GSE252265空间转录组原始输入文件目录，用于脚本12–14。

| 文件 | 作用 |
|---|---|
| `data/raw/GSE252265/GSE252265_RAW.tar` | GSE252265补充原始文件压缩包。 |
| `data/raw/GSE252265/GSE252265_aggr_tissue_positions.csv.gz` | GSE252265空间坐标表，用于将表达矩阵和组织位置对应。 |
| `data/raw/GSE252265/GSE252265_aggregation.csv.gz` | GSE252265聚合信息表。 |
| `data/raw/GSE252265/GSE252265_barcodes.tsv.gz` | GSE252265空间spot/cell barcode列表。 |
| `data/raw/GSE252265/GSE252265_features.tsv.gz` | GSE252265基因/feature注释表。 |
| `data/raw/GSE252265/GSE252265_filtered_feature_bc_matrix.h5` | GSE252265过滤后的feature-barcode表达矩阵，HDF5格式。 |
| `data/raw/GSE252265/GSE252265_matrix.mtx.gz` | GSE252265表达矩阵，Matrix Market压缩格式。 |

### scripts/

主分析脚本目录，01–16按顺序构成当前复现流程。

| 文件 | 作用 |
|---|---|
| `scripts/01_download_and_prepare_scRNA.R` | 下载并整理GSE215403单细胞RNA-seq原始数据，建立后续分析所需的输入目录和基础文件。 |
| `scripts/02_read_and_QC_scRNA.R` | 读取单细胞表达矩阵，创建Seurat对象，进行基础质控、过滤和标准化前处理。 |
| `scripts/03_QC_reproduction_candidate.R` | 复现并检查论文单细胞质控指标，输出候选QC阈值、细胞数统计和相关图表。 |
| `scripts/04_standard_Seurat_PCA_UMAP_resolution_scan.R` | 执行标准Seurat流程，包括归一化、高变基因、PCA、UMAP和不同resolution下的聚类扫描。 |
| `scripts/05_manual_annotation_and_target_gene_summary.R` | 基于marker和聚类结果进行人工细胞类型注释，并总结SASH1等目标基因在各细胞类型中的表达。 |
| `scripts/06_malignant_candidate_diagnostic.R` | 对上皮/肿瘤相关cluster进行诊断，筛选可能的恶性细胞候选群体。 |
| `scripts/07_CopyKAT_malignant_call.R` | 运行CopyKAT推断拷贝数变化，用于辅助判断恶性细胞。 |
| `scripts/07_CopyKAT_malignant_call_batch.sh` | CopyKAT分析的批处理脚本，适合长时间运行或分批运行。 |
| `scripts/08_finalize_malignant_call.R` | 整合人工注释、cluster信息和CopyKAT结果，生成最终恶性细胞标签。 |
| `scripts/08_finalize_malignant_call_batch.sh` | 最终恶性细胞注释流程的批处理脚本。 |
| `scripts/09_malignant_cell_composition_check.R` | 检查恶性细胞在样本、cluster和细胞类型中的组成分布。 |
| `scripts/09_malignant_cell_composition_check_batch.sh` | 恶性细胞组成检查的批处理脚本。 |
| `scripts/10_manual_review_epithelial_core.R` | 对上皮细胞和核心肿瘤相关cluster进行人工复核。 |
| `scripts/11a_malignant_focus_cluster_audit.R` | 审查核心恶性focus cluster和候选恶性相关cluster，输出审查表和诊断图。 |
| `scripts/11b_core_extended_malignant_overlay_rebuild.R` | 重建核心/扩展恶性细胞状态覆盖图，并生成后续轨迹分析需要的输入文件。 |
| `scripts/12_spatial_download_QC_gene_maps.R` | 下载并读取GSE252265空间转录组数据，完成空间spot质控和SASH1、COL1A1、EMP1、MYH11空间表达图。 |
| `scripts/13_spatial_SASH1_COL1A1_neighborhood_analysis.R` | 分析SASH1-high和COL1A1-high spots的共定位、最近邻距离和空间置换检验。 |
| `scripts/14_spatial_domain_annotation_and_core_gene_maps.R` | 基于marker score近似注释空间结构域，并生成Fig.6-like空间结构域和核心基因表达图。 |
| `scripts/15_scRNA_core_gene_expression_Figure4_like.R` | 生成Fig.4-like单细胞核心基因表达图，重点展示主要细胞类型层面的SASH1、COL1A1、EMP1和MYH11表达。 |
| `scripts/15b_scRNA_core_gene_expression_Figure4_cluster_like.R` | 生成cluster/cell type级别的Fig.4-like核心基因表达图，更接近论文中按cluster展示表达的形式。 |
| `scripts/16_scRNA_malignant_trajectory_Figure5_like.R` | 生成Fig.5-like恶性细胞轨迹、核心基因表达和gene-high状态图。 |

### scripts_utils/

辅助脚本目录，目前用于split文件重建等非主流程工具。

| 文件 | 作用 |
|---|---|
| `scripts_utils/reconstruct_split_files.sh` | 根据split_file_manifest和.parts文件夹重建被切分的大型RDS文件。 |

### split_file_manifest/

大文件切分后的manifest目录，用于记录分片文件和重建信息。

| 文件 | 作用 |
|---|---|
| `split_file_manifest/oversized_files.txt` | 大文件切分和重建所需的manifest记录。 |

### results/figures/

生成的图像结果，包括Fig.4-like、Fig.5-like、Fig.6-like和空间邻近分析图。

| 文件 | 作用 |
|---|---|
| `results/figures/02_QC_scatter_nCount_vs_feature.pdf` | 脚本02生成的结果图。 |
| `results/figures/02_QC_scatter_nCount_vs_mt.pdf` | 脚本02生成的结果图。 |
| `results/figures/02_QC_violin_by_sample.pdf` | 脚本02生成的结果图。 |
| `results/figures/02_QC_violin_by_sample.png` | 脚本02生成的结果图。 |
| `results/figures/03_QC_before_after_reproduction_candidate.pdf` | 脚本03生成的结果图。 |
| `results/figures/03_QC_before_after_reproduction_candidate.png` | 脚本03生成的结果图。 |
| `results/figures/03_QC_violin_after_filtering.pdf` | 脚本03生成的结果图。 |
| `results/figures/03_cell_number_before_after_QC.pdf` | 脚本03生成的结果图。 |
| `results/figures/03_cell_number_before_after_QC.png` | 脚本03生成的结果图。 |
| `results/figures/03_cell_number_old03_vs_reproduction_candidate.pdf` | 脚本03生成的结果图。 |
| `results/figures/04_Harmony_UMAP_sample_and_cluster.pdf` | 脚本04生成的结果图。 |
| `results/figures/04_Harmony_UMAP_sample_and_cluster.png` | 脚本04生成的结果图。 |
| `results/figures/04_PCA_elbow_plot.pdf` | 脚本04生成的结果图。 |
| `results/figures/04_PCA_elbow_plot.png` | 脚本04生成的结果图。 |
| `results/figures/04_UMAP_cluster_resolution_0.2.pdf` | 脚本04生成的结果图。 |
| `results/figures/04_UMAP_cluster_resolution_0.3.pdf` | 脚本04生成的结果图。 |
| `results/figures/04_UMAP_cluster_resolution_0.5.pdf` | 脚本04生成的结果图。 |
| `results/figures/04_UMAP_sample_and_primary_cluster.pdf` | 脚本04生成的结果图。 |
| `results/figures/04_UMAP_sample_and_primary_cluster.png` | 脚本04生成的结果图。 |
| `results/figures/05_UMAP_cluster_and_manual_celltype.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_UMAP_cluster_and_manual_celltype.png` | 脚本05生成的结果图。 |
| `results/figures/05_UMAP_cluster_and_manual_celltype_polished.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_UMAP_cluster_and_manual_celltype_polished.png` | 脚本05生成的结果图。 |
| `results/figures/05_UMAP_cluster_and_preliminary_annotation.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_UMAP_cluster_and_preliminary_annotation.png` | 脚本05生成的结果图。 |
| `results/figures/05_canonical_marker_DotPlot.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_canonical_marker_DotPlot.png` | 脚本05生成的结果图。 |
| `results/figures/05_diagnostic_FeaturePlot_Epithelial_markers.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_diagnostic_FeaturePlot_Immune_markers.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_diagnostic_FeaturePlot_Pericyte_mural_markers.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_diagnostic_FeaturePlot_Stromal_vascular_markers.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_diagnostic_cluster_by_sample_heatmap_resolution_0.2.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_diagnostic_major_lineage_DotPlot_resolution_0.2.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_diagnostic_major_lineage_DotPlot_resolution_0.2.png` | 脚本05生成的结果图。 |
| `results/figures/05_diagnostic_target_genes_UMAP_resolution_0.2.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_major_cell_populations_UMAP_final.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_major_cell_populations_UMAP_final.png` | 脚本05生成的结果图。 |
| `results/figures/05_manual_celltype_UMAP_paper_style.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_manual_celltype_UMAP_paper_style.png` | 脚本05生成的结果图。 |
| `results/figures/05_target_genes_DotPlot_by_manual_celltype.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_target_genes_UMAP.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_target_genes_UMAP.png` | 脚本05生成的结果图。 |
| `results/figures/05_target_genes_UMAP_quantile_scaled.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_target_genes_VlnPlot_by_manual_celltype.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_target_genes_by_preliminary_celltype.pdf` | 脚本05生成的结果图。 |
| `results/figures/05_top5_markers_heatmap.pdf` | 脚本05生成的结果图。 |
| `results/figures/06_epithelial_pseudotime_UMAP.pdf` | 脚本06生成的结果图。 |
| `results/figures/06_epithelial_pseudotime_UMAP.png` | 脚本06生成的结果图。 |
| `results/figures/06_epithelial_salivary_tumor_marker_DotPlot.pdf` | 脚本06生成的结果图。 |
| `results/figures/06_epithelial_state_UMAP.pdf` | 脚本06生成的结果图。 |
| `results/figures/06_epithelial_state_scores_along_pseudotime.pdf` | 脚本06生成的结果图。 |
| `results/figures/06_epithelial_subset_UMAP_cluster_and_sample.pdf` | 脚本06生成的结果图。 |
| `results/figures/06_epithelial_subset_UMAP_cluster_and_sample.png` | 脚本06生成的结果图。 |
| `results/figures/06_gene_expression_trends_along_pseudotime.pdf` | 脚本06生成的结果图。 |
| `results/figures/06_target_genes_DotPlot_by_diagnostic_status.pdf` | 脚本06生成的结果图。 |
| `results/figures/06_target_genes_VlnPlot_by_diagnostic_status.pdf` | 脚本06生成的结果图。 |
| `results/figures/06_target_genes_by_epithelial_cluster.pdf` | 脚本06生成的结果图。 |
| `results/figures/06_target_genes_in_epithelial_like_UMAP.pdf` | 脚本06生成的结果图。 |
| `results/figures/06_tumor_candidate_and_salivary_reference_UMAP.pdf` | 脚本06生成的结果图。 |
| `results/figures/06_tumor_candidate_and_salivary_reference_UMAP.png` | 脚本06生成的结果图。 |
| `results/figures/06c_chromosome_aware_CNV_expression_heatmap.pdf` | 脚本06生成的结果图。 |
| `results/figures/06c_chromosome_aware_CNV_expression_scores.pdf` | 脚本06生成的结果图。 |
| `results/figures/06c_chromosome_level_deviation_heatmap.pdf` | 脚本06生成的结果图。 |
| `results/figures/07_CopyKAT_aneuploid_tumor_candidate_UMAP.pdf` | 脚本07生成的结果图。 |
| `results/figures/07_CopyKAT_prediction_UMAP.pdf` | 脚本07生成的结果图。 |
| `results/figures/07a_cluster6_strict_malignant_UMAP_sample_and_subcluster.pdf` | 脚本07生成的结果图。 |
| `results/figures/07a_cluster6_strict_malignant_core_gene_UMAP.pdf` | 脚本07生成的结果图。 |
| `results/figures/07a_cluster6_strict_malignant_core_gene_by_sample.pdf` | 脚本07生成的结果图。 |
| `results/figures/07b_refined_epithelial_malignant_UMAP_sample_and_subcluster.pdf` | 脚本07生成的结果图。 |
| `results/figures/07b_refined_epithelial_malignant_core_gene_UMAP.pdf` | 脚本07生成的结果图。 |
| `results/figures/07b_strict_malignant_epithelial_refinement_UMAP.pdf` | 脚本07生成的结果图。 |
| `results/figures/07b_strict_malignant_lineage_module_scores_UMAP.pdf` | 脚本07生成的结果图。 |
| `results/figures/07c_core_gene_expression_by_malignant_state.pdf` | 脚本07生成的结果图。 |
| `results/figures/07c_malignant_state_cluster_and_annotation_UMAP.pdf` | 脚本07生成的结果图。 |
| `results/figures/07c_malignant_state_module_scores_UMAP.pdf` | 脚本07生成的结果图。 |
| `results/figures/08_CopyKAT_aneuploid_fraction_by_cluster.pdf` | 脚本08生成的结果图。 |
| `results/figures/08_CopyKAT_aneuploid_fraction_by_sample_cluster.pdf` | 脚本08生成的结果图。 |
| `results/figures/08_core_gene_expression_by_final_status.pdf` | 脚本08生成的结果图。 |
| `results/figures/08_final_malignant_status_UMAP.pdf` | 脚本08生成的结果图。 |
| `results/figures/08a_epithelial_core_sample_cluster_heatmap.pdf` | 脚本08生成的结果图。 |
| `results/figures/08a_within_sample_trajectory_candidate_UMAP_panel.pdf` | 脚本08生成的结果图。 |
| `results/figures/09_strict_malignant_cell_composition_by_cluster.pdf` | 脚本09生成的结果图。 |
| `results/figures/09_strict_malignant_cell_composition_by_sample.pdf` | 脚本09生成的结果图。 |
| `results/figures/09_strict_malignant_cells_by_sample_cluster.pdf` | 脚本09生成的结果图。 |
| `results/figures/10_core_gene_expression_by_epithelial_core_state.pdf` | 脚本10生成的结果图。 |
| `results/figures/10_epithelial_core_UMAP_sample_cluster_program.pdf` | 脚本10生成的结果图。 |
| `results/figures/10_epithelial_core_relative_program_scores_UMAP.pdf` | 脚本10生成的结果图。 |
| `results/figures/10_manual_review_cluster_decision_UMAP.pdf` | 脚本10生成的结果图。 |
| `results/figures/11_global_cellular_trajectory_vertex_bins_paper_style.pdf` | 脚本11生成的结果图。 |
| `results/figures/11_global_cellular_trajectory_vertex_groups.pdf` | 脚本11生成的结果图。 |
| `results/figures/11_global_trajectory_major_cell_type.pdf` | 脚本11生成的结果图。 |
| `results/figures/11_global_trajectory_malignant_focused_cells.pdf` | 脚本11生成的结果图。 |
| `results/figures/11_malignant_focused_gene_trajectory_overlays.pdf` | 脚本11生成的结果图。 |
| `results/figures/11_malignant_focused_gene_trajectory_overlays_paper_style.pdf` | 脚本11生成的结果图。 |
| `results/figures/11a_candidate_cluster_CopyKAT_support.pdf` | 脚本11生成的结果图。 |
| `results/figures/11a_candidate_cluster_marker_DotPlot.pdf` | 脚本11生成的结果图。 |
| `results/figures/11a_candidate_cluster_sample_composition_heatmap.pdf` | 脚本11生成的结果图。 |
| `results/figures/11a_candidate_cluster_vertex_bin_distribution.pdf` | 脚本11生成的结果图。 |
| `results/figures/11a_candidate_malignant_clusters_UMAP.pdf` | 脚本11生成的结果图。 |
| `results/figures/11b_core_malignant_focused_gene_overlays.pdf` | 脚本11生成的结果图。 |
| `results/figures/11b_core_vs_extended_malignant_focus_comparison.pdf` | 脚本11生成的结果图。 |
| `results/figures/11b_extended_malignant_focused_gene_overlays.pdf` | 脚本11生成的结果图。 |
| `results/figures/12_SASH1_COL1A1_EMP1_MYH11_spatial_expression.pdf` | 脚本12生成的结果图。 |
| `results/figures/12_spatial_QC_UMI_distribution.pdf` | 脚本12生成的结果图。 |
| `results/figures/12_spatial_QC_detected_gene_distribution.pdf` | 脚本12生成的结果图。 |
| `results/figures/12_spatial_QC_percent_mt_distribution.pdf` | 脚本12生成的结果图。 |
| `results/figures/13_SASH1_COL1A1_high_colocalization_barplot.pdf` | 脚本13生成的结果图。 |
| `results/figures/13_SASH1_COL1A1_high_neighborhood_permutation_test.pdf` | 脚本13生成的结果图。 |
| `results/figures/13_SASH1_high_COL1A1_high_spatial_overlay.pdf` | 脚本13生成的结果图。 |
| `results/figures/13_core_gene_spatial_expression_panel.pdf` | 脚本13生成的结果图。 |
| `results/figures/14_spatial_domain_COL1A1_expression.pdf` | 脚本14生成的结果图。 |
| `results/figures/14_spatial_domain_EMP1_expression.pdf` | 脚本14生成的结果图。 |
| `results/figures/14_spatial_domain_Figure6_like_panel.pdf` | 脚本14生成的结果图。 |
| `results/figures/14_spatial_domain_MYH11_expression.pdf` | 脚本14生成的结果图。 |
| `results/figures/14_spatial_domain_SASH1_expression.pdf` | 脚本14生成的结果图。 |
| `results/figures/14_spatial_domain_annotation.pdf` | 脚本14生成的结果图。 |
| `results/figures/14_spatial_domain_tissue_layout.pdf` | 脚本14生成的结果图。 |
| `results/figures/15_Figure4_core_gene_FeaturePlot_panel.pdf` | 脚本15生成的结果图。 |
| `results/figures/15_Figure4_core_gene_VlnPlot_panel.pdf` | 脚本15生成的结果图。 |
| `results/figures/15_Figure4_like_scRNA_core_gene_expression_panel.pdf` | 脚本15生成的结果图。 |
| `results/figures/15_Figure4a_scRNA_celltype_UMAP.pdf` | 脚本15生成的结果图。 |
| `results/figures/15_Figure4a_scRNA_celltype_UMAP.png` | 脚本15生成的结果图。 |
| `results/figures/15_Figure4b_core_gene_DotPlot_by_celltype.pdf` | 脚本15生成的结果图。 |
| `results/figures/15b_Figure4_FeaturePlot_COL1A1.pdf` | 脚本15b生成的结果图。 |
| `results/figures/15b_Figure4_FeaturePlot_EMP1.pdf` | 脚本15b生成的结果图。 |
| `results/figures/15b_Figure4_FeaturePlot_MYH11.pdf` | 脚本15b生成的结果图。 |
| `results/figures/15b_Figure4_FeaturePlot_SASH1.pdf` | 脚本15b生成的结果图。 |
| `results/figures/15b_Figure4_VlnPlot_COL1A1_by_cluster_celltype.pdf` | 脚本15b生成的结果图。 |
| `results/figures/15b_Figure4_VlnPlot_EMP1_by_cluster_celltype.pdf` | 脚本15b生成的结果图。 |
| `results/figures/15b_Figure4_VlnPlot_MYH11_by_cluster_celltype.pdf` | 脚本15b生成的结果图。 |
| `results/figures/15b_Figure4_VlnPlot_SASH1_by_cluster_celltype.pdf` | 脚本15b生成的结果图。 |
| `results/figures/15b_Figure4_cluster_celltype_UMAP.pdf` | 脚本15b生成的结果图。 |
| `results/figures/15b_Figure4_core_gene_FeaturePlot_panel.pdf` | 脚本15b生成的结果图。 |
| `results/figures/15b_Figure4_core_gene_VlnPlot_by_cluster_celltype_panel.pdf` | 脚本15b生成的结果图。 |
| `results/figures/15b_Figure4_like_scRNA_core_gene_expression_cluster_panel.pdf` | 脚本15b生成的结果图。 |
| `results/figures/15b_Figure4a_scRNA_celltype_UMAP.pdf` | 脚本15b生成的结果图。 |
| `results/figures/15b_Figure4a_scRNA_celltype_UMAP.png` | 脚本15b生成的结果图。 |
| `results/figures/15b_Figure4b_core_gene_DotPlot_by_cluster_celltype.pdf` | 脚本15b生成的结果图。 |
| `results/figures/16_Figure5_core_gene_high_status_trajectory_panel.pdf` | 脚本16生成的结果图。 |
| `results/figures/16_Figure5_core_gene_trajectory_expression_panel.pdf` | 脚本16生成的结果图。 |
| `results/figures/16_Figure5_like_malignant_trajectory_core_gene_panel.pdf` | 脚本16生成的结果图。 |
| `results/figures/16_Figure5_supplementary_core_gene_pseudotime_or_proxy_trend.pdf` | 脚本16生成的结果图。 |
| `results/figures/16_Figure5_trajectory_expression_COL1A1.pdf` | 脚本16生成的结果图。 |
| `results/figures/16_Figure5_trajectory_expression_EMP1.pdf` | 脚本16生成的结果图。 |
| `results/figures/16_Figure5_trajectory_expression_MYH11.pdf` | 脚本16生成的结果图。 |
| `results/figures/16_Figure5_trajectory_expression_SASH1.pdf` | 脚本16生成的结果图。 |
| `results/figures/16_Figure5_trajectory_high_status_COL1A1.pdf` | 脚本16生成的结果图。 |
| `results/figures/16_Figure5_trajectory_high_status_EMP1.pdf` | 脚本16生成的结果图。 |
| `results/figures/16_Figure5_trajectory_high_status_MYH11.pdf` | 脚本16生成的结果图。 |
| `results/figures/16_Figure5_trajectory_high_status_SASH1.pdf` | 脚本16生成的结果图。 |
| `results/figures/16_Figure5a_trajectory_status.pdf` | 脚本16生成的结果图。 |
| `results/figures/16_Figure5b_trajectory_celltype.pdf` | 脚本16生成的结果图。 |
| `results/figures/16_Figure5c_supplementary_pseudotime_or_proxy.pdf` | 脚本16生成的结果图。 |

### results/tables/

生成的表格结果，包括QC统计、表达汇总、空间邻近检验、轨迹metadata和sessionInfo。

| 文件 | 作用 |
|---|---|
| `results/tables/01_10x_folder_check.csv` | 脚本01生成的结果表或统计汇总。 |
| `results/tables/01_core_package_check.csv` | 脚本01生成的结果表或统计汇总。 |
| `results/tables/01_raw_file_list.csv` | 脚本01生成的结果表或统计汇总。 |
| `results/tables/01_raw_sample_file_check.csv` | 脚本01生成的结果表或统计汇总。 |
| `results/tables/02_QC_threshold_suggestions.csv` | 脚本02生成的结果表或统计汇总。 |
| `results/tables/02_sample_QC_summary.csv` | 脚本02生成的结果表或统计汇总。 |
| `results/tables/02_sessionInfo.txt` | 脚本02运行时的R环境和包版本记录。 |
| `results/tables/02_target_gene_check.csv` | 脚本02生成的结果表或统计汇总。 |
| `results/tables/03_QC_comparison_old03_vs_reproduction_candidate.csv` | 脚本03生成的结果表或统计汇总。 |
| `results/tables/03_QC_filter_summary_by_sample.csv` | 脚本03生成的结果表或统计汇总。 |
| `results/tables/03_QC_filter_thresholds_by_sample.csv` | 脚本03生成的结果表或统计汇总。 |
| `results/tables/03_QC_reproduction_candidate_parameters.csv` | 脚本03生成的结果表或统计汇总。 |
| `results/tables/03_QC_reproduction_candidate_summary.csv` | 脚本03生成的结果表或统计汇总。 |
| `results/tables/03_sessionInfo.txt` | 脚本03运行时的R环境和包版本记录。 |
| `results/tables/04_cluster_by_sample_cell_numbers.csv` | 脚本04生成的结果表或统计汇总。 |
| `results/tables/04_cluster_summary_resolution_0.2.csv` | 脚本04生成的结果表或统计汇总。 |
| `results/tables/04_cluster_summary_resolution_0.3.csv` | 脚本04生成的结果表或统计汇总。 |
| `results/tables/04_cluster_summary_resolution_0.5.csv` | 脚本04生成的结果表或统计汇总。 |
| `results/tables/04_sessionInfo.txt` | 脚本04运行时的R环境和包版本记录。 |
| `results/tables/04_variable_features.csv` | 脚本04生成的结果表或统计汇总。 |
| `results/tables/04_variable_features_2000.csv` | 脚本04生成的结果表或统计汇总。 |
| `results/tables/05_DotPlot_marker_genes_used.csv` | 脚本05生成的结果表或统计汇总。 |
| `results/tables/05_all_cluster_markers.csv` | 脚本05生成的结果表或统计汇总。 |
| `results/tables/05_cluster_preliminary_annotation_template.csv` | 脚本05生成的结果表或统计汇总。 |
| `results/tables/05_diagnostic_all_markers_resolution_0.2.csv` | 脚本05生成的结果表或统计汇总。 |
| `results/tables/05_diagnostic_cluster_QC_summary_resolution_0.2.csv` | 脚本05生成的结果表或统计汇总。 |
| `results/tables/05_diagnostic_cluster_by_sample_resolution_0.2.csv` | 脚本05生成的结果表或统计汇总。 |
| `results/tables/05_diagnostic_major_lineage_markers_used.csv` | 脚本05生成的结果表或统计汇总。 |
| `results/tables/05_diagnostic_manual_annotation_template_resolution_0.2.csv` | 脚本05生成的结果表或统计汇总。 |
| `results/tables/05_diagnostic_sessionInfo.txt` | 脚本05运行时的R环境和包版本记录。 |
| `results/tables/05_diagnostic_top30_markers_resolution_0.2.csv` | 脚本05生成的结果表或统计汇总。 |
| `results/tables/05_manual_annotation_table_resolution_0.2.csv` | 脚本05生成的结果表或统计汇总。 |
| `results/tables/05_marker_panel_gene_check.csv` | 脚本05生成的结果表或统计汇总。 |
| `results/tables/05_sessionInfo.txt` | 脚本05运行时的R环境和包版本记录。 |
| `results/tables/05_target_gene_expression_by_cluster.csv` | 脚本05生成的结果表或统计汇总。 |
| `results/tables/05_target_gene_expression_by_manual_celltype.csv` | 脚本05生成的结果表或统计汇总。 |
| `results/tables/05_top15_markers_by_cluster.csv` | 脚本05生成的结果表或统计汇总。 |
| `results/tables/06_diagnostic_gene_check.csv` | 脚本06生成的结果表或统计汇总。 |
| `results/tables/06_epithelial_cell_level_pseudotime.csv` | 脚本06生成的结果表或统计汇总。 |
| `results/tables/06_epithelial_cluster_cell_numbers.csv` | 脚本06生成的结果表或统计汇总。 |
| `results/tables/06_epithelial_cluster_state_scores.csv` | 脚本06生成的结果表或统计汇总。 |
| `results/tables/06_epithelial_state_marker_genes_used.csv` | 脚本06生成的结果表或统计汇总。 |
| `results/tables/06_malignant_candidate_by_sample.csv` | 脚本06生成的结果表或统计汇总。 |
| `results/tables/06_malignant_candidate_cluster_summary.csv` | 脚本06生成的结果表或统计汇总。 |
| `results/tables/06_pseudotime_root_cluster.csv` | 脚本06生成的结果表或统计汇总。 |
| `results/tables/06_sessionInfo.txt` | 脚本06运行时的R环境和包版本记录。 |
| `results/tables/06_target_gene_expression_by_diagnostic_status.csv` | 脚本06生成的结果表或统计汇总。 |
| `results/tables/06_target_gene_pseudobulk_by_sample.csv` | 脚本06生成的结果表或统计汇总。 |
| `results/tables/06c_chromosome_aware_CNV_expression_scores.csv` | 脚本06生成的结果表或统计汇总。 |
| `results/tables/07_CopyKAT_run_summary.csv` | 脚本07生成的结果表或统计汇总。 |
| `results/tables/07_CopyKAT_tumor_candidate_prediction_by_sample_cluster.csv` | 脚本07生成的结果表或统计汇总。 |
| `results/tables/07_CopyKAT_tumor_candidate_prediction_overall.csv` | 脚本07生成的结果表或统计汇总。 |
| `results/tables/07a_cluster6_internal_subcluster_top20_markers.csv` | 脚本07生成的结果表或统计汇总。 |
| `results/tables/07a_cluster6_shared_state_selected_samples.csv` | 脚本07生成的结果表或统计汇总。 |
| `results/tables/07a_cluster6_top20_sample_diagnostic_markers.csv` | 脚本07生成的结果表或统计汇总。 |
| `results/tables/07b_refined_epithelial_malignant_by_sample_cluster.csv` | 脚本07生成的结果表或统计汇总。 |
| `results/tables/07b_refinement_summary.csv` | 脚本07生成的结果表或统计汇总。 |
| `results/tables/07c_core_gene_expression_by_sample_and_malignant_state.csv` | 脚本07生成的结果表或统计汇总。 |
| `results/tables/07c_internal_malignant_state_markers_limited_features.csv` | 脚本07生成的结果表或统计汇总。 |
| `results/tables/07c_internal_malignant_state_top20_markers.csv` | 脚本07生成的结果表或统计汇总。 |
| `results/tables/07c_malignant_state_cluster_summary.csv` | 脚本07生成的结果表或统计汇总。 |
| `results/tables/07c_sessionInfo.txt` | 脚本07运行时的R环境和包版本记录。 |
| `results/tables/07c_state_gene_sets_available.csv` | 脚本07生成的结果表或统计汇总。 |
| `results/tables/08_aneuploid_fraction_by_cluster.csv` | 脚本08生成的结果表或统计汇总。 |
| `results/tables/08_aneuploid_fraction_by_sample_cluster.csv` | 脚本08生成的结果表或统计汇总。 |
| `results/tables/08_core_gene_expression_by_final_status.csv` | 脚本08生成的结果表或统计汇总。 |
| `results/tables/08a_epithelial_core_cell_metadata_for_trajectory.csv` | 脚本08生成的结果表或统计汇总。 |
| `results/tables/08a_epithelial_core_sample_cluster_composition.csv` | 脚本08生成的结果表或统计汇总。 |
| `results/tables/08a_sessionInfo.txt` | 脚本08运行时的R环境和包版本记录。 |
| `results/tables/08a_within_sample_trajectory_feasibility_summary.csv` | 脚本08生成的结果表或统计汇总。 |
| `results/tables/08c_output_file_check.csv` | 脚本08生成的结果表或统计汇总。 |
| `results/tables/08c_sessionInfo.txt` | 脚本08运行时的R环境和包版本记录。 |
| `results/tables/09_pseudotime_sample_recommendation.csv` | 脚本09生成的结果表或统计汇总。 |
| `results/tables/09_strict_malignant_cells_per_sample_summary.csv` | 脚本09生成的结果表或统计汇总。 |
| `results/tables/10_core_gene_expression_by_sample_cluster_and_program.csv` | 脚本10生成的结果表或统计汇总。 |
| `results/tables/10_epithelial_core_by_sample_cluster_and_program.csv` | 脚本10生成的结果表或统计汇总。 |
| `results/tables/10_epithelial_core_relative_program_summary.csv` | 脚本10生成的结果表或统计汇总。 |
| `results/tables/10_manual_review_cell_summary.csv` | 脚本10生成的结果表或统计汇总。 |
| `results/tables/10_manual_review_cluster_decision.csv` | 脚本10生成的结果表或统计汇总。 |
| `results/tables/10_relative_program_gene_sets_available.csv` | 脚本10生成的结果表或统计汇总。 |
| `results/tables/10_sessionInfo.txt` | 脚本10运行时的R环境和包版本记录。 |
| `results/tables/11_global_trajectory_cell_metadata.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11_global_trajectory_cell_metadata_with_vertex_bins.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11_malignant_focused_gene_overlay_summary.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11_output_file_check.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11_paper_faithful_cell_type_and_malignant_focus_summary.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11_paper_style_malignant_focused_gene_overlay_summary.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11_paper_style_vertex_bin_cell_type_summary.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11_paper_style_vertex_bin_definition.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11_principal_graph_edge_coordinates.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11_principal_graph_vertex_coordinates.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11_sessionInfo.txt` | 脚本11运行时的R环境和包版本记录。 |
| `results/tables/11_vertex_group_cell_type_summary.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11a_candidate_cluster_CopyKAT_support_summary.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11a_candidate_cluster_integrated_review_summary.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11a_candidate_cluster_sample_composition.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11a_candidate_cluster_size_summary.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11a_candidate_cluster_vertex_bin_distribution.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11a_candidate_malignant_cluster_cell_metadata.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11a_sessionInfo.txt` | 脚本11运行时的R环境和包版本记录。 |
| `results/tables/11b_core_extended_gene_overlay_summary.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11b_malignant_focus_definition_summary.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11b_monocle3_cell_umap_coordinates_and_focus_labels.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11b_output_file_check.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11b_principal_graph_edge_coordinates.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11b_principal_graph_vertex_coordinates.csv` | 脚本11生成的结果表或统计汇总。 |
| `results/tables/11b_sessionInfo.txt` | 脚本11运行时的R环境和包版本记录。 |
| `results/tables/12_SASH1_COL1A1_spot_detection_summary.csv` | 脚本12生成的结果表或统计汇总。 |
| `results/tables/12_barcode_coordinate_match_summary.csv` | 脚本12生成的结果表或统计汇总。 |
| `results/tables/12_input_file_inventory.csv` | 脚本12生成的结果表或统计汇总。 |
| `results/tables/12_output_file_check.csv` | 脚本12生成的结果表或统计汇总。 |
| `results/tables/12_raw_tissue_positions_table.csv` | 脚本12生成的结果表或统计汇总。 |
| `results/tables/12_sessionInfo.txt` | 脚本12运行时的R环境和包版本记录。 |
| `results/tables/12_spatial_spot_QC_cell_metadata.csv` | 脚本12生成的结果表或统计汇总。 |
| `results/tables/12_spatial_spot_QC_summary.csv` | 脚本12生成的结果表或统计汇总。 |
| `results/tables/12_spatial_tissue_spot_expression_metadata.csv` | 脚本12生成的结果表或统计汇总。 |
| `results/tables/13_SASH1_COL1A1_high_colocalization_summary.csv` | 脚本13生成的结果表或统计汇总。 |
| `results/tables/13_SASH1_COL1A1_high_neighborhood_permutation_summary.csv` | 脚本13生成的结果表或统计汇总。 |
| `results/tables/13_SASH1_COL1A1_high_thresholds_by_sample.csv` | 脚本13生成的结果表或统计汇总。 |
| `results/tables/13_output_file_check.csv` | 脚本13生成的结果表或统计汇总。 |
| `results/tables/13_permutation_distribution_All_spots.csv` | 脚本13生成的结果表或统计汇总。 |
| `results/tables/13_sessionInfo.txt` | 脚本13运行时的R环境和包版本记录。 |
| `results/tables/13_spatial_gene_expression_with_coordinates.csv` | 脚本13生成的结果表或统计汇总。 |
| `results/tables/13_spatial_high_status_metadata.csv` | 脚本13生成的结果表或统计汇总。 |
| `results/tables/14_cluster_domain_marker_score_summary.csv` | 脚本14生成的结果表或统计汇总。 |
| `results/tables/14_cluster_paper_style_domain_annotation.csv` | 脚本14生成的结果表或统计汇总。 |
| `results/tables/14_domain_core_gene_expression_summary.csv` | 脚本14生成的结果表或统计汇总。 |
| `results/tables/14_domain_marker_genes_found.csv` | 脚本14生成的结果表或统计汇总。 |
| `results/tables/14_output_file_check.csv` | 脚本14生成的结果表或统计汇总。 |
| `results/tables/14_possible_spatial_image_files.csv` | 脚本14生成的结果表或统计汇总。 |
| `results/tables/14_sessionInfo.txt` | 脚本14运行时的R环境和包版本记录。 |
| `results/tables/14_spatial_domain_and_core_gene_metadata.csv` | 脚本14生成的结果表或统计汇总。 |
| `results/tables/15_Figure4_cell_metadata_with_core_gene_expression.csv` | 脚本15生成的结果表或统计汇总。 |
| `results/tables/15_Figure4_celltype_core_gene_expression_summary.csv` | 脚本15生成的结果表或统计汇总。 |
| `results/tables/15_output_file_check.csv` | 脚本15生成的结果表或统计汇总。 |
| `results/tables/15_sessionInfo.txt` | 脚本15运行时的R环境和包版本记录。 |
| `results/tables/15b_Figure4_cell_metadata_with_core_gene_expression.csv` | 脚本15b生成的结果表或统计汇总。 |
| `results/tables/15b_Figure4_celltype_core_gene_expression_summary.csv` | 脚本15b生成的结果表或统计汇总。 |
| `results/tables/15b_Figure4_cluster_celltype_core_gene_expression_summary.csv` | 脚本15b生成的结果表或统计汇总。 |
| `results/tables/15b_output_file_check.csv` | 脚本15b生成的结果表或统计汇总。 |
| `results/tables/15b_sessionInfo.txt` | 脚本15b运行时的R环境和包版本记录。 |
| `results/tables/16_Figure5_core_gene_high_thresholds.csv` | 脚本16生成的结果表或统计汇总。 |
| `results/tables/16_Figure5_gene_high_status_distribution_summary.csv` | 脚本16生成的结果表或统计汇总。 |
| `results/tables/16_Figure5_status_source_summary.csv` | 脚本16生成的结果表或统计汇总。 |
| `results/tables/16_Figure5_trajectory_metadata_with_core_gene_expression.csv` | 脚本16生成的结果表或统计汇总。 |
| `results/tables/16_Figure5_trajectory_status_core_gene_summary.csv` | 脚本16生成的结果表或统计汇总。 |
| `results/tables/16_available_trajectory_annotation_columns.csv` | 脚本16生成的结果表或统计汇总。 |
| `results/tables/16_output_file_check.csv` | 脚本16生成的结果表或统计汇总。 |
| `results/tables/16_sessionInfo.txt` | 脚本16运行时的R环境和包版本记录。 |

### results/objects/

生成的中间对象目录，主要保存Seurat对象、Monocle轨迹对象或大文件分片。

| 文件 | 作用 |
|---|---|
| `results/objects/02_raw_before_QC_filtering.rds` | 分析流程生成的Seurat/轨迹中间对象，供后续脚本读取。 |
| `results/objects/03_QC_reproduction_candidate.rds` | 分析流程生成的Seurat/轨迹中间对象，供后续脚本读取。 |
| `results/objects/04_standard_Seurat_multi_resolution.rds` | 分析流程生成的Seurat/轨迹中间对象，供后续脚本读取。 |
| `results/objects/05_diagnostic_manual_annotation_diagnostic.rds` | 分析流程生成的Seurat/轨迹中间对象，供后续脚本读取。 |
| `results/objects/05_manual_annotated_before_malignant_call.rds` | 分析流程生成的Seurat/轨迹中间对象，供后续脚本读取。 |
| `results/objects/05_manual_annotated_plot_ready.rds` | 分析流程生成的Seurat/轨迹中间对象，供后续脚本读取。 |
| `results/objects/06_malignant_candidate_diagnostic.rds` | 分析流程生成的Seurat/轨迹中间对象，供后续脚本读取。 |
| `results/objects/07_CopyKAT_malignant_call.rds` | 分析流程生成的Seurat/轨迹中间对象，供后续脚本读取。 |
| `results/objects/08_final_malignant_call.rds` | 分析流程生成的Seurat/轨迹中间对象，供后续脚本读取。 |
| `results/objects/10_malignant_epithelial_state_characterization.rds` | 分析流程生成的Seurat/轨迹中间对象，供后续脚本读取。 |
| `results/objects/10_manual_review_epithelial_core.rds` | 分析流程生成的Seurat/轨迹中间对象，供后续脚本读取。 |
| `results/objects/11_global_trajectory_Seurat.rds.parts/11_global_trajectory_Seurat.rds.part_aa` | 大型RDS对象的分片文件，用于在本地重建完整对象。 |
| `results/objects/11_global_trajectory_Seurat.rds.parts/11_global_trajectory_Seurat.rds.part_ab` | 大型RDS对象的分片文件，用于在本地重建完整对象。 |

## Git LFS说明

仓库中的大文件使用Git LFS管理。对于特别大的、可重新生成的中间对象，可以先不上传，后续根据需要逐个补充。
