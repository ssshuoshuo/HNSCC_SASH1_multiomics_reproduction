# 05_manual_annotation_and_target_gene_summary.R

# 本脚本功能：
# 1. 读取04标准Seurat分析后的multi-resolution对象
# 2. 基于cluster_res_0.2写入人工major cell-type annotation
# 3. 写入初步epithelial status标签
# 4. 生成cluster和manual cell-type的UMAP图
# 5. 计算SASH1、MYH11、EMP1、COL1A1在cluster层面的表达汇总
# 6. 计算SASH1、MYH11、EMP1、COL1A1在manual cell-type层面的表达汇总
# 7. 输出目标基因DotPlot、VlnPlot和UMAP表达图
# 8. 保存供后续malignant-call和trajectory分析使用的对象

# 本项目专用数据：
# GSE215403
# 12个OSCC单细胞样本：
# OSCC, scB1, scB2, scB5, scB7, scB8,
# scB9, scB10, scB12, scB13, scB14, scB15
#
# 本脚本只写入major cell-type annotation和初步epithelial status。
# 最终malignant/non-malignant判定不在本脚本完成。
# 后续步骤会结合tumor-related epithelial identity和CopyKAT CNV结果
# 进一步定义恶性细胞候选群。
#
# 通用代码修改位置：
# 1. 换数据集时：
#    修改input_object_file和cluster_column
#
# 2. 换聚类分辨率时：
#    修改cluster_column和manual_annotation
#
# 3. 换细胞类型注释时：
#    修改manual_annotation
#
# 4. 换关注基因时：
#    修改target_genes


# ============================================================
# A. 加载包
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 3600)

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tidyr",
  "ggplot2",
  "patchwork",
  "ggrastr"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(Seurat)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ============================================================
# B. 项目路径与文件夹
# ============================================================

project_dir <- getwd()

object_dir <- file.path(
  project_dir,
  "results",
  "objects"
)

table_dir <- file.path(
  project_dir,
  "results",
  "tables"
)

figure_dir <- file.path(
  project_dir,
  "results",
  "figures"
)

dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# C. 读取04对象
# ============================================================

input_object_file <- file.path(
  object_dir,
  "04_standard_Seurat_multi_resolution.rds"
)

if (!file.exists(input_object_file)) {
  
  stop(
    paste0(
      "找不到04对象：\n",
      input_object_file,
      "\n请先运行04_standard_Seurat_PCA_UMAP_resolution_scan.R"
    )
  )
}

sc <- readRDS(input_object_file)

DefaultAssay(sc) <- "RNA"

cluster_column <- "cluster_res_0.2"

if (!cluster_column %in% colnames(sc@meta.data)) {
  stop("找不到cluster_res_0.2。")
}

message("读取完成。")
message("细胞数：", ncol(sc))

# ============================================================
# D. 设置人工annotation map
# ============================================================

# manual_annotation将04得到的Seurat cluster映射为major cell type。
#
# 当前map综合参考：
# 1. cluster-level top markers
# 2. canonical lineage marker DotPlot
# 3. marker FeaturePlot
# 4. cluster-by-sample composition
#
# 注意：
# Tumor_Epithelial、Differentiated_Tumor等标签在本步骤表示
# tumor-related epithelial candidate clusters。
# 它们不是最终malignant-cell calls。
#
# 最终恶性细胞判定会在后续步骤中结合CopyKAT inferred CNV、
# epithelial identity和marker evidence进一步确认。

manual_annotation <- c(
  
  "0"  = "Cytotoxic_T_NKT",
  "1"  = "Macrophage",
  "2"  = "Differentiated_epithelial_tumor_candidate",
  "3"  = "Cycling_epithelial_tumor_candidate",
  "4"  = "Cancer_testis_epithelial_tumor_candidate",
  "5"  = "Treg",
  "6"  = "Epithelial_tumor_candidate",
  "7"  = "Fibroblast_CAF",
  "8"  = "B_cell",
  "9"  = "Plasma_cell",
  "10" = "Blood_endothelial",
  "11" = "Epithelial_tumor_candidate",
  "12" = "Mast_cell",
  "13" = "Plasmacytoid_DC",
  "14" = "Lymphatic_endothelial",
  "15" = "Salivary_epithelial_normal_like"
)

# ============================================================
# E. 检查annotation是否覆盖全部cluster
# ============================================================

observed_clusters <- sort(
  unique(
    as.character(
      sc[[cluster_column, drop = TRUE]]
    )
  ),
  decreasing = FALSE
)

missing_clusters <- setdiff(
  observed_clusters,
  names(manual_annotation)
)

extra_clusters <- setdiff(
  names(manual_annotation),
  observed_clusters
)

if (length(missing_clusters) > 0) {
  
  stop(
    paste0(
      "annotation map缺少cluster：",
      paste(missing_clusters, collapse = ", ")
    )
  )
}

if (length(extra_clusters) > 0) {
  
  warning(
    paste0(
      "annotation map中存在当前对象没有的cluster：",
      paste(extra_clusters, collapse = ", ")
    )
  )
}

# ============================================================
# F. 写入人工annotation
# ============================================================

cell_cluster <- as.character(
  sc[[cluster_column, drop = TRUE]]
)

celltype_manual <- unname(
  manual_annotation[cell_cluster]
)

names(celltype_manual) <- colnames(sc)

sc[["celltype_manual"]] <- celltype_manual

# ============================================================
# G. 写入初步epithelial status标签
# ============================================================

# epithelial_status_initial用于后续选择tumor-related epithelial candidates
# 和normal-like epithelial reference。
#
# malignant_candidate：
# tumor-related epithelial candidate clusters，尚未经过CNV确认。
#
# nonmalignant_epithelial_candidate：
# normal-like/salivary epithelial candidate。
#
# non_epithelial：
# 免疫、基质、血管和其他非上皮细胞。

epithelial_tumor_labels <- c(
  "Differentiated_epithelial_tumor_candidate",
  "Cycling_epithelial_tumor_candidate",
  "Cancer_testis_epithelial_tumor_candidate",
  "Epithelial_tumor_candidate"
)

celltype_vector <- sc$celltype_manual

sc$epithelial_status_initial <- ifelse(
  celltype_vector %in% epithelial_tumor_labels,
  "malignant_candidate",
  ifelse(
    celltype_vector == "Salivary_epithelial_normal_like",
    "nonmalignant_epithelial_candidate",
    "non_epithelial"
  )
)

# ============================================================
# H. 输出annotation表
# ============================================================

annotation_table <- sc@meta.data %>%
  mutate(
    cluster = as.character(
      .data[[cluster_column]]
    )
  ) %>%
  count(
    cluster,
    celltype_manual,
    epithelial_status_initial,
    name = "cell_number"
  ) %>%
  mutate(
    percent_of_all_cells = round(
      100 * cell_number / sum(cell_number),
      2
    )
  ) %>%
  arrange(
    suppressWarnings(as.numeric(cluster))
  )

write.csv(
  annotation_table,
  file.path(
    table_dir,
    "05_manual_annotation_table_resolution_0.2.csv"
  ),
  row.names = FALSE
)

print(annotation_table)

# ============================================================
# I. 输出major cell-type UMAP
# ============================================================

p_umap_cluster <- DimPlot(
  object = sc,
  reduction = "umap_pca",
  group.by = cluster_column,
  label = TRUE,
  repel = TRUE,
  pt.size = 0.15,
  raster = TRUE
) +
  ggtitle("PCA-UMAP by cluster")

p_umap_celltype <- DimPlot(
  object = sc,
  reduction = "umap_pca",
  group.by = "celltype_manual",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.15,
  raster = TRUE
) +
  ggtitle("PCA-UMAP by manual cell-type annotation")

ggsave(
  filename = file.path(
    figure_dir,
    "05_UMAP_cluster_and_manual_celltype.pdf"
  ),
  plot = p_umap_cluster + p_umap_celltype,
  width = 18,
  height = 8
)

ggsave(
  filename = file.path(
    figure_dir,
    "05_UMAP_cluster_and_manual_celltype.png"
  ),
  plot = p_umap_cluster + p_umap_celltype,
  width = 18,
  height = 8,
  dpi = 300
)

# ============================================================
# J. 目标基因cluster/celltype表达汇总
# ============================================================

target_genes <- c(
  "SASH1",
  "MYH11",
  "EMP1",
  "COL1A1"
)

target_genes_present <- intersect(
  target_genes,
  rownames(sc)
)

if (length(target_genes_present) != length(target_genes)) {
  
  warning(
    paste0(
      "部分目标基因不在矩阵中：",
      paste(
        setdiff(target_genes, target_genes_present),
        collapse = ", "
      )
    )
  )
}

target_expression <- FetchData(
  object = sc,
  vars = c(
    target_genes_present,
    cluster_column,
    "celltype_manual",
    "epithelial_status_initial"
  )
)

target_by_cluster <- target_expression %>%
  mutate(
    cluster = as.character(
      .data[[cluster_column]]
    )
  ) %>%
  select(
    -all_of(cluster_column)
  ) %>%
  pivot_longer(
    cols = all_of(target_genes_present),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  group_by(cluster, gene) %>%
  summarise(
    cell_number = n(),
    percent_expressed = round(
      100 * mean(expression > 0),
      2
    ),
    mean_expression = mean(expression),
    median_expression = median(expression),
    .groups = "drop"
  ) %>%
  arrange(
    gene,
    suppressWarnings(as.numeric(cluster))
  )

write.csv(
  target_by_cluster,
  file.path(
    table_dir,
    "05_target_gene_expression_by_cluster.csv"
  ),
  row.names = FALSE
)

target_by_celltype <- target_expression %>%
  pivot_longer(
    cols = all_of(target_genes_present),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  group_by(celltype_manual, epithelial_status_initial, gene) %>%
  summarise(
    cell_number = n(),
    percent_expressed = round(
      100 * mean(expression > 0),
      2
    ),
    mean_expression = mean(expression),
    median_expression = median(expression),
    .groups = "drop"
  ) %>%
  arrange(gene, desc(mean_expression))

write.csv(
  target_by_celltype,
  file.path(
    table_dir,
    "05_target_gene_expression_by_manual_celltype.csv"
  ),
  row.names = FALSE
)

# ============================================================
# K. 目标基因DotPlot
# ============================================================

p_target_dotplot <- DotPlot(
  object = sc,
  features = target_genes_present,
  group.by = "celltype_manual",
  assay = "RNA",
  dot.scale = 7
) +
  RotatedAxis() +
  ggtitle("Core genes across manually annotated cell types")

ggsave(
  filename = file.path(
    figure_dir,
    "05_target_genes_DotPlot_by_manual_celltype.pdf"
  ),
  plot = p_target_dotplot,
  width = 14,
  height = 8
)

# ============================================================
# L. 目标基因VlnPlot
# ============================================================

p_target_violin <- VlnPlot(
  object = sc,
  features = target_genes_present,
  group.by = "celltype_manual",
  pt.size = 0,
  ncol = 2
)

ggsave(
  filename = file.path(
    figure_dir,
    "05_target_genes_VlnPlot_by_manual_celltype.pdf"
  ),
  plot = p_target_violin,
  width = 18,
  height = 12
)

# ============================================================
# M. 目标基因UMAP表达图
# ============================================================

# 全局FeaturePlot容易让低表达信号显得很淡。
# 这里使用min.cutoff="q05"和max.cutoff="q95"改善展示效果。
# 该设置只影响可视化，不改变原始表达矩阵和统计结果。

p_target_feature_adjusted <- FeaturePlot(
  object = sc,
  features = target_genes_present,
  reduction = "umap_pca",
  order = TRUE,
  pt.size = 0.15,
  raster = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q95",
  ncol = 2
)

ggsave(
  filename = file.path(
    figure_dir,
    "05_target_genes_UMAP_quantile_scaled.pdf"
  ),
  plot = p_target_feature_adjusted,
  width = 12,
  height = 10
)

# ============================================================
# N. 保存对象和环境信息
# ============================================================

sc$analysis_stage <- "manual_major_annotation_before_malignant_call"

saveRDS(
  sc,
  file.path(
    object_dir,
    "05_manual_annotated_before_malignant_call.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "05_sessionInfo.txt"
  )
)

# ============================================================
# O. 最终提示
# ============================================================

message("\n============================================================")
message("05_manual_annotation_and_target_gene_summary.R 运行完成。")
message("")
message("已保存对象：")
message("results/objects/05_manual_annotated_before_malignant_call.rds")
message("")
message("请重点查看：")
message("1. results/figures/05_UMAP_cluster_and_manual_celltype.pdf")
message("2. results/figures/05_target_genes_DotPlot_by_manual_celltype.pdf")
message("3. results/figures/05_target_genes_VlnPlot_by_manual_celltype.pdf")
message("4. results/figures/05_target_genes_UMAP_quantile_scaled.pdf")
message("5. results/tables/05_target_gene_expression_by_cluster.csv")
message("6. results/tables/05_target_gene_expression_by_manual_celltype.csv")
message("============================================================\n")