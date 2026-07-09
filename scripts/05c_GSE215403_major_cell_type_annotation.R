# ============================================================
# 05c_manual_annotation_and_target_gene_summary.R
#
# 目标：
# 1. 在 05b 诊断基础上写入人工 major cell-type annotation
# 2. 生成论文式 major-cell-type UMAP
# 3. 计算 SASH1 / MYH11 / EMP1 / COL1A1
#    在 cluster 与 cell type 层面的表达汇总
# 4. 生成 DotPlot、VlnPlot 和目标基因表达表
# 5. 保存供 malignant-call 和 Monocle3 使用的对象
#
# - inferCNV / CopyKAT
# - malignant 最终判定
# - 拟时序
#
# ============================================================

# ============================================================
# 用户配置说明
# ============================================================
# 运行前请检查以下设置：
# 1. project_dir：项目根目录。
# 2. raw_dir：原始数据目录。
# 3. object_dir：RDS对象输出目录。
# 4. table_dir：CSV和TXT结果输出目录。
# 5. figure_dir：PDF图输出目录。
# 6. 输入文件名：若本地文件名不同，请在对应input_file处修改。
# 7. 线程数、内存和运行位置：CopyKAT、Seurat聚类和Monocle3建议在服务器或高内存本地环境运行。
# ============================================================


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
# B. 项目路径
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
# C. 读取 05b 对象
# ============================================================

input_object_file <- file.path(
  object_dir,
  "05b_GSE215403_manual_annotation_diagnostic.rds"
)

if (!file.exists(input_object_file)) {
  
  stop(
    paste0(
      "找不到 05b 对象：\n",
      input_object_file,
      "\n请先运行 05b_manual_annotation_diagnostic.R"
    )
  )
}

sc <- readRDS(input_object_file)

DefaultAssay(sc) <- "RNA"

cluster_column <- "cluster_res_0.2"

if (!cluster_column %in% colnames(sc@meta.data)) {
  stop("找不到 cluster_res_0.2。")
}

message("读取完成。")
message("细胞数：", ncol(sc))

# ============================================================
# D. 人工 annotation map
# ============================================================
#
# 当前 map 基于：
# 1. 05b top markers
# 2. major lineage DotPlot
# 3. FeaturePlot
# 4. cluster-by-sample composition
#
# 通用项目：
# 最常替换的是这一段。
#
# 注意：
# Epithelial_tumor_candidate 只是“候选恶性上皮”。
# 最终 malignant / non-malignant 判定会在 06b 用 CNV
# 和 malignant-marker evidence 再进一步确认。
# ============================================================

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
# E. 检查 annotation 是否覆盖全部 cluster
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
      "annotation map 缺少 cluster：",
      paste(missing_clusters, collapse = ", ")
    )
  )
}

if (length(extra_clusters) > 0) {
  
  warning(
    paste0(
      "annotation map 中存在当前对象没有的 cluster：",
      paste(extra_clusters, collapse = ", ")
    )
  )
}

# ============================================================
# F. 写入人工 annotation
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
# G. 写入 epithelial status 初步标签
# ============================================================
#
# 这里只是为了后续 06b 选择测试细胞与参考细胞。
#
# malignant_candidate：
# 可能恶性，但尚未被 CNV 确认。
#
# nonmalignant_epithelial_candidate：
# 主要用于标识 salivary / normal-like epithelial。
#
# non_epithelial：
# 免疫、基质、血管等。
# ============================================================

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
# H. 生成 annotation 表
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
    "05c_manual_annotation_table_resolution_0.2.csv"
  ),
  row.names = FALSE
)

print(annotation_table)

# ============================================================
# I. Major cell type UMAP
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
    "05c_UMAP_cluster_and_manual_celltype.pdf"
  ),
  plot = p_umap_cluster + p_umap_celltype,
  width = 18,
  height = 8
)

ggsave(
  filename = file.path(
    figure_dir,
    "05c_UMAP_cluster_and_manual_celltype.png"
  ),
  plot = p_umap_cluster + p_umap_celltype,
  width = 18,
  height = 8,
  dpi = 300
)

# ============================================================
# J. 目标基因 cluster / celltype 表达汇总
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
    "05c_target_gene_expression_by_cluster.csv"
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
    "05c_target_gene_expression_by_manual_celltype.csv"
  ),
  row.names = FALSE
)

# ============================================================
# K. 目标基因 DotPlot
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
    "05c_target_genes_DotPlot_by_manual_celltype.pdf"
  ),
  plot = p_target_dotplot,
  width = 14,
  height = 8
)

# ============================================================
# L. 目标基因 VlnPlot
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
    "05c_target_genes_VlnPlot_by_manual_celltype.pdf"
  ),
  plot = p_target_violin,
  width = 18,
  height = 12
)

# ============================================================
# M. 高表达细胞可视化版本
# ============================================================
#
# 全局 FeaturePlot 容易让低表达信号显得很淡。
# 这里通过 min.cutoff = "q05"、max.cutoff = "q95"
# 改善展示，但不改变原始数据和统计结果。
#
# 通用项目：
# 可改为 q01/q99、q10/q90，或指定数值。
# ============================================================

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
    "05c_target_genes_UMAP_quantile_scaled.pdf"
  ),
  plot = p_target_feature_adjusted,
  width = 12,
  height = 10
)

# ============================================================
# N. 保存对象
# ============================================================

sc$analysis_stage <- "manual_major_annotation_before_malignant_call"

saveRDS(
  sc,
  file.path(
    object_dir,
    "05c_GSE215403_manual_annotated_before_malignant_call.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "05c_sessionInfo.txt"
  )
)

# ============================================================
# O. 完成提示
# ============================================================

message("\n============================================================")
message("05c_manual_annotation_and_target_gene_summary.R 运行完成。")
message("")
message("已保存对象：")
message("results/objects/05c_GSE215403_manual_annotated_before_malignant_call.rds")
message("")
message("重点查看：")
message("1. results/figures/05c_UMAP_cluster_and_manual_celltype.pdf")
message("2. results/figures/05c_target_genes_DotPlot_by_manual_celltype.pdf")
message("3. results/figures/05c_target_genes_VlnPlot_by_manual_celltype.pdf")
message("4. results/figures/05c_target_genes_UMAP_quantile_scaled.pdf")
message("5. results/tables/05c_target_gene_expression_by_cluster.csv")
message("6. results/tables/05c_target_gene_expression_by_manual_celltype.csv")
message("============================================================\n")