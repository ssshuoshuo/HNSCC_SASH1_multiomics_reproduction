# ============================================================
# 06e_finalize_malignant_call.R
#
# 功能：
# 1. 读取06d CopyKAT结果对象
# 3. 定义严格malignant epithelial cells
# 4. 输出cluster和sample层面的aneuploid比例
# 5. 输出SASH1等核心基因在最终分组中的描述性比较
# 6. 保存后续Monocle3拟时序所需对象
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
# A. R library与包
# ============================================================

options(timeout = 3600)

user_r_library <- Sys.getenv("R_LIBS_USER")

if (nzchar(user_r_library)) {
  
  dir.create(
    user_r_library,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  .libPaths(
    c(
      user_r_library,
      .libPaths()
    )
  )
}

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tidyr",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  
  stop(
    paste0(
      "缺少R包：",
      paste(missing_packages, collapse = ", ")
    )
  )
}

library(Seurat)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# B. 路径
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
# C. 读取06d对象
# ============================================================

input_object_file <- file.path(
  object_dir,
  "06d_GSE215403_CopyKAT_malignant_call.rds"
)

if (!file.exists(input_object_file)) {
  
  stop(
    paste0(
      "找不到06d对象：\n",
      input_object_file
    )
  )
}

sc <- readRDS(input_object_file)

DefaultAssay(sc) <- "RNA"

cluster_column <- "cluster_res_0.2"

required_metadata <- c(
  cluster_column,
  "sample_id",
  "celltype_manual",
  "copykat_prediction"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(sc@meta.data)
)

if (length(missing_metadata) > 0) {
  
  stop(
    paste0(
      "缺少metadata：",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

# ============================================================
# D. 固定参数
# ============================================================

tumor_candidate_clusters <- c(
  "2",
  "3",
  "4",
  "6",
  "11"
)

salivary_reference_cluster <- "15"

target_genes <- c(
  "SASH1",
  "MYH11",
  "EMP1",
  "COL1A1"
)

# ============================================================
# E. 清理CopyKAT prediction标签
# ============================================================

cluster_vector <- as.character(
  sc[[cluster_column, drop = TRUE]]
)

copykat_prediction_clean <- as.character(
  sc$copykat_prediction
)

copykat_prediction_clean[
  is.na(copykat_prediction_clean) |
    copykat_prediction_clean == ""
] <- "not.run"

valid_copykat_labels <- c(
  "aneuploid",
  "diploid",
  "not.defined",
  "not.run"
)

copykat_prediction_clean[
  !copykat_prediction_clean %in%
    valid_copykat_labels
] <- "not.run"

sc$copykat_prediction_clean <- factor(
  copykat_prediction_clean,
  levels = valid_copykat_labels
)

# ============================================================
# F. 定义最终分析分组
# ============================================================

is_tumor_candidate <- cluster_vector %in%
  tumor_candidate_clusters

is_strict_malignant <- is_tumor_candidate &
  copykat_prediction_clean == "aneuploid"

is_salivary_reference <- cluster_vector ==
  salivary_reference_cluster

malignant_status_final <- rep(
  "Other_cells",
  ncol(sc)
)

malignant_status_final[
  is_tumor_candidate
] <- "Tumor_candidate_not_confirmed"

malignant_status_final[
  is_strict_malignant
] <- "Strict_malignant_CopyKAT_aneuploid"

malignant_status_final[
  is_salivary_reference
] <- "Salivary_epithelial_reference"

sc$malignant_status_final <- factor(
  malignant_status_final,
  levels = c(
    "Strict_malignant_CopyKAT_aneuploid",
    "Tumor_candidate_not_confirmed",
    "Salivary_epithelial_reference",
    "Other_cells"
  )
)

sc$is_strict_malignant_copykat <- is_strict_malignant

# ============================================================
# G. 整理metadata
# ============================================================

meta <- sc@meta.data %>%
  mutate(
    cell_barcode = rownames(sc@meta.data),
    cluster = as.character(
      .data[[cluster_column]]
    ),
    copykat_prediction_clean = as.character(
      copykat_prediction_clean
    ),
    malignant_status_final = as.character(
      malignant_status_final
    )
  )

# ============================================================
# H. 总体汇总
# ============================================================

final_status_summary <- meta %>%
  count(
    malignant_status_final,
    name = "cell_number"
  ) %>%
  mutate(
    percent_of_all_cells = round(
      100 * cell_number / sum(cell_number),
      2
    )
  )

write.csv(
  final_status_summary,
  file.path(
    table_dir,
    "06e_final_malignant_status_summary.csv"
  ),
  row.names = FALSE
)

candidate_prediction_overall <- meta %>%
  filter(
    cluster %in% tumor_candidate_clusters
  ) %>%
  count(
    cluster,
    celltype_manual,
    copykat_prediction_clean,
    name = "cell_number"
  ) %>%
  group_by(
    cluster
  ) %>%
  mutate(
    percent_within_cluster = round(
      100 * cell_number / sum(cell_number),
      2
    )
  ) %>%
  ungroup() %>%
  arrange(
    suppressWarnings(as.numeric(cluster)),
    desc(cell_number)
  )

write.csv(
  candidate_prediction_overall,
  file.path(
    table_dir,
    "06e_CopyKAT_prediction_overall_clean.csv"
  ),
  row.names = FALSE
)

# ============================================================
# I. 按cluster计算aneuploid比例
# ============================================================

aneuploid_by_cluster <- meta %>%
  filter(
    cluster %in% tumor_candidate_clusters
  ) %>%
  group_by(
    cluster,
    celltype_manual
  ) %>%
  summarise(
    total_candidate_cells = n(),
    aneuploid_cells = sum(
      copykat_prediction_clean == "aneuploid"
    ),
    aneuploid_percent = round(
      100 * aneuploid_cells /
        total_candidate_cells,
      2
    ),
    .groups = "drop"
  ) %>%
  arrange(
    suppressWarnings(as.numeric(cluster))
  )

write.csv(
  aneuploid_by_cluster,
  file.path(
    table_dir,
    "06e_aneuploid_fraction_by_cluster.csv"
  ),
  row.names = FALSE
)

# ============================================================
# J. 按sample和cluster计算aneuploid比例
# ============================================================

aneuploid_by_sample_cluster <- meta %>%
  filter(
    cluster %in% tumor_candidate_clusters
  ) %>%
  group_by(
    sample_id,
    cluster,
    celltype_manual
  ) %>%
  summarise(
    total_candidate_cells = n(),
    aneuploid_cells = sum(
      copykat_prediction_clean == "aneuploid"
    ),
    aneuploid_percent = round(
      100 * aneuploid_cells /
        total_candidate_cells,
      2
    ),
    .groups = "drop"
  ) %>%
  arrange(
    sample_id,
    suppressWarnings(as.numeric(cluster))
  )

write.csv(
  aneuploid_by_sample_cluster,
  file.path(
    table_dir,
    "06e_aneuploid_fraction_by_sample_cluster.csv"
  ),
  row.names = FALSE
)

# ============================================================
# K. UMAP：最终malignant status
# ============================================================

malignant_status_colors <- c(
  "Strict_malignant_CopyKAT_aneuploid" = "#D73027",
  "Tumor_candidate_not_confirmed" = "#FDAE61",
  "Salivary_epithelial_reference" = "#7570B3",
  "Other_cells" = "#D9D9D9"
)

p_final_status_umap <- DimPlot(
  object = sc,
  reduction = "umap_pca",
  group.by = "malignant_status_final",
  cols = malignant_status_colors,
  pt.size = 0.22,
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  ggtitle(
    "CopyKAT-supported malignant epithelial cells"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(
    figure_dir,
    "06e_final_malignant_status_UMAP.pdf"
  ),
  plot = p_final_status_umap,
  width = 12,
  height = 8
)

# ============================================================
# L. 柱状图：每个cluster的aneuploid比例
# ============================================================

p_aneuploid_cluster <- ggplot(
  aneuploid_by_cluster,
  aes(
    x = factor(
      cluster,
      levels = tumor_candidate_clusters
    ),
    y = aneuploid_percent,
    fill = factor(
      cluster,
      levels = tumor_candidate_clusters
    )
  )
) +
  geom_col(
    width = 0.72,
    show.legend = FALSE
  ) +
  geom_text(
    aes(
      label = paste0(
        aneuploid_percent,
        "%"
      )
    ),
    vjust = -0.4,
    size = 4
  ) +
  labs(
    title = "CopyKAT aneuploid fraction across tumor epithelial candidates",
    x = "Cluster",
    y = "Aneuploid cells (%)"
  ) +
  ylim(
    0,
    110
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    figure_dir,
    "06e_CopyKAT_aneuploid_fraction_by_cluster.pdf"
  ),
  plot = p_aneuploid_cluster,
  width = 10,
  height = 7
)

# ============================================================
# M. Heatmap：sample内各cluster的aneuploid比例
# ============================================================

p_aneuploid_sample_cluster <- ggplot(
  aneuploid_by_sample_cluster,
  aes(
    x = factor(
      cluster,
      levels = tumor_candidate_clusters
    ),
    y = sample_id,
    fill = aneuploid_percent
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.35
  ) +
  geom_text(
    aes(
      label = paste0(
        aneuploid_percent,
        "%"
      )
    ),
    size = 3
  ) +
  scale_fill_gradient(
    low = "#FFF7EC",
    high = "#D7301F",
    limits = c(0, 100),
    name = "Aneuploid\ncells (%)"
  ) +
  labs(
    title = "Sample-wise CopyKAT aneuploid fraction",
    x = "Tumor epithelial candidate cluster",
    y = "Sample"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    figure_dir,
    "06e_CopyKAT_aneuploid_fraction_by_sample_cluster.pdf"
  ),
  plot = p_aneuploid_sample_cluster,
  width = 10,
  height = 8
)

# ============================================================
# N. 核心基因表达：最终细胞分组
# ============================================================

comparison_group <- rep(
  NA_character_,
  ncol(sc)
)

comparison_group[
  is_strict_malignant
] <- "Strict_malignant"

comparison_group[
  is_tumor_candidate &
    !is_strict_malignant
] <- "Candidate_not_confirmed"

comparison_group[
  is_salivary_reference
] <- "Salivary_reference"

sc$core_gene_comparison_group <- factor(
  comparison_group,
  levels = c(
    "Strict_malignant",
    "Candidate_not_confirmed",
    "Salivary_reference"
  )
)

target_genes_present <- intersect(
  target_genes,
  rownames(sc)
)

expression_data <- FetchData(
  object = sc,
  vars = c(
    "sample_id",
    "core_gene_comparison_group",
    target_genes_present
  )
)

expression_long <- expression_data %>%
  filter(
    !is.na(core_gene_comparison_group)
  ) %>%
  pivot_longer(
    cols = all_of(target_genes_present),
    names_to = "gene",
    values_to = "expression"
  )

target_expression_summary <- expression_long %>%
  group_by(
    core_gene_comparison_group,
    gene
  ) %>%
  summarise(
    cell_number = n(),
    percent_expressed = round(
      100 * mean(expression > 0),
      2
    ),
    mean_expression = mean(expression),
    median_expression = median(expression),
    .groups = "drop"
  )

write.csv(
  target_expression_summary,
  file.path(
    table_dir,
    "06e_core_gene_expression_by_final_status.csv"
  ),
  row.names = FALSE
)

target_expression_pseudobulk <- expression_long %>%
  group_by(
    sample_id,
    core_gene_comparison_group,
    gene
  ) %>%
  summarise(
    cell_number = n(),
    percent_expressed = round(
      100 * mean(expression > 0),
      2
    ),
    mean_expression = mean(expression),
    median_expression = median(expression),
    .groups = "drop"
  )

write.csv(
  target_expression_pseudobulk,
  file.path(
    table_dir,
    "06e_core_gene_pseudobulk_by_sample.csv"
  ),
  row.names = FALSE
)

p_target_expression <- ggplot(
  target_expression_summary,
  aes(
    x = gene,
    y = core_gene_comparison_group,
    size = percent_expressed,
    color = mean_expression
  )
) +
  geom_point() +
  scale_size(
    range = c(1, 12),
    name = "Percent expressed"
  ) +
  scale_color_gradient(
    low = "#F7FBFF",
    high = "#08306B",
    name = "Mean expression"
  ) +
  labs(
    title = "Core genes across final epithelial-cell groups",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    figure_dir,
    "06e_core_gene_expression_by_final_status.pdf"
  ),
  plot = p_target_expression,
  width = 10,
  height = 6
)

# ============================================================
# O. 保存对象与session信息
# ============================================================

sc$analysis_stage <- "final_CopyKAT_supported_malignant_call"

saveRDS(
  sc,
  file.path(
    object_dir,
    "06e_GSE215403_final_malignant_call.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "06e_sessionInfo.txt"
  )
)

# ============================================================
# P. 完成提示
# ============================================================

message("\n============================================================")
message("06e_finalize_malignant_call.R运行完成。")
message("")
message("重点查看：")
message("1. results/figures/06e_final_malignant_status_UMAP.pdf")
message("2. results/figures/06e_CopyKAT_aneuploid_fraction_by_cluster.pdf")
message("3. results/figures/06e_CopyKAT_aneuploid_fraction_by_sample_cluster.pdf")
message("4. results/figures/06e_core_gene_expression_by_final_status.pdf")
message("5. results/tables/06e_aneuploid_fraction_by_cluster.csv")
message("6. results/tables/06e_aneuploid_fraction_by_sample_cluster.csv")
message("7. results/tables/06e_core_gene_expression_by_final_status.csv")
message("============================================================\n")
