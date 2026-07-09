# ============================================================
# 06f_malignant_cell_composition_check.R
#
# 功能：
# 1. 统计Strict_malignant_CopyKAT_aneuploid细胞的sample组成
# 2. 统计其cluster组成
# 3. 检查是否存在明显sample / cluster主导
# 4. 为后续Monocle3拟时序决定主分析与敏感性分析策略
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
  "ggplot2",
  "scales"
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
library(scales)

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
# C. 读取06e最终对象
# ============================================================

input_object_file <- file.path(
  object_dir,
  "06e_GSE215403_final_malignant_call.rds"
)

if (!file.exists(input_object_file)) {
  
  stop(
    paste0(
      "找不到06e对象：\n",
      input_object_file
    )
  )
}

sc <- readRDS(input_object_file)

cluster_column <- "cluster_res_0.2"
status_column <- "malignant_status_final"

required_metadata <- c(
  "sample_id",
  cluster_column,
  "celltype_manual",
  status_column,
  "copykat_prediction_clean"
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
# D. 提取严格malignant细胞metadata
# ============================================================

strict_label <- "Strict_malignant_CopyKAT_aneuploid"

strict_meta <- sc@meta.data %>%
  mutate(
    cell_barcode = rownames(sc@meta.data),
    cluster = as.character(
      .data[[cluster_column]]
    ),
    malignant_status_final = as.character(
      .data[[status_column]]
    ),
    copykat_prediction_clean = as.character(
      copykat_prediction_clean
    )
  ) %>%
  filter(
    malignant_status_final == strict_label
  )

if (nrow(strict_meta) == 0) {
  
  stop(
    "未找到Strict_malignant_CopyKAT_aneuploid细胞。"
  )
}

message(
  "Strict malignant cell number: ",
  nrow(strict_meta)
)

# ============================================================
# E. Sample组成
# ============================================================

strict_by_sample <- strict_meta %>%
  count(
    sample_id,
    name = "strict_malignant_cells"
  ) %>%
  mutate(
    percent_of_all_strict_malignant = round(
      100 * strict_malignant_cells /
        sum(strict_malignant_cells),
      2
    )
  ) %>%
  arrange(
    desc(strict_malignant_cells)
  )

write.csv(
  strict_by_sample,
  file.path(
    table_dir,
    "06f_strict_malignant_cells_by_sample.csv"
  ),
  row.names = FALSE
)

# ============================================================
# F. Cluster组成
# ============================================================

strict_by_cluster <- strict_meta %>%
  count(
    cluster,
    celltype_manual,
    name = "strict_malignant_cells"
  ) %>%
  mutate(
    percent_of_all_strict_malignant = round(
      100 * strict_malignant_cells /
        sum(strict_malignant_cells),
      2
    )
  ) %>%
  arrange(
    suppressWarnings(as.numeric(cluster))
  )

write.csv(
  strict_by_cluster,
  file.path(
    table_dir,
    "06f_strict_malignant_cells_by_cluster.csv"
  ),
  row.names = FALSE
)

# ============================================================
# G. Sample × cluster组成
# ============================================================

strict_by_sample_cluster <- strict_meta %>%
  count(
    sample_id,
    cluster,
    celltype_manual,
    name = "strict_malignant_cells"
  ) %>%
  group_by(
    sample_id
  ) %>%
  mutate(
    percent_within_sample = round(
      100 * strict_malignant_cells /
        sum(strict_malignant_cells),
      2
    )
  ) %>%
  ungroup() %>%
  arrange(
    sample_id,
    suppressWarnings(as.numeric(cluster))
  )

write.csv(
  strict_by_sample_cluster,
  file.path(
    table_dir,
    "06f_strict_malignant_cells_by_sample_cluster.csv"
  ),
  row.names = FALSE
)

# ============================================================
# H. 每个sample的肿瘤细胞量与总体比例
# ============================================================

sample_total_summary <- sc@meta.data %>%
  mutate(
    malignant_status_final = as.character(
      .data[[status_column]]
    )
  ) %>%
  count(
    sample_id,
    name = "total_cells"
  ) %>%
  left_join(
    strict_by_sample,
    by = "sample_id"
  ) %>%
  mutate(
    strict_malignant_cells = replace_na(
      strict_malignant_cells,
      0
    ),
    percent_of_sample_cells_strict_malignant = round(
      100 * strict_malignant_cells /
        total_cells,
      2
    )
  ) %>%
  arrange(
    desc(strict_malignant_cells)
  )

write.csv(
  sample_total_summary,
  file.path(
    table_dir,
    "06f_strict_malignant_cells_per_sample_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# I. 拟时序候选sample建议
# ============================================================
#
# 后续拟时序优先使用strict malignant cell数>=100的sample。
# ============================================================

minimum_cells_for_sample_level_analysis <- 100

pseudotime_sample_recommendation <- sample_total_summary %>%
  mutate(
    recommended_for_sample_aware_pseudotime =
      strict_malignant_cells >=
      minimum_cells_for_sample_level_analysis
  )

write.csv(
  pseudotime_sample_recommendation,
  file.path(
    table_dir,
    "06f_pseudotime_sample_recommendation.csv"
  ),
  row.names = FALSE
)

# ============================================================
# J. 图1：严格malignant细胞按sample组成
# ============================================================

p_sample_composition <- ggplot(
  strict_by_sample,
  aes(
    x = reorder(
      sample_id,
      strict_malignant_cells
    ),
    y = strict_malignant_cells
  )
) +
  geom_col(
    fill = "#D73027",
    width = 0.75
  ) +
  geom_text(
    aes(
      label = paste0(
        percent_of_all_strict_malignant,
        "%"
      )
    ),
    hjust = -0.15,
    size = 3.7
  ) +
  coord_flip() +
  labs(
    title = "Sample composition of strict malignant epithelial cells",
    x = NULL,
    y = "Strict malignant cell number"
  ) +
  expand_limits(
    y = max(
      strict_by_sample$strict_malignant_cells
    ) * 1.18
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
    "06f_strict_malignant_cell_composition_by_sample.pdf"
  ),
  plot = p_sample_composition,
  width = 10,
  height = 7
)

# ============================================================
# K. 图2：严格malignant细胞按cluster组成
# ============================================================

p_cluster_composition <- ggplot(
  strict_by_cluster,
  aes(
    x = factor(
      cluster,
      levels = sort(
        unique(
          cluster
        )
      )
    ),
    y = strict_malignant_cells,
    fill = celltype_manual
  )
) +
  geom_col(
    width = 0.72,
    show.legend = FALSE
  ) +
  geom_text(
    aes(
      label = paste0(
        percent_of_all_strict_malignant,
        "%"
      )
    ),
    vjust = -0.4,
    size = 3.7
  ) +
  labs(
    title = "Cluster composition of strict malignant epithelial cells",
    x = "Cluster",
    y = "Strict malignant cell number"
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
    "06f_strict_malignant_cell_composition_by_cluster.pdf"
  ),
  plot = p_cluster_composition,
  width = 10,
  height = 7
)

# ============================================================
# L. 图3：sample × cluster热图
# ============================================================

all_samples <- sort(
  unique(
    as.character(
      sc$sample_id
    )
  )
)

all_clusters <- c(
  "2",
  "3",
  "4",
  "6",
  "11"
)

strict_heatmap_data <- expand.grid(
  sample_id = all_samples,
  cluster = all_clusters,
  stringsAsFactors = FALSE
) %>%
  left_join(
    strict_by_sample_cluster %>%
      select(
        sample_id,
        cluster,
        strict_malignant_cells,
        percent_within_sample
      ),
    by = c(
      "sample_id",
      "cluster"
    )
  ) %>%
  mutate(
    strict_malignant_cells = replace_na(
      strict_malignant_cells,
      0
    ),
    percent_within_sample = replace_na(
      percent_within_sample,
      0
    )
  )

p_sample_cluster_heatmap <- ggplot(
  strict_heatmap_data,
  aes(
    x = factor(
      cluster,
      levels = all_clusters
    ),
    y = factor(
      sample_id,
      levels = rev(all_samples)
    ),
    fill = strict_malignant_cells
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.35
  ) +
  geom_text(
    aes(
      label = strict_malignant_cells
    ),
    size = 3
  ) +
  scale_fill_gradient(
    low = "#FFF7EC",
    high = "#D7301F",
    name = "Strict malignant\ncell number"
  ) +
  labs(
    title = "Strict malignant cells across samples and clusters",
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
    "06f_strict_malignant_cells_by_sample_cluster.pdf"
  ),
  plot = p_sample_cluster_heatmap,
  width = 10,
  height = 8
)

# ============================================================
# M. 保存轻量metadata对象
# ============================================================

saveRDS(
  strict_meta,
  file.path(
    object_dir,
    "06f_strict_malignant_cell_metadata.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "06f_sessionInfo.txt"
  )
)

# ============================================================
# N. 完成提示
# ============================================================

message("\n============================================================")
message("06f_malignant_cell_composition_check.R运行完成。")
message("")
message("重点查看：")
message("1. results/figures/06f_strict_malignant_cell_composition_by_sample.pdf")
message("2. results/figures/06f_strict_malignant_cell_composition_by_cluster.pdf")
message("3. results/figures/06f_strict_malignant_cells_by_sample_cluster.pdf")
message("4. results/tables/06f_strict_malignant_cells_per_sample_summary.csv")
message("5. results/tables/06f_strict_malignant_cells_by_sample_cluster.csv")
message("6. results/tables/06f_pseudotime_sample_recommendation.csv")
message("============================================================\n")
