# 09_malignant_cell_composition_check.R

# 本脚本功能：
# 1. 读取08 final malignant call对象
# 2. 提取Strict_malignant_CopyKAT_aneuploid细胞
# 3. 统计strict malignant细胞的sample组成
# 4. 统计strict malignant细胞的cluster组成
# 5. 检查是否存在明显sample或cluster主导
# 6. 输出sample、cluster和sample×cluster组成图
# 7. 为后续trajectory分析制定主分析和敏感性分析策略
# 8. 保存strict malignant细胞metadata和session信息

# 本项目专用数据：
# GSE215403
# 12个OSCC单细胞样本：
# OSCC, scB1, scB2, scB5, scB7, scB8,
# scB9, scB10, scB12, scB13, scB14, scB15
#
# 本脚本关注08中定义的严格恶性细胞：
# Strict_malignant_CopyKAT_aneuploid
#
# 该步骤用于判断后续trajectory分析是否可能被单一样本
# 或单一cluster主导。
#
# 通用代码修改位置：
# 1. 换数据集时：
#    修改input_object_file、cluster_column和status_column
#
# 2. 换final malignant标签时：
#    修改strict_label
#
# 3. 换candidate cluster范围时：
#    修改all_clusters
#
# 4. 调整sample-level分析门槛时：
#    修改minimum_cells_for_sample_level_analysis


# ============================================================
# A. 加载包
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
# C. 读取08最终对象
# ============================================================

input_object_file <- file.path(
  object_dir,
  "08_final_malignant_call.rds"
)

if (!file.exists(input_object_file)) {
  
  stop(
    paste0(
      "找不到08对象：\n",
      input_object_file,
      "\n请先运行08_finalize_malignant_call.R"
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
# D. 提取strict malignant细胞metadata
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
    "09_strict_malignant_cells_by_sample.csv"
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
    "09_strict_malignant_cells_by_cluster.csv"
  ),
  row.names = FALSE
)

# ============================================================
# G. Sample×cluster组成
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
    "09_strict_malignant_cells_by_sample_cluster.csv"
  ),
  row.names = FALSE
)

# ============================================================
# H. 每个sample的strict malignant细胞量与总体比例
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
    "09_strict_malignant_cells_per_sample_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# I. Trajectory候选sample建议
# ============================================================

# 后续sample-aware trajectory或敏感性分析优先考虑
# strict malignant cell数>=100的sample。

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
    "09_pseudotime_sample_recommendation.csv"
  ),
  row.names = FALSE
)

# ============================================================
# J. 图1：strict malignant细胞按sample组成
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
    "09_strict_malignant_cell_composition_by_sample.pdf"
  ),
  plot = p_sample_composition,
  width = 10,
  height = 7
)

# ============================================================
# K. 图2：strict malignant细胞按cluster组成
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
    "09_strict_malignant_cell_composition_by_cluster.pdf"
  ),
  plot = p_cluster_composition,
  width = 10,
  height = 7
)

# ============================================================
# L. 图3：sample×cluster热图
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
    "09_strict_malignant_cells_by_sample_cluster.pdf"
  ),
  plot = p_sample_cluster_heatmap,
  width = 10,
  height = 8
)

# ============================================================
# M. 保存轻量metadata对象和环境信息
# ============================================================

saveRDS(
  strict_meta,
  file.path(
    object_dir,
    "09_strict_malignant_cell_metadata.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "09_sessionInfo.txt"
  )
)

# ============================================================
# N. 最终提示
# ============================================================

message("\n============================================================")
message("09_malignant_cell_composition_check.R运行完成。")
message("")
message("已保存对象：")
message("results/objects/09_strict_malignant_cell_metadata.rds")
message("")
message("请重点查看：")
message("1. results/figures/09_strict_malignant_cell_composition_by_sample.pdf")
message("2. results/figures/09_strict_malignant_cell_composition_by_cluster.pdf")
message("3. results/figures/09_strict_malignant_cells_by_sample_cluster.pdf")
message("4. results/tables/09_strict_malignant_cells_per_sample_summary.csv")
message("5. results/tables/09_strict_malignant_cells_by_sample_cluster.csv")
message("6. results/tables/09_pseudotime_sample_recommendation.csv")
message("============================================================\n")