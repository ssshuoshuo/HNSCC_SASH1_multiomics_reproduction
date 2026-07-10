# 13_spatial_SASH1_COL1A1_neighborhood_analysis.R

# 本脚本功能：
# 1. 读取12号脚本生成的空间Seurat对象
# 2. 提取SASH1、COL1A1、EMP1、MYH11空间表达
# 3. 用每个样本内分位数定义SASH1-high和COL1A1-high spot
# 4. 绘制SASH1-high与COL1A1-high空间叠加图
# 5. 计算SASH1-high与COL1A1-high的共定位比例
# 6. 计算SASH1-high到最近COL1A1-high spot的空间距离
# 7. 通过随机置换COL1A1-high标签进行空间邻近性检验
# 8. 输出论文式空间验证统计图和表格

# 本项目专用数据：
# GSE252265空间转录组数据
# 输入对象来自：
# results/objects/12_spatial_tissue_spots_Seurat.rds
#
# 注意：
# 当前12号结果中spatial_sample_id为All_spots。
# 因此本脚本先完成聚合空间层面的验证。
# 后续若能从RAW文件中恢复8个切片标签，可直接按spatial_sample_id分组重跑。

# 通用代码修改位置：
# 1. 修改target_gene_1和target_gene_2可替换空间验证基因
# 2. 修改high_quantile可改变high spot阈值
# 3. 修改permutation_number可改变置换检验次数


# ============================================================
# A. 加载包
# ============================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "ggplot2",
  "patchwork"
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
      paste(
        missing_packages,
        collapse = ", "
      ),
      "\n请先安装缺少包后再运行。"
    )
  )
}

library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(patchwork)

# ============================================================
# B. 路径与参数
# ============================================================

project_dir <- normalizePath(
  "~/Desktop/HNSCC_SASH1_reproduction"
)

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

dir.create(
  object_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

input_object_file <- file.path(
  object_dir,
  "12_spatial_tissue_spots_Seurat.rds"
)

target_gene_1 <- "SASH1"
target_gene_2 <- "COL1A1"

support_genes <- c(
  "EMP1",
  "MYH11"
)

target_genes <- unique(
  c(
    target_gene_1,
    target_gene_2,
    support_genes
  )
)

high_quantile <- 0.75

permutation_number <- 1000

set.seed(
  20260710
)

# ============================================================
# C. 读取空间对象
# ============================================================

if (!file.exists(input_object_file)) {
  stop(
    paste0(
      "未找到输入对象：",
      input_object_file,
      "\n请先成功运行12_spatial_download_QC_gene_maps.R。"
    )
  )
}

spatial_object <- readRDS(
  input_object_file
)

DefaultAssay(spatial_object) <- "SpatialRNA"

required_metadata_columns <- c(
  "spatial_x",
  "spatial_y",
  "spatial_sample_id",
  "in_tissue"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(
    spatial_object@meta.data
  )
)

if (length(missing_metadata_columns) > 0) {
  stop(
    paste0(
      "空间对象缺少必要metadata列：",
      paste(
        missing_metadata_columns,
        collapse = ", "
      )
    )
  )
}

target_genes_found <- intersect(
  target_genes,
  rownames(spatial_object)
)

if (
  !all(
    c(
      target_gene_1,
      target_gene_2
    ) %in%
    target_genes_found
  )
) {
  stop(
    paste0(
      "空间对象中未同时找到",
      target_gene_1,
      "和",
      target_gene_2,
      "。"
    )
  )
}

message(
  "空间对象spot数：",
  ncol(spatial_object)
)

message(
  "空间对象样本标签：",
  paste(
    unique(
      spatial_object$spatial_sample_id
    ),
    collapse = ", "
  )
)

message(
  "找到目标基因：",
  paste(
    target_genes_found,
    collapse = ", "
  )
)

# ============================================================
# D. 提取表达和坐标
# ============================================================

normalized_expression <- LayerData(
  object = spatial_object,
  assay = "SpatialRNA",
  layer = "data"
)

spatial_plot_data <- spatial_object@meta.data

spatial_plot_data$barcode <- rownames(
  spatial_plot_data
)

for (current_gene in target_genes_found) {
  
  spatial_plot_data[
    ,
    current_gene
  ] <- as.numeric(
    normalized_expression[
      current_gene,
      spatial_plot_data$barcode
    ]
  )
}

write.csv(
  spatial_plot_data,
  file.path(
    table_dir,
    "13_spatial_gene_expression_with_coordinates.csv"
  ),
  row.names = FALSE
)

# ============================================================
# E. 样本内分位数定义high spot
# ============================================================

calculate_positive_quantile_threshold <- function(
    expression_vector,
    quantile_value
) {
  
  positive_values <- expression_vector[
    expression_vector > 0
  ]
  
  if (length(positive_values) < 10) {
    return(
      Inf
    )
  }
  
  as.numeric(
    stats::quantile(
      positive_values,
      probs = quantile_value,
      na.rm = TRUE
    )
  )
}

threshold_table <- spatial_plot_data %>%
  dplyr::group_by(
    spatial_sample_id
  ) %>%
  dplyr::summarise(
    SASH1_positive_spot_number = sum(
      .data[[target_gene_1]] > 0,
      na.rm = TRUE
    ),
    COL1A1_positive_spot_number = sum(
      .data[[target_gene_2]] > 0,
      na.rm = TRUE
    ),
    SASH1_high_threshold = calculate_positive_quantile_threshold(
      .data[[target_gene_1]],
      high_quantile
    ),
    COL1A1_high_threshold = calculate_positive_quantile_threshold(
      .data[[target_gene_2]],
      high_quantile
    ),
    .groups = "drop"
  )

spatial_plot_data <- spatial_plot_data %>%
  dplyr::left_join(
    threshold_table,
    by = "spatial_sample_id"
  )

spatial_plot_data$SASH1_high <- spatial_plot_data[
  ,
  target_gene_1
] >= spatial_plot_data$SASH1_high_threshold

spatial_plot_data$COL1A1_high <- spatial_plot_data[
  ,
  target_gene_2
] >= spatial_plot_data$COL1A1_high_threshold

spatial_plot_data$SASH1_high[
  is.na(
    spatial_plot_data$SASH1_high
  )
] <- FALSE

spatial_plot_data$COL1A1_high[
  is.na(
    spatial_plot_data$COL1A1_high
  )
] <- FALSE

spatial_plot_data$spatial_high_status <- "Other"

spatial_plot_data$spatial_high_status[
  spatial_plot_data$SASH1_high &
    !spatial_plot_data$COL1A1_high
] <- "SASH1-high only"

spatial_plot_data$spatial_high_status[
  !spatial_plot_data$SASH1_high &
    spatial_plot_data$COL1A1_high
] <- "COL1A1-high only"

spatial_plot_data$spatial_high_status[
  spatial_plot_data$SASH1_high &
    spatial_plot_data$COL1A1_high
] <- "SASH1-high and COL1A1-high"

spatial_plot_data$spatial_high_status <- factor(
  spatial_plot_data$spatial_high_status,
  levels = c(
    "Other",
    "SASH1-high only",
    "COL1A1-high only",
    "SASH1-high and COL1A1-high"
  )
)

write.csv(
  threshold_table,
  file.path(
    table_dir,
    "13_SASH1_COL1A1_high_thresholds_by_sample.csv"
  ),
  row.names = FALSE
)

write.csv(
  spatial_plot_data,
  file.path(
    table_dir,
    "13_spatial_high_status_metadata.csv"
  ),
  row.names = FALSE
)

# ============================================================
# F. 共定位统计
# ============================================================

colocalization_summary <- spatial_plot_data %>%
  dplyr::group_by(
    spatial_sample_id
  ) %>%
  dplyr::summarise(
    tissue_spot_number = dplyr::n(),
    SASH1_detected_spot_number = sum(
      .data[[target_gene_1]] > 0,
      na.rm = TRUE
    ),
    COL1A1_detected_spot_number = sum(
      .data[[target_gene_2]] > 0,
      na.rm = TRUE
    ),
    SASH1_high_spot_number = sum(
      SASH1_high,
      na.rm = TRUE
    ),
    COL1A1_high_spot_number = sum(
      COL1A1_high,
      na.rm = TRUE
    ),
    SASH1_COL1A1_high_overlap_spot_number = sum(
      SASH1_high &
        COL1A1_high,
      na.rm = TRUE
    ),
    SASH1_high_overlap_COL1A1_high_percent = round(
      100 *
        SASH1_COL1A1_high_overlap_spot_number /
        SASH1_high_spot_number,
      2
    ),
    COL1A1_high_overlap_SASH1_high_percent = round(
      100 *
        SASH1_COL1A1_high_overlap_spot_number /
        COL1A1_high_spot_number,
      2
    ),
    expected_overlap_if_independent = round(
      SASH1_high_spot_number *
        COL1A1_high_spot_number /
        tissue_spot_number,
      2
    ),
    observed_to_expected_overlap_ratio = round(
      SASH1_COL1A1_high_overlap_spot_number /
        expected_overlap_if_independent,
      4
    ),
    .groups = "drop"
  )

write.csv(
  colocalization_summary,
  file.path(
    table_dir,
    "13_SASH1_COL1A1_high_colocalization_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# G. 最近邻距离函数
# ============================================================

calculate_nearest_distance <- function(
    query_coordinates,
    target_coordinates
) {
  
  if (
    nrow(query_coordinates) == 0 ||
    nrow(target_coordinates) == 0
  ) {
    return(
      numeric(0)
    )
  }
  
  nearest_distances <- vapply(
    seq_len(
      nrow(query_coordinates)
    ),
    function(current_index) {
      
      current_x <- query_coordinates$spatial_x[
        current_index
      ]
      
      current_y <- query_coordinates$spatial_y[
        current_index
      ]
      
      min(
        sqrt(
          (
            target_coordinates$spatial_x -
              current_x
          )^2 +
            (
              target_coordinates$spatial_y -
                current_y
            )^2
        ),
        na.rm = TRUE
      )
    },
    numeric(1)
  )
  
  nearest_distances
}

# ============================================================
# H. 计算真实最近邻距离
# ============================================================

distance_result_list <- list()

for (current_sample in unique(spatial_plot_data$spatial_sample_id)) {
  
  current_data <- spatial_plot_data[
    spatial_plot_data$spatial_sample_id ==
      current_sample,
    ,
    drop = FALSE
  ]
  
  SASH1_high_data <- current_data[
    current_data$SASH1_high,
    ,
    drop = FALSE
  ]
  
  COL1A1_high_data <- current_data[
    current_data$COL1A1_high,
    ,
    drop = FALSE
  ]
  
  nearest_distances <- calculate_nearest_distance(
    query_coordinates = SASH1_high_data,
    target_coordinates = COL1A1_high_data
  )
  
  distance_result_list[[current_sample]] <- data.frame(
    spatial_sample_id = current_sample,
    barcode = SASH1_high_data$barcode,
    nearest_COL1A1_high_distance = nearest_distances,
    stringsAsFactors = FALSE
  )

nearest_distance_table <- dplyr::bind_rows(
  distance_result_list
)

write.csv(
  nearest_distance_table,
  file.path(
    table_dir,
    "13_SASH1_high_to_nearest_COL1A1_high_distance.csv"
  ),
  row.names = FALSE
)

nearest_distance_summary <- nearest_distance_table %>%
  dplyr::group_by(
    spatial_sample_id
  ) %>%
  dplyr::summarise(
    SASH1_high_spot_number = dplyr::n(),
    median_nearest_COL1A1_high_distance = round(
      median(
        nearest_COL1A1_high_distance,
        na.rm = TRUE
      ),
      3
    ),
    mean_nearest_COL1A1_high_distance = round(
      mean(
        nearest_COL1A1_high_distance,
        na.rm = TRUE
      ),
      3
    ),
    .groups = "drop"
  )

write.csv(
  nearest_distance_summary,
  file.path(
    table_dir,
    "13_SASH1_high_to_nearest_COL1A1_high_distance_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# I. 随机置换COL1A1-high标签检验空间邻近性
# ============================================================

permutation_result_list <- list()

for (current_sample in unique(spatial_plot_data$spatial_sample_id)) {
  
  current_data <- spatial_plot_data[
    spatial_plot_data$spatial_sample_id ==
      current_sample,
    ,
    drop = FALSE
  ]
  
  SASH1_high_data <- current_data[
    current_data$SASH1_high,
    ,
    drop = FALSE
  ]
  
  COL1A1_high_number <- sum(
    current_data$COL1A1_high,
    na.rm = TRUE
  )
  
  if (
    nrow(SASH1_high_data) == 0 ||
    COL1A1_high_number == 0
  ) {
    next
  }
  
  observed_distance <- median(
    calculate_nearest_distance(
      query_coordinates = SASH1_high_data,
      target_coordinates = current_data[
        current_data$COL1A1_high,
        ,
        drop = FALSE
      ]
    ),
    na.rm = TRUE
  )
  
  permuted_distances <- vapply(
    seq_len(
      permutation_number
    ),
    function(current_permutation) {
      
      permuted_COL1A1_high <- rep(
        FALSE,
        nrow(current_data)
      )
      
      permuted_high_index <- sample(
        seq_len(
          nrow(current_data)
        ),
        size = COL1A1_high_number,
        replace = FALSE
      )
      
      permuted_COL1A1_high[
        permuted_high_index
      ] <- TRUE
      
      permuted_target_data <- current_data[
        permuted_COL1A1_high,
        ,
        drop = FALSE
      ]
      
      median(
        calculate_nearest_distance(
          query_coordinates = SASH1_high_data,
          target_coordinates = permuted_target_data
        ),
        na.rm = TRUE
      )
    },
    numeric(1)
  )
  
  empirical_p_closer <- (
    sum(
      permuted_distances <= observed_distance,
      na.rm = TRUE
    ) + 1
  ) / (
    permutation_number + 1
  )
  
  empirical_p_farther <- (
    sum(
      permuted_distances >= observed_distance,
      na.rm = TRUE
    ) + 1
  ) / (
    permutation_number + 1
  )
  
  permutation_result_list[[current_sample]] <- data.frame(
    spatial_sample_id = current_sample,
    observed_median_distance = observed_distance,
    permutation_median_distance_mean = mean(
      permuted_distances,
      na.rm = TRUE
    ),
    permutation_median_distance_median = median(
      permuted_distances,
      na.rm = TRUE
    ),
    empirical_p_for_closer_than_random = empirical_p_closer,
    empirical_p_for_farther_than_random = empirical_p_farther,
    permutation_number = permutation_number,
    stringsAsFactors = FALSE
  )
  
  permutation_distribution_table <- data.frame(
    spatial_sample_id = current_sample,
    permutation_id = seq_len(
      permutation_number
    ),
    permuted_median_distance = permuted_distances,
    observed_median_distance = observed_distance,
    stringsAsFactors = FALSE
  )
  
  write.csv(
    permutation_distribution_table,
    file.path(
      table_dir,
      paste0(
        "13_permutation_distribution_",
        current_sample,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}

permutation_summary <- dplyr::bind_rows(
  permutation_result_list
)

write.csv(
  permutation_summary,
  file.path(
    table_dir,
    "13_SASH1_COL1A1_high_neighborhood_permutation_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# J. 空间叠加图
# ============================================================

p_high_overlay <- ggplot(
  spatial_plot_data,
  aes(
    x = spatial_x,
    y = spatial_y,
    color = spatial_high_status
  )
) +
  geom_point(
    size = 0.7,
    alpha = 0.9
  ) +
  scale_y_reverse() +
  coord_equal() +
  facet_wrap(
    ~spatial_sample_id,
    scales = "fixed"
  ) +
  labs(
    title = "Spatial Overlay of SASH1-high and COL1A1-high Spots",
    x = NULL,
    y = NULL,
    color = "Spatial status"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    strip.text = element_text(
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    figure_dir,
    "13_SASH1_high_COL1A1_high_spatial_overlay.pdf"
  ),
  plot = p_high_overlay,
  width = 12,
  height = 9
)

# ============================================================
# K. 连续表达组合图
# ============================================================

make_single_gene_plot <- function(
    plot_data,
    gene_name
) {
  
  ggplot(
    plot_data,
    aes(
      x = spatial_x,
      y = spatial_y,
      color = .data[[gene_name]]
    )
  ) +
    geom_point(
      size = 0.65,
      alpha = 0.85
    ) +
    scale_y_reverse() +
    coord_equal() +
    facet_wrap(
      ~spatial_sample_id,
      scales = "fixed"
    ) +
    labs(
      title = paste0(
        gene_name,
        " spatial expression"
      ),
      x = NULL,
      y = NULL,
      color = gene_name
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      strip.text = element_text(
        face = "bold"
      )
    )
}

p_SASH1 <- make_single_gene_plot(
  spatial_plot_data,
  "SASH1"
)

p_COL1A1 <- make_single_gene_plot(
  spatial_plot_data,
  "COL1A1"
)

p_EMP1 <- make_single_gene_plot(
  spatial_plot_data,
  "EMP1"
)

p_MYH11 <- make_single_gene_plot(
  spatial_plot_data,
  "MYH11"
)

p_gene_panel <- (
  p_SASH1 +
    p_COL1A1
) / (
  p_EMP1 +
    p_MYH11
) +
  patchwork::plot_annotation(
    title = "Spatial Expression of SASH1, COL1A1, EMP1 and MYH11"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "13_core_gene_spatial_expression_panel.pdf"
  ),
  plot = p_gene_panel,
  width = 16,
  height = 13
)

# ============================================================
# L. 共定位柱状图
# ============================================================

colocalization_long <- data.frame(
  spatial_sample_id = rep(
    colocalization_summary$spatial_sample_id,
    times = 3
  ),
  metric = rep(
    c(
      "SASH1-high",
      "COL1A1-high",
      "Overlap"
    ),
    each = nrow(
      colocalization_summary
    )
  ),
  spot_number = c(
    colocalization_summary$SASH1_high_spot_number,
    colocalization_summary$COL1A1_high_spot_number,
    colocalization_summary$SASH1_COL1A1_high_overlap_spot_number
  ),
  stringsAsFactors = FALSE
)

p_colocalization_bar <- ggplot(
  colocalization_long,
  aes(
    x = metric,
    y = spot_number,
    fill = metric
  )
) +
  geom_col(
    width = 0.7
  ) +
  facet_wrap(
    ~spatial_sample_id,
    scales = "free_y"
  ) +
  labs(
    title = "SASH1-high and COL1A1-high Spot Counts",
    x = NULL,
    y = "Spot number"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "13_SASH1_COL1A1_high_colocalization_barplot.pdf"
  ),
  plot = p_colocalization_bar,
  width = 8,
  height = 5
)

# ============================================================
# M. 置换检验图
# ============================================================

if (
  exists(
    "permutation_distribution_table"
  ) &&
  nrow(
    permutation_summary
  ) > 0
) {
  
  all_permutation_files <- list.files(
    table_dir,
    pattern = "^13_permutation_distribution_.*\\.csv$",
    full.names = TRUE
  )
  
  all_permutation_data <- dplyr::bind_rows(
    lapply(
      all_permutation_files,
      read.csv,
      stringsAsFactors = FALSE
    )
  )
  
  p_permutation <- ggplot(
    all_permutation_data,
    aes(
      x = permuted_median_distance
    )
  ) +
    geom_histogram(
      bins = 40
    ) +
    geom_vline(
      aes(
        xintercept = observed_median_distance
      ),
      linewidth = 0.8
    ) +
    facet_wrap(
      ~spatial_sample_id,
      scales = "free"
    ) +
    labs(
      title = "Permutation Test for SASH1-high to COL1A1-high Distance",
      x = "Permuted median nearest distance",
      y = "Permutation count"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
  
  ggsave(
    filename = file.path(
      figure_dir,
      "13_SASH1_COL1A1_high_neighborhood_permutation_test.pdf"
    ),
    plot = p_permutation,
    width = 9,
    height = 6
  )
}

# ============================================================
# N. 保存对象和运行信息
# ============================================================

spatial_object$SASH1_high <- spatial_plot_data[
  colnames(spatial_object),
  "SASH1_high"
]

spatial_object$COL1A1_high <- spatial_plot_data[
  colnames(spatial_object),
  "COL1A1_high"
]

spatial_object$spatial_high_status <- spatial_plot_data[
  colnames(spatial_object),
  "spatial_high_status"
]

saveRDS(
  spatial_object,
  file.path(
    object_dir,
    "13_spatial_SASH1_COL1A1_high_annotated_Seurat.rds"
  ),
  compress = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    table_dir,
    "13_sessionInfo.txt"
  )
)

# ============================================================
# O. 输出检查
# ============================================================

required_output_files <- c(
  file.path(
    object_dir,
    "13_spatial_SASH1_COL1A1_high_annotated_Seurat.rds"
  ),
  file.path(
    table_dir,
    "13_spatial_gene_expression_with_coordinates.csv"
  ),
  file.path(
    table_dir,
    "13_SASH1_COL1A1_high_thresholds_by_sample.csv"
  ),
  file.path(
    table_dir,
    "13_spatial_high_status_metadata.csv"
  ),
  file.path(
    table_dir,
    "13_SASH1_COL1A1_high_colocalization_summary.csv"
  ),
  file.path(
    table_dir,
    "13_SASH1_high_to_nearest_COL1A1_high_distance_summary.csv"
  ),
  file.path(
    table_dir,
    "13_SASH1_COL1A1_high_neighborhood_permutation_summary.csv"
  ),
  file.path(
    figure_dir,
    "13_SASH1_high_COL1A1_high_spatial_overlay.pdf"
  ),
  file.path(
    figure_dir,
    "13_core_gene_spatial_expression_panel.pdf"
  ),
  file.path(
    figure_dir,
    "13_SASH1_COL1A1_high_colocalization_barplot.pdf"
  ),
  file.path(
    figure_dir,
    "13_SASH1_COL1A1_high_neighborhood_permutation_test.pdf"
  )
)

output_status <- data.frame(
  file = required_output_files,
  exists = file.exists(
    required_output_files
  ),
  stringsAsFactors = FALSE
)

write.csv(
  output_status,
  file.path(
    table_dir,
    "13_output_file_check.csv"
  ),
  row.names = FALSE
)

print(
  output_status
)

message("\n============================================================")
message("13 SASH1与COL1A1空间邻近/共定位分析完成。")
message("")
message("重点查看：")
message("1. 13_SASH1_high_COL1A1_high_spatial_overlay.pdf")
message("2. 13_SASH1_COL1A1_high_colocalization_summary.csv")
message("3. 13_SASH1_high_to_nearest_COL1A1_high_distance_summary.csv")
message("4. 13_SASH1_COL1A1_high_neighborhood_permutation_summary.csv")
message("5. 13_SASH1_COL1A1_high_neighborhood_permutation_test.pdf")
message("")
message("解释方向：")
message("若observed_median_distance小于置换距离，并且closer p值较小，说明SASH1-high更靠近COL1A1-high。")
message("若observed_median_distance大于置换距离，并且farther p值较小，说明SASH1-high远离COL1A1-high。")
message("若两者都不显著，说明当前聚合空间数据不足以支持明确空间邻近或排他。")
message("============================================================\n")