# 11a_malignant_focus_cluster_audit.R

# 本脚本功能：
# 1. 审查当前用于Malignant_Focused的五个候选cluster：2、3、4、6、11
# 2. 整合04标准Seurat对象、08 final malignant call对象和trajectory metadata
# 3. 汇总候选cluster的细胞数、样本组成、CopyKAT支持比例和vertex-bin分布
# 4. 输出候选cluster UMAP、CopyKAT支持比例图、样本组成热图和marker DotPlot
# 5. 输出综合审查表，用于人工判断哪些cluster更接近论文中的Malignant_Cells核心区域
# 6. 本步骤不自动删除任何cluster，也不自动重定义恶性细胞
# 7. 记录本步骤session信息

# 本项目专用数据：
# GSE215403
# 12个OSCC单细胞样本：
# OSCC, scB1, scB2, scB5, scB7, scB8,
# scB9, scB10, scB12, scB13, scB14, scB15
#
# 本脚本关注五个tumor-related epithelial candidate clusters：
# 2=Differentiated_Tumor
# 3=Cycling_Tumor
# 4=CT_Antigen_Tumor
# 6=Tumor_Epithelial
# 11=Tumor_Epithelial
#
# 本步骤是人工审查步骤。
# 它只输出辅助决策图表，不直接改变final malignant call。
#
# 通用代码修改位置：
# 1. 换数据集时：
#    修改input_04_file、input_08_file和input_trajectory_metadata_file
#
# 2. 换cluster分辨率时：
#    修改cluster_column
#
# 3. 换candidate cluster定义时：
#    修改candidate_clusters和candidate_cluster_labels
#
# 4. 换marker审查基因时：
#    修改marker_genes


# ============================================================
# A. 加载包
# ============================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tidyr",
  "ggplot2",
  "patchwork",
  "pheatmap"
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
      )
    )
  )
}

library(Seurat)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(pheatmap)

# ============================================================
# B. 项目路径与文件夹
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
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# C. 设置输入文件
# ============================================================

input_04_file <- file.path(
  object_dir,
  "04_standard_Seurat_multi_resolution.rds"
)

input_08_file <- file.path(
  object_dir,
  "08_final_malignant_call.rds"
)

input_trajectory_metadata_file <- file.path(
  table_dir,
  "11_global_trajectory_cell_metadata_with_vertex_bins.csv"
)

if (!file.exists(input_04_file)) {
  stop(
    paste0(
      "找不到04对象：\n",
      input_04_file,
      "\n请先确认04结果文件是否已完成重命名。"
    )
  )
}

if (!file.exists(input_08_file)) {
  stop(
    paste0(
      "找不到08对象：\n",
      input_08_file,
      "\n请先确认08结果文件是否已完成重命名。"
    )
  )
}

if (!file.exists(input_trajectory_metadata_file)) {
  stop(
    paste0(
      "找不到trajectory metadata文件：\n",
      input_trajectory_metadata_file,
      "\n请先确认原08c结果是否已重命名为11前缀。"
    )
  )
}

# ============================================================
# D. 读取04、08与trajectory metadata
# ============================================================

sc_04 <- readRDS(
  input_04_file
)

sc_08 <- readRDS(
  input_08_file
)

DefaultAssay(sc_04) <- "RNA"
DefaultAssay(sc_08) <- "RNA"

trajectory_metadata <- read.csv(
  input_trajectory_metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

message(
  "04细胞数：",
  ncol(sc_04)
)

message(
  "08细胞数：",
  ncol(sc_08)
)

message(
  "trajectory metadata行数：",
  nrow(trajectory_metadata)
)

# ============================================================
# E. 恢复主cluster与候选恶性cluster标签
# ============================================================

cluster_column <- "cluster_res_0.2"

if (!cluster_column %in% colnames(sc_04@meta.data)) {
  
  alternative_cluster_column <- "RNA_snn_res.0.2"
  
  if (alternative_cluster_column %in% colnames(sc_04@meta.data)) {
    
    cluster_column <- alternative_cluster_column
    
  } else {
    
    stop(
      paste0(
        "04对象中未找到主cluster列。",
        "\n已尝试：cluster_res_0.2和RNA_snn_res.0.2"
      )
    )
  }
}

sc_04$paper_cluster_res_0_2 <- as.character(
  sc_04@meta.data[
    ,
    cluster_column
  ]
)

candidate_clusters <- c(
  "2",
  "3",
  "4",
  "6",
  "11"
)

candidate_cluster_labels <- c(
  "2" = "2_Differentiated_Tumor",
  "3" = "3_Cycling_Tumor",
  "4" = "4_CT_Antigen_Tumor",
  "6" = "6_Tumor_Epithelial",
  "11" = "11_Tumor_Epithelial"
)

sc_04$candidate_malignant_cluster <- ifelse(
  sc_04$paper_cluster_res_0_2 %in%
    candidate_clusters,
  candidate_cluster_labels[
    sc_04$paper_cluster_res_0_2
  ],
  "Other_cells"
)

sc_04$candidate_malignant_cluster <- factor(
  sc_04$candidate_malignant_cluster,
  levels = c(
    "Other_cells",
    unname(
      candidate_cluster_labels
    )
  )
)

# ============================================================
# F. 自动识别08中的CopyKAT严格恶性状态列
# ============================================================

# 不强制假设metadata列名。
# 自动搜索包含Strict_malignant_CopyKAT_aneuploid的列。

metadata_08 <- sc_08@meta.data

strict_status_column_candidates <- colnames(
  metadata_08
)[
  vapply(
    metadata_08,
    function(current_column) {
      
      any(
        as.character(current_column) ==
          "Strict_malignant_CopyKAT_aneuploid",
        na.rm = TRUE
      )
    },
    logical(1)
  )
]

if (length(strict_status_column_candidates) == 0) {
  stop(
    paste0(
      "08对象中未找到包含Strict_malignant_CopyKAT_aneuploid的状态列。",
      "\n请检查08对象metadata中的malignant_status_final等列。"
    )
  )
}

strict_status_column <- strict_status_column_candidates[
  1
]

message(
  "识别到08严格恶性状态列：",
  strict_status_column
)

copykat_status <- as.character(
  metadata_08[
    ,
    strict_status_column
  ]
)

names(copykat_status) <- rownames(
  metadata_08
)

# ============================================================
# G. 建立统一cell-level审查表
# ============================================================

audit_metadata <- sc_04@meta.data

audit_metadata$cell_barcode <- rownames(
  audit_metadata
)

audit_metadata$paper_cluster_res_0_2 <- as.character(
  audit_metadata$paper_cluster_res_0_2
)

audit_metadata$candidate_cluster <- ifelse(
  audit_metadata$paper_cluster_res_0_2 %in%
    candidate_clusters,
  candidate_cluster_labels[
    audit_metadata$paper_cluster_res_0_2
  ],
  "Other_cells"
)

audit_metadata$copykat_strict_status <- copykat_status[
  audit_metadata$cell_barcode
]

trajectory_metadata <- trajectory_metadata[
  ,
  intersect(
    c(
      "cell_barcode",
      "closest_vertex",
      "vertex_group",
      "vertex_bin",
      "paper_malignant_focus"
    ),
    colnames(trajectory_metadata)
  ),
  drop = FALSE
]

audit_metadata <- merge(
  x = audit_metadata,
  y = trajectory_metadata,
  by = "cell_barcode",
  all.x = TRUE,
  sort = FALSE
)

audit_metadata$candidate_cluster <- factor(
  audit_metadata$candidate_cluster,
  levels = c(
    "Other_cells",
    unname(
      candidate_cluster_labels
    )
  )
)

write.csv(
  audit_metadata,
  file.path(
    table_dir,
    "11a_candidate_malignant_cluster_cell_metadata.csv"
  ),
  row.names = FALSE
)

# ============================================================
# H. 候选cluster基础组成与样本来源
# ============================================================

candidate_cell_metadata <- audit_metadata[
  audit_metadata$paper_cluster_res_0_2 %in%
    candidate_clusters,
  ,
  drop = FALSE
]

cluster_size_summary <- candidate_cell_metadata %>%
  dplyr::count(
    paper_cluster_res_0_2,
    candidate_cluster,
    name = "cell_number"
  ) %>%
  dplyr::mutate(
    percent_of_candidate_cells = round(
      100 *
        cell_number /
        sum(cell_number),
      2
    )
  ) %>%
  dplyr::arrange(
    suppressWarnings(
      as.numeric(
        paper_cluster_res_0_2
      )
    )
  )

write.csv(
  cluster_size_summary,
  file.path(
    table_dir,
    "11a_candidate_cluster_size_summary.csv"
  ),
  row.names = FALSE
)

if (!"sample_id" %in% colnames(candidate_cell_metadata)) {
  stop(
    "04对象metadata中未找到sample_id列。"
  )
}

candidate_sample_summary <- candidate_cell_metadata %>%
  dplyr::count(
    paper_cluster_res_0_2,
    candidate_cluster,
    sample_id,
    name = "cell_number"
  ) %>%
  dplyr::group_by(
    paper_cluster_res_0_2
  ) %>%
  dplyr::mutate(
    percent_within_cluster = round(
      100 *
        cell_number /
        sum(cell_number),
      2
    )
  ) %>%
  dplyr::ungroup()

write.csv(
  candidate_sample_summary,
  file.path(
    table_dir,
    "11a_candidate_cluster_sample_composition.csv"
  ),
  row.names = FALSE
)

# ============================================================
# I. CopyKAT严格恶性支持比例
# ============================================================

copykat_support_summary <- candidate_cell_metadata %>%
  dplyr::mutate(
    strict_copykat_support = copykat_strict_status ==
      "Strict_malignant_CopyKAT_aneuploid"
  ) %>%
  dplyr::group_by(
    paper_cluster_res_0_2,
    candidate_cluster
  ) %>%
  dplyr::summarise(
    total_cell_number = dplyr::n(),
    strict_copykat_cell_number = sum(
      strict_copykat_support,
      na.rm = TRUE
    ),
    strict_copykat_percent = round(
      100 *
        strict_copykat_cell_number /
        total_cell_number,
      2
    ),
    copykat_missing_cell_number = sum(
      is.na(copykat_strict_status)
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      strict_copykat_percent
    )
  )

write.csv(
  copykat_support_summary,
  file.path(
    table_dir,
    "11a_candidate_cluster_CopyKAT_support_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# J. Trajectory vertex-bin分布
# ============================================================

vertex_bin_summary <- candidate_cell_metadata %>%
  dplyr::filter(
    !is.na(vertex_bin)
  ) %>%
  dplyr::count(
    paper_cluster_res_0_2,
    candidate_cluster,
    vertex_bin,
    name = "cell_number"
  ) %>%
  dplyr::group_by(
    paper_cluster_res_0_2
  ) %>%
  dplyr::mutate(
    percent_within_cluster = round(
      100 *
        cell_number /
        sum(cell_number),
      2
    )
  ) %>%
  dplyr::ungroup()

write.csv(
  vertex_bin_summary,
  file.path(
    table_dir,
    "11a_candidate_cluster_vertex_bin_distribution.csv"
  ),
  row.names = FALSE
)

# ============================================================
# K. 图1：五个候选cluster的UMAP位置
# ============================================================

find_umap_reduction <- function(
    seurat_object
) {
  
  reduction_names <- names(
    seurat_object@reductions
  )
  
  if (length(reduction_names) == 0) {
    return(NA_character_)
  }
  
  umap_candidates <- reduction_names[
    grepl(
      "umap",
      reduction_names,
      ignore.case = TRUE
    )
  ]
  
  if (length(umap_candidates) > 0) {
    return(
      umap_candidates[
        1
      ]
    )
  }
  
  two_dimensional_candidates <- reduction_names[
    vapply(
      reduction_names,
      function(current_reduction) {
        
        current_embedding <- Embeddings(
          seurat_object,
          reduction = current_reduction
        )
        
        ncol(current_embedding) >= 2
      },
      logical(1)
    )
  ]
  
  if (length(two_dimensional_candidates) > 0) {
    return(
      two_dimensional_candidates[
        1
      ]
    )
  }
  
  return(NA_character_)
}

plot_object <- sc_04

plot_reduction <- find_umap_reduction(
  sc_04
)

plot_object_label <- "04"

if (is.na(plot_reduction)) {
  
  message(
    "04对象未找到UMAP或其他二维降维，改用08对象。"
  )
  
  plot_object <- sc_08
  
  plot_reduction <- find_umap_reduction(
    sc_08
  )
  
  plot_object_label <- "08"
}

if (is.na(plot_reduction)) {
  
  stop(
    paste0(
      "04和08对象中均未找到可用二维降维。",
      "\n\n",
      "04可用reductions：",
      paste(
        names(sc_04@reductions),
        collapse = ", "
      ),
      "\n",
      "08可用reductions：",
      paste(
        names(sc_08@reductions),
        collapse = ", "
      )
    )
  )
}

message(
  "用于候选cluster UMAP的对象：",
  plot_object_label
)

message(
  "使用的reduction：",
  plot_reduction
)

audit_metadata_for_plot <- audit_metadata
rownames(audit_metadata_for_plot) <- audit_metadata_for_plot$cell_barcode

plot_object$candidate_malignant_cluster <- audit_metadata_for_plot[
  colnames(plot_object),
  "candidate_cluster"
]

plot_object$candidate_malignant_cluster <- as.character(
  plot_object$candidate_malignant_cluster
)

plot_object$candidate_malignant_cluster[
  is.na(
    plot_object$candidate_malignant_cluster
  )
] <- "Other_cells"

plot_object$candidate_malignant_cluster <- factor(
  plot_object$candidate_malignant_cluster,
  levels = c(
    "Other_cells",
    unname(
      candidate_cluster_labels
    )
  )
)

p_candidate_cluster_umap <- DimPlot(
  plot_object,
  reduction = plot_reduction,
  group.by = "candidate_malignant_cluster",
  pt.size = 0.25,
  label = FALSE,
  shuffle = TRUE
) +
  ggtitle(
    paste0(
      "Candidate Tumor-Related Clusters in ",
      plot_object_label,
      " ",
      plot_reduction,
      " Space"
    )
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
    "11a_candidate_malignant_clusters_UMAP.pdf"
  ),
  plot = p_candidate_cluster_umap,
  width = 11,
  height = 8
)

# ============================================================
# L. 图2：CopyKAT严格支持比例
# ============================================================

p_copykat_support <- ggplot(
  copykat_support_summary,
  aes(
    x = reorder(
      candidate_cluster,
      strict_copykat_percent
    ),
    y = strict_copykat_percent
  )
) +
  geom_col() +
  geom_text(
    aes(
      label = paste0(
        strict_copykat_percent,
        "%"
      )
    ),
    hjust = -0.08,
    size = 3.5
  ) +
  coord_flip(
    clip = "off"
  ) +
  expand_limits(
    y = max(
      copykat_support_summary$strict_copykat_percent,
      na.rm = TRUE
    ) + 12
  ) +
  labs(
    title = "CopyKAT Strict Aneuploid Support by Candidate Cluster",
    x = NULL,
    y = "Strict CopyKAT-supported cells (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.margin = margin(
      5.5,
      35,
      5.5,
      5.5
    )
  )

ggsave(
  filename = file.path(
    figure_dir,
    "11a_candidate_cluster_CopyKAT_support.pdf"
  ),
  plot = p_copykat_support,
  width = 10,
  height = 6
)

# ============================================================
# M. 图3：候选cluster样本组成热图
# ============================================================

sample_heatmap_matrix <- candidate_sample_summary %>%
  dplyr::select(
    candidate_cluster,
    sample_id,
    percent_within_cluster
  ) %>%
  tidyr::pivot_wider(
    names_from = sample_id,
    values_from = percent_within_cluster,
    values_fill = 0
  )

sample_heatmap_matrix <- as.data.frame(
  sample_heatmap_matrix
)

rownames(sample_heatmap_matrix) <- sample_heatmap_matrix$candidate_cluster

sample_heatmap_matrix$candidate_cluster <- NULL

pheatmap::pheatmap(
  mat = as.matrix(
    sample_heatmap_matrix
  ),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  main = "Sample Composition of Candidate Tumor-Related Clusters",
  filename = file.path(
    figure_dir,
    "11a_candidate_cluster_sample_composition_heatmap.pdf"
  ),
  width = 11,
  height = 5
)

# ============================================================
# N. 图4：候选cluster marker DotPlot
# ============================================================

marker_genes <- c(
  "EPCAM",
  "KRT8",
  "KRT18",
  "KRT19",
  "KRT5",
  "KRT14",
  "KRT17",
  "PTPRC",
  "LST1",
  "TYROBP",
  "DCN",
  "LUM",
  "COL3A1",
  "SASH1",
  "MYH11",
  "EMP1",
  "COL1A1"
)

marker_genes <- intersect(
  marker_genes,
  rownames(sc_04)
)

sc_marker_plot <- subset(
  sc_04,
  cells = candidate_cell_metadata$cell_barcode
)

candidate_cluster_for_plot <- unname(
  candidate_cluster_labels[
    as.character(
      sc_marker_plot$paper_cluster_res_0_2
    )
  ]
)

if (any(is.na(candidate_cluster_for_plot))) {
  
  unknown_cluster_values <- unique(
    as.character(
      sc_marker_plot$paper_cluster_res_0_2[
        is.na(candidate_cluster_for_plot)
      ]
    )
  )
  
  stop(
    paste0(
      "候选cluster标签映射失败：",
      paste(
        unknown_cluster_values,
        collapse = ", "
      )
    )
  )
}

sc_marker_plot$candidate_cluster <- factor(
  candidate_cluster_for_plot,
  levels = unname(
    candidate_cluster_labels
  )
)

p_marker_dotplot <- DotPlot(
  sc_marker_plot,
  features = marker_genes,
  group.by = "candidate_cluster",
  dot.scale = 7
) +
  RotatedAxis() +
  ggtitle(
    "Marker Review of Candidate Tumor-Related Clusters"
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
    "11a_candidate_cluster_marker_DotPlot.pdf"
  ),
  plot = p_marker_dotplot,
  width = 15,
  height = 6
)

# ============================================================
# O. 图5：候选cluster在trajectory vertex-bin中的分布
# ============================================================

p_vertex_bin_distribution <- ggplot(
  vertex_bin_summary,
  aes(
    x = vertex_bin,
    y = percent_within_cluster,
    fill = candidate_cluster
  )
) +
  geom_col(
    position = "dodge"
  ) +
  labs(
    title = "Vertex-bin Distribution of Candidate Tumor-Related Clusters",
    x = "Trajectory vertex bin",
    y = "Cells within candidate cluster (%)",
    fill = "Candidate cluster"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

ggsave(
  filename = file.path(
    figure_dir,
    "11a_candidate_cluster_vertex_bin_distribution.pdf"
  ),
  plot = p_vertex_bin_distribution,
  width = 13,
  height = 7
)

# ============================================================
# P. 输出综合审查表
# ============================================================

sample_diversity_summary <- candidate_sample_summary %>%
  dplyr::group_by(
    paper_cluster_res_0_2,
    candidate_cluster
  ) %>%
  dplyr::summarise(
    sample_number_with_cells = dplyr::n_distinct(
      sample_id
    ),
    sample_number_with_at_least_5_percent = sum(
      percent_within_cluster >= 5
    ),
    dominant_sample_percent = round(
      max(
        percent_within_cluster
      ),
      2
    ),
    .groups = "drop"
  )

vertex_diversity_summary <- vertex_bin_summary %>%
  dplyr::group_by(
    paper_cluster_res_0_2,
    candidate_cluster
  ) %>%
  dplyr::summarise(
    vertex_bin_number_present = dplyr::n_distinct(
      vertex_bin
    ),
    dominant_vertex_bin_percent = round(
      max(
        percent_within_cluster
      ),
      2
    ),
    .groups = "drop"
  )

candidate_cluster_review_summary <- cluster_size_summary %>%
  dplyr::left_join(
    copykat_support_summary,
    by = c(
      "paper_cluster_res_0_2",
      "candidate_cluster"
    )
  ) %>%
  dplyr::left_join(
    sample_diversity_summary,
    by = c(
      "paper_cluster_res_0_2",
      "candidate_cluster"
    )
  ) %>%
  dplyr::left_join(
    vertex_diversity_summary,
    by = c(
      "paper_cluster_res_0_2",
      "candidate_cluster"
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      strict_copykat_percent
    )
  )

write.csv(
  candidate_cluster_review_summary,
  file.path(
    table_dir,
    "11a_candidate_cluster_integrated_review_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# Q. 保存运行信息
# ============================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    table_dir,
    "11a_sessionInfo.txt"
  )
)

# ============================================================
# R. 最终提示
# ============================================================

message("\n============================================================")
message("11a_malignant_focus_cluster_audit.R运行完成。")
message("")
message("本步骤没有删除cluster，也没有重跑trajectory。")
message("请结合CopyKAT支持比例、marker表达、样本来源、vertex-bin分布和UMAP位置，")
message("再决定论文式Malignant_Focused的最终cluster定义。")
message("")
message("请重点查看：")
message("1. results/figures/11a_candidate_malignant_clusters_UMAP.pdf")
message("2. results/figures/11a_candidate_cluster_CopyKAT_support.pdf")
message("3. results/figures/11a_candidate_cluster_sample_composition_heatmap.pdf")
message("4. results/figures/11a_candidate_cluster_marker_DotPlot.pdf")
message("5. results/figures/11a_candidate_cluster_vertex_bin_distribution.pdf")
message("6. results/tables/11a_candidate_cluster_integrated_review_summary.csv")
message("============================================================\n")