# 11b_core_extended_malignant_overlay_rebuild.R

# 本脚本功能：
# 1. 基于04全细胞Seurat对象重新构建Monocle3 global graph
# 2. 比较两套论文式Malignant_Focused定义
# 3. Core_Malignant_Focused定义为cluster 6+11
# 4. Extended_Malignant_Focused定义为cluster 4+6+11
# 5. 分别输出Core和Extended的Malignant-Focused binary overlay
# 6. 分别输出SASH1、MYH11、EMP1、COL1A1 gene overlays
# 7. 输出gene overlay摘要表
# 8. 保存cell-level Monocle3 UMAP坐标、principal graph vertex和edge坐标
# 9. 记录本步骤session信息，避免后续再次重建trajectory

# 本项目专用数据：
# GSE215403
# 12个OSCC单细胞样本：
# OSCC, scB1, scB2, scB5, scB7, scB8,
# scB9, scB10, scB12, scB13, scB14, scB15
#
# 本步骤的global graph是全细胞state graph。
# 它用于展示Malignant_Focused细胞在全局细胞状态图中的位置，
# 不等同于严格意义上的pseudotime排序。
#
# 通用代码修改位置：
# 1. 换数据集时：
#    修改input_file和cluster_column
#
# 2. 换Core/Extended定义时：
#    修改core_malignant_focus和extended_malignant_focus对应的cluster
#
# 3. 换关注基因时：
#    修改target_genes
#
# 4. 调整Monocle3图构建参数时：
#    修改preprocess_cds、reduce_dimension、cluster_cells和learn_graph参数


# ============================================================
# A. 加载包
# ============================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "monocle3",
  "Matrix",
  "dplyr",
  "ggplot2",
  "patchwork",
  "igraph"
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
      "\n\n请不要整体更新R包，把缺少包名称发回。"
    )
  )
}

library(Seurat)
library(SeuratObject)
library(monocle3)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)
library(igraph)

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
# C. 读取04全细胞对象
# ============================================================

input_file <- file.path(
  object_dir,
  "04_standard_Seurat_multi_resolution.rds"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "找不到04对象：\n",
      input_file
    )
  )
}

sc <- readRDS(input_file)

DefaultAssay(sc) <- "RNA"

message(
  "11b输入细胞数：",
  ncol(sc)
)

message(
  "11b输入基因数：",
  nrow(sc)
)

# ============================================================
# D. 恢复04主cluster标签
# ============================================================

cluster_column <- "RNA_snn_res.0.2"

if (!cluster_column %in% colnames(sc@meta.data)) {
  stop(
    paste0(
      "04对象中未找到：",
      cluster_column
    )
  )
}

sc$paper_cluster_res_0_2 <- as.character(
  sc@meta.data[
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

sc$tumor_related_state <- ifelse(
  sc$paper_cluster_res_0_2 %in%
    candidate_clusters,
  "Tumor_related_state",
  "Background"
)

sc$core_malignant_focus <- ifelse(
  sc$paper_cluster_res_0_2 %in%
    c(
      "6",
      "11"
    ),
  "Core_Malignant_Focused",
  "Background"
)

sc$extended_malignant_focus <- ifelse(
  sc$paper_cluster_res_0_2 %in%
    c(
      "4",
      "6",
      "11"
    ),
  "Extended_Malignant_Focused",
  "Background"
)

sc$core_malignant_focus <- factor(
  sc$core_malignant_focus,
  levels = c(
    "Background",
    "Core_Malignant_Focused"
  )
)

sc$extended_malignant_focus <- factor(
  sc$extended_malignant_focus,
  levels = c(
    "Background",
    "Extended_Malignant_Focused"
  )
)

focus_definition_summary <- data.frame(
  focus_definition = c(
    "Core_Malignant_Focused",
    "Extended_Malignant_Focused",
    "All_tumor_related_states"
  ),
  included_clusters = c(
    "6, 11",
    "4, 6, 11",
    "2, 3, 4, 6, 11"
  ),
  cell_number = c(
    sum(
      sc$core_malignant_focus ==
        "Core_Malignant_Focused"
    ),
    sum(
      sc$extended_malignant_focus ==
        "Extended_Malignant_Focused"
    ),
    sum(
      sc$tumor_related_state ==
        "Tumor_related_state"
    )
  ),
  stringsAsFactors = FALSE
)

focus_definition_summary$percent_of_all_cells <- round(
  100 *
    focus_definition_summary$cell_number /
    ncol(sc),
  2
)

write.csv(
  focus_definition_summary,
  file.path(
    table_dir,
    "11b_malignant_focus_definition_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# E. 创建Monocle3 cell_data_set
# ============================================================

raw_counts <- LayerData(
  object = sc,
  assay = "RNA",
  layer = "counts"
)

gene_metadata <- data.frame(
  gene_short_name = rownames(raw_counts),
  row.names = rownames(raw_counts),
  stringsAsFactors = FALSE
)

cell_metadata <- sc@meta.data[
  colnames(raw_counts),
  ,
  drop = FALSE
]

cds <- monocle3::new_cell_data_set(
  expression_data = raw_counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

rm(raw_counts)
gc()

# ============================================================
# F. 重建global Monocle3 UMAP与principal graph
# ============================================================

# 参数与原global trajectory步骤保持一致。
# 该步骤只重建Monocle3图形空间，
# 不影响已有Seurat、CopyKAT和11a审查结果。

set.seed(1234)

cds <- monocle3::preprocess_cds(
  cds,
  num_dim = 50,
  method = "PCA",
  norm_method = "log",
  verbose = TRUE
)

cds <- monocle3::reduce_dimension(
  cds,
  reduction_method = "UMAP",
  preprocess_method = "PCA",
  umap.metric = "cosine",
  umap.min_dist = 0.3,
  umap.n_neighbors = 30,
  verbose = TRUE
)

cds <- monocle3::cluster_cells(
  cds,
  reduction_method = "UMAP",
  k = 20,
  cluster_method = "leiden",
  verbose = TRUE
)

cds <- monocle3::learn_graph(
  cds,
  use_partition = FALSE,
  close_loop = FALSE,
  learn_graph_control = list(
    prune_graph = TRUE
  ),
  verbose = TRUE
)

# ============================================================
# G. 提取每个细胞的Monocle3 UMAP坐标
# ============================================================

monocle_umap <- reducedDims(
  cds
)[["UMAP"]]

if (
  !all(
    colnames(cds) %in%
    rownames(monocle_umap)
  )
) {
  stop(
    "Monocle3 UMAP坐标与cds细胞名不一致。"
  )
}

umap_plot_data <- data.frame(
  cell_barcode = colnames(cds),
  UMAP_1 = monocle_umap[
    colnames(cds),
    1
  ],
  UMAP_2 = monocle_umap[
    colnames(cds),
    2
  ],
  paper_cluster_res_0_2 = as.character(
    colData(cds)$paper_cluster_res_0_2
  ),
  core_malignant_focus = as.character(
    colData(cds)$core_malignant_focus
  ),
  extended_malignant_focus = as.character(
    colData(cds)$extended_malignant_focus
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# H. 提取closest vertex与principal graph坐标
# ============================================================

principal_auxiliary <- monocle3::principal_graph_aux(
  cds
)[["UMAP"]]

closest_vertex_index <- as.character(
  as.data.frame(
    principal_auxiliary$pr_graph_cell_proj_closest_vertex
  )[
    colnames(cds),
    1
  ]
)

names(closest_vertex_index) <- colnames(cds)

trajectory_graph <- monocle3::principal_graph(
  cds
)[["UMAP"]]

graph_vertex_ids <- igraph::V(
  trajectory_graph
)$name

candidate_vertex_id_y <- paste0(
  "Y_",
  closest_vertex_index
)

if (
  sum(
    candidate_vertex_id_y %in%
    graph_vertex_ids
  ) < 2
) {
  stop(
    "closest vertex数字索引无法映射为Y_x格式principal graph ID。"
  )
}

closest_vertex <- candidate_vertex_id_y

graph_coordinate_matrix <- principal_auxiliary$dp_mst

if (
  is.null(graph_coordinate_matrix) ||
  !is.matrix(graph_coordinate_matrix)
) {
  stop(
    "无法提取dp_mst principal graph坐标。"
  )
}

if (
  ncol(graph_coordinate_matrix) ==
  length(graph_vertex_ids) &&
  nrow(graph_coordinate_matrix) >= 2
) {
  
  vertex_coordinates <- data.frame(
    vertex_id = colnames(
      graph_coordinate_matrix
    ),
    UMAP_1 = as.numeric(
      graph_coordinate_matrix[
        1,
        
      ]
    ),
    UMAP_2 = as.numeric(
      graph_coordinate_matrix[
        2,
        
      ]
    ),
    stringsAsFactors = FALSE
  )
  
} else if (
  nrow(graph_coordinate_matrix) ==
  length(graph_vertex_ids) &&
  ncol(graph_coordinate_matrix) >= 2
) {
  
  vertex_coordinates <- data.frame(
    vertex_id = rownames(
      graph_coordinate_matrix
    ),
    UMAP_1 = as.numeric(
      graph_coordinate_matrix[
        ,
        1
      ]
    ),
    UMAP_2 = as.numeric(
      graph_coordinate_matrix[
        ,
        2
      ]
    ),
    stringsAsFactors = FALSE
  )
  
} else {
  
  stop(
    paste0(
      "dp_mst维度异常：",
      paste(
        dim(graph_coordinate_matrix),
        collapse = "×"
      )
    )
  )
}

vertex_coordinates <- vertex_coordinates[
  vertex_coordinates$vertex_id %in%
    graph_vertex_ids,
  ,
  drop = FALSE
]

if (
  sum(
    unique(closest_vertex) %in%
    vertex_coordinates$vertex_id
  ) < 2
) {
  stop(
    "closest vertex与principal graph坐标未成功匹配。"
  )
}

# ============================================================
# I. 提取principal graph edge坐标
# ============================================================

edge_table <- as.data.frame(
  igraph::as_data_frame(
    trajectory_graph,
    what = "edges"
  ),
  stringsAsFactors = FALSE
)

if (
  !all(
    c(
      "from",
      "to"
    ) %in%
    colnames(edge_table)
  )
) {
  stop(
    "principal graph edge表中未找到from/to列。"
  )
}

colnames(edge_table) <- c(
  "from_vertex",
  "to_vertex"
)

vertex_coordinates_from <- vertex_coordinates[
  ,
  c(
    "vertex_id",
    "UMAP_1",
    "UMAP_2"
  ),
  drop = FALSE
]

colnames(vertex_coordinates_from) <- c(
  "from_vertex",
  "from_UMAP_1",
  "from_UMAP_2"
)

vertex_coordinates_to <- vertex_coordinates[
  ,
  c(
    "vertex_id",
    "UMAP_1",
    "UMAP_2"
  ),
  drop = FALSE
]

colnames(vertex_coordinates_to) <- c(
  "to_vertex",
  "to_UMAP_1",
  "to_UMAP_2"
)

edge_coordinates <- merge(
  edge_table,
  vertex_coordinates_from,
  by = "from_vertex",
  all.x = TRUE,
  sort = FALSE
)

edge_coordinates <- merge(
  edge_coordinates,
  vertex_coordinates_to,
  by = "to_vertex",
  all.x = TRUE,
  sort = FALSE
)

edge_coordinates <- edge_coordinates[
  complete.cases(
    edge_coordinates[
      ,
      c(
        "from_UMAP_1",
        "from_UMAP_2",
        "to_UMAP_1",
        "to_UMAP_2"
      )
    ]
  ),
  ,
  drop = FALSE
]

message(
  "11b principal graph edge数量：",
  nrow(edge_coordinates)
)

# ============================================================
# J. 提取四个目标基因的归一化表达
# ============================================================

target_genes <- c(
  "SASH1",
  "MYH11",
  "EMP1",
  "COL1A1"
)

cds_gene_symbols <- as.character(
  rowData(cds)$gene_short_name
)

names(cds_gene_symbols) <- rownames(cds)

target_gene_rows <- names(
  cds_gene_symbols[
    cds_gene_symbols %in%
      target_genes
  ]
)

target_gene_symbols_found <- cds_gene_symbols[
  target_gene_rows
]

if (length(target_gene_rows) == 0) {
  stop(
    "四个目标基因均未在cds中找到。"
  )
}

size_factor_values <- as.numeric(
  colData(cds)$Size_Factor
)

names(size_factor_values) <- colnames(cds)

invalid_size_factor <- is.na(
  size_factor_values
) |
  !is.finite(
    size_factor_values
  ) |
  size_factor_values <= 0

if (any(invalid_size_factor)) {
  
  count_matrix_for_sf <- SummarizedExperiment::assay(
    cds,
    "counts"
  )
  
  total_umi <- Matrix::colSums(
    count_matrix_for_sf
  )
  
  size_factor_values <- total_umi /
    median(
      total_umi[
        total_umi > 0
      ]
    )
  
  names(size_factor_values) <- colnames(cds)
  
  rm(count_matrix_for_sf)
  gc()
}

count_matrix <- SummarizedExperiment::assay(
  cds,
  "counts"
)

target_raw_counts <- as.matrix(
  count_matrix[
    target_gene_rows,
    colnames(cds),
    drop = FALSE
  ]
)

target_expression <- sweep(
  target_raw_counts,
  MARGIN = 2,
  STATS = size_factor_values[
    colnames(cds)
  ],
  FUN = "/"
)

target_expression <- log1p(
  target_expression
)

rownames(target_expression) <- target_gene_rows

colnames(target_expression) <- colnames(cds)

rm(
  count_matrix,
  target_raw_counts
)

gc()

# ============================================================
# K. 绘制binary focus overlay
# ============================================================

make_focus_overlay <- function(
    base_plot_data,
    focus_column,
    focus_label,
    graph_edge_data
) {
  
  plot_data <- base_plot_data
  
  plot_data$focus_status <- ifelse(
    plot_data[
      ,
      focus_column
    ] == focus_label,
    focus_label,
    "Background"
  )
  
  plot_data$focus_status <- factor(
    plot_data$focus_status,
    levels = c(
      "Background",
      focus_label
    )
  )
  
  ggplot() +
    geom_segment(
      data = graph_edge_data,
      aes(
        x = from_UMAP_1,
        y = from_UMAP_2,
        xend = to_UMAP_1,
        yend = to_UMAP_2
      ),
      color = "grey45",
      linewidth = 0.30,
      alpha = 0.85
    ) +
    geom_point(
      data = plot_data,
      aes(
        x = UMAP_1,
        y = UMAP_2,
        color = focus_status
      ),
      size = 0.34,
      alpha = 0.82
    ) +
    scale_color_manual(
      values = c(
        "Background" = "grey82",
        focus_label = "#00A6A6"
      )
    ) +
    coord_equal() +
    labs(
      title = paste0(
        "Global HNSCC Trajectory with ",
        focus_label
      ),
      x = "UMAP_1",
      y = "UMAP_2",
      color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      legend.position = "top"
    )
}

# ============================================================
# L. 绘制gene overlay
# ============================================================

# Background：
# 非当前focus定义中的细胞。
#
# Low：
# 当前focus中除High以外的所有细胞，
# 包括该基因表达为0的细胞。
#
# High：
# 当前focus内表达>0，并位于阳性细胞上四分位数的细胞。

make_gene_overlay <- function(
    gene_symbol,
    focus_column,
    focus_label,
    base_plot_data,
    expression_matrix,
    gene_symbol_vector,
    graph_edge_data
) {
  
  current_gene_row <- names(
    gene_symbol_vector[
      gene_symbol_vector ==
        gene_symbol
    ]
  )
  
  if (length(current_gene_row) == 0) {
    stop(
      paste0(
        "未找到目标基因：",
        gene_symbol
      )
    )
  }
  
  current_gene_row <- current_gene_row[
    1
  ]
  
  plot_data <- base_plot_data
  
  plot_data$gene_expression <- as.numeric(
    expression_matrix[
      current_gene_row,
      plot_data$cell_barcode
    ]
  )
  
  focus_index <- plot_data[
    ,
    focus_column
  ] == focus_label
  
  plot_data$expression_class <- "Background"
  
  plot_data$expression_class[
    focus_index
  ] <- "Low"
  
  positive_values <- plot_data$gene_expression[
    focus_index &
      plot_data$gene_expression > 0
  ]
  
  if (length(positive_values) >= 20) {
    
    high_cutoff <- as.numeric(
      quantile(
        positive_values,
        probs = 0.75,
        na.rm = TRUE
      )
    )
    
  } else {
    
    high_cutoff <- Inf
  }
  
  plot_data$expression_class[
    focus_index &
      plot_data$gene_expression > 0 &
      plot_data$gene_expression >= high_cutoff
  ] <- "High"
  
  plot_data$expression_class <- factor(
    plot_data$expression_class,
    levels = c(
      "Background",
      "Low",
      "High"
    )
  )
  
  p <- ggplot() +
    geom_segment(
      data = graph_edge_data,
      aes(
        x = from_UMAP_1,
        y = from_UMAP_2,
        xend = to_UMAP_1,
        yend = to_UMAP_2
      ),
      color = "grey45",
      linewidth = 0.30,
      alpha = 0.85
    ) +
    geom_point(
      data = plot_data,
      aes(
        x = UMAP_1,
        y = UMAP_2,
        color = expression_class
      ),
      size = 0.34,
      alpha = 0.85
    ) +
    scale_color_manual(
      values = c(
        "Background" = "grey84",
        "Low" = "#2C7FB8",
        "High" = "#D7301F"
      )
    ) +
    coord_equal() +
    labs(
      title = paste0(
        gene_symbol,
        " in ",
        focus_label
      ),
      x = "UMAP_1",
      y = "UMAP_2",
      color = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 10
      ),
      legend.position = "top"
    )
  
  focus_class <- plot_data$expression_class[
    focus_index
  ]
  
  summary_table <- as.data.frame(
    table(
      focus_class
    ),
    stringsAsFactors = FALSE
  )
  
  colnames(summary_table) <- c(
    "expression_class",
    "cell_number"
  )
  
  summary_table$gene <- gene_symbol
  
  summary_table$focus_definition <- focus_label
  
  summary_table$percent_within_focus <- round(
    100 *
      summary_table$cell_number /
      sum(
        summary_table$cell_number
      ),
    2
  )
  
  summary_table$high_cutoff_log_normalized_expression <- high_cutoff
  
  summary_table$detected_cell_number <- sum(
    plot_data$gene_expression[
      focus_index
    ] > 0
  )
  
  summary_table$detected_percent_within_focus <- round(
    100 *
      summary_table$detected_cell_number[
        1
      ] /
      sum(focus_index),
    2
  )
  
  return(
    list(
      plot = p,
      summary = summary_table
    )
  )
}

# ============================================================
# M. 输出Core与Extended两套focus binary图
# ============================================================

p_core_focus <- make_focus_overlay(
  base_plot_data = umap_plot_data,
  focus_column = "core_malignant_focus",
  focus_label = "Core_Malignant_Focused",
  graph_edge_data = edge_coordinates
)

p_extended_focus <- make_focus_overlay(
  base_plot_data = umap_plot_data,
  focus_column = "extended_malignant_focus",
  focus_label = "Extended_Malignant_Focused",
  graph_edge_data = edge_coordinates
)

p_focus_comparison <- p_core_focus +
  p_extended_focus +
  patchwork::plot_annotation(
    title = "Core versus Extended Malignant-Focused Definitions"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "11b_core_vs_extended_malignant_focus_comparison.pdf"
  ),
  plot = p_focus_comparison,
  width = 16,
  height = 8
)

# ============================================================
# N. 输出Core与Extended两套四基因overlay
# ============================================================

target_genes_for_plot <- intersect(
  target_genes,
  target_gene_symbols_found
)

run_gene_overlay_set <- function(
    focus_column,
    focus_label
) {
  
  overlay_results <- lapply(
    target_genes_for_plot,
    function(current_gene) {
      
      make_gene_overlay(
        gene_symbol = current_gene,
        focus_column = focus_column,
        focus_label = focus_label,
        base_plot_data = umap_plot_data,
        expression_matrix = target_expression,
        gene_symbol_vector = cds_gene_symbols,
        graph_edge_data = edge_coordinates
      )
    }
  )
  
  names(overlay_results) <- target_genes_for_plot
  
  overlay_plots <- lapply(
    overlay_results,
    function(x) {
      x$plot
    }
  )
  
  overlay_summary <- do.call(
    rbind,
    lapply(
      overlay_results,
      function(x) {
        x$summary
      }
    )
  )
  
  rownames(overlay_summary) <- NULL
  
  return(
    list(
      plots = overlay_plots,
      summary = overlay_summary
    )
  )
}

core_overlay_results <- run_gene_overlay_set(
  focus_column = "core_malignant_focus",
  focus_label = "Core_Malignant_Focused"
)

extended_overlay_results <- run_gene_overlay_set(
  focus_column = "extended_malignant_focus",
  focus_label = "Extended_Malignant_Focused"
)

p_core_gene_overlay <- patchwork::wrap_plots(
  core_overlay_results$plots,
  ncol = 2
) +
  patchwork::plot_annotation(
    title = "Core Malignant-Focused Gene Overlays"
  )

p_extended_gene_overlay <- patchwork::wrap_plots(
  extended_overlay_results$plots,
  ncol = 2
) +
  patchwork::plot_annotation(
    title = "Extended Malignant-Focused Gene Overlays"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "11b_core_malignant_focused_gene_overlays.pdf"
  ),
  plot = p_core_gene_overlay,
  width = 16,
  height = 13
)

ggsave(
  filename = file.path(
    figure_dir,
    "11b_extended_malignant_focused_gene_overlays.pdf"
  ),
  plot = p_extended_gene_overlay,
  width = 16,
  height = 13
)

gene_overlay_summary <- rbind(
  core_overlay_results$summary,
  extended_overlay_results$summary
)

write.csv(
  gene_overlay_summary,
  file.path(
    table_dir,
    "11b_core_extended_gene_overlay_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# O. 保存未来重画所需轻量化数据
# ============================================================

umap_plot_data$closest_vertex <- closest_vertex[
  umap_plot_data$cell_barcode
]

write.csv(
  umap_plot_data,
  file.path(
    table_dir,
    "11b_monocle3_cell_umap_coordinates_and_focus_labels.csv"
  ),
  row.names = FALSE
)

write.csv(
  vertex_coordinates,
  file.path(
    table_dir,
    "11b_principal_graph_vertex_coordinates.csv"
  ),
  row.names = FALSE
)

write.csv(
  edge_coordinates,
  file.path(
    table_dir,
    "11b_principal_graph_edge_coordinates.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    table_dir,
    "11b_sessionInfo.txt"
  )
)

# ============================================================
# P. 输出检查与完成提示
# ============================================================

required_output_files <- c(
  file.path(
    figure_dir,
    "11b_core_vs_extended_malignant_focus_comparison.pdf"
  ),
  file.path(
    figure_dir,
    "11b_core_malignant_focused_gene_overlays.pdf"
  ),
  file.path(
    figure_dir,
    "11b_extended_malignant_focused_gene_overlays.pdf"
  ),
  file.path(
    table_dir,
    "11b_malignant_focus_definition_summary.csv"
  ),
  file.path(
    table_dir,
    "11b_core_extended_gene_overlay_summary.csv"
  ),
  file.path(
    table_dir,
    "11b_monocle3_cell_umap_coordinates_and_focus_labels.csv"
  ),
  file.path(
    table_dir,
    "11b_principal_graph_vertex_coordinates.csv"
  ),
  file.path(
    table_dir,
    "11b_principal_graph_edge_coordinates.csv"
  ),
  file.path(
    table_dir,
    "11b_sessionInfo.txt"
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
    "11b_output_file_check.csv"
  ),
  row.names = FALSE
)

print(
  output_status
)

message("\n============================================================")
message("11b Core/Extended Malignant-Focused overlay完成。")
message("")
message("关键图：")
message("1. results/figures/11b_core_vs_extended_malignant_focus_comparison.pdf")
message("2. results/figures/11b_core_malignant_focused_gene_overlays.pdf")
message("3. results/figures/11b_extended_malignant_focused_gene_overlays.pdf")
message("")
message("后续优先用两套focus定义与论文Figure 5对照，")
message("再决定主图使用Core还是Extended。")
message("============================================================\n")