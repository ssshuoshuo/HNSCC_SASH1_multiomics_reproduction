# 16_scRNA_malignant_trajectory_Figure5_like.R

# 本脚本功能：
# 1. 读取前面11/11b生成的trajectory坐标表
# 2. 从多个Seurat对象和metadata表中补充细胞状态、cluster、cell type和malignant focus信息
# 3. 绘制论文Fig.5风格的malignant-focused trajectory status图
# 4. 绘制SASH1、MYH11、EMP1、COL1A1在trajectory上的表达图
# 5. 绘制4个核心基因gene-high状态在trajectory上的分布图
# 6. 输出接近论文Fig.5结构的组合panel
# 7. 若没有真实pseudotime，则只输出pseudotime_proxy作为supplementary/QC，不放入主图

# 本项目专用数据：
# GSE215403 scRNA-seq
#
# 输入优先使用：
# results/tables/11b_monocle3_cell_umap_coordinates_and_focus_labels.csv
# results/tables/11b_principal_graph_edge_coordinates.csv
# results/tables/11b_principal_graph_vertex_coordinates.csv
#
# 表达对象优先使用：
# results/objects/11_global_trajectory_Seurat.rds
# 若不存在，则回退到10、08、05对象。
#
# 项目兜底规则：
# 如果metadata中没有明确core/extended malignant focus标签，
# 使用前面人工审查确定的cluster规则：
# cluster6和11=Core malignant focus
# cluster4=Extended malignant focus
# cluster2和3=Candidate malignant-related tumor
#
# 通用代码修改位置：
# 1. 修改core_genes可替换核心基因
# 2. 修改high_quantile可调整gene-high阈值
# 3. 修改core_focus_clusters和extended_focus_clusters可换成其他项目的cluster规则


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
# B. 路径和参数
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

core_genes <- c(
  "SASH1",
  "MYH11",
  "EMP1",
  "COL1A1"
)

high_quantile <- 0.75

core_focus_clusters <- c(
  "6",
  "11"
)

extended_focus_clusters <- c(
  "4"
)

candidate_malignant_related_clusters <- c(
  "2",
  "3"
)

set.seed(
  20260710
)

candidate_trajectory_metadata_files <- c(
  file.path(
    table_dir,
    "11b_monocle3_cell_umap_coordinates_and_focus_labels.csv"
  ),
  file.path(
    table_dir,
    "11_global_trajectory_cell_metadata_with_vertex_bins.csv"
  ),
  file.path(
    table_dir,
    "11_global_trajectory_cell_metadata.csv"
  ),
  file.path(
    table_dir,
    "08e_monocle3_cell_umap_coordinates_and_focus_labels.csv"
  ),
  file.path(
    table_dir,
    "08c_global_trajectory_cell_metadata_with_vertex_bins.csv"
  ),
  file.path(
    table_dir,
    "08b_global_trajectory_cell_metadata.csv"
  )
)

candidate_graph_edge_files <- c(
  file.path(
    table_dir,
    "11b_principal_graph_edge_coordinates.csv"
  ),
  file.path(
    table_dir,
    "11_principal_graph_edge_coordinates.csv"
  ),
  file.path(
    table_dir,
    "08e_principal_graph_edge_coordinates.csv"
  ),
  file.path(
    table_dir,
    "08b_principal_graph_edge_coordinates.csv"
  )
)

candidate_graph_vertex_files <- c(
  file.path(
    table_dir,
    "11b_principal_graph_vertex_coordinates.csv"
  ),
  file.path(
    table_dir,
    "11_principal_graph_vertex_coordinates.csv"
  ),
  file.path(
    table_dir,
    "08e_principal_graph_vertex_coordinates.csv"
  ),
  file.path(
    table_dir,
    "08b_principal_graph_vertex_coordinates.csv"
  )
)

candidate_seurat_object_files <- c(
  file.path(
    object_dir,
    "11_global_trajectory_Seurat.rds"
  ),
  file.path(
    object_dir,
    "10_manual_review_epithelial_core.rds"
  ),
  file.path(
    object_dir,
    "08_final_malignant_call.rds"
  ),
  file.path(
    object_dir,
    "05_manual_annotated_plot_ready.rds"
  ),
  file.path(
    object_dir,
    "05_manual_annotated_before_malignant_call.rds"
  )
)

candidate_annotation_csv_files <- c(
  file.path(
    table_dir,
    "11b_monocle3_cell_umap_coordinates_and_focus_labels.csv"
  ),
  file.path(
    table_dir,
    "11a_candidate_malignant_cluster_cell_metadata.csv"
  ),
  file.path(
    table_dir,
    "11_global_trajectory_cell_metadata_with_vertex_bins.csv"
  ),
  file.path(
    table_dir,
    "11_global_trajectory_cell_metadata.csv"
  ),
  file.path(
    table_dir,
    "10_manual_review_cell_summary.csv"
  ),
  file.path(
    table_dir,
    "08e_monocle3_cell_umap_coordinates_and_focus_labels.csv"
  ),
  file.path(
    table_dir,
    "08d_candidate_malignant_cluster_cell_metadata.csv"
  )
)

# ============================================================
# C. 工具函数
# ============================================================

find_existing_file <- function(
    candidate_files,
    label
) {
  
  existing_files <- candidate_files[
    file.exists(
      candidate_files
    )
  ]
  
  if (length(existing_files) == 0) {
    stop(
      paste0(
        "未找到",
        label,
        "。已尝试：\n",
        paste(
          candidate_files,
          collapse = "\n"
        )
      )
    )
  }
  
  existing_files[1]
}

find_optional_file <- function(
    candidate_files
) {
  
  existing_files <- candidate_files[
    file.exists(
      candidate_files
    )
  ]
  
  if (length(existing_files) == 0) {
    return(
      NA_character_
    )
  }
  
  existing_files[1]
}

find_column_by_patterns <- function(
    data_frame,
    patterns,
    required = TRUE
) {
  
  column_names <- colnames(
    data_frame
  )
  
  lower_names <- tolower(
    column_names
  )
  
  matched_index <- integer(0)
  
  for (current_pattern in patterns) {
    
    current_match <- grep(
      current_pattern,
      lower_names,
      perl = TRUE
    )
    
    matched_index <- unique(
      c(
        matched_index,
        current_match
      )
    )
  }
  
  if (length(matched_index) == 0) {
    
    if (required) {
      stop(
        paste0(
          "未找到列。搜索模式：",
          paste(
            patterns,
            collapse = ", "
          ),
          "\n当前列名：",
          paste(
            column_names,
            collapse = ", "
          )
        )
      )
    }
    
    return(
      NA_character_
    )
  }
  
  column_names[
    matched_index[1]
  ]
}

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

make_unique_nonempty <- function(
    x
) {
  
  x <- as.character(
    x
  )
  
  x[
    is.na(
      x
    ) |
      x == ""
  ] <- NA_character_
  
  x
}

extract_barcode_column <- function(
    data_frame
) {
  
  barcode_column <- find_column_by_patterns(
    data_frame,
    patterns = c(
      "^cell$",
      "^barcode$",
      "^cell_barcode$",
      "cell.*barcode",
      "^cells$",
      "^barcodes$"
    ),
    required = FALSE
  )
  
  if (is.na(barcode_column)) {
    barcode_column <- colnames(
      data_frame
    )[1]
  }
  
  barcode_column
}

standardize_cluster <- function(
    cluster_vector
) {
  
  cluster_vector <- as.character(
    cluster_vector
  )
  
  cluster_vector <- gsub(
    "^c-",
    "",
    cluster_vector
  )
  
  cluster_vector <- gsub(
    "_.*$",
    "",
    cluster_vector
  )
  
  cluster_vector
}

# ============================================================
# D. 读取trajectory坐标和主要Seurat对象
# ============================================================

trajectory_metadata_file <- find_existing_file(
  candidate_trajectory_metadata_files,
  "trajectory metadata表"
)

seurat_object_file <- find_existing_file(
  candidate_seurat_object_files,
  "Seurat对象"
)

graph_edge_file <- find_optional_file(
  candidate_graph_edge_files
)

graph_vertex_file <- find_optional_file(
  candidate_graph_vertex_files
)

message(
  "使用trajectory metadata：",
  trajectory_metadata_file
)

message(
  "使用表达Seurat对象：",
  seurat_object_file
)

if (!is.na(graph_edge_file)) {
  message(
    "使用principal graph edge：",
    graph_edge_file
  )
}

if (!is.na(graph_vertex_file)) {
  message(
    "使用principal graph vertex：",
    graph_vertex_file
  )
}

trajectory_metadata <- read.csv(
  trajectory_metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

seurat_object <- readRDS(
  seurat_object_file
)

if ("RNA" %in% Assays(seurat_object)) {
  DefaultAssay(seurat_object) <- "RNA"
}

# ============================================================
# E. 自动识别trajectory坐标和barcode
# ============================================================

barcode_column <- extract_barcode_column(
  trajectory_metadata
)

x_column <- find_column_by_patterns(
  trajectory_metadata,
  patterns = c(
    "^umap_1$",
    "^umap1$",
    "umap.*1",
    "monocle.*1",
    "component_1",
    "dim_1",
    "^x$"
  )
)

y_column <- find_column_by_patterns(
  trajectory_metadata,
  patterns = c(
    "^umap_2$",
    "^umap2$",
    "umap.*2",
    "monocle.*2",
    "component_2",
    "dim_2",
    "^y$"
  )
)

trajectory_plot_data <- trajectory_metadata

trajectory_plot_data$cell_barcode <- as.character(
  trajectory_plot_data[
    ,
    barcode_column
  ]
)

trajectory_plot_data$trajectory_x <- as.numeric(
  trajectory_plot_data[
    ,
    x_column
  ]
)

trajectory_plot_data$trajectory_y <- as.numeric(
  trajectory_plot_data[
    ,
    y_column
  ]
)

trajectory_plot_data <- trajectory_plot_data[
  !is.na(
    trajectory_plot_data$trajectory_x
  ) &
    !is.na(
      trajectory_plot_data$trajectory_y
    ),
  ,
  drop = FALSE
]

trajectory_plot_data <- trajectory_plot_data[
  !duplicated(
    trajectory_plot_data$cell_barcode
  ),
  ,
  drop = FALSE
]

message(
  "trajectory细胞数：",
  nrow(
    trajectory_plot_data
  )
)

message(
  "barcode列：",
  barcode_column
)

message(
  "trajectory坐标列：",
  x_column,
  ", ",
  y_column
)

# ============================================================
# F. 合并额外Seurat对象metadata
# ============================================================

merge_metadata_by_barcode <- function(
    main_data,
    metadata_data,
    source_label
) {
  
  if (!"cell_barcode" %in% colnames(metadata_data)) {
    barcode_column_local <- extract_barcode_column(
      metadata_data
    )
    
    metadata_data$cell_barcode <- as.character(
      metadata_data[
        ,
        barcode_column_local
      ]
    )
  }
  
  metadata_data <- metadata_data[
    !duplicated(
      metadata_data$cell_barcode
    ),
    ,
    drop = FALSE
  ]
  
  rownames(metadata_data) <- metadata_data$cell_barcode
  
  shared_cells <- intersect(
    main_data$cell_barcode,
    metadata_data$cell_barcode
  )
  
  if (length(shared_cells) == 0) {
    return(
      main_data
    )
  }
  
  candidate_columns <- setdiff(
    colnames(metadata_data),
    c(
      "cell_barcode"
    )
  )
  
  for (current_column in candidate_columns) {
    
    if (
      current_column %in%
      c(
        "trajectory_x",
        "trajectory_y"
      )
    ) {
      next
    }
    
    new_column_name <- current_column
    
    if (new_column_name %in% colnames(main_data)) {
      
      existing_values <- main_data[
        ,
        new_column_name
      ]
      
      existing_empty <- all(
        is.na(
          existing_values
        )
      )
      
      if (!existing_empty) {
        new_column_name <- paste0(
          current_column,
          "_from_",
          source_label
        )
      }
    }
    
    value_vector <- rep(
      NA,
      nrow(
        main_data
      )
    )
    
    names(value_vector) <- main_data$cell_barcode
    
    value_vector[
      shared_cells
    ] <- metadata_data[
      shared_cells,
      current_column
    ]
    
    main_data[
      ,
      new_column_name
    ] <- value_vector[
      main_data$cell_barcode
    ]
  }
  
  main_data
}

for (current_object_file in candidate_seurat_object_files) {
  
  if (!file.exists(current_object_file)) {
    next
  }
  
  message(
    "合并metadata对象：",
    basename(
      current_object_file
    )
  )
  
  current_object <- readRDS(
    current_object_file
  )
  
  current_metadata <- current_object@meta.data
  
  current_metadata$cell_barcode <- rownames(
    current_metadata
  )
  
  source_label <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    tools::file_path_sans_ext(
      basename(
        current_object_file
      )
    )
  )
  
  trajectory_plot_data <- merge_metadata_by_barcode(
    main_data = trajectory_plot_data,
    metadata_data = current_metadata,
    source_label = source_label
  )
}

for (current_csv_file in candidate_annotation_csv_files) {
  
  if (!file.exists(current_csv_file)) {
    next
  }
  
  message(
    "合并metadata表：",
    basename(
      current_csv_file
    )
  )
  
  current_metadata <- read.csv(
    current_csv_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  source_label <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    tools::file_path_sans_ext(
      basename(
        current_csv_file
      )
    )
  )
  
  trajectory_plot_data <- merge_metadata_by_barcode(
    main_data = trajectory_plot_data,
    metadata_data = current_metadata,
    source_label = source_label
  )
}

interesting_columns <- grep(
  "malignant|focus|core|extended|final|status|celltype|cell_type|manual|cluster|pseudotime|vertex|bin",
  colnames(
    trajectory_plot_data
  ),
  value = TRUE,
  ignore.case = TRUE
)

write.csv(
  data.frame(
    column_name = interesting_columns,
    stringsAsFactors = FALSE
  ),
  file.path(
    table_dir,
    "16_available_trajectory_annotation_columns.csv"
  ),
  row.names = FALSE
)

# ============================================================
# G. 匹配Seurat表达
# ============================================================

matched_cells <- intersect(
  trajectory_plot_data$cell_barcode,
  colnames(
    seurat_object
  )
)

if (length(matched_cells) < 100) {
  stop(
    paste0(
      "trajectory表与表达Seurat对象匹配细胞数过少：",
      length(
        matched_cells
      ),
      "\n请检查barcode列。"
    )
  )
}

trajectory_plot_data <- trajectory_plot_data[
  trajectory_plot_data$cell_barcode %in%
    matched_cells,
  ,
  drop = FALSE
]

trajectory_plot_data <- trajectory_plot_data[
  match(
    matched_cells,
    trajectory_plot_data$cell_barcode
  ),
  ,
  drop = FALSE
]

seurat_object <- subset(
  seurat_object,
  cells = matched_cells
)

seurat_object <- seurat_object[
  ,
  trajectory_plot_data$cell_barcode
]

message(
  "匹配表达对象细胞数：",
  ncol(
    seurat_object
  )
)

data_layer_available <- TRUE

expression_matrix <- tryCatch(
  {
    LayerData(
      object = seurat_object,
      assay = DefaultAssay(
        seurat_object
      ),
      layer = "data"
    )
  },
  error = function(e) {
    data_layer_available <<- FALSE
    NULL
  }
)

if (!data_layer_available) {
  
  seurat_object <- NormalizeData(
    seurat_object,
    assay = DefaultAssay(
      seurat_object
    ),
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  
  expression_matrix <- LayerData(
    object = seurat_object,
    assay = DefaultAssay(
      seurat_object
    ),
    layer = "data"
  )
}

core_genes_found <- intersect(
  core_genes,
  rownames(
    seurat_object
  )
)

if (!all(core_genes %in% core_genes_found)) {
  stop(
    paste0(
      "表达对象未找到全部核心基因。找到：",
      paste(
        core_genes_found,
        collapse = ", "
      )
    )
  )
}

for (current_gene in core_genes_found) {
  
  trajectory_plot_data[
    ,
    current_gene
  ] <- as.numeric(
    expression_matrix[
      current_gene,
      trajectory_plot_data$cell_barcode
    ]
  )
}

# ============================================================
# G2. 强制补充稳定的celltype和cluster注释
# ============================================================

annotation_object_file <- file.path(
  object_dir,
  "05_manual_annotated_plot_ready.rds"
)

if (file.exists(annotation_object_file)) {
  
  annotation_object <- readRDS(
    annotation_object_file
  )
  
  annotation_metadata <- annotation_object@meta.data
  
  annotation_metadata$cell_barcode <- rownames(
    annotation_metadata
  )
  
  shared_annotation_cells <- intersect(
    trajectory_plot_data$cell_barcode,
    annotation_metadata$cell_barcode
  )
  
  message(
    "从05对象补充注释，匹配细胞数：",
    length(shared_annotation_cells)
  )
  
  rownames(annotation_metadata) <- annotation_metadata$cell_barcode
  
  trajectory_plot_data$Figure5_forced_celltype <- NA_character_
  trajectory_plot_data$Figure5_forced_cluster <- NA_character_
  
  trajectory_plot_data[
    trajectory_plot_data$cell_barcode %in% shared_annotation_cells,
    "Figure5_forced_celltype"
  ] <- annotation_metadata[
    trajectory_plot_data[
      trajectory_plot_data$cell_barcode %in% shared_annotation_cells,
      "cell_barcode"
    ],
    "celltype_plot"
  ]
  
  trajectory_plot_data[
    trajectory_plot_data$cell_barcode %in% shared_annotation_cells,
    "Figure5_forced_cluster"
  ] <- annotation_metadata[
    trajectory_plot_data[
      trajectory_plot_data$cell_barcode %in% shared_annotation_cells,
      "cell_barcode"
    ],
    "seurat_clusters"
  ]
  
} else {
  
  message(
    "未找到05_manual_annotated_plot_ready.rds，跳过强制注释补充。"
  )
  
  trajectory_plot_data$Figure5_forced_celltype <- NA_character_
  trajectory_plot_data$Figure5_forced_cluster <- NA_character_
}

if (file.exists(annotation_object_file)) {
  
  annotation_object <- readRDS(
    annotation_object_file
  )
  
  annotation_metadata <- annotation_object@meta.data
  annotation_metadata$cell_barcode <- rownames(annotation_metadata)
  rownames(annotation_metadata) <- annotation_metadata$cell_barcode
  
  shared_annotation_cells <- intersect(
    trajectory_plot_data$cell_barcode,
    annotation_metadata$cell_barcode
  )
  
  trajectory_plot_data$Figure5_forced_celltype <- NA_character_
  trajectory_plot_data$Figure5_forced_cluster <- NA_character_
  
  trajectory_plot_data[
    trajectory_plot_data$cell_barcode %in% shared_annotation_cells,
    "Figure5_forced_celltype"
  ] <- as.character(
    annotation_metadata[
      trajectory_plot_data[
        trajectory_plot_data$cell_barcode %in% shared_annotation_cells,
        "cell_barcode"
      ],
      "celltype_plot"
    ]
  )
  
  trajectory_plot_data[
    trajectory_plot_data$cell_barcode %in% shared_annotation_cells,
    "Figure5_forced_cluster"
  ] <- as.character(
    annotation_metadata[
      trajectory_plot_data[
        trajectory_plot_data$cell_barcode %in% shared_annotation_cells,
        "cell_barcode"
      ],
      "seurat_clusters"
    ]
  )
  
  message("Figure5_forced_celltype检查：")
  print(table(trajectory_plot_data$Figure5_forced_celltype, useNA = "ifany"))
  
  message("Figure5_forced_cluster检查：")
  print(table(trajectory_plot_data$Figure5_forced_cluster, useNA = "ifany"))
}

# ============================================================
# H. 构建论文Fig.5风格状态标签
# ============================================================

trajectory_plot_data$Figure5_cluster <- standardize_cluster(
  trajectory_plot_data$Figure5_forced_cluster
)

trajectory_plot_data$Figure5_celltype <- as.character(
  trajectory_plot_data$Figure5_forced_celltype
)

trajectory_plot_data$Figure5_celltype[
  is.na(
    trajectory_plot_data$Figure5_celltype
  ) |
    trajectory_plot_data$Figure5_celltype == ""
] <- "Unknown"

trajectory_plot_data$Figure5_status <- "Other non-tumor"
trajectory_plot_data$Figure5_status_source <- "forced_cluster_celltype_rule"

trajectory_plot_data$Figure5_status[
  trajectory_plot_data$Figure5_cluster %in%
    core_focus_clusters
] <- "Core malignant focus"

trajectory_plot_data$Figure5_status[
  trajectory_plot_data$Figure5_cluster %in%
    extended_focus_clusters
] <- "Extended malignant focus"

trajectory_plot_data$Figure5_status[
  trajectory_plot_data$Figure5_cluster %in%
    candidate_malignant_related_clusters
] <- "Candidate malignant-related tumor"

trajectory_plot_data$Figure5_status[
  !trajectory_plot_data$Figure5_cluster %in%
    c(
      core_focus_clusters,
      extended_focus_clusters,
      candidate_malignant_related_clusters
    ) &
    grepl(
      "tumor|epithelial|ct-antigen|cycling|differentiated",
      trajectory_plot_data$Figure5_celltype,
      ignore.case = TRUE
    )
] <- "Other tumor/epithelial"

trajectory_plot_data$Figure5_status <- factor(
  trajectory_plot_data$Figure5_status,
  levels = c(
    "Core malignant focus",
    "Extended malignant focus",
    "Candidate malignant-related tumor",
    "Other tumor/epithelial",
    "Other non-tumor"
  )
)

message(
  "Fig.5状态使用强制规则：cluster6/11 core，cluster4 extended，cluster2/3 candidate。"
)

print(
  table(
    trajectory_plot_data$Figure5_status,
    useNA = "ifany"
  )
)

print(
  table(
    trajectory_plot_data$Figure5_celltype,
    useNA = "ifany"
  )
)

# ============================================================
# I. pseudotime和vertex bin
# ============================================================

candidate_pseudotime_columns <- c(
  "pseudotime",
  "monocle3_pseudotime",
  "principal_graph_pseudotime",
  "trajectory_pseudotime",
  "pseudotime_order"
)

pseudotime_column <- candidate_pseudotime_columns[
  candidate_pseudotime_columns %in%
    colnames(
      trajectory_plot_data
    )
][1]

real_pseudotime_available <- TRUE

if (is.na(pseudotime_column)) {
  
  trajectory_plot_data$pseudotime_proxy <- rank(
    trajectory_plot_data$trajectory_x,
    ties.method = "average"
  )
  
  pseudotime_column <- "pseudotime_proxy"
  
  real_pseudotime_available <- FALSE
  
  message(
    "未找到真实pseudotime列，使用trajectory_x rank生成pseudotime_proxy。"
  )
}

candidate_vertex_bin_columns <- c(
  "vertex_bin",
  "paper_style_vertex_bin",
  "principal_graph_vertex_bin",
  "trajectory_vertex_bin",
  "vertex_group"
)

vertex_bin_column <- candidate_vertex_bin_columns[
  candidate_vertex_bin_columns %in%
    colnames(
      trajectory_plot_data
    )
][1]

if (is.na(vertex_bin_column)) {
  
  vertex_bin_column <- grep(
    "vertex|bin",
    colnames(
      trajectory_plot_data
    ),
    value = TRUE,
    ignore.case = TRUE
  )[1]
}

if (is.na(vertex_bin_column)) {
  trajectory_plot_data$Figure5_vertex_bin <- "Vertex bin not available"
  vertex_bin_column <- "Figure5_vertex_bin"
}

# ============================================================
# J. 定义核心基因high状态
# ============================================================

gene_threshold_table <- data.frame(
  gene = character(0),
  positive_cell_number = integer(0),
  high_threshold = numeric(0),
  high_cell_number = integer(0),
  high_cell_percent = numeric(0),
  stringsAsFactors = FALSE
)

for (current_gene in core_genes) {
  
  current_threshold <- calculate_positive_quantile_threshold(
    trajectory_plot_data[
      ,
      current_gene
    ],
    high_quantile
  )
  
  current_high_column <- paste0(
    current_gene,
    "_high"
  )
  
  trajectory_plot_data[
    ,
    current_high_column
  ] <- trajectory_plot_data[
    ,
    current_gene
  ] >= current_threshold
  
  trajectory_plot_data[
    is.na(
      trajectory_plot_data[
        ,
        current_high_column
      ]
    ),
    current_high_column
  ] <- FALSE
  
  gene_threshold_table <- rbind(
    gene_threshold_table,
    data.frame(
      gene = current_gene,
      positive_cell_number = sum(
        trajectory_plot_data[
          ,
          current_gene
        ] > 0,
        na.rm = TRUE
      ),
      high_threshold = current_threshold,
      high_cell_number = sum(
        trajectory_plot_data[
          ,
          current_high_column
        ],
        na.rm = TRUE
      ),
      high_cell_percent = round(
        100 *
          mean(
            trajectory_plot_data[
              ,
              current_high_column
            ],
            na.rm = TRUE
          ),
        2
      ),
      stringsAsFactors = FALSE
    )
  )
}

# ============================================================
# K. 输出metadata和统计表
# ============================================================

write.csv(
  gene_threshold_table,
  file.path(
    table_dir,
    "16_Figure5_core_gene_high_thresholds.csv"
  ),
  row.names = FALSE
)

write.csv(
  trajectory_plot_data,
  file.path(
    table_dir,
    "16_Figure5_trajectory_metadata_with_core_gene_expression.csv"
  ),
  row.names = FALSE
)

trajectory_status_summary <- trajectory_plot_data %>%
  dplyr::group_by(
    Figure5_status
  ) %>%
  dplyr::summarise(
    cell_number = dplyr::n(),
    median_pseudotime_value = median(
      .data[[pseudotime_column]],
      na.rm = TRUE
    ),
    SASH1_mean = mean(
      SASH1,
      na.rm = TRUE
    ),
    MYH11_mean = mean(
      MYH11,
      na.rm = TRUE
    ),
    EMP1_mean = mean(
      EMP1,
      na.rm = TRUE
    ),
    COL1A1_mean = mean(
      COL1A1,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  trajectory_status_summary,
  file.path(
    table_dir,
    "16_Figure5_trajectory_status_core_gene_summary.csv"
  ),
  row.names = FALSE
)

gene_high_status_summary <- data.frame(
  gene = character(0),
  Figure5_status = character(0),
  cell_number = integer(0),
  percent_of_gene_high_cells = numeric(0),
  stringsAsFactors = FALSE
)

for (current_gene in core_genes) {
  
  current_high_column <- paste0(
    current_gene,
    "_high"
  )
  
  current_table <- trajectory_plot_data[
    trajectory_plot_data[
      ,
      current_high_column
    ],
    ,
    drop = FALSE
  ] %>%
    dplyr::group_by(
      Figure5_status
    ) %>%
    dplyr::summarise(
      cell_number = dplyr::n(),
      .groups = "drop"
    )
  
  if (nrow(current_table) > 0) {
    
    current_table$gene <- current_gene
    
    current_table$percent_of_gene_high_cells <- round(
      100 *
        current_table$cell_number /
        sum(
          current_table$cell_number
        ),
      2
    )
    
    gene_high_status_summary <- rbind(
      gene_high_status_summary,
      current_table[
        ,
        c(
          "gene",
          "Figure5_status",
          "cell_number",
          "percent_of_gene_high_cells"
        )
      ]
    )
  }
}

write.csv(
  gene_high_status_summary,
  file.path(
    table_dir,
    "16_Figure5_gene_high_status_distribution_summary.csv"
  ),
  row.names = FALSE
)

status_source_summary <- trajectory_plot_data %>%
  dplyr::group_by(
    Figure5_status_source,
    Figure5_status
  ) %>%
  dplyr::summarise(
    cell_number = dplyr::n(),
    .groups = "drop"
  )

write.csv(
  status_source_summary,
  file.path(
    table_dir,
    "16_Figure5_status_source_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# L. 读取principal graph
# ============================================================

graph_edge_data <- NULL

if (!is.na(graph_edge_file)) {
  
  raw_graph_edge_data <- read.csv(
    graph_edge_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  edge_x_column <- find_column_by_patterns(
    raw_graph_edge_data,
    patterns = c(
      "^x$",
      "x_start",
      "source.*x",
      "from.*x",
      "start.*x",
      "x1"
    ),
    required = FALSE
  )
  
  edge_y_column <- find_column_by_patterns(
    raw_graph_edge_data,
    patterns = c(
      "^y$",
      "y_start",
      "source.*y",
      "from.*y",
      "start.*y",
      "y1"
    ),
    required = FALSE
  )
  
  edge_xend_column <- find_column_by_patterns(
    raw_graph_edge_data,
    patterns = c(
      "^xend$",
      "x_end",
      "target.*x",
      "to.*x",
      "end.*x",
      "x2"
    ),
    required = FALSE
  )
  
  edge_yend_column <- find_column_by_patterns(
    raw_graph_edge_data,
    patterns = c(
      "^yend$",
      "y_end",
      "target.*y",
      "to.*y",
      "end.*y",
      "y2"
    ),
    required = FALSE
  )
  
  if (
    !any(
      is.na(
        c(
          edge_x_column,
          edge_y_column,
          edge_xend_column,
          edge_yend_column
        )
      )
    )
  ) {
    
    graph_edge_data <- data.frame(
      x = as.numeric(
        raw_graph_edge_data[
          ,
          edge_x_column
        ]
      ),
      y = as.numeric(
        raw_graph_edge_data[
          ,
          edge_y_column
        ]
      ),
      xend = as.numeric(
        raw_graph_edge_data[
          ,
          edge_xend_column
        ]
      ),
      yend = as.numeric(
        raw_graph_edge_data[
          ,
          edge_yend_column
        ]
      ),
      stringsAsFactors = FALSE
    )
    
    graph_edge_data <- graph_edge_data[
      stats::complete.cases(
        graph_edge_data
      ),
      ,
      drop = FALSE
    ]
  }
}

graph_vertex_data <- NULL

if (!is.na(graph_vertex_file)) {
  
  raw_graph_vertex_data <- read.csv(
    graph_vertex_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  vertex_x_column <- find_column_by_patterns(
    raw_graph_vertex_data,
    patterns = c(
      "^x$",
      "umap.*1",
      "component_1",
      "dim_1"
    ),
    required = FALSE
  )
  
  vertex_y_column <- find_column_by_patterns(
    raw_graph_vertex_data,
    patterns = c(
      "^y$",
      "umap.*2",
      "component_2",
      "dim_2"
    ),
    required = FALSE
  )
  
  if (
    !is.na(
      vertex_x_column
    ) &&
    !is.na(
      vertex_y_column
    )
  ) {
    
    graph_vertex_data <- data.frame(
      x = as.numeric(
        raw_graph_vertex_data[
          ,
          vertex_x_column
        ]
      ),
      y = as.numeric(
        raw_graph_vertex_data[
          ,
          vertex_y_column
        ]
      ),
      stringsAsFactors = FALSE
    )
    
    graph_vertex_data <- graph_vertex_data[
      stats::complete.cases(
        graph_vertex_data
      ),
      ,
      drop = FALSE
    ]
  }
}

# ============================================================
# M. 绘图函数
# ============================================================

add_graph_layers <- function(
    plot_object
) {
  
  if (
    !is.null(
      graph_edge_data
    ) &&
    nrow(
      graph_edge_data
    ) > 0
  ) {
    
    plot_object <- plot_object +
      geom_segment(
        data = graph_edge_data,
        aes(
          x = x,
          y = y,
          xend = xend,
          yend = yend
        ),
        inherit.aes = FALSE,
        linewidth = 0.28,
        color = "black",
        alpha = 0.75
      )
  }
  
  if (
    !is.null(
      graph_vertex_data
    ) &&
    nrow(
      graph_vertex_data
    ) > 0
  ) {
    
    plot_object <- plot_object +
      geom_point(
        data = graph_vertex_data,
        aes(
          x = x,
          y = y
        ),
        inherit.aes = FALSE,
        size = 0.7,
        color = "black",
        alpha = 0.8
      )
  }
  
  plot_object
}

base_trajectory_theme <- function() {
  
  theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      legend.position = "right"
    )
}

status_colors <- c(
  "Core malignant focus" = "#D73027",
  "Extended malignant focus" = "#FC8D59",
  "Candidate malignant-related tumor" = "#FEE08B",
  "Other tumor/epithelial" = "#91BFDB",
  "Other non-tumor" = "#D9D9D9"
)

status_levels_present <- levels(
  trajectory_plot_data$Figure5_status
)

status_colors_present <- status_colors[
  names(
    status_colors
  ) %in%
    status_levels_present
]

plot_trajectory_status <- function() {
  
  p <- ggplot(
    trajectory_plot_data,
    aes(
      x = trajectory_x,
      y = trajectory_y
    )
  )
  
  p <- add_graph_layers(
    p
  )
  
  p +
    geom_point(
      aes(
        color = Figure5_status
      ),
      size = 0.55,
      alpha = 0.88
    ) +
    scale_color_manual(
      values = status_colors_present,
      na.value = "grey80"
    ) +
    labs(
      title = "(a) Malignant-focused trajectory status",
      x = "Trajectory UMAP_1",
      y = "Trajectory UMAP_2",
      color = "Status"
    ) +
    base_trajectory_theme()
}

plot_trajectory_celltype <- function() {
  
  p <- ggplot(
    trajectory_plot_data,
    aes(
      x = trajectory_x,
      y = trajectory_y
    )
  )
  
  p <- add_graph_layers(
    p
  )
  
  p +
    geom_point(
      aes(
        color = Figure5_celltype
      ),
      size = 0.5,
      alpha = 0.85
    ) +
    labs(
      title = "(b) Cell type along trajectory",
      x = "Trajectory UMAP_1",
      y = "Trajectory UMAP_2",
      color = "Cell type"
    ) +
    base_trajectory_theme()
}

plot_trajectory_continuous <- function(
    color_column,
    plot_title
) {
  
  current_plot_data <- trajectory_plot_data
  
  current_plot_data$positive_expression <- current_plot_data[
    ,
    color_column
  ] > 0
  
  p <- ggplot(
    current_plot_data,
    aes(
      x = trajectory_x,
      y = trajectory_y
    )
  )
  
  p <- add_graph_layers(
    p
  )
  
  p +
    geom_point(
      color = "grey82",
      size = 0.45,
      alpha = 0.75
    ) +
    geom_point(
      data = current_plot_data[
        current_plot_data$positive_expression,
        ,
        drop = FALSE
      ],
      aes(
        color = .data[[color_column]]
      ),
      size = 0.55,
      alpha = 0.9
    ) +
    scale_color_gradientn(
      colors = c(
        "#D8C7FF",
        "#7B61FF",
        "#1D4ED8"
      )
    ) +
    labs(
      title = plot_title,
      x = "Trajectory UMAP_1",
      y = "Trajectory UMAP_2",
      color = color_column
    ) +
    base_trajectory_theme()
}

plot_gene_high_status <- function(
    gene_name
) {
  
  current_high_column <- paste0(
    gene_name,
    "_high"
  )
  
  current_plot_data <- trajectory_plot_data
  
  current_plot_data$gene_high_status <- ifelse(
    current_plot_data[
      ,
      current_high_column
    ],
    paste0(
      gene_name,
      "-high"
    ),
    "Other"
  )
  
  current_plot_data$gene_high_status <- factor(
    current_plot_data$gene_high_status,
    levels = c(
      "Other",
      paste0(
        gene_name,
        "-high"
      )
    )
  )
  
  p <- ggplot(
    current_plot_data,
    aes(
      x = trajectory_x,
      y = trajectory_y
    )
  )
  
  p <- add_graph_layers(
    p
  )
  
  p +
    geom_point(
      data = current_plot_data[
        current_plot_data$gene_high_status == "Other",
        ,
        drop = FALSE
      ],
      color = "grey82",
      size = 0.45,
      alpha = 0.7
    ) +
    geom_point(
      data = current_plot_data[
        current_plot_data$gene_high_status != "Other",
        ,
        drop = FALSE
      ],
      color = "#D73027",
      size = 0.65,
      alpha = 0.95
    ) +
    labs(
      title = paste0(
        gene_name,
        "-high cells on trajectory"
      ),
      x = "Trajectory UMAP_1",
      y = "Trajectory UMAP_2",
      color = NULL
    ) +
    base_trajectory_theme() +
    theme(
      legend.position = "none"
    )
}

plot_pseudotime_proxy <- function() {
  
  p <- ggplot(
    trajectory_plot_data,
    aes(
      x = trajectory_x,
      y = trajectory_y
    )
  )
  
  p <- add_graph_layers(
    p
  )
  
  p +
    geom_point(
      aes(
        color = .data[[pseudotime_column]]
      ),
      size = 0.55,
      alpha = 0.9
    ) +
    scale_color_gradientn(
      colors = c(
        "grey85",
        "#7B61FF",
        "#1D4ED8"
      )
    ) +
    labs(
      title = "Supplementary pseudotime/proxy",
      x = "Trajectory UMAP_1",
      y = "Trajectory UMAP_2",
      color = pseudotime_column
    ) +
    base_trajectory_theme()
}

# ============================================================
# N. 输出trajectory状态图
# ============================================================

p_status <- plot_trajectory_status()

p_celltype <- plot_trajectory_celltype()

p_pseudotime <- plot_pseudotime_proxy()

ggsave(
  filename = file.path(
    figure_dir,
    "16_Figure5a_trajectory_status.pdf"
  ),
  plot = p_status,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(
    figure_dir,
    "16_Figure5b_trajectory_celltype.pdf"
  ),
  plot = p_celltype,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(
    figure_dir,
    "16_Figure5c_supplementary_pseudotime_or_proxy.pdf"
  ),
  plot = p_pseudotime,
  width = 8,
  height = 6
)

# ============================================================
# O. 核心基因表达trajectory overlay
# ============================================================

gene_expression_plot_list <- list()

for (current_gene in core_genes) {
  
  gene_expression_plot_list[[current_gene]] <- plot_trajectory_continuous(
    color_column = current_gene,
    plot_title = paste0(
      current_gene,
      " expression on trajectory"
    )
  )
  
  ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        "16_Figure5_trajectory_expression_",
        current_gene,
        ".pdf"
      )
    ),
    plot = gene_expression_plot_list[[current_gene]],
    width = 8,
    height = 6
  )
}

p_gene_expression_panel <- patchwork::wrap_plots(
  gene_expression_plot_list,
  ncol = 2
) +
  patchwork::plot_annotation(
    title = "Core Gene Expression Along Malignant-focused Trajectory"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "16_Figure5_core_gene_trajectory_expression_panel.pdf"
  ),
  plot = p_gene_expression_panel,
  width = 14,
  height = 11
)

# ============================================================
# P. 核心基因high状态trajectory overlay
# ============================================================

gene_high_plot_list <- list()

for (current_gene in core_genes) {
  
  gene_high_plot_list[[current_gene]] <- plot_gene_high_status(
    current_gene
  )
  
  ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        "16_Figure5_trajectory_high_status_",
        current_gene,
        ".pdf"
      )
    ),
    plot = gene_high_plot_list[[current_gene]],
    width = 8,
    height = 6
  )
}

p_gene_high_panel <- patchwork::wrap_plots(
  gene_high_plot_list,
  ncol = 2
) +
  patchwork::plot_annotation(
    title = "Core Gene-high Cells Along Malignant-focused Trajectory"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "16_Figure5_core_gene_high_status_trajectory_panel.pdf"
  ),
  plot = p_gene_high_panel,
  width = 14,
  height = 11
)

# ============================================================
# Q. supplementary pseudotime趋势图
# ============================================================

trend_plot_data <- trajectory_plot_data[
  ,
  c(
    "cell_barcode",
    pseudotime_column,
    core_genes
  )
]

colnames(
  trend_plot_data
)[
  colnames(
    trend_plot_data
  ) ==
    pseudotime_column
] <- "pseudotime_value"

trend_long_data <- data.frame(
  cell_barcode = rep(
    trend_plot_data$cell_barcode,
    times = length(
      core_genes
    )
  ),
  pseudotime_value = rep(
    trend_plot_data$pseudotime_value,
    times = length(
      core_genes
    )
  ),
  gene = rep(
    core_genes,
    each = nrow(
      trend_plot_data
    )
  ),
  expression = unlist(
    trend_plot_data[
      ,
      core_genes
    ],
    use.names = FALSE
  ),
  stringsAsFactors = FALSE
)

p_gene_trend <- ggplot(
  trend_long_data,
  aes(
    x = pseudotime_value,
    y = expression
  )
) +
  geom_point(
    size = 0.18,
    alpha = 0.12
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = FALSE,
    linewidth = 0.7
  ) +
  facet_wrap(
    ~gene,
    scales = "free_y",
    ncol = 2
  ) +
  labs(
    title = ifelse(
      real_pseudotime_available,
      "Core Gene Expression Trend Along Pseudotime",
      "Supplementary Core Gene Trend Along Pseudotime Proxy"
    ),
    x = ifelse(
      real_pseudotime_available,
      "Pseudotime",
      "Pseudotime proxy"
    ),
    y = "Normalized expression"
  ) +
  theme_classic(base_size = 10) +
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
    "16_Figure5_supplementary_core_gene_pseudotime_or_proxy_trend.pdf"
  ),
  plot = p_gene_trend,
  width = 11,
  height = 8
)

# 兼容旧文件名，方便不破坏README或旧引用
ggsave(
  filename = file.path(
    figure_dir,
    "16_Figure5_core_gene_pseudotime_trend.pdf"
  ),
  plot = p_gene_trend,
  width = 11,
  height = 8
)

# ============================================================
# R. 合并论文Fig.5-like主图
# ============================================================

p_figure5_like <- (
  p_status +
    p_celltype
) / (
  p_gene_expression_panel
) / (
  p_gene_high_panel
) +
  patchwork::plot_annotation(
    title = "Figure5-like Malignant-focused Trajectory Analysis",
    subtitle = "Trajectory status, core gene expression, and gene-high cells"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "16_Figure5_like_malignant_trajectory_core_gene_panel.pdf"
  ),
  plot = p_figure5_like,
  width = 16,
  height = 22
)

# ============================================================
# S. 保存对象和运行信息
# ============================================================

rownames(
  trajectory_plot_data
) <- trajectory_plot_data$cell_barcode

for (current_gene in core_genes) {
  
  current_high_column <- paste0(
    current_gene,
    "_high"
  )
  
  seurat_object[[current_high_column]] <- trajectory_plot_data[
    colnames(
      seurat_object
    ),
    current_high_column
  ]
}

seurat_object$Figure5_status <- trajectory_plot_data[
  colnames(
    seurat_object
  ),
  "Figure5_status"
]

saveRDS(
  seurat_object,
  file.path(
    object_dir,
    "16_Figure5_like_malignant_trajectory_core_gene_Seurat.rds"
  ),
  compress = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    table_dir,
    "16_sessionInfo.txt"
  )
)

# ============================================================
# T. 输出检查
# ============================================================

required_output_files <- c(
  file.path(
    object_dir,
    "16_Figure5_like_malignant_trajectory_core_gene_Seurat.rds"
  ),
  file.path(
    table_dir,
    "16_available_trajectory_annotation_columns.csv"
  ),
  file.path(
    table_dir,
    "16_Figure5_core_gene_high_thresholds.csv"
  ),
  file.path(
    table_dir,
    "16_Figure5_trajectory_metadata_with_core_gene_expression.csv"
  ),
  file.path(
    table_dir,
    "16_Figure5_trajectory_status_core_gene_summary.csv"
  ),
  file.path(
    table_dir,
    "16_Figure5_gene_high_status_distribution_summary.csv"
  ),
  file.path(
    table_dir,
    "16_Figure5_status_source_summary.csv"
  ),
  file.path(
    figure_dir,
    "16_Figure5a_trajectory_status.pdf"
  ),
  file.path(
    figure_dir,
    "16_Figure5b_trajectory_celltype.pdf"
  ),
  file.path(
    figure_dir,
    "16_Figure5c_supplementary_pseudotime_or_proxy.pdf"
  ),
  file.path(
    figure_dir,
    "16_Figure5_core_gene_trajectory_expression_panel.pdf"
  ),
  file.path(
    figure_dir,
    "16_Figure5_core_gene_high_status_trajectory_panel.pdf"
  ),
  file.path(
    figure_dir,
    "16_Figure5_supplementary_core_gene_pseudotime_or_proxy_trend.pdf"
  ),
  file.path(
    figure_dir,
    "16_Figure5_like_malignant_trajectory_core_gene_panel.pdf"
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
    "16_output_file_check.csv"
  ),
  row.names = FALSE
)

print(
  output_status
)

message("\n============================================================")
message("16 论文Fig.5风格恶性细胞trajectory核心基因图完成。")
message("")
message("重点查看：")
message("1. 16_Figure5_like_malignant_trajectory_core_gene_panel.pdf")
message("2. 16_Figure5a_trajectory_status.pdf")
message("3. 16_Figure5_core_gene_trajectory_expression_panel.pdf")
message("4. 16_Figure5_core_gene_high_status_trajectory_panel.pdf")
message("5. 16_Figure5_trajectory_status_core_gene_summary.csv")
message("6. 16_Figure5_gene_high_status_distribution_summary.csv")
message("7. 16_Figure5_status_source_summary.csv")
message("")
message("说明：")
message("主图不再使用pseudotime_proxy趋势作为核心结论。")
message("如果没有真实pseudotime，proxy趋势只作为supplementary/QC输出。")
message("============================================================\n")