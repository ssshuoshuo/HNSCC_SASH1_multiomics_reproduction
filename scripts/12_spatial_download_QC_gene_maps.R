# 12_spatial_download_QC_gene_maps.R

# 本脚本功能：
# 1. 下载GSE252265公开Visium空间转录组数据
# 2. 读取聚合后的filtered_feature_bc_matrix.h5
# 3. 读取组织spot坐标文件
# 4. 检查表达矩阵barcode与空间坐标barcode的匹配情况
# 5. 建立不含H&E图像的空间Seurat对象
# 6. 完成spot-level QC
# 7. 输出SASH1、COL1A1、EMP1、MYH11及基础QC指标的空间表达图
# 8. 为下一步空间排他性或邻近关系分析筛选合格样本
# 9. 保存空间Seurat对象和运行信息

# 本项目专用数据：
# GSE252265
# 8例舌癌患者的10x Genomics Visium空间转录组数据
#
# 本步骤先确保表达矩阵、spot坐标和样本标签可靠。
# 当前聚合坐标文件未提供明确样本列时，本步骤先标记为All_spots。
#
# 通用代码修改位置：
# 1. 换数据集时：
#    修改raw_dir和download_geo_file中的accession及filename
#
# 2. 换目标基因时：
#    修改target_genes
#
# 3. 后续获得可靠切片或样本标签后：
#    修改spatial_sample_id的赋值逻辑


# ============================================================
# A. 加载包
# ============================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "hdf5r",
  "dplyr",
  "tidyr",
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
      "\n\n请不要一次性更新全部R包。",
      "\n将缺少包名称和报错发回。"
    )
  )
}

library(Seurat)
library(SeuratObject)
library(Matrix)
library(hdf5r)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ============================================================
# B. 项目路径与文件夹
# ============================================================

project_dir <- normalizePath(
  "~/Desktop/HNSCC_SASH1_reproduction"
)

raw_dir <- file.path(
  project_dir,
  "data",
  "raw",
  "GSE252265"
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
  raw_dir,
  recursive = TRUE,
  showWarnings = FALSE
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

options(
  timeout = 1800
)

# ============================================================
# C. GEO官方FTP固定目录下载函数
# ============================================================

make_geo_ftp_url <- function(
    accession,
    filename
) {
  
  accession_bucket <- sub(
    "[0-9]{3}$",
    "nnn",
    accession
  )
  
  paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/",
    accession_bucket,
    "/",
    accession,
    "/suppl/",
    filename
  )
}

download_geo_file <- function(
    accession,
    filename,
    output_directory,
    minimum_size_mb
) {
  
  output_file <- file.path(
    output_directory,
    filename
  )
  
  temporary_file <- paste0(
    output_file,
    ".download"
  )
  
  minimum_size_bytes <- minimum_size_mb *
    1024^2
  
  if (file.exists(output_file)) {
    
    existing_size <- file.info(
      output_file
    )$size
    
    if (
      !is.na(existing_size) &&
      existing_size >= minimum_size_bytes
    ) {
      
      message(
        "文件已完整存在，跳过下载：",
        filename,
        " (",
        round(
          existing_size / 1024^2,
          2
        ),
        " MB)"
      )
      
      return(output_file)
    }
    
    message(
      "发现过小或残缺文件，删除后重新下载：",
      filename,
      " (当前大小 ",
      round(
        existing_size / 1024^2,
        3
      ),
      " MB)"
    )
    
    unlink(
      output_file,
      force = TRUE
    )
  }
  
  if (file.exists(temporary_file)) {
    
    message(
      "删除上次中断留下的临时文件：",
      basename(temporary_file)
    )
    
    unlink(
      temporary_file,
      force = TRUE
    )
  }
  
  download_url <- make_geo_ftp_url(
    accession = accession,
    filename = filename
  )
  
  message(
    "下载：",
    filename
  )
  
  message(
    "来源：",
    download_url
  )
  
  download_status <- tryCatch(
    {
      utils::download.file(
        url = download_url,
        destfile = temporary_file,
        method = "libcurl",
        mode = "wb",
        quiet = FALSE
      )
    },
    error = function(e) {
      message(
        "下载错误：",
        conditionMessage(e)
      )
      return(1)
    }
  )
  
  if (
    !identical(
      download_status,
      0L
    )
  ) {
    
    if (file.exists(temporary_file)) {
      unlink(
        temporary_file,
        force = TRUE
      )
    }
    
    stop(
      paste0(
        "下载失败：",
        filename,
        "\n请将本次下载输出完整发回。"
      )
    )
  }
  
  if (!file.exists(temporary_file)) {
    stop(
      paste0(
        "下载后未找到临时文件：",
        filename
      )
    )
  }
  
  downloaded_size <- file.info(
    temporary_file
  )$size
  
  if (
    is.na(downloaded_size) ||
    downloaded_size < minimum_size_bytes
  ) {
    
    unlink(
      temporary_file,
      force = TRUE
    )
    
    stop(
      paste0(
        "下载文件大小异常，疑似不完整：",
        filename,
        "\n实际大小：",
        round(
          downloaded_size / 1024^2,
          3
        ),
        " MB",
        "\n最低期望大小：",
        minimum_size_mb,
        " MB"
      )
    )
  }
  
  renamed_successfully <- file.rename(
    temporary_file,
    output_file
  )
  
  if (!renamed_successfully) {
    stop(
      paste0(
        "下载完成但无法将临时文件改名为正式文件：",
        filename
      )
    )
  }
  
  message(
    "下载完成：",
    filename,
    " (",
    round(
      downloaded_size / 1024^2,
      2
    ),
    " MB)"
  )
  
  return(output_file)
}

# ============================================================
# D. 下载GSE252265所需输入文件
# ============================================================

spatial_h5_file <- download_geo_file(
  accession = "GSE252265",
  filename = "GSE252265_filtered_feature_bc_matrix.h5",
  output_directory = raw_dir,
  minimum_size_mb = 60
)

tissue_positions_file_gz <- download_geo_file(
  accession = "GSE252265",
  filename = "GSE252265_aggr_tissue_positions.csv.gz",
  output_directory = raw_dir,
  minimum_size_mb = 0.4
)

aggregation_file_gz <- download_geo_file(
  accession = "GSE252265",
  filename = "GSE252265_aggregation.csv.gz",
  output_directory = raw_dir,
  minimum_size_mb = 0.0001
)

barcodes_file_gz <- download_geo_file(
  accession = "GSE252265",
  filename = "GSE252265_barcodes.tsv.gz",
  output_directory = raw_dir,
  minimum_size_mb = 0.05
)

features_file_gz <- download_geo_file(
  accession = "GSE252265",
  filename = "GSE252265_features.tsv.gz",
  output_directory = raw_dir,
  minimum_size_mb = 0.2
)

# ============================================================
# E. 保存输入文件清单
# ============================================================

input_file_paths <- c(
  spatial_h5_file,
  tissue_positions_file_gz,
  aggregation_file_gz,
  barcodes_file_gz,
  features_file_gz
)

input_file_inventory <- data.frame(
  file_name = basename(
    input_file_paths
  ),
  file_path = input_file_paths,
  file_size_mb = round(
    file.info(
      input_file_paths
    )$size / 1024^2,
    3
  ),
  stringsAsFactors = FALSE
)

write.csv(
  input_file_inventory,
  file.path(
    table_dir,
    "12_input_file_inventory.csv"
  ),
  row.names = FALSE
)

# ============================================================
# F. 读取10x H5表达矩阵
# ============================================================

raw_expression <- Seurat::Read10X_h5(
  filename = spatial_h5_file,
  use.names = TRUE,
  unique.features = TRUE
)

if (is.list(raw_expression)) {
  
  expression_layer_names <- names(
    raw_expression
  )
  
  if (
    "Gene Expression" %in%
    expression_layer_names
  ) {
    
    raw_expression <- raw_expression[["Gene Expression"]]
    
  } else {
    
    raw_expression <- raw_expression[[1]]
    
    message(
      "H5中未找到Gene Expression层，默认使用第一个层：",
      expression_layer_names[
        1
      ]
    )
  }
}

if (!inherits(raw_expression, "dgCMatrix")) {
  raw_expression <- as(
    raw_expression,
    "dgCMatrix"
  )
}

message(
  "表达矩阵基因数：",
  nrow(raw_expression)
)

message(
  "表达矩阵barcode数：",
  ncol(raw_expression)
)

# ============================================================
# G. 读取组织spot坐标
# ============================================================

read_spatial_positions <- function(
    file_path
) {
  
  positions_try_header <- tryCatch(
    read.csv(
      gzfile(file_path),
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) {
      NULL
    }
  )
  
  if (
    !is.null(positions_try_header) &&
    ncol(positions_try_header) >= 5
  ) {
    
    current_names <- tolower(
      colnames(positions_try_header)
    )
    
    has_barcode_name <- any(
      grepl(
        "barcode",
        current_names
      )
    )
    
    if (has_barcode_name) {
      return(positions_try_header)
    }
  }
  
  positions_no_header <- read.csv(
    gzfile(file_path),
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (ncol(positions_no_header) < 6) {
    stop(
      paste0(
        "组织spot坐标列数少于6，当前列数：",
        ncol(positions_no_header)
      )
    )
  }
  
  default_position_names <- c(
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres"
  )
  
  colnames(positions_no_header)[
    seq_len(
      length(default_position_names)
    )
  ] <- default_position_names
  
  return(positions_no_header)
}

spatial_positions <- read_spatial_positions(
  tissue_positions_file_gz
)

message(
  "组织spot坐标行数：",
  nrow(spatial_positions)
)

message(
  "组织spot坐标列：",
  paste(
    colnames(spatial_positions),
    collapse = ", "
  )
)

write.csv(
  spatial_positions,
  file.path(
    table_dir,
    "12_raw_tissue_positions_table.csv"
  ),
  row.names = FALSE
)

# ============================================================
# H. 统一barcode列与坐标列名称
# ============================================================

find_column_by_pattern <- function(
    data_frame,
    patterns,
    required = TRUE
) {
  
  lower_names <- tolower(
    colnames(data_frame)
  )
  
  matched_columns <- unlist(
    lapply(
      patterns,
      function(current_pattern) {
        which(
          grepl(
            current_pattern,
            lower_names
          )
        )
      }
    )
  )
  
  matched_columns <- unique(
    matched_columns
  )
  
  if (length(matched_columns) == 0) {
    
    if (required) {
      stop(
        paste0(
          "未找到所需列。搜索模式：",
          paste(
            patterns,
            collapse = ", "
          ),
          "\n当前列名：",
          paste(
            colnames(data_frame),
            collapse = ", "
          )
        )
      )
    }
    
    return(NA_character_)
  }
  
  return(
    colnames(data_frame)[
      matched_columns[
        1
      ]
    ]
  )
}

barcode_column <- find_column_by_pattern(
  spatial_positions,
  patterns = c(
    "^barcode$",
    "barcode"
  )
)

in_tissue_column <- find_column_by_pattern(
  spatial_positions,
  patterns = c(
    "^in_tissue$",
    "in.tissue"
  ),
  required = FALSE
)

row_coordinate_column <- find_column_by_pattern(
  spatial_positions,
  patterns = c(
    "pxl_row",
    "pixel_row",
    "image_row"
  ),
  required = FALSE
)

column_coordinate_column <- find_column_by_pattern(
  spatial_positions,
  patterns = c(
    "pxl_col",
    "pixel_col",
    "image_col"
  ),
  required = FALSE
)

if (
  is.na(row_coordinate_column) ||
  is.na(column_coordinate_column)
) {
  stop(
    paste0(
      "未能识别像素坐标列。",
      "\n当前列名：",
      paste(
        colnames(spatial_positions),
        collapse = ", "
      )
    )
  )
}

spatial_coordinates <- data.frame(
  barcode = as.character(
    spatial_positions[
      ,
      barcode_column
    ]
  ),
  spatial_x = as.numeric(
    spatial_positions[
      ,
      column_coordinate_column
    ]
  ),
  spatial_y = as.numeric(
    spatial_positions[
      ,
      row_coordinate_column
    ]
  ),
  stringsAsFactors = FALSE
)

if (!is.na(in_tissue_column)) {
  
  spatial_coordinates$in_tissue <- as.numeric(
    spatial_positions[
      ,
      in_tissue_column
    ]
  )
  
} else {
  
  spatial_coordinates$in_tissue <- 1
  
  message(
    "未识别in_tissue列，暂将所有spot视为组织spot。"
  )
}

# ============================================================
# I. 识别或设置空间样本标签
# ============================================================

sample_column <- find_column_by_pattern(
  spatial_positions,
  patterns = c(
    "^sample_id$",
    "^sample$",
    "^library_id$",
    "^library$",
    "orig.ident",
    "^slide$",
    "section"
  ),
  required = FALSE
)

if (!is.na(sample_column)) {
  
  spatial_coordinates$spatial_sample_id <- as.character(
    spatial_positions[
      ,
      sample_column
    ]
  )
  
  message(
    "识别到空间样本列：",
    sample_column
  )
  
} else {
  
  spatial_coordinates$spatial_sample_id <- "All_spots"
  
  message(
    "坐标文件中未识别到明确样本列。",
    "\n本步骤先保留All_spots；",
    "\n下一步将结合aggregation文件和原始JSON确认切片标签。"
  )
}

# ============================================================
# J. barcode匹配与安全检查
# ============================================================

expression_barcodes <- colnames(
  raw_expression
)

coordinate_barcodes <- spatial_coordinates$barcode

exact_match_number <- sum(
  expression_barcodes %in%
    coordinate_barcodes
)

exact_match_percent <- round(
  100 *
    exact_match_number /
    length(expression_barcodes),
  2
)

barcode_match_summary <- data.frame(
  expression_barcode_number = length(
    expression_barcodes
  ),
  coordinate_barcode_number = length(
    coordinate_barcodes
  ),
  exact_match_number = exact_match_number,
  exact_match_percent = exact_match_percent,
  stringsAsFactors = FALSE
)

write.csv(
  barcode_match_summary,
  file.path(
    table_dir,
    "12_barcode_coordinate_match_summary.csv"
  ),
  row.names = FALSE
)

message(
  "表达矩阵与空间坐标精确barcode匹配比例：",
  exact_match_percent,
  "%"
)

if (exact_match_percent < 80) {
  stop(
    paste0(
      "表达矩阵与坐标文件barcode精确匹配比例不足80%。",
      "\n当前为：",
      exact_match_percent,
      "%。",
      "\n\n请将以下两个文件发回：",
      "\n1. 12_barcode_coordinate_match_summary.csv",
      "\n2. 12_raw_tissue_positions_table.csv",
      "\n不要继续做空间表达图。"
    )
  )
}

# ============================================================
# K. 建立空间Seurat对象并加入坐标信息
# ============================================================

spatial_seurat <- CreateSeuratObject(
  counts = raw_expression,
  assay = "SpatialRNA",
  project = "GSE252265",
  min.cells = 0,
  min.features = 0
)

DefaultAssay(spatial_seurat) <- "SpatialRNA"

spatial_coordinates <- spatial_coordinates[
  match(
    colnames(spatial_seurat),
    spatial_coordinates$barcode
  ),
  ,
  drop = FALSE
]

if (
  any(
    is.na(
      spatial_coordinates$barcode
    )
  )
) {
  stop(
    "barcode匹配通过后仍出现未匹配坐标，停止。"
  )
}

rownames(spatial_coordinates) <- spatial_coordinates$barcode

spatial_seurat$spatial_x <- spatial_coordinates[
  colnames(spatial_seurat),
  "spatial_x"
]

spatial_seurat$spatial_y <- spatial_coordinates[
  colnames(spatial_seurat),
  "spatial_y"
]

spatial_seurat$in_tissue <- spatial_coordinates[
  colnames(spatial_seurat),
  "in_tissue"
]

spatial_seurat$spatial_sample_id <- spatial_coordinates[
  colnames(spatial_seurat),
  "spatial_sample_id"
]

# ============================================================
# L. 计算spot-level QC指标
# ============================================================

spatial_seurat[["percent.mt"]] <- PercentageFeatureSet(
  spatial_seurat,
  pattern = "^MT-"
)

spatial_seurat <- NormalizeData(
  spatial_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

spot_qc_table <- spatial_seurat@meta.data

spot_qc_table$barcode <- rownames(
  spot_qc_table
)

spot_qc_summary <- spot_qc_table %>%
  dplyr::group_by(
    spatial_sample_id
  ) %>%
  dplyr::summarise(
    spot_number = dplyr::n(),
    tissue_spot_number = sum(
      in_tissue == 1,
      na.rm = TRUE
    ),
    median_UMI = round(
      median(
        nCount_SpatialRNA
      ),
      2
    ),
    median_detected_gene = round(
      median(
        nFeature_SpatialRNA
      ),
      2
    ),
    median_percent_mt = round(
      median(
        percent.mt
      ),
      2
    ),
    .groups = "drop"
  )

write.csv(
  spot_qc_summary,
  file.path(
    table_dir,
    "12_spatial_spot_QC_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  spot_qc_table,
  file.path(
    table_dir,
    "12_spatial_spot_QC_cell_metadata.csv"
  ),
  row.names = FALSE
)

# ============================================================
# M. 组织spot筛选
# ============================================================

tissue_spot_barcodes <- rownames(
  spatial_seurat@meta.data
)[
  spatial_seurat$in_tissue == 1
]

if (length(tissue_spot_barcodes) < 100) {
  stop(
    paste0(
      "组织spot数过少：",
      length(tissue_spot_barcodes),
      "\n请检查in_tissue列及坐标文件。"
    )
  )
}

spatial_tissue <- subset(
  spatial_seurat,
  cells = tissue_spot_barcodes
)

message(
  "保留组织spot数：",
  ncol(spatial_tissue)
)

# ============================================================
# N. 构建空间绘图数据
# ============================================================

target_genes <- c(
  "SASH1",
  "COL1A1",
  "EMP1",
  "MYH11"
)

target_genes_found <- intersect(
  target_genes,
  rownames(spatial_tissue)
)

message(
  "空间数据中找到目标基因：",
  paste(
    target_genes_found,
    collapse = ", "
  )
)

if (
  !all(
    c(
      "SASH1",
      "COL1A1"
    ) %in%
    target_genes_found
  )
) {
  stop(
    "未同时找到SASH1和COL1A1，无法进入论文空间验证主线。"
  )
}

normalized_expression <- LayerData(
  object = spatial_tissue,
  assay = "SpatialRNA",
  layer = "data"
)

spatial_plot_data <- spatial_tissue@meta.data

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

spatial_plot_data$log10_UMI <- log10(
  spatial_plot_data$nCount_SpatialRNA + 1
)

spatial_plot_data$log10_detected_gene <- log10(
  spatial_plot_data$nFeature_SpatialRNA + 1
)

# ============================================================
# O. 空间绘图函数
# ============================================================

make_spatial_feature_plot <- function(
    plot_data,
    feature_name,
    plot_title
) {
  
  ggplot(
    plot_data,
    aes(
      x = spatial_x,
      y = spatial_y,
      color = .data[[feature_name]]
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
      title = plot_title,
      x = NULL,
      y = NULL,
      color = feature_name
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

make_spatial_qc_plot <- function(
    plot_data,
    feature_name,
    plot_title,
    legend_title
) {
  
  ggplot(
    plot_data,
    aes(
      x = spatial_x,
      y = spatial_y,
      color = .data[[feature_name]]
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
      title = plot_title,
      x = NULL,
      y = NULL,
      color = legend_title
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

# ============================================================
# P. 输出基础QC空间图
# ============================================================

p_spatial_umi <- make_spatial_qc_plot(
  plot_data = spatial_plot_data,
  feature_name = "log10_UMI",
  plot_title = "GSE252265 Spatial UMI Distribution",
  legend_title = "log10(UMI+1)"
)

p_spatial_feature <- make_spatial_qc_plot(
  plot_data = spatial_plot_data,
  feature_name = "log10_detected_gene",
  plot_title = "GSE252265 Spatial Detected Gene Distribution",
  legend_title = "log10(Detected genes+1)"
)

p_spatial_mt <- make_spatial_qc_plot(
  plot_data = spatial_plot_data,
  feature_name = "percent.mt",
  plot_title = "GSE252265 Spatial Mitochondrial Percentage",
  legend_title = "Percent mitochondrial"
)

ggsave(
  filename = file.path(
    figure_dir,
    "12_spatial_QC_UMI_distribution.pdf"
  ),
  plot = p_spatial_umi,
  width = 14,
  height = 9
)

ggsave(
  filename = file.path(
    figure_dir,
    "12_spatial_QC_detected_gene_distribution.pdf"
  ),
  plot = p_spatial_feature,
  width = 14,
  height = 9
)

ggsave(
  filename = file.path(
    figure_dir,
    "12_spatial_QC_percent_mt_distribution.pdf"
  ),
  plot = p_spatial_mt,
  width = 14,
  height = 9
)

# ============================================================
# Q. 输出SASH1、COL1A1、EMP1、MYH11空间图
# ============================================================

gene_plot_list <- list()

for (current_gene in target_genes_found) {
  
  gene_plot_list[[current_gene]] <- make_spatial_feature_plot(
    plot_data = spatial_plot_data,
    feature_name = current_gene,
    plot_title = paste0(
      current_gene,
      " Spatial Expression"
    )
  )
}

p_target_gene_spatial <- patchwork::wrap_plots(
  gene_plot_list,
  ncol = 2
) +
  patchwork::plot_annotation(
    title = "GSE252265 Spatial Expression of Core Candidate Genes"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "12_SASH1_COL1A1_EMP1_MYH11_spatial_expression.pdf"
  ),
  plot = p_target_gene_spatial,
  width = 16,
  height = 13
)

# ============================================================
# R. 初步SASH1与COL1A1共同表达审查表
# ============================================================

spatial_gene_summary <- spatial_plot_data %>%
  dplyr::group_by(
    spatial_sample_id
  ) %>%
  dplyr::summarise(
    tissue_spot_number = dplyr::n(),
    SASH1_detected_spot_number = sum(
      SASH1 > 0
    ),
    SASH1_detected_percent = round(
      100 *
        SASH1_detected_spot_number /
        tissue_spot_number,
      2
    ),
    COL1A1_detected_spot_number = sum(
      COL1A1 > 0
    ),
    COL1A1_detected_percent = round(
      100 *
        COL1A1_detected_spot_number /
        tissue_spot_number,
      2
    ),
    SASH1_COL1A1_co_detected_spot_number = sum(
      SASH1 > 0 &
        COL1A1 > 0
    ),
    SASH1_COL1A1_co_detected_percent = round(
      100 *
        SASH1_COL1A1_co_detected_spot_number /
        tissue_spot_number,
      2
    ),
    .groups = "drop"
  )

write.csv(
  spatial_gene_summary,
  file.path(
    table_dir,
    "12_SASH1_COL1A1_spot_detection_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  spatial_plot_data,
  file.path(
    table_dir,
    "12_spatial_tissue_spot_expression_metadata.csv"
  ),
  row.names = FALSE
)

# ============================================================
# S. 保存对象与运行信息
# ============================================================

saveRDS(
  spatial_tissue,
  file.path(
    object_dir,
    "12_spatial_tissue_spots_Seurat.rds"
  ),
  compress = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    table_dir,
    "12_sessionInfo.txt"
  )
)

# ============================================================
# T. 输出检查与完成提示
# ============================================================

required_output_files <- c(
  file.path(
    object_dir,
    "12_spatial_tissue_spots_Seurat.rds"
  ),
  file.path(
    table_dir,
    "12_input_file_inventory.csv"
  ),
  file.path(
    table_dir,
    "12_raw_tissue_positions_table.csv"
  ),
  file.path(
    table_dir,
    "12_barcode_coordinate_match_summary.csv"
  ),
  file.path(
    table_dir,
    "12_spatial_spot_QC_summary.csv"
  ),
  file.path(
    table_dir,
    "12_SASH1_COL1A1_spot_detection_summary.csv"
  ),
  file.path(
    figure_dir,
    "12_spatial_QC_UMI_distribution.pdf"
  ),
  file.path(
    figure_dir,
    "12_spatial_QC_detected_gene_distribution.pdf"
  ),
  file.path(
    figure_dir,
    "12_spatial_QC_percent_mt_distribution.pdf"
  ),
  file.path(
    figure_dir,
    "12_SASH1_COL1A1_EMP1_MYH11_spatial_expression.pdf"
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
    "12_output_file_check.csv"
  ),
  row.names = FALSE
)

print(
  output_status
)

message("\n============================================================")
message("12 GSE252265空间转录组下载、QC和核心基因空间图完成。")
message("")
message("重点查看：")
message("1. 12_spatial_spot_QC_summary.csv")
message("2. 12_SASH1_COL1A1_spot_detection_summary.csv")
message("3. 12_spatial_QC_UMI_distribution.pdf")
message("4. 12_SASH1_COL1A1_EMP1_MYH11_spatial_expression.pdf")
message("")
message("下一步将根据合格样本，正式进行：")
message("SASH1与COL1A1-high纤维化区域的空间排他或邻近统计分析。")
message("============================================================\n")