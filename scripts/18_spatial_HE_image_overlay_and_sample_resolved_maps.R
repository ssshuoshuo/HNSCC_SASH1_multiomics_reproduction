# 18_spatial_HE_image_overlay_and_sample_resolved_maps.R

# 本脚本功能：
# 1. 读取GSE252265聚合空间转录组表达矩阵
# 2. 根据barcode后缀识别8个真实空间样本
# 3. 将gem_id映射到UH8、UH12、UH17、UH19、UH20、UH21、UH22和UH24
# 4. 解压并读取每个样本对应的H&E组织图像
# 5. 读取每个样本对应的Space Ranger scalefactor
# 6. 将full-resolution空间坐标转换为low-resolution图像坐标
# 7. 生成逐样本H&E+组织内spot覆盖图
# 8. 生成SASH1、COL1A1、EMP1和MYH11逐样本H&E表达叠加图
# 9. 保存带真实spatial_sample_id的Seurat对象
# 10. 输出barcode、表达矩阵、坐标、样本和图像的完整匹配检查表
#
# 本项目专用空间数据：
# GSE252265
#
# 8个空间转录组样本：
# UH8、UH12、UH17、UH19、
# UH20、UH21、UH22、UH24
#
# gem_id与样本映射：
# 1 -> UH8
# 2 -> UH12
# 3 -> UH17
# 4 -> UH19
# 5 -> UH20
# 6 -> UH21
# 7 -> UH22
# 8 -> UH24
#
# 重要说明：
# 1. 本脚本使用low-resolution H&E图进行绘图，以降低内存消耗
# 2. 表达矩阵和组织内坐标已经逐样本验证完全匹配
# 3. barcode后缀-1至-8对应聚合数据中的gem_id
# 4. 本脚本不会修改12–14阶段的已有结果
# 5. 新对象和新图统一使用18开头命名
#
# 通用代码修改位置：
# 1. 换项目路径时：
#    修改project_dir
#
# 2. 换空间数据集时：
#    修改spatial_h5_file、position_file和raw_extracted_dir
#
# 3. 换样本与gem_id映射时：
#    修改sample_mapping
#
# 4. 换目标基因时：
#    修改target_genes
#
# 5. 调整spot大小时：
#    修改plot_params$spot_size
#
# 6. 调整表达图透明度时：
#    修改plot_params$expression_alpha
#
# ============================================================
# A. 项目路径与输出文件夹
# ============================================================

# getwd()应当是R Project根目录。
# 为避免从错误目录运行，本项目默认使用绝对路径。
# 换电脑或移动项目文件夹时，只需要修改这里。

project_dir <- "/Users/yaoshuo/Desktop/HNSCC_SASH1_reproduction"

data_dir <- file.path(
  project_dir,
  "data"
)

spatial_raw_dir <- file.path(
  data_dir,
  "raw",
  "GSE252265"
)

raw_extracted_dir <- file.path(
  spatial_raw_dir,
  "raw_extracted"
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

object_dir <- file.path(
  project_dir,
  "results",
  "objects"
)

he_figure_dir <- file.path(
  figure_dir,
  "18_spatial_HE"
)

he_spot_figure_dir <- file.path(
  he_figure_dir,
  "spot_overlays"
)

he_gene_figure_dir <- file.path(
  he_figure_dir,
  "core_gene_overlays"
)

processed_image_dir <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE252265_HE_images"
)

required_dirs <- c(
  table_dir,
  figure_dir,
  object_dir,
  he_figure_dir,
  he_spot_figure_dir,
  he_gene_figure_dir,
  processed_image_dir
)

for (path_use in required_dirs) {
  dir.create(
    path_use,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ============================================================
# B. 输入文件与分析参数
# ============================================================

# -----------------------------
# B1. 空间表达矩阵和坐标文件
# -----------------------------

spatial_h5_file <- file.path(
  spatial_raw_dir,
  "GSE252265_filtered_feature_bc_matrix.h5"
)

position_file <- file.path(
  spatial_raw_dir,
  "GSE252265_aggr_tissue_positions.csv.gz"
)

aggregation_file <- file.path(
  spatial_raw_dir,
  "GSE252265_aggregation.csv.gz"
)

input_object_candidates <- c(
  file.path(
    object_dir,
    "14_spatial_domain_annotated_paper_style_Seurat.rds"
  ),
  file.path(
    object_dir,
    "13_spatial_SASH1_COL1A1_high_annotated_Seurat.rds"
  ),
  file.path(
    object_dir,
    "12_spatial_tissue_spots_Seurat.rds"
  )
)

# -----------------------------
# B2. 样本与gem_id映射
# -----------------------------

sample_mapping <- data.frame(
  gem_id = as.character(1:8),
  spatial_sample_id = c(
    "UH8",
    "UH12",
    "UH17",
    "UH19",
    "UH20",
    "UH21",
    "UH22",
    "UH24"
  ),
  gsm_id = c(
    "GSM7998252",
    "GSM7998253",
    "GSM7998254",
    "GSM7998255",
    "GSM7998256",
    "GSM7998257",
    "GSM7998258",
    "GSM7998259"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# B3. 目标基因
# -----------------------------

target_genes <- c(
  "SASH1",
  "COL1A1",
  "EMP1",
  "MYH11"
)

# -----------------------------
# B4. 绘图参数
# -----------------------------

plot_params <- list(
  spot_size = 0.55,
  spot_alpha = 0.70,
  expression_alpha = 0.85,
  background_spot_size = 0.25,
  background_spot_alpha = 0.15,
  panel_width = 5.5,
  panel_height = 5.5,
  dpi = 200L
)

# ============================================================
# C. 安装与加载R包
# ============================================================

options(
  repos = c(
    CRAN = "https://cloud.r-project.org"
  ),
  timeout = 3600
)

cran_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "patchwork",
  "png",
  "jsonlite",
  "scales"
)

for (pkg in cran_packages) {
  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {
    message(
      "正在安装R包：",
      pkg
    )
    
    install.packages(
      pkg
    )
  } else {
    message(
      "已安装，跳过：",
      pkg
    )
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(png)
  library(jsonlite)
  library(scales)
})

# ============================================================
# D. 通用辅助函数
# ============================================================

# -----------------------------
# D1. 文件检查
# -----------------------------

assert_file_exists <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(
      "缺少输入文件：",
      file_path
    )
  }
  
  invisible(TRUE)
}

assert_directory_exists <- function(directory_path) {
  if (!dir.exists(directory_path)) {
    stop(
      "缺少输入文件夹：",
      directory_path
    )
  }
  
  invisible(TRUE)
}

first_existing_file <- function(paths) {
  existing <- paths[
    file.exists(paths)
  ]
  
  if (length(existing) == 0) {
    return(NA_character_)
  }
  
  existing[1]
}

# -----------------------------
# D2. UTF-8输出函数
# -----------------------------

write_utf8_csv <- function(
    x,
    file_name
) {
  write.csv(
    x,
    file = file_name,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

write_utf8_text <- function(
    lines,
    file_name
) {
  con <- file(
    file_name,
    open = "w",
    encoding = "UTF-8"
  )
  
  on.exit(
    close(con),
    add = TRUE
  )
  
  writeLines(
    lines,
    con = con
  )
}

# -----------------------------
# D3. 解压单个gzip文件
# -----------------------------

gunzip_to_file <- function(
    input_file,
    output_file
) {
  if (
    file.exists(output_file) &&
    file.info(output_file)$size > 0
  ) {
    return(output_file)
  }
  
  input_connection <- gzfile(
    input_file,
    open = "rb"
  )
  
  output_connection <- file(
    output_file,
    open = "wb"
  )
  
  on.exit(
    close(input_connection),
    add = TRUE
  )
  
  on.exit(
    close(output_connection),
    add = TRUE
  )
  
  repeat {
    buffer <- readBin(
      input_connection,
      what = "raw",
      n = 1024 * 1024
    )
    
    if (length(buffer) == 0) {
      break
    }
    
    writeBin(
      buffer,
      output_connection
    )
  }
  
  output_file
}

# -----------------------------
# D4. 提取barcode的gem_id
# -----------------------------

extract_gem_id <- function(barcode_vector) {
  stringr::str_extract(
    barcode_vector,
    "(?<=-)[0-9]+$"
  )
}

# -----------------------------
# D5. 读取Seurat counts矩阵
# -----------------------------

get_count_matrix <- function(
    seurat_object,
    assay_name = "Spatial"
) {
  tryCatch(
    SeuratObject::GetAssayData(
      seurat_object,
      assay = assay_name,
      layer = "counts"
    ),
    error = function(e) {
      SeuratObject::GetAssayData(
        seurat_object,
        assay = assay_name,
        slot = "counts"
      )
    }
  )
}

# -----------------------------
# D6. 读取Seurat normalized data
# -----------------------------

get_normalized_matrix <- function(
    seurat_object,
    assay_name
) {
  tryCatch(
    SeuratObject::GetAssayData(
      seurat_object,
      assay = assay_name,
      layer = "data"
    ),
    error = function(e) {
      SeuratObject::GetAssayData(
        seurat_object,
        assay = assay_name,
        slot = "data"
      )
    }
  )
}

# -----------------------------
# D7. 读取样本H&E和scalefactor
# -----------------------------

read_sample_image_information <- function(
    sample_id,
    gsm_id,
    raw_directory,
    processed_directory
) {
  lowres_pattern <- paste0(
    "^",
    gsm_id,
    "_",
    sample_id,
    "_tissue_lowres_image\\.png\\.gz$"
  )
  
  hires_pattern <- paste0(
    "^",
    gsm_id,
    "_",
    sample_id,
    "_tissue_hires_image\\.png\\.gz$"
  )
  
  scale_pattern <- paste0(
    "^",
    gsm_id,
    "_",
    sample_id,
    "_scalefactors_json\\.json\\.gz$"
  )
  
  lowres_gz <- list.files(
    raw_directory,
    pattern = lowres_pattern,
    full.names = TRUE
  )
  
  hires_gz <- list.files(
    raw_directory,
    pattern = hires_pattern,
    full.names = TRUE
  )
  
  scale_gz <- list.files(
    raw_directory,
    pattern = scale_pattern,
    full.names = TRUE
  )
  
  if (length(lowres_gz) != 1) {
    stop(
      sample_id,
      "没有找到唯一lowres图像。检测数量：",
      length(lowres_gz)
    )
  }
  
  if (length(hires_gz) != 1) {
    stop(
      sample_id,
      "没有找到唯一hires图像。检测数量：",
      length(hires_gz)
    )
  }
  
  if (length(scale_gz) != 1) {
    stop(
      sample_id,
      "没有找到唯一scalefactor文件。检测数量：",
      length(scale_gz)
    )
  }
  
  sample_output_dir <- file.path(
    processed_directory,
    sample_id
  )
  
  dir.create(
    sample_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  lowres_png <- file.path(
    sample_output_dir,
    paste0(
      sample_id,
      "_tissue_lowres_image.png"
    )
  )
  
  hires_png <- file.path(
    sample_output_dir,
    paste0(
      sample_id,
      "_tissue_hires_image.png"
    )
  )
  
  scale_json <- file.path(
    sample_output_dir,
    paste0(
      sample_id,
      "_scalefactors_json.json"
    )
  )
  
  gunzip_to_file(
    lowres_gz,
    lowres_png
  )
  
  gunzip_to_file(
    hires_gz,
    hires_png
  )
  
  gunzip_to_file(
    scale_gz,
    scale_json
  )
  
  scale_information <- jsonlite::fromJSON(
    scale_json
  )
  
  lowres_image <- png::readPNG(
    lowres_png
  )
  
  image_height <- dim(
    lowres_image
  )[1]
  
  image_width <- dim(
    lowres_image
  )[2]
  
  list(
    sample_id = sample_id,
    gsm_id = gsm_id,
    lowres_gz = lowres_gz,
    hires_gz = hires_gz,
    scale_gz = scale_gz,
    lowres_png = lowres_png,
    hires_png = hires_png,
    scale_json = scale_json,
    lowres_image = lowres_image,
    lowres_scale = as.numeric(
      scale_information$tissue_lowres_scalef
    ),
    hires_scale = as.numeric(
      scale_information$tissue_hires_scalef
    ),
    spot_diameter_fullres = as.numeric(
      scale_information$spot_diameter_fullres
    ),
    fiducial_diameter_fullres = as.numeric(
      scale_information$fiducial_diameter_fullres
    ),
    json_gem_id = as.character(
      scale_information$gem_id
    ),
    image_width = image_width,
    image_height = image_height
  )
}

# -----------------------------
# D8. 绘制H&E背景
# -----------------------------

create_he_base_plot <- function(
    sample_table,
    image_information,
    plot_title
) {
  ggplot2::ggplot() +
    ggplot2::annotation_raster(
      raster = image_information$lowres_image,
      xmin = 0,
      xmax = image_information$image_width,
      ymin = image_information$image_height,
      ymax = 0
    ) +
    ggplot2::coord_fixed(
      xlim = c(
        0,
        image_information$image_width
      ),
      ylim = c(
        image_information$image_height,
        0
      ),
      expand = FALSE,
      clip = "on"
    ) +
    ggplot2::labs(
      title = plot_title,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 12
      ),
      plot.margin = ggplot2::margin(
        2,
        2,
        2,
        2
      )
    )
}

# -----------------------------
# D9. 绘制H&E+spot覆盖图
# -----------------------------

create_spot_overlay_plot <- function(
    sample_table,
    image_information,
    sample_id
) {
  base_plot <- create_he_base_plot(
    sample_table = sample_table,
    image_information = image_information,
    plot_title = paste0(
      sample_id,
      " H&E + tissue spots"
    )
  )
  
  base_plot +
    ggplot2::geom_point(
      data = sample_table,
      mapping = ggplot2::aes(
        x = .data$image_x_lowres,
        y = .data$image_y_lowres
      ),
      shape = 21,
      size = plot_params$spot_size,
      stroke = 0,
      fill = "#1F78B4",
      alpha = plot_params$spot_alpha
    )
}

# -----------------------------
# D10. 绘制H&E+基因表达图
# -----------------------------

create_gene_overlay_plot <- function(
    sample_table,
    image_information,
    sample_id,
    gene_name
) {
  expression_column <- paste0(
    gene_name,
    "_expression"
  )
  
  expression_values <- sample_table[[expression_column]
  ]
  
  positive_values <- expression_values[
    is.finite(expression_values) &
      expression_values > 0
  ]
  
  if (length(positive_values) > 0) {
    upper_limit <- as.numeric(
      stats::quantile(
        positive_values,
        probs = 0.99,
        na.rm = TRUE
      )
    )
  } else {
    upper_limit <- 1
  }
  
  if (
    !is.finite(upper_limit) ||
    upper_limit <= 0
  ) {
    upper_limit <- 1
  }
  
  plotting_table <- sample_table %>%
    dplyr::mutate(
      expression_plot = pmin(
        .data[[expression_column]
        ],
        upper_limit
      )
    )
  
  base_plot <- create_he_base_plot(
    sample_table = sample_table,
    image_information = image_information,
    plot_title = paste0(
      sample_id,
      " ",
      gene_name
    )
  )
  
  base_plot +
    ggplot2::geom_point(
      data = plotting_table,
      mapping = ggplot2::aes(
        x = .data$image_x_lowres,
        y = .data$image_y_lowres
      ),
      size = plot_params$background_spot_size,
      color = "grey85",
      alpha = plot_params$background_spot_alpha
    ) +
    ggplot2::geom_point(
      data = plotting_table %>%
        dplyr::filter(
          .data$expression_plot > 0
        ),
      mapping = ggplot2::aes(
        x = .data$image_x_lowres,
        y = .data$image_y_lowres,
        color = .data$expression_plot
      ),
      size = plot_params$spot_size,
      alpha = plot_params$expression_alpha
    ) +
    ggplot2::scale_color_viridis_c(
      option = "magma",
      name = gene_name,
      limits = c(
        0,
        upper_limit
      ),
      oob = scales::squish
    ) +
    ggplot2::theme(
      legend.position = "right",
      legend.key.height = grid::unit(
        0.7,
        "cm"
      )
    )
}

# ============================================================
# E. 检查输入文件
# ============================================================

assert_file_exists(
  spatial_h5_file
)

assert_file_exists(
  position_file
)

assert_file_exists(
  aggregation_file
)

assert_directory_exists(
  raw_extracted_dir
)

write_utf8_csv(
  sample_mapping,
  file.path(
    table_dir,
    "18_gem_id_sample_mapping.csv"
  )
)

# ============================================================
# F. 读取坐标表和表达矩阵
# ============================================================

# -----------------------------
# F1. 读取聚合坐标表
# -----------------------------

position_table <- read.csv(
  position_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_position_columns <- c(
  "barcode",
  "in_tissue",
  "array_row",
  "array_col",
  "pxl_row_in_fullres",
  "pxl_col_in_fullres"
)

missing_position_columns <- setdiff(
  required_position_columns,
  colnames(position_table)
)

if (length(missing_position_columns) > 0) {
  stop(
    "坐标表缺少列：",
    paste(
      missing_position_columns,
      collapse = ", "
    )
  )
}

position_table <- position_table %>%
  dplyr::mutate(
    barcode = as.character(
      .data$barcode
    ),
    gem_id = extract_gem_id(
      .data$barcode
    )
  ) %>%
  dplyr::left_join(
    sample_mapping,
    by = "gem_id"
  )

if (
  any(
    is.na(
      position_table$spatial_sample_id
    )
  )
) {
  stop(
    "部分坐标barcode无法映射到空间样本。"
  )
}

tissue_position_table <- position_table %>%
  dplyr::filter(
    .data$in_tissue == 1
  )

# -----------------------------
# F2. 读取H5表达矩阵
# -----------------------------

raw_expression <- Seurat::Read10X_h5(
  spatial_h5_file
)

if (is.list(raw_expression)) {
  if (
    "Gene Expression" %in%
    names(raw_expression)
  ) {
    raw_expression <- raw_expression[["Gene Expression"]
    ]
  } else {
    raw_expression <- raw_expression[[1]
    ]
  }
}

if (
  !inherits(
    raw_expression,
    "Matrix"
  )
) {
  raw_expression <- Matrix::Matrix(
    raw_expression,
    sparse = TRUE
  )
}

expression_barcodes <- colnames(
  raw_expression
)

# -----------------------------
# F3. 验证表达矩阵和坐标匹配
# -----------------------------

expression_only_barcodes <- setdiff(
  expression_barcodes,
  tissue_position_table$barcode
)

coordinate_only_barcodes <- setdiff(
  tissue_position_table$barcode,
  expression_barcodes
)

matched_barcodes <- intersect(
  expression_barcodes,
  tissue_position_table$barcode
)

barcode_match_summary <- data.frame(
  metric = c(
    "Expression_spots",
    "In_tissue_coordinate_spots",
    "Matched_spots",
    "Expression_only_spots",
    "Coordinate_only_spots"
  ),
  value = c(
    length(expression_barcodes),
    nrow(tissue_position_table),
    length(matched_barcodes),
    length(expression_only_barcodes),
    length(coordinate_only_barcodes)
  ),
  stringsAsFactors = FALSE
)

print(
  barcode_match_summary
)

write_utf8_csv(
  barcode_match_summary,
  file.path(
    table_dir,
    "18_barcode_expression_coordinate_match_summary.csv"
  )
)

if (
  length(expression_only_barcodes) > 0 ||
  length(coordinate_only_barcodes) > 0
) {
  stop(
    "表达矩阵和组织内坐标未完全匹配。"
  )
}

# ============================================================
# G. 读取或建立Seurat对象
# ============================================================

existing_object_file <- first_existing_file(
  input_object_candidates
)

if (!is.na(existing_object_file)) {
  message(
    "读取已有空间Seurat对象：",
    existing_object_file
  )
  
  spatial_object <- readRDS(
    existing_object_file
  )
  
  common_cells <- intersect(
    colnames(spatial_object),
    matched_barcodes
  )
  
  if (
    length(common_cells) !=
    length(matched_barcodes)
  ) {
    warning(
      "已有Seurat对象的barcode与H5对象不完全一致，",
      "将从H5重新建立空间对象。"
    )
    
    spatial_object <- Seurat::CreateSeuratObject(
      counts = raw_expression[
        ,
        matched_barcodes,
        drop = FALSE
      ],
      assay = "Spatial",
      project = "GSE252265"
    )
  } else {
    spatial_object <- subset(
      spatial_object,
      cells = matched_barcodes
    )
  }
} else {
  message(
    "没有检测到12–14阶段空间对象，",
    "从H5重新建立Seurat对象。"
  )
  
  spatial_object <- Seurat::CreateSeuratObject(
    counts = raw_expression[
      ,
      matched_barcodes,
      drop = FALSE
    ],
    assay = "Spatial",
    project = "GSE252265"
  )
}

if (
  "Spatial" %in%
  names(spatial_object@assays)
) {
  assay_to_use <- "Spatial"
} else if (
  "RNA" %in%
  names(spatial_object@assays)
) {
  assay_to_use <- "RNA"
} else {
  assay_to_use <- Seurat::DefaultAssay(
    spatial_object
  )
}

Seurat::DefaultAssay(
  spatial_object
) <- assay_to_use

# -----------------------------
# G1. 加入真实样本metadata
# -----------------------------

metadata_table <- tissue_position_table %>%
  dplyr::filter(
    .data$barcode %in%
      colnames(spatial_object)
  ) %>%
  dplyr::arrange(
    match(
      .data$barcode,
      colnames(spatial_object)
    )
  )

if (
  !identical(
    metadata_table$barcode,
    colnames(spatial_object)
  )
) {
  stop(
    "metadata排序与Seurat对象barcode顺序不一致。"
  )
}

spatial_object$gem_id <- metadata_table$gem_id
spatial_object$spatial_sample_id <-
  metadata_table$spatial_sample_id
spatial_object$spatial_gsm_id <-
  metadata_table$gsm_id
spatial_object$array_row <-
  metadata_table$array_row
spatial_object$array_col <-
  metadata_table$array_col
spatial_object$pxl_row_in_fullres <-
  metadata_table$pxl_row_in_fullres
spatial_object$pxl_col_in_fullres <-
  metadata_table$pxl_col_in_fullres

# -----------------------------
# G2. 确保存在normalized data
# -----------------------------

normalized_matrix <- tryCatch(
  get_normalized_matrix(
    spatial_object,
    assay_to_use
  ),
  error = function(e) NULL
)

if (
  is.null(normalized_matrix) ||
  nrow(normalized_matrix) == 0 ||
  ncol(normalized_matrix) == 0
) {
  spatial_object <- Seurat::NormalizeData(
    spatial_object,
    assay = assay_to_use,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
}

normalized_matrix <- get_normalized_matrix(
  spatial_object,
  assay_to_use
)

# -----------------------------
# G3. 检查目标基因
# -----------------------------

genes_found <- target_genes[
  target_genes %in%
    rownames(normalized_matrix)
]

genes_missing <- setdiff(
  target_genes,
  genes_found
)

gene_check <- data.frame(
  gene = target_genes,
  found = target_genes %in%
    genes_found,
  stringsAsFactors = FALSE
)

write_utf8_csv(
  gene_check,
  file.path(
    table_dir,
    "18_core_gene_check.csv"
  )
)

if (length(genes_found) == 0) {
  stop(
    "4个目标基因均未在表达矩阵中检测到。"
  )
}

if (length(genes_missing) > 0) {
  warning(
    "以下基因未检测到：",
    paste(
      genes_missing,
      collapse = ", "
    )
  )
}

# ============================================================
# H. 读取8个样本的H&E图像
# ============================================================

image_information_list <- list()
image_inventory_rows <- list()

for (
  sample_index in
  seq_len(
    nrow(sample_mapping)
  )
) {
  current_sample <- sample_mapping[
    sample_index,
    ,
    drop = FALSE
  ]
  
  sample_id <- current_sample$spatial_sample_id
  gsm_id <- current_sample$gsm_id
  gem_id <- current_sample$gem_id
  
  image_information <- read_sample_image_information(
    sample_id = sample_id,
    gsm_id = gsm_id,
    raw_directory = raw_extracted_dir,
    processed_directory = processed_image_dir
  )
  
  if (
    image_information$json_gem_id !=
    gem_id
  ) {
    stop(
      sample_id,
      "的JSON gem_id与样本映射不一致。JSON=",
      image_information$json_gem_id,
      "，映射表=",
      gem_id
    )
  }
  
  image_information_list[[sample_id]
  ] <- image_information
  
  image_inventory_rows[[sample_id]
  ] <- data.frame(
    spatial_sample_id = sample_id,
    gsm_id = gsm_id,
    gem_id = gem_id,
    lowres_image_file =
      image_information$lowres_png,
    hires_image_file =
      image_information$hires_png,
    scalefactor_file =
      image_information$scale_json,
    lowres_scale =
      image_information$lowres_scale,
    hires_scale =
      image_information$hires_scale,
    image_width =
      image_information$image_width,
    image_height =
      image_information$image_height,
    spot_diameter_fullres =
      image_information$spot_diameter_fullres,
    json_gem_id =
      image_information$json_gem_id,
    stringsAsFactors = FALSE
  )
}

image_inventory <- dplyr::bind_rows(
  image_inventory_rows
)

write_utf8_csv(
  image_inventory,
  file.path(
    table_dir,
    "18_spatial_HE_image_inventory.csv"
  )
)

# ============================================================
# I. 建立逐spot表达和图像坐标表
# ============================================================

plotting_table <- metadata_table

for (gene_name in genes_found) {
  plotting_table[[paste0(
      gene_name,
      "_expression"
    )]
  ] <- as.numeric(
    normalized_matrix[
      gene_name,
      plotting_table$barcode
    ]
  )
}

plotting_table$image_x_lowres <- NA_real_
plotting_table$image_y_lowres <- NA_real_
plotting_table$image_x_hires <- NA_real_
plotting_table$image_y_hires <- NA_real_

for (sample_id in sample_mapping$spatial_sample_id) {
  current_rows <- which(
    plotting_table$spatial_sample_id ==
      sample_id
  )
  
  current_image_information <-
    image_information_list[[sample_id]
    ]
  
  plotting_table$image_x_lowres[
    current_rows
  ] <- plotting_table$pxl_col_in_fullres[
    current_rows
  ] *
    current_image_information$lowres_scale
  
  plotting_table$image_y_lowres[
    current_rows
  ] <- plotting_table$pxl_row_in_fullres[
    current_rows
  ] *
    current_image_information$lowres_scale
  
  plotting_table$image_x_hires[
    current_rows
  ] <- plotting_table$pxl_col_in_fullres[
    current_rows
  ] *
    current_image_information$hires_scale
  
  plotting_table$image_y_hires[
    current_rows
  ] <- plotting_table$pxl_row_in_fullres[
    current_rows
  ] *
    current_image_information$hires_scale
}

write_utf8_csv(
  plotting_table,
  file.path(
    table_dir,
    "18_spatial_sample_resolved_expression_coordinates.csv"
  )
)

# ============================================================
# J. 坐标边界和样本spot数检查
# ============================================================

coordinate_check_rows <- list()

for (sample_id in sample_mapping$spatial_sample_id) {
  sample_table <- plotting_table %>%
    dplyr::filter(
      .data$spatial_sample_id ==
        sample_id
    )
  
  image_information <-
    image_information_list[[sample_id]
    ]
  
  coordinate_check_rows[[sample_id]
  ] <- data.frame(
    spatial_sample_id = sample_id,
    gem_id = unique(
      sample_table$gem_id
    ),
    tissue_spots = nrow(
      sample_table
    ),
    min_image_x = min(
      sample_table$image_x_lowres,
      na.rm = TRUE
    ),
    max_image_x = max(
      sample_table$image_x_lowres,
      na.rm = TRUE
    ),
    min_image_y = min(
      sample_table$image_y_lowres,
      na.rm = TRUE
    ),
    max_image_y = max(
      sample_table$image_y_lowres,
      na.rm = TRUE
    ),
    image_width =
      image_information$image_width,
    image_height =
      image_information$image_height,
    all_x_inside_image = all(
      sample_table$image_x_lowres >= 0 &
        sample_table$image_x_lowres <=
        image_information$image_width
    ),
    all_y_inside_image = all(
      sample_table$image_y_lowres >= 0 &
        sample_table$image_y_lowres <=
        image_information$image_height
    ),
    stringsAsFactors = FALSE
  )
}

coordinate_check <- dplyr::bind_rows(
  coordinate_check_rows
)

print(
  coordinate_check
)

write_utf8_csv(
  coordinate_check,
  file.path(
    table_dir,
    "18_sample_spot_and_image_coordinate_check.csv"
  )
)

if (
  any(
    !coordinate_check$all_x_inside_image
  ) ||
  any(
    !coordinate_check$all_y_inside_image
  )
) {
  warning(
    "部分spot坐标超出lowres图像边界，",
    "请重点检查18_sample_spot_and_image_coordinate_check.csv。"
  )
}

# ============================================================
# K. 生成逐样本H&E+spot覆盖图
# ============================================================

spot_plot_list <- list()

for (sample_id in sample_mapping$spatial_sample_id) {
  sample_table <- plotting_table %>%
    dplyr::filter(
      .data$spatial_sample_id ==
        sample_id
    )
  
  image_information <-
    image_information_list[[sample_id]
    ]
  
  current_plot <- create_spot_overlay_plot(
    sample_table = sample_table,
    image_information = image_information,
    sample_id = sample_id
  )
  
  spot_plot_list[[sample_id]
  ] <- current_plot
  
  ggplot2::ggsave(
    filename = file.path(
      he_spot_figure_dir,
      paste0(
        "18_",
        sample_id,
        "_HE_tissue_spot_overlay.pdf"
      )
    ),
    plot = current_plot,
    width = plot_params$panel_width,
    height = plot_params$panel_height
  )
  
  ggplot2::ggsave(
    filename = file.path(
      he_spot_figure_dir,
      paste0(
        "18_",
        sample_id,
        "_HE_tissue_spot_overlay.png"
      )
    ),
    plot = current_plot,
    width = plot_params$panel_width,
    height = plot_params$panel_height,
    dpi = plot_params$dpi
  )
}

all_spot_panel <- patchwork::wrap_plots(
  spot_plot_list,
  ncol = 4
)

ggplot2::ggsave(
  filename = file.path(
    he_figure_dir,
    "18_all_samples_HE_tissue_spot_overlay_panel.pdf"
  ),
  plot = all_spot_panel,
  width = 20,
  height = 10
)

ggplot2::ggsave(
  filename = file.path(
    he_figure_dir,
    "18_all_samples_HE_tissue_spot_overlay_panel.png"
  ),
  plot = all_spot_panel,
  width = 20,
  height = 10,
  dpi = plot_params$dpi
)

# ============================================================
# L. 生成逐样本核心基因H&E叠加图
# ============================================================

gene_sample_plot_list <- list()

for (gene_name in genes_found) {
  gene_sample_plot_list[[gene_name]
  ] <- list()
  
  for (sample_id in sample_mapping$spatial_sample_id) {
    sample_table <- plotting_table %>%
      dplyr::filter(
        .data$spatial_sample_id ==
          sample_id
      )
    
    image_information <-
      image_information_list[[sample_id]
      ]
    
    current_plot <- create_gene_overlay_plot(
      sample_table = sample_table,
      image_information = image_information,
      sample_id = sample_id,
      gene_name = gene_name
    )
    
    gene_sample_plot_list[[gene_name]
    ][[sample_id]
    ] <- current_plot
    
    ggplot2::ggsave(
      filename = file.path(
        he_gene_figure_dir,
        paste0(
          "18_",
          sample_id,
          "_",
          gene_name,
          "_HE_expression_overlay.pdf"
        )
      ),
      plot = current_plot,
      width = plot_params$panel_width,
      height = plot_params$panel_height
    )
  }
  
  current_gene_panel <- patchwork::wrap_plots(
    gene_sample_plot_list[[gene_name]
    ],
    ncol = 4
  )
  
  ggplot2::ggsave(
    filename = file.path(
      he_figure_dir,
      paste0(
        "18_all_samples_",
        gene_name,
        "_HE_expression_panel.pdf"
      )
    ),
    plot = current_gene_panel,
    width = 20,
    height = 10
  )
}

# -----------------------------
# L1. 生成8样本×4基因总图
# -----------------------------

all_gene_plots <- list()

for (sample_id in sample_mapping$spatial_sample_id) {
  for (gene_name in genes_found) {
    plot_id <- paste(
      sample_id,
      gene_name,
      sep = "_"
    )
    
    all_gene_plots[[plot_id]
    ] <- gene_sample_plot_list[[gene_name]
    ][[sample_id]
    ]
  }
}

all_gene_panel <- patchwork::wrap_plots(
  all_gene_plots,
  ncol = length(
    genes_found
  )
)

ggplot2::ggsave(
  filename = file.path(
    he_figure_dir,
    "18_all_samples_core_gene_HE_expression_panel.pdf"
  ),
  plot = all_gene_panel,
  width = 5 * length(
    genes_found
  ),
  height = 4.8 * nrow(
    sample_mapping
  ),
  limitsize = FALSE
)

# ============================================================
# M. 样本分辨率表达统计
# ============================================================

sample_gene_summary_rows <- list()

for (sample_id in sample_mapping$spatial_sample_id) {
  sample_table <- plotting_table %>%
    dplyr::filter(
      .data$spatial_sample_id ==
        sample_id
    )
  
  for (gene_name in genes_found) {
    expression_column <- paste0(
      gene_name,
      "_expression"
    )
    
    expression_values <- sample_table[[expression_column]
    ]
    
    sample_gene_summary_rows[[paste0(
        sample_id,
        "_",
        gene_name
      )]
    ] <- data.frame(
      spatial_sample_id = sample_id,
      gene = gene_name,
      total_spots = length(
        expression_values
      ),
      expressing_spots = sum(
        expression_values > 0,
        na.rm = TRUE
      ),
      expression_pct = mean(
        expression_values > 0,
        na.rm = TRUE
      ) * 100,
      mean_expression = mean(
        expression_values,
        na.rm = TRUE
      ),
      median_expression = stats::median(
        expression_values,
        na.rm = TRUE
      ),
      maximum_expression = max(
        expression_values,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }
}

sample_gene_summary <- dplyr::bind_rows(
  sample_gene_summary_rows
)

write_utf8_csv(
  sample_gene_summary,
  file.path(
    table_dir,
    "18_sample_resolved_core_gene_expression_summary.csv"
  )
)

# ============================================================
# N. 保存更新后的Seurat对象
# ============================================================

output_rds_file <- file.path(
  object_dir,
  "18_spatial_sample_resolved_HE_Seurat.rds"
)

if (file.exists(output_rds_file)) {
  file.remove(output_rds_file)
}

gc()

saveRDS(
  spatial_object,
  file = output_rds_file,
  compress = "gzip"
)

message(
  "Seurat对象保存完成：",
  output_rds_file
)

message(
  "文件大小：",
  round(
    file.info(output_rds_file)$size / 1024^3,
    3
  ),
  "GB"
)

# ============================================================
# O. 输出文件检查
# ============================================================

required_output_files <- c(
  file.path(
    object_dir,
    "18_spatial_sample_resolved_HE_Seurat.rds"
  ),
  file.path(
    table_dir,
    "18_gem_id_sample_mapping.csv"
  ),
  file.path(
    table_dir,
    "18_barcode_expression_coordinate_match_summary.csv"
  ),
  file.path(
    table_dir,
    "18_spatial_HE_image_inventory.csv"
  ),
  file.path(
    table_dir,
    "18_spatial_sample_resolved_expression_coordinates.csv"
  ),
  file.path(
    table_dir,
    "18_sample_spot_and_image_coordinate_check.csv"
  ),
  file.path(
    table_dir,
    "18_sample_resolved_core_gene_expression_summary.csv"
  ),
  file.path(
    he_figure_dir,
    "18_all_samples_HE_tissue_spot_overlay_panel.pdf"
  ),
  file.path(
    he_figure_dir,
    "18_all_samples_core_gene_HE_expression_panel.pdf"
  ),
  file.path(
    table_dir,
    "18_sessionInfo.txt"
  )
)

for (gene_name in genes_found) {
  required_output_files <- c(
    required_output_files,
    file.path(
      he_figure_dir,
      paste0(
        "18_all_samples_",
        gene_name,
        "_HE_expression_panel.pdf"
      )
    )
  )
}

output_status <- data.frame(
  output = required_output_files,
  exists = file.exists(
    required_output_files
  ),
  stringsAsFactors = FALSE
)

write_utf8_csv(
  output_status,
  file.path(
    table_dir,
    "18_output_file_check.csv"
  )
)

print(
  output_status
)

# ============================================================
# P. 最终提示
# ============================================================

message(
  "\n============================================================"
)

message(
  "18_spatial_HE_image_overlay_and_sample_resolved_maps.R运行完成。"
)

message(
  "表达矩阵与组织内坐标匹配spot数：",
  length(
    matched_barcodes
  )
)

message(
  "真实空间样本数：",
  length(
    unique(
      plotting_table$spatial_sample_id
    )
  )
)

message(
  "样本列表：",
  paste(
    unique(
      plotting_table$spatial_sample_id
    ),
    collapse = ", "
  )
)

message(
  "实际绘制基因：",
  paste(
    genes_found,
    collapse = ", "
  )
)

message(
  "请重点检查以下文件："
)

message(
  "1. results/figures/18_spatial_HE/18_all_samples_HE_tissue_spot_overlay_panel.pdf"
)

message(
  "2. results/figures/18_spatial_HE/18_all_samples_core_gene_HE_expression_panel.pdf"
)

message(
  "3. results/figures/18_spatial_HE/18_all_samples_SASH1_HE_expression_panel.pdf"
)

message(
  "4. results/tables/18_sample_spot_and_image_coordinate_check.csv"
)

message(
  "5. results/tables/18_sample_resolved_core_gene_expression_summary.csv"
)

message(
  "6. results/tables/18_output_file_check.csv"
)

message(
  "7. results/objects/18_spatial_sample_resolved_HE_Seurat.rds"
)

message(
  "============================================================\n"
)