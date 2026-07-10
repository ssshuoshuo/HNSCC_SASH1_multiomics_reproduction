# 14_spatial_domain_annotation_paper_style.R

# 本脚本功能：
# 1. 读取12号空间Seurat对象
# 2. 对GSE252265空间spot进行标准化、降维、聚类
# 3. 用marker module score近似注释论文Fig.6中的空间功能domain
# 4. 绘制接近论文Fig.6b的空间domain annotation图
# 5. 绘制接近论文Fig.6c-f的SASH1、COL1A1、EMP1、MYH11空间表达图
# 6. 合并输出paper-style空间6联图
# 7. 保留通用marker score表和domain annotation表，方便换数据集复用

# 本项目专用数据：
# GSE252265 Visium空间转录组
# 输入对象：
# results/objects/12_spatial_tissue_spots_Seurat.rds
#
# 注意：
# 论文Fig.6a使用H&E图像。
# 当前12号对象是基于聚合表达矩阵和坐标构建的无图像Seurat对象。
# 因此本脚本用spot tissue layout近似替代Fig.6a。
# 若后续能定位每张切片的H&E图和scale factor，可再补SpatialFeaturePlot原生H&E叠加图。

# 通用代码修改位置：
# 1. 换数据集时：
#    修改input_object_file
#
# 2. 换论文domain时：
#    修改domain_marker_list
#
# 3. 换目标基因时：
#    修改core_genes
#
# 4. 调整空间聚类精细度时：
#    修改cluster_resolution


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

raw_dir <- file.path(
  project_dir,
  "data",
  "raw",
  "GSE252265"
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

cluster_resolution <- 0.45

core_genes <- c(
  "SASH1",
  "COL1A1",
  "EMP1",
  "MYH11"
)

# ============================================================
# C. 读取空间对象
# ============================================================

if (!file.exists(input_object_file)) {
  stop(
    paste0(
      "未找到输入对象：",
      input_object_file,
      "\n请先运行12_spatial_download_QC_gene_maps.R。"
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

message(
  "输入空间spot数：",
  ncol(spatial_object)
)

message(
  "空间样本标签：",
  paste(
    unique(
      spatial_object$spatial_sample_id
    ),
    collapse = ", "
  )
)

# ============================================================
# D. 检查RAW里是否存在图像文件
# ============================================================

raw_extracted_dir <- file.path(
  raw_dir,
  "raw_extracted"
)

possible_image_files <- character(0)

if (dir.exists(raw_extracted_dir)) {
  
  possible_image_files <- list.files(
    raw_extracted_dir,
    pattern = "tissue_hires_image\\.png$|tissue_lowres_image\\.png$|\\.jpg$|\\.jpeg$|\\.png$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
}

image_inventory <- data.frame(
  image_file = possible_image_files,
  stringsAsFactors = FALSE
)

write.csv(
  image_inventory,
  file.path(
    table_dir,
    "14_possible_spatial_image_files.csv"
  ),
  row.names = FALSE
)

if (length(possible_image_files) == 0) {
  message(
    "未在raw_extracted中识别到可直接使用的H&E/PNG图像。",
    "\n本脚本将用spot tissue layout近似替代论文Fig.6a。"
  )
} else {
  message(
    "识别到可能的空间图像文件，已写入14_possible_spatial_image_files.csv。"
  )
}

# ============================================================
# E. 标准化、降维和聚类
# ============================================================

spatial_object <- NormalizeData(
  spatial_object,
  assay = "SpatialRNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

spatial_object <- FindVariableFeatures(
  spatial_object,
  assay = "SpatialRNA",
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

spatial_object <- ScaleData(
  spatial_object,
  assay = "SpatialRNA",
  features = VariableFeatures(
    spatial_object
  ),
  verbose = FALSE
)

spatial_object <- RunPCA(
  spatial_object,
  assay = "SpatialRNA",
  features = VariableFeatures(
    spatial_object
  ),
  npcs = 30,
  verbose = FALSE
)

spatial_object <- FindNeighbors(
  spatial_object,
  reduction = "pca",
  dims = 1:20,
  verbose = FALSE
)

spatial_object <- FindClusters(
  spatial_object,
  resolution = cluster_resolution,
  verbose = FALSE
)

spatial_object <- RunUMAP(
  spatial_object,
  reduction = "pca",
  dims = 1:20,
  verbose = FALSE
)

spatial_cluster_column <- paste0(
  "SpatialRNA_snn_res.",
  cluster_resolution
)

if (!spatial_cluster_column %in% colnames(spatial_object@meta.data)) {
  spatial_cluster_column <- "seurat_clusters"
}

spatial_object$paper_style_spatial_cluster <- as.character(
  spatial_object@meta.data[
    ,
    spatial_cluster_column
  ]
)

# ============================================================
# F. 定义论文Fig.6风格的空间domain marker
# ============================================================

domain_marker_list <- list(
  "Fibrotic Stroma" = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "COL6A1",
    "DCN",
    "LUM",
    "FAP",
    "ACTA2",
    "TAGLN"
  ),
  "Myeloid-rich Inflammatory Zone" = c(
    "LYZ",
    "LST1",
    "TYROBP",
    "FCER1G",
    "AIF1",
    "C1QA",
    "C1QB",
    "C1QC",
    "CD68"
  ),
  "Differentiated/Keratinizing Tumor" = c(
    "KRT1",
    "KRT4",
    "KRT10",
    "KRT13",
    "IVL",
    "SPRR1B",
    "SPRR2A",
    "KRTDAP"
  ),
  "Metabolic/Secretory Tumor" = c(
    "KRT8",
    "KRT18",
    "KRT19",
    "MUC1",
    "SLPI",
    "LCN2",
    "AGR2",
    "TFF3"
  ),
  "Basal-like Proliferative Tumor" = c(
    "KRT5",
    "KRT14",
    "KRT15",
    "TP63",
    "MKI67",
    "TOP2A",
    "PCNA",
    "UBE2C"
  ),
  "CSC-like Niche" = c(
    "SOX2",
    "ALDH1A1",
    "PROM1",
    "EPCAM",
    "KRT17",
    "ITGA6",
    "LGR5",
    "CD44"
  )
)

domain_marker_list_found <- lapply(
  domain_marker_list,
  function(current_genes) {
    intersect(
      current_genes,
      rownames(spatial_object)
    )
  }
)

domain_marker_inventory <- data.frame(
  domain = rep(
    names(domain_marker_list_found),
    lengths(domain_marker_list_found)
  ),
  marker_gene = unlist(
    domain_marker_list_found,
    use.names = FALSE
  ),
  stringsAsFactors = FALSE
)

write.csv(
  domain_marker_inventory,
  file.path(
    table_dir,
    "14_domain_marker_genes_found.csv"
  ),
  row.names = FALSE
)

valid_domain_marker_list <- domain_marker_list_found[
  lengths(domain_marker_list_found) >= 2
]

if (length(valid_domain_marker_list) < 2) {
  stop(
    "可用domain marker gene set少于2组，无法进行空间domain注释。"
  )
}

message(
  "可用空间domain marker set：",
  paste(
    names(valid_domain_marker_list),
    collapse = ", "
  )
)

# ============================================================
# G. 计算domain module score
# ============================================================

for (current_domain in names(valid_domain_marker_list)) {
  
  current_clean_name <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    current_domain
  )
  
  current_score_prefix <- paste0(
    "DomainScore_",
    current_clean_name
  )
  
  spatial_object <- AddModuleScore(
    object = spatial_object,
    features = list(
      valid_domain_marker_list[[current_domain]]
    ),
    assay = "SpatialRNA",
    name = current_score_prefix,
    search = FALSE
  )
  
  generated_score_column <- paste0(
    current_score_prefix,
    "1"
  )
  
  final_score_column <- paste0(
    "score_",
    current_clean_name
  )
  
  colnames(spatial_object@meta.data)[
    colnames(spatial_object@meta.data) ==
      generated_score_column
  ] <- final_score_column
}

domain_score_columns <- grep(
  "^score_",
  colnames(spatial_object@meta.data),
  value = TRUE
)

domain_score_to_label <- data.frame(
  score_column = domain_score_columns,
  domain_label = gsub(
    "^score_",
    "",
    domain_score_columns
  ),
  stringsAsFactors = FALSE
)

domain_score_to_label$domain_label <- gsub(
  "_",
  " ",
  domain_score_to_label$domain_label
)

# ============================================================
# H. 按cluster平均module score注释空间domain
# ============================================================

metadata_for_domain <- spatial_object@meta.data

cluster_domain_score <- metadata_for_domain %>%
  dplyr::group_by(
    paper_style_spatial_cluster
  ) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(
        domain_score_columns
      ),
      mean,
      na.rm = TRUE
    ),
    spot_number = dplyr::n(),
    .groups = "drop"
  )

cluster_domain_annotation <- cluster_domain_score

cluster_domain_annotation$best_score_column <- apply(
  cluster_domain_annotation[
    ,
    domain_score_columns,
    drop = FALSE
  ],
  1,
  function(current_scores) {
    domain_score_columns[
      which.max(
        current_scores
      )
    ]
  }
)

cluster_domain_annotation$paper_style_spatial_domain <- domain_score_to_label$domain_label[
  match(
    cluster_domain_annotation$best_score_column,
    domain_score_to_label$score_column
  )
]

write.csv(
  cluster_domain_score,
  file.path(
    table_dir,
    "14_cluster_domain_marker_score_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  cluster_domain_annotation,
  file.path(
    table_dir,
    "14_cluster_paper_style_domain_annotation.csv"
  ),
  row.names = FALSE
)

spatial_object$paper_style_spatial_domain <- cluster_domain_annotation$paper_style_spatial_domain[
  match(
    spatial_object$paper_style_spatial_cluster,
    cluster_domain_annotation$paper_style_spatial_cluster
  )
]

spatial_object$paper_style_spatial_domain <- factor(
  spatial_object$paper_style_spatial_domain,
  levels = c(
    "Fibrotic Stroma",
    "Myeloid rich Inflammatory Zone",
    "Differentiated Keratinizing Tumor",
    "Metabolic Secretory Tumor",
    "Basal like Proliferative Tumor",
    "CSC like Niche"
  )
)

# ============================================================
# I. 提取核心基因表达和绘图数据
# ============================================================

core_genes_found <- intersect(
  core_genes,
  rownames(spatial_object)
)

if (!all(core_genes %in% core_genes_found)) {
  stop(
    paste0(
      "空间对象未找到全部核心基因。找到：",
      paste(
        core_genes_found,
        collapse = ", "
      )
    )
  )
}

normalized_expression <- LayerData(
  object = spatial_object,
  assay = "SpatialRNA",
  layer = "data"
)

plot_data <- spatial_object@meta.data

plot_data$barcode <- rownames(
  plot_data
)

for (current_gene in core_genes_found) {
  
  plot_data[
    ,
    current_gene
  ] <- as.numeric(
    normalized_expression[
      current_gene,
      plot_data$barcode
    ]
  )
}

write.csv(
  plot_data,
  file.path(
    table_dir,
    "14_spatial_domain_and_core_gene_metadata.csv"
  ),
  row.names = FALSE
)

domain_composition_summary <- plot_data %>%
  dplyr::group_by(
    spatial_sample_id,
    paper_style_spatial_domain
  ) %>%
  dplyr::summarise(
    spot_number = dplyr::n(),
    median_SASH1 = median(
      SASH1,
      na.rm = TRUE
    ),
    median_COL1A1 = median(
      COL1A1,
      na.rm = TRUE
    ),
    median_EMP1 = median(
      EMP1,
      na.rm = TRUE
    ),
    median_MYH11 = median(
      MYH11,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  domain_composition_summary,
  file.path(
    table_dir,
    "14_domain_core_gene_expression_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# J. 论文Fig.6风格绘图函数
# ============================================================

plot_spatial_layout <- function(
    current_plot_data
) {
  
  ggplot(
    current_plot_data,
    aes(
      x = spatial_x,
      y = spatial_y
    )
  ) +
    geom_point(
      size = 0.55,
      alpha = 0.85,
      color = "grey35"
    ) +
    scale_y_reverse() +
    coord_equal() +
    facet_wrap(
      ~spatial_sample_id,
      scales = "fixed"
    ) +
    labs(
      title = "(a) Tissue spot layout",
      subtitle = "H&E image not embedded in current Seurat object",
      x = NULL,
      y = NULL
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 8
      ),
      strip.text = element_text(
        face = "bold"
      )
    )
}

plot_spatial_domain <- function(
    current_plot_data
) {
  
  ggplot(
    current_plot_data,
    aes(
      x = spatial_x,
      y = spatial_y,
      color = paper_style_spatial_domain
    )
  ) +
    geom_point(
      size = 0.65,
      alpha = 0.9
    ) +
    scale_y_reverse() +
    coord_equal() +
    facet_wrap(
      ~spatial_sample_id,
      scales = "fixed"
    ) +
    labs(
      title = "(b) Spatial architecture of HNSCC functional domains",
      x = NULL,
      y = NULL,
      color = "Domain"
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      strip.text = element_text(
        face = "bold"
      ),
      legend.position = "right"
    )
}

plot_spatial_gene <- function(
    current_plot_data,
    gene_name,
    panel_label
) {
  
  ggplot(
    current_plot_data,
    aes(
      x = spatial_x,
      y = spatial_y,
      color = .data[[gene_name]]
    )
  ) +
    geom_point(
      size = 0.65,
      alpha = 0.9
    ) +
    scale_color_gradientn(
      colors = c(
        "#1b0c41",
        "#482878",
        "#b5367a",
        "#f98e52",
        "#fcfdbf"
      )
    ) +
    scale_y_reverse() +
    coord_equal() +
    facet_wrap(
      ~spatial_sample_id,
      scales = "fixed"
    ) +
    labs(
      title = paste0(
        panel_label,
        " ",
        gene_name
      ),
      x = NULL,
      y = NULL,
      color = gene_name
    ) +
    theme_void(base_size = 10) +
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
# K. 输出单独图
# ============================================================

p_layout <- plot_spatial_layout(
  plot_data
)

p_domain <- plot_spatial_domain(
  plot_data
)

p_SASH1 <- plot_spatial_gene(
  current_plot_data = plot_data,
  gene_name = "SASH1",
  panel_label = "(c)"
)

p_COL1A1 <- plot_spatial_gene(
  current_plot_data = plot_data,
  gene_name = "COL1A1",
  panel_label = "(d)"
)

p_EMP1 <- plot_spatial_gene(
  current_plot_data = plot_data,
  gene_name = "EMP1",
  panel_label = "(e)"
)

p_MYH11 <- plot_spatial_gene(
  current_plot_data = plot_data,
  gene_name = "MYH11",
  panel_label = "(f)"
)

ggsave(
  filename = file.path(
    figure_dir,
    "14_spatial_domain_tissue_layout.pdf"
  ),
  plot = p_layout,
  width = 7,
  height = 6
)

ggsave(
  filename = file.path(
    figure_dir,
    "14_spatial_domain_annotation.pdf"
  ),
  plot = p_domain,
  width = 9,
  height = 7
)

ggsave(
  filename = file.path(
    figure_dir,
    "14_spatial_domain_SASH1_expression.pdf"
  ),
  plot = p_SASH1,
  width = 7,
  height = 6
)

ggsave(
  filename = file.path(
    figure_dir,
    "14_spatial_domain_COL1A1_expression.pdf"
  ),
  plot = p_COL1A1,
  width = 7,
  height = 6
)

ggsave(
  filename = file.path(
    figure_dir,
    "14_spatial_domain_EMP1_expression.pdf"
  ),
  plot = p_EMP1,
  width = 7,
  height = 6
)

ggsave(
  filename = file.path(
    figure_dir,
    "14_spatial_domain_MYH11_expression.pdf"
  ),
  plot = p_MYH11,
  width = 7,
  height = 6
)

# ============================================================
# L. 合并论文Fig.6风格6联图
# ============================================================

p_figure6_like <- (
  p_layout +
    p_domain
) / (
  p_SASH1 +
    p_COL1A1
) / (
  p_EMP1 +
    p_MYH11
) +
  patchwork::plot_annotation(
    title = "Paper-style Spatial Transcriptomics Analysis of Core Genes in GSE252265",
    subtitle = "Approximation of Fig.6 using pooled spatial coordinates without embedded H&E image"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "14_spatial_domain_Figure6_like_panel.pdf"
  ),
  plot = p_figure6_like,
  width = 16,
  height = 20
)

# ============================================================
# M. 保存对象和运行信息
# ============================================================

saveRDS(
  spatial_object,
  file.path(
    object_dir,
    "14_spatial_domain_annotated_paper_style_Seurat.rds"
  ),
  compress = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    table_dir,
    "14_sessionInfo.txt"
  )
)

# ============================================================
# N. 输出检查
# ============================================================

required_output_files <- c(
  file.path(
    object_dir,
    "14_spatial_domain_annotated_paper_style_Seurat.rds"
  ),
  file.path(
    table_dir,
    "14_possible_spatial_image_files.csv"
  ),
  file.path(
    table_dir,
    "14_domain_marker_genes_found.csv"
  ),
  file.path(
    table_dir,
    "14_cluster_domain_marker_score_summary.csv"
  ),
  file.path(
    table_dir,
    "14_cluster_paper_style_domain_annotation.csv"
  ),
  file.path(
    table_dir,
    "14_spatial_domain_and_core_gene_metadata.csv"
  ),
  file.path(
    table_dir,
    "14_domain_core_gene_expression_summary.csv"
  ),
  file.path(
    figure_dir,
    "14_spatial_domain_tissue_layout.pdf"
  ),
  file.path(
    figure_dir,
    "14_spatial_domain_annotation.pdf"
  ),
  file.path(
    figure_dir,
    "14_spatial_domain_SASH1_expression.pdf"
  ),
  file.path(
    figure_dir,
    "14_spatial_domain_COL1A1_expression.pdf"
  ),
  file.path(
    figure_dir,
    "14_spatial_domain_EMP1_expression.pdf"
  ),
  file.path(
    figure_dir,
    "14_spatial_domain_MYH11_expression.pdf"
  ),
  file.path(
    figure_dir,
    "14_spatial_domain_Figure6_like_panel.pdf"
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
    "14_output_file_check.csv"
  ),
  row.names = FALSE
)

print(
  output_status
)

message("\n============================================================")
message("14 论文Fig.6风格空间domain和核心基因图完成。")
message("")
message("重点查看：")
message("1. 14_spatial_domain_Figure6_like_panel.pdf")
message("2. 14_spatial_domain_annotation.pdf")
message("3. 14_domain_core_gene_expression_summary.csv")
message("4. 14_cluster_paper_style_domain_annotation.csv")
message("")
message("解释限制：")
message("当前结果基于All_spots聚合坐标，尚未恢复论文Fig.6中的代表性单切片H&E背景。")
message("如果14_possible_spatial_image_files.csv中发现H&E/PNG图像，后续可继续补H&E背景叠加。")
message("============================================================\n")