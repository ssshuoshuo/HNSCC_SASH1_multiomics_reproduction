# 15b_scRNA_core_gene_expression_Figure4_cluster_like.R

# 本脚本功能：
# 1. 读取已注释的GSE215403单细胞Seurat对象
# 2. 自动识别UMAP、cluster列和细胞类型注释列
# 3. 构建论文Fig.4风格的cluster+cell type标签，例如c-6_Fibroblast / CAF
# 4. 绘制major cell type UMAP
# 5. 绘制cluster+cell type标签UMAP
# 6. 绘制SASH1、MYH11、EMP1、COL1A1的cluster级DotPlot
# 7. 绘制4个核心基因FeaturePlot
# 8. 绘制4个核心基因cluster级VlnPlot
# 9. 合并输出更接近论文Fig.4的panel图

# 本项目专用数据：
# GSE215403 scRNA-seq
# 优先使用：
# results/objects/05_manual_annotated_plot_ready.rds
#
# 与15号脚本区别：
# 15号按major cell type分组，图更清楚。
# 15b按cluster+cell type分组，图更接近论文Fig.4。


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

candidate_input_files <- c(
  file.path(
    object_dir,
    "05_manual_annotated_plot_ready.rds"
  ),
  file.path(
    object_dir,
    "05_manual_annotated_before_malignant_call.rds"
  ),
  file.path(
    object_dir,
    "08_final_malignant_call.rds"
  ),
  file.path(
    object_dir,
    "04_standard_Seurat_multi_resolution.rds"
  )
)

input_object_file <- candidate_input_files[
  file.exists(
    candidate_input_files
  )
][1]

if (
  is.na(
    input_object_file
  )
) {
  stop(
    paste0(
      "未找到可用输入对象。已尝试：\n",
      paste(
        candidate_input_files,
        collapse = "\n"
      )
    )
  )
}

core_genes <- c(
  "SASH1",
  "MYH11",
  "EMP1",
  "COL1A1"
)

candidate_celltype_columns <- c(
  "celltype_plot",
  "celltype_manual",
  "manual_celltype",
  "manual_cell_type",
  "manual_annotation",
  "major_celltype",
  "major_cell_type",
  "celltype",
  "cell_type",
  "CellType",
  "main_celltype",
  "main_cell_type",
  "annotated_celltype",
  "paper_celltype"
)

candidate_cluster_columns <- c(
  "seurat_clusters",
  "RNA_snn_res.0.2",
  "RNA_snn_res.0.3",
  "RNA_snn_res.0.5",
  "SCT_snn_res.0.2",
  "SCT_snn_res.0.3",
  "SCT_snn_res.0.5",
  "cluster_res_0.2",
  "cluster_res_0.3",
  "cluster_res_0.5"
)

# ============================================================
# C. 读取对象并识别metadata
# ============================================================

seurat_object <- readRDS(
  input_object_file
)

message(
  "读取对象：",
  input_object_file
)

if ("RNA" %in% Assays(seurat_object)) {
  DefaultAssay(seurat_object) <- "RNA"
}

available_reductions <- Reductions(
  seurat_object
)

if ("umap" %in% available_reductions) {
  
  umap_reduction <- "umap"
  
} else {
  
  umap_candidates <- grep(
    "umap",
    available_reductions,
    value = TRUE,
    ignore.case = TRUE
  )
  
  if (length(umap_candidates) == 0) {
    stop(
      paste0(
        "对象中没有UMAP reduction。当前reductions：",
        paste(
          available_reductions,
          collapse = ", "
        )
      )
    )
  }
  
  umap_reduction <- umap_candidates[1]
}

celltype_column <- candidate_celltype_columns[
  candidate_celltype_columns %in%
    colnames(
      seurat_object@meta.data
    )
][1]

if (
  is.na(
    celltype_column
  )
) {
  stop(
    paste0(
      "未自动识别细胞类型注释列。",
      "\n当前metadata列名包括：\n",
      paste(
        colnames(
          seurat_object@meta.data
        ),
        collapse = ", "
      )
    )
  )
}

cluster_column <- candidate_cluster_columns[
  candidate_cluster_columns %in%
    colnames(
      seurat_object@meta.data
    )
][1]

if (
  is.na(
    cluster_column
  )
) {
  stop(
    paste0(
      "未自动识别cluster列。",
      "\n当前metadata列名包括：\n",
      paste(
        colnames(
          seurat_object@meta.data
        ),
        collapse = ", "
      )
    )
  )
}

seurat_object$Figure4_celltype <- as.character(
  seurat_object@meta.data[
    ,
    celltype_column
  ]
)

seurat_object$Figure4_cluster <- as.character(
  seurat_object@meta.data[
    ,
    cluster_column
  ]
)

seurat_object$Figure4_cluster <- gsub(
  "^c-",
  "",
  seurat_object$Figure4_cluster
)

seurat_object$Figure4_cluster_celltype <- paste0(
  "c-",
  seurat_object$Figure4_cluster,
  "_",
  seurat_object$Figure4_celltype
)

cluster_order_table <- seurat_object@meta.data %>%
  dplyr::group_by(
    Figure4_cluster,
    Figure4_cluster_celltype
  ) %>%
  dplyr::summarise(
    cell_number = dplyr::n(),
    .groups = "drop"
  )

cluster_order_table$cluster_numeric <- suppressWarnings(
  as.numeric(
    as.character(
      cluster_order_table$Figure4_cluster
    )
  )
)

cluster_order_table <- cluster_order_table %>%
  dplyr::arrange(
    cluster_numeric,
    Figure4_cluster_celltype
  )

seurat_object$Figure4_cluster_celltype <- factor(
  seurat_object$Figure4_cluster_celltype,
  levels = unique(
    cluster_order_table$Figure4_cluster_celltype
  )
)

Idents(
  seurat_object
) <- seurat_object$Figure4_cluster_celltype

core_genes_found <- intersect(
  core_genes,
  rownames(
    seurat_object
  )
)

if (!all(core_genes %in% core_genes_found)) {
  stop(
    paste0(
      "未找到全部核心基因。找到：",
      paste(
        core_genes_found,
        collapse = ", "
      )
    )
  )
}

message(
  "使用UMAP reduction：",
  umap_reduction
)

message(
  "使用细胞类型列：",
  celltype_column
)

message(
  "使用cluster列：",
  cluster_column
)

message(
  "核心基因：",
  paste(
    core_genes_found,
    collapse = ", "
  )
)

message(
  "cluster+celltype标签数：",
  length(
    levels(
      seurat_object$Figure4_cluster_celltype
    )
  )
)

# ============================================================
# D. 处理表达矩阵
# ============================================================

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
  
  message(
    "未找到data layer，重新NormalizeData。"
  )
  
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

metadata_table <- seurat_object@meta.data

metadata_table$cell_barcode <- rownames(
  metadata_table
)

for (current_gene in core_genes_found) {
  
  metadata_table[
    ,
    current_gene
  ] <- as.numeric(
    expression_matrix[
      current_gene,
      metadata_table$cell_barcode
    ]
  )
}

# ============================================================
# E. 统计表
# ============================================================

celltype_gene_summary <- metadata_table %>%
  dplyr::group_by(
    Figure4_celltype
  ) %>%
  dplyr::summarise(
    cell_number = dplyr::n(),
    SASH1_average_expression = mean(
      SASH1,
      na.rm = TRUE
    ),
    SASH1_positive_percent = round(
      100 *
        mean(
          SASH1 > 0,
          na.rm = TRUE
        ),
      2
    ),
    MYH11_average_expression = mean(
      MYH11,
      na.rm = TRUE
    ),
    MYH11_positive_percent = round(
      100 *
        mean(
          MYH11 > 0,
          na.rm = TRUE
        ),
      2
    ),
    EMP1_average_expression = mean(
      EMP1,
      na.rm = TRUE
    ),
    EMP1_positive_percent = round(
      100 *
        mean(
          EMP1 > 0,
          na.rm = TRUE
        ),
      2
    ),
    COL1A1_average_expression = mean(
      COL1A1,
      na.rm = TRUE
    ),
    COL1A1_positive_percent = round(
      100 *
        mean(
          COL1A1 > 0,
          na.rm = TRUE
        ),
      2
    ),
    .groups = "drop"
  )

cluster_celltype_gene_summary <- metadata_table %>%
  dplyr::group_by(
    Figure4_cluster,
    Figure4_cluster_celltype,
    Figure4_celltype
  ) %>%
  dplyr::summarise(
    cell_number = dplyr::n(),
    SASH1_average_expression = mean(
      SASH1,
      na.rm = TRUE
    ),
    SASH1_positive_percent = round(
      100 *
        mean(
          SASH1 > 0,
          na.rm = TRUE
        ),
      2
    ),
    MYH11_average_expression = mean(
      MYH11,
      na.rm = TRUE
    ),
    MYH11_positive_percent = round(
      100 *
        mean(
          MYH11 > 0,
          na.rm = TRUE
        ),
      2
    ),
    EMP1_average_expression = mean(
      EMP1,
      na.rm = TRUE
    ),
    EMP1_positive_percent = round(
      100 *
        mean(
          EMP1 > 0,
          na.rm = TRUE
        ),
      2
    ),
    COL1A1_average_expression = mean(
      COL1A1,
      na.rm = TRUE
    ),
    COL1A1_positive_percent = round(
      100 *
        mean(
          COL1A1 > 0,
          na.rm = TRUE
        ),
      2
    ),
    .groups = "drop"
  )

write.csv(
  celltype_gene_summary,
  file.path(
    table_dir,
    "15b_Figure4_celltype_core_gene_expression_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  cluster_celltype_gene_summary,
  file.path(
    table_dir,
    "15b_Figure4_cluster_celltype_core_gene_expression_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  metadata_table,
  file.path(
    table_dir,
    "15b_Figure4_cell_metadata_with_core_gene_expression.csv"
  ),
  row.names = FALSE
)

# ============================================================
# F. Fig.4a major cell type UMAP
# ============================================================

p_umap_celltype <- DimPlot(
  object = seurat_object,
  reduction = umap_reduction,
  group.by = "Figure4_celltype",
  label = TRUE,
  repel = TRUE,
  label.size = 3.0,
  raster = FALSE
) +
  ggtitle(
    "(a) Manual Cell Type Annotation"
  ) +
  xlab(
    "UMAP_1"
  ) +
  ylab(
    "UMAP_2"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "right"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "15b_Figure4a_scRNA_celltype_UMAP.pdf"
  ),
  plot = p_umap_celltype,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(
    figure_dir,
    "15b_Figure4a_scRNA_celltype_UMAP.png"
  ),
  plot = p_umap_celltype,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# G. cluster+celltype UMAP
# ============================================================

p_umap_cluster_celltype <- DimPlot(
  object = seurat_object,
  reduction = umap_reduction,
  group.by = "Figure4_cluster_celltype",
  label = TRUE,
  repel = TRUE,
  label.size = 2.4,
  raster = FALSE
) +
  ggtitle(
    "Cluster-level Cell Type Annotation"
  ) +
  xlab(
    "UMAP_1"
  ) +
  ylab(
    "UMAP_2"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "15b_Figure4_cluster_celltype_UMAP.pdf"
  ),
  plot = p_umap_cluster_celltype,
  width = 8,
  height = 6
)

# ============================================================
# H. Fig.4b cluster级DotPlot
# ============================================================

p_dotplot_cluster <- DotPlot(
  object = seurat_object,
  features = core_genes,
  group.by = "Figure4_cluster_celltype"
) +
  coord_flip() +
  ggtitle(
    "(b) Gene Expression DotPlot in Clusters"
  ) +
  xlab(
    "Cluster and cell type"
  ) +
  ylab(
    "Core genes"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 7
    ),
    axis.text.y = element_text(
      size = 9
    )
  )

ggsave(
  filename = file.path(
    figure_dir,
    "15b_Figure4b_core_gene_DotPlot_by_cluster_celltype.pdf"
  ),
  plot = p_dotplot_cluster,
  width = 16,
  height = 6
)

# ============================================================
# I. Fig.4c/e/g/i FeaturePlot
# ============================================================

feature_plot_list <- list()

for (current_gene in core_genes) {
  
  feature_plot_list[[current_gene]] <- FeaturePlot(
    object = seurat_object,
    features = current_gene,
    reduction = umap_reduction,
    order = TRUE,
    raster = FALSE
  ) +
    ggtitle(
      paste0(
        current_gene,
        " FeaturePlot"
      )
    ) +
    xlab(
      "UMAP_1"
    ) +
    ylab(
      "UMAP_2"
    ) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
}

p_feature_panel <- patchwork::wrap_plots(
  feature_plot_list,
  ncol = 2
) +
  patchwork::plot_annotation(
    title = "Core Gene FeaturePlots on scRNA UMAP"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "15b_Figure4_core_gene_FeaturePlot_panel.pdf"
  ),
  plot = p_feature_panel,
  width = 11,
  height = 9
)

# ============================================================
# J. Fig.4d/f/h/j cluster级VlnPlot
# ============================================================

vln_plot_list <- list()

for (current_gene in core_genes) {
  
  vln_plot_list[[current_gene]] <- VlnPlot(
    object = seurat_object,
    features = current_gene,
    group.by = "Figure4_cluster_celltype",
    pt.size = 0,
    raster = FALSE
  ) +
    ggtitle(
      paste0(
        current_gene,
        " Expression by Cluster"
      )
    ) +
    xlab(
      "Cluster and cell type"
    ) +
    ylab(
      "Expression"
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 60,
        hjust = 1,
        vjust = 1,
        size = 6
      ),
      legend.position = "none"
    )
}

p_vln_panel <- patchwork::wrap_plots(
  vln_plot_list,
  ncol = 2
) +
  patchwork::plot_annotation(
    title = "Core Gene Violin Plots by Cluster"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "15b_Figure4_core_gene_VlnPlot_by_cluster_celltype_panel.pdf"
  ),
  plot = p_vln_panel,
  width = 18,
  height = 11
)

# ============================================================
# K. 单基因单独输出，方便替换论文panel
# ============================================================

for (current_gene in core_genes) {
  
  ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        "15b_Figure4_FeaturePlot_",
        current_gene,
        ".pdf"
      )
    ),
    plot = feature_plot_list[[current_gene]],
    width = 7,
    height = 5.5
  )
  
  ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        "15b_Figure4_VlnPlot_",
        current_gene,
        "_by_cluster_celltype.pdf"
      )
    ),
    plot = vln_plot_list[[current_gene]],
    width = 14,
    height = 5.5
  )
}

# ============================================================
# L. 合并论文Fig.4-like总图
# ============================================================

p_figure4_like <- (
  p_umap_celltype +
    p_dotplot_cluster
) / (
  p_feature_panel
) / (
  p_vln_panel
) +
  patchwork::plot_annotation(
    title = "Figure4-like scRNA-seq Core Gene Expression Analysis",
    subtitle = "GSE215403 cluster-level expression of SASH1/MYH11/EMP1/COL1A1"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "15b_Figure4_like_scRNA_core_gene_expression_cluster_panel.pdf"
  ),
  plot = p_figure4_like,
  width = 18,
  height = 25
)

# ============================================================
# M. 保存对象和运行信息
# ============================================================

saveRDS(
  seurat_object,
  file.path(
    object_dir,
    "15b_Figure4_like_scRNA_cluster_core_gene_expression_Seurat.rds"
  ),
  compress = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    table_dir,
    "15b_sessionInfo.txt"
  )
)

# ============================================================
# N. 输出检查
# ============================================================

required_output_files <- c(
  file.path(
    object_dir,
    "15b_Figure4_like_scRNA_cluster_core_gene_expression_Seurat.rds"
  ),
  file.path(
    table_dir,
    "15b_Figure4_celltype_core_gene_expression_summary.csv"
  ),
  file.path(
    table_dir,
    "15b_Figure4_cluster_celltype_core_gene_expression_summary.csv"
  ),
  file.path(
    table_dir,
    "15b_Figure4_cell_metadata_with_core_gene_expression.csv"
  ),
  file.path(
    figure_dir,
    "15b_Figure4a_scRNA_celltype_UMAP.pdf"
  ),
  file.path(
    figure_dir,
    "15b_Figure4_cluster_celltype_UMAP.pdf"
  ),
  file.path(
    figure_dir,
    "15b_Figure4b_core_gene_DotPlot_by_cluster_celltype.pdf"
  ),
  file.path(
    figure_dir,
    "15b_Figure4_core_gene_FeaturePlot_panel.pdf"
  ),
  file.path(
    figure_dir,
    "15b_Figure4_core_gene_VlnPlot_by_cluster_celltype_panel.pdf"
  ),
  file.path(
    figure_dir,
    "15b_Figure4_like_scRNA_core_gene_expression_cluster_panel.pdf"
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
    "15b_output_file_check.csv"
  ),
  row.names = FALSE
)

print(
  output_status
)

message("\n============================================================")
message("15b 论文Fig.4 cluster级单细胞核心基因表达图完成。")
message("")
message("重点查看：")
message("1. 15b_Figure4_like_scRNA_core_gene_expression_cluster_panel.pdf")
message("2. 15b_Figure4b_core_gene_DotPlot_by_cluster_celltype.pdf")
message("3. 15b_Figure4_core_gene_VlnPlot_by_cluster_celltype_panel.pdf")
message("4. 15b_Figure4_cluster_celltype_UMAP.pdf")
message("5. 15b_Figure4_cluster_celltype_core_gene_expression_summary.csv")
message("")
message("解释方向：")
message("这版比15号更接近论文Fig.4，因为DotPlot和VlnPlot使用cluster+cell type细分标签。")
message("检查COL1A1是否集中在CAF/fibroblast cluster。")
message("检查SASH1是否在肿瘤相关cluster中整体较低。")
message("检查EMP1是否在部分肿瘤和基质/内皮cluster中表达。")
message("检查MYH11是否整体低表达。")
message("============================================================\n")