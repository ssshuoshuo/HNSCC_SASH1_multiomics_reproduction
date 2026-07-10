# 15_scRNA_core_gene_expression_Figure4_like.R

# 本脚本功能：
# 1. 读取已注释的GSE215403单细胞Seurat对象
# 2. 自动识别UMAP降维和细胞类型注释列
# 3. 绘制论文Fig.4风格的细胞类型UMAP
# 4. 绘制SASH1、MYH11、EMP1、COL1A1的DotPlot
# 5. 绘制4个核心基因的FeaturePlot
# 6. 绘制4个核心基因的VlnPlot
# 7. 合并输出接近论文Fig.4结构的panel图

# 本项目专用数据：
# GSE215403 scRNA-seq
# 优先使用：
# results/objects/05_manual_annotated_plot_ready.rds
# 若不存在，则依次尝试其他已注释对象

# 通用代码修改位置：
# 1. 换输入对象时：
#    修改candidate_input_files
#
# 2. 换核心基因时：
#    修改core_genes
#
# 3. 换细胞类型列时：
#    修改candidate_celltype_columns


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
  "SCT_snn_res.0.5"
)

# ============================================================
# C. 读取对象并识别关键信息
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
      ),
      "\n请把真实细胞类型列名填入candidate_celltype_columns。"
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
  cluster_column <- "seurat_clusters"
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

seurat_object$Figure4_cluster_celltype <- paste0(
  "c-",
  seurat_object$Figure4_cluster,
  "_",
  seurat_object$Figure4_celltype
)

cluster_celltype_order <- seurat_object@meta.data %>%
  dplyr::group_by(
    Figure4_cluster,
    Figure4_cluster_celltype
  ) %>%
  dplyr::summarise(
    cell_number = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    as.numeric(
      Figure4_cluster
    )
  ) %>%
  dplyr::pull(
    Figure4_cluster_celltype
  )

seurat_object$Figure4_cluster_celltype <- factor(
  seurat_object$Figure4_cluster_celltype,
  levels = unique(
    cluster_celltype_order
  )
)

Idents(
  seurat_object
) <- seurat_object$Figure4_celltype

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

# ============================================================
# D. 细胞类型和核心基因表达统计
# ============================================================

expression_matrix <- LayerData(
  object = seurat_object,
  assay = DefaultAssay(
    seurat_object
  ),
  layer = "data"
)

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

write.csv(
  celltype_gene_summary,
  file.path(
    table_dir,
    "15_Figure4_celltype_core_gene_expression_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  metadata_table,
  file.path(
    table_dir,
    "15_Figure4_cell_metadata_with_core_gene_expression.csv"
  ),
  row.names = FALSE
)

# ============================================================
# E. 绘制Fig.4a细胞类型UMAP
# ============================================================

p_umap_celltype <- DimPlot(
  object = seurat_object,
  reduction = umap_reduction,
  group.by = "Figure4_celltype",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggtitle(
    "(a) Manual Cell Type Annotation"
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
    "15_Figure4a_scRNA_celltype_UMAP.pdf"
  ),
  plot = p_umap_celltype,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(
    figure_dir,
    "15_Figure4a_scRNA_celltype_UMAP.png"
  ),
  plot = p_umap_celltype,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# F. 绘制Fig.4b DotPlot
# ============================================================

p_dotplot_celltype <- DotPlot(
  object = seurat_object,
  features = core_genes,
  group.by = "Figure4_celltype"
) +
  coord_flip() +
  ggtitle(
    "(b) Core Gene Expression DotPlot by Cell Type"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  filename = file.path(
    figure_dir,
    "15_Figure4b_core_gene_DotPlot_by_celltype.pdf"
  ),
  plot = p_dotplot_celltype,
  width = 9,
  height = 6
)

# ============================================================
# G. 绘制Fig.4c/e/g/i FeaturePlot
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
    "15_Figure4_core_gene_FeaturePlot_panel.pdf"
  ),
  plot = p_feature_panel,
  width = 11,
  height = 9
)

# ============================================================
# H. 绘制Fig.4d/f/h/j VlnPlot
# ============================================================

vln_plot_list <- list()

for (current_gene in core_genes) {
  
  vln_plot_list[[current_gene]] <- VlnPlot(
    object = seurat_object,
    features = current_gene,
    group.by = "Figure4_celltype",
    pt.size = 0,
    raster = FALSE
  ) +
    ggtitle(
      paste0(
        current_gene,
        " Expression by Cell Type"
      )
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
        size = 7
      )
    )
}

p_vln_panel <- patchwork::wrap_plots(
  vln_plot_list,
  ncol = 2
) +
  patchwork::plot_annotation(
    title = "Core Gene Violin Plots by Cell Type"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "15_Figure4_core_gene_VlnPlot_panel.pdf"
  ),
  plot = p_vln_panel,
  width = 14,
  height = 10
)

# ============================================================
# I. 合并论文Fig.4风格总图
# ============================================================

p_figure4_like <- (
  p_umap_celltype +
    p_dotplot_celltype
) / (
  p_feature_panel
) / (
  p_vln_panel
) +
  patchwork::plot_annotation(
    title = "Figure4-like scRNA-seq Core Gene Expression Analysis",
    subtitle = "GSE215403 manual cell type annotation and SASH1/MYH11/EMP1/COL1A1 expression"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "15_Figure4_like_scRNA_core_gene_expression_panel.pdf"
  ),
  plot = p_figure4_like,
  width = 16,
  height = 24
)

# ============================================================
# J. 保存对象和运行信息
# ============================================================

saveRDS(
  seurat_object,
  file.path(
    object_dir,
    "15_Figure4_like_scRNA_core_gene_expression_Seurat.rds"
  ),
  compress = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    table_dir,
    "15_sessionInfo.txt"
  )
)

# ============================================================
# K. 输出检查
# ============================================================

required_output_files <- c(
  file.path(
    object_dir,
    "15_Figure4_like_scRNA_core_gene_expression_Seurat.rds"
  ),
  file.path(
    table_dir,
    "15_Figure4_celltype_core_gene_expression_summary.csv"
  ),
  file.path(
    table_dir,
    "15_Figure4_cell_metadata_with_core_gene_expression.csv"
  ),
  file.path(
    figure_dir,
    "15_Figure4a_scRNA_celltype_UMAP.pdf"
  ),
  file.path(
    figure_dir,
    "15_Figure4b_core_gene_DotPlot_by_celltype.pdf"
  ),
  file.path(
    figure_dir,
    "15_Figure4_core_gene_FeaturePlot_panel.pdf"
  ),
  file.path(
    figure_dir,
    "15_Figure4_core_gene_VlnPlot_panel.pdf"
  ),
  file.path(
    figure_dir,
    "15_Figure4_like_scRNA_core_gene_expression_panel.pdf"
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
    "15_output_file_check.csv"
  ),
  row.names = FALSE
)

print(
  output_status
)

message("\n============================================================")
message("15 论文Fig.4风格单细胞核心基因表达图完成。")
message("")
message("重点查看：")
message("1. 15_Figure4_like_scRNA_core_gene_expression_panel.pdf")
message("2. 15_Figure4a_scRNA_celltype_UMAP.pdf")
message("3. 15_Figure4b_core_gene_DotPlot_by_celltype.pdf")
message("4. 15_Figure4_core_gene_FeaturePlot_panel.pdf")
message("5. 15_Figure4_core_gene_VlnPlot_panel.pdf")
message("6. 15_Figure4_celltype_core_gene_expression_summary.csv")
message("")
message("解释方向：")
message("检查COL1A1是否主要集中在Fibroblasts。")
message("检查SASH1是否在恶性细胞中整体偏低。")
message("检查EMP1是否在部分恶性细胞和基质/内皮细胞中有表达。")
message("检查MYH11是否整体低表达。")
message("============================================================\n")