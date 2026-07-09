# 04_standard_Seurat_PCA_UMAP_resolution_scan.R

# 本脚本功能：
# 1. 读取03中QC过滤后的Seurat object
# 2. 采用标准Seurat主流程进行降维和聚类
# 3. 执行LogNormalize、HVG筛选、ScaleData、PCA、UMAP和聚类
# 4. 同时计算多个clustering resolution
# 5. 输出PCA elbow plot、UMAP图和每个resolution的cluster统计表
# 6. 保存供后续marker分析和人工细胞类型注释使用的Seurat object
# 7. 记录本步骤使用的参数和session信息

# 本项目专用数据：
# GSE215403
# 12个OSCC单细胞样本：
# OSCC, scB1, scB2, scB5, scB7, scB8,
# scB9, scB10, scB12, scB13, scB14, scB15
#
# 本脚本采用标准Seurat流程：
# LogNormalize → FindVariableFeatures → ScaleData → PCA
# → UMAP → FindNeighbors → FindClusters
#
# 本脚本同时保存多个resolution，便于后续根据marker coherence、
# 样本组成和主要细胞类型分离情况选择主分辨率。
#
# 通用代码修改位置：
# 1. 换数据集时：
#    修改input_object_file
#
# 2. 调整高变基因数量时：
#    修改variable_feature_number
#
# 3. 调整PCA和UMAP维度时：
#    修改n_pcs_total和n_dims_for_umap
#
# 4. 调整聚类分辨率时：
#    修改clustering_resolutions和primary_resolution


# ============================================================
# A. 加载包
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 3600)

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tidyr",
  "ggplot2",
  "patchwork"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(Seurat)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

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
# C. 读取03 QC对象
# ============================================================

input_object_file <- file.path(
  object_dir,
  "03_QC_reproduction_candidate.rds"
)

if (!file.exists(input_object_file)) {
  
  stop(
    paste0(
      "找不到03 QC对象：\n",
      input_object_file,
      "\n请先运行03_QC_reproduction_candidate.R"
    )
  )
}

sc <- readRDS(input_object_file)

DefaultAssay(sc) <- "RNA"

message("读取完成。")
message("当前细胞数：", ncol(sc))
message("当前基因数：", nrow(sc))

# ============================================================
# D. 设置分析参数
# ============================================================

# variable_feature_number：
# 参与PCA的高变基因数量。
#
# n_pcs_total：
# 计算PCA时保留的PC数量。
#
# n_dims_for_umap：
# 用于UMAP、邻居图和聚类的PC数量。
#
# clustering_resolutions：
# 同时扫描的聚类分辨率。
#
# primary_resolution：
# 临时主分辨率，用于生成综合诊断图。
# 最终主分辨率会结合后续marker分析确定。

variable_feature_number <- 2000

n_pcs_total <- 50
n_dims_for_umap <- 30

clustering_resolutions <- c(
  0.2,
  0.3,
  0.5
)

primary_resolution <- 0.3

set.seed(1234)

# ============================================================
# E. 合并RNA layers
# ============================================================

# Seurat v5合并多个样本后，RNA assay中可能存在多个counts layer。
# 标准化前先合并为统一layer，便于后续NormalizeData和ScaleData。

sc <- JoinLayers(
  object = sc,
  assay = "RNA"
)

# ============================================================
# F. LogNormalize
# ============================================================

sc <- NormalizeData(
  object = sc,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

# ============================================================
# G. 筛选高变基因
# ============================================================

sc <- FindVariableFeatures(
  object = sc,
  selection.method = "vst",
  nfeatures = variable_feature_number,
  verbose = TRUE
)

variable_features <- VariableFeatures(sc)

write.csv(
  data.frame(
    rank = seq_along(variable_features),
    gene = variable_features
  ),
  file.path(
    table_dir,
    "04_variable_features_2000.csv"
  ),
  row.names = FALSE
)

message("前20个高变基因：")
print(head(variable_features, 20))

# ============================================================
# H. ScaleData
# ============================================================

# 这里只缩放高变基因。
# 这样可以减少内存占用，并与后续PCA使用的features保持一致。

sc <- ScaleData(
  object = sc,
  features = variable_features,
  vars.to.regress = NULL,
  verbose = TRUE
)

# ============================================================
# I. PCA
# ============================================================

sc <- RunPCA(
  object = sc,
  features = variable_features,
  npcs = n_pcs_total,
  verbose = TRUE
)

p_elbow <- ElbowPlot(
  object = sc,
  reduction = "pca",
  ndims = n_pcs_total
) +
  ggtitle("04 PCA Elbow Plot")

ggsave(
  filename = file.path(
    figure_dir,
    "04_PCA_elbow_plot.pdf"
  ),
  plot = p_elbow,
  width = 8,
  height = 5
)

ggsave(
  filename = file.path(
    figure_dir,
    "04_PCA_elbow_plot.png"
  ),
  plot = p_elbow,
  width = 8,
  height = 5,
  dpi = 300
)

# ============================================================
# J. UMAP与邻居图
# ============================================================

# UMAP和邻居图都基于前n_dims_for_umap个PC。
# reduction.name设为umap_pca，避免和后续其他UMAP结果混淆。

sc <- RunUMAP(
  object = sc,
  reduction = "pca",
  dims = 1:n_dims_for_umap,
  reduction.name = "umap_pca",
  reduction.key = "UMAPPCA_",
  seed.use = 1234,
  verbose = TRUE
)

sc <- FindNeighbors(
  object = sc,
  reduction = "pca",
  dims = 1:n_dims_for_umap,
  verbose = TRUE
)

# ============================================================
# K. 多分辨率聚类
# ============================================================

# 每个resolution都会生成一个metadata列：
# cluster_res_0.2
# cluster_res_0.3
# cluster_res_0.5
#
# 后续可以对比不同resolution下的marker清晰度和细胞类型分离效果。

for (res in clustering_resolutions) {
  
  res_label <- format(
    res,
    nsmall = 1,
    trim = TRUE
  )
  
  message(
    "\n正在计算resolution=",
    res_label
  )
  
  sc <- FindClusters(
    object = sc,
    resolution = res,
    verbose = TRUE
  )
  
  cluster_column <- paste0(
    "cluster_res_",
    res_label
  )
  
  sc[[cluster_column]] <- as.character(
    Idents(sc)
  )
}

# ============================================================
# L. 设置临时主分辨率
# ============================================================

# 这里只是为了生成统一诊断图。
# 最终使用哪个resolution，会在后续marker分析和人工注释后确定。

primary_cluster_column <- paste0(
  "cluster_res_",
  format(primary_resolution, nsmall = 1, trim = TRUE)
)

Idents(sc) <- sc[[primary_cluster_column, drop = TRUE]]

# ============================================================
# M. 输出UMAP诊断图
# ============================================================

p_umap_sample <- DimPlot(
  object = sc,
  reduction = "umap_pca",
  group.by = "sample_id",
  pt.size = 0.15,
  raster = TRUE,
  shuffle = TRUE
) +
  ggtitle("04 PCA-UMAP colored by sample")

p_umap_primary_cluster <- DimPlot(
  object = sc,
  reduction = "umap_pca",
  group.by = primary_cluster_column,
  label = TRUE,
  repel = TRUE,
  pt.size = 0.15,
  raster = TRUE
) +
  ggtitle(
    paste0(
      "04 PCA-UMAP: resolution=",
      primary_resolution
    )
  )

ggsave(
  filename = file.path(
    figure_dir,
    "04_UMAP_sample_and_primary_cluster.pdf"
  ),
  plot = p_umap_sample + p_umap_primary_cluster,
  width = 16,
  height = 7
)

ggsave(
  filename = file.path(
    figure_dir,
    "04_UMAP_sample_and_primary_cluster.png"
  ),
  plot = p_umap_sample + p_umap_primary_cluster,
  width = 16,
  height = 7,
  dpi = 300
)

# ============================================================
# N. 输出每个resolution的UMAP和统计表
# ============================================================

for (res in clustering_resolutions) {
  
  res_label <- format(
    res,
    nsmall = 1,
    trim = TRUE
  )
  
  cluster_column <- paste0(
    "cluster_res_",
    res_label
  )
  
  p_cluster <- DimPlot(
    object = sc,
    reduction = "umap_pca",
    group.by = cluster_column,
    label = TRUE,
    repel = TRUE,
    pt.size = 0.15,
    raster = TRUE
  ) +
    ggtitle(
      paste0(
        "04 PCA-UMAP: clustering resolution=",
        res_label
      )
    )
  
  ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        "04_UMAP_cluster_resolution_",
        res_label,
        ".pdf"
      )
    ),
    plot = p_cluster,
    width = 10,
    height = 8
  )
  
  cluster_summary <- sc@meta.data %>%
    count(
      .data[[cluster_column]],
      name = "cell_number"
    ) %>%
    mutate(
      percent_of_all_cells = round(
        100 * cell_number / sum(cell_number),
        2
      )
    ) %>%
    rename(
      cluster = 1
    ) %>%
    arrange(
      suppressWarnings(as.numeric(cluster))
    )
  
  write.csv(
    cluster_summary,
    file.path(
      table_dir,
      paste0(
        "04_cluster_summary_resolution_",
        res_label,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}

# ============================================================
# O. 保存对象和环境信息
# ============================================================

sc$analysis_stage <- "standard_Seurat_PCA_UMAP_multi_resolution"

saveRDS(
  sc,
  file.path(
    object_dir,
    "04_standard_Seurat_multi_resolution.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "04_sessionInfo.txt"
  )
)

# ============================================================
# P. 最终提示
# ============================================================

message("\n============================================================")
message("04_standard_Seurat_PCA_UMAP_resolution_scan.R 运行完成。")
message("")
message("已保存对象：")
message("results/objects/04_standard_Seurat_multi_resolution.rds")
message("")
message("请重点查看：")
message("1. results/figures/04_PCA_elbow_plot.pdf")
message("2. results/figures/04_UMAP_sample_and_primary_cluster.pdf")
message("3. results/figures/04_UMAP_cluster_resolution_0.2.pdf")
message("4. results/figures/04_UMAP_cluster_resolution_0.3.pdf")
message("5. results/figures/04_UMAP_cluster_resolution_0.5.pdf")
message("6. results/tables/04_cluster_summary_resolution_0.2.csv")
message("7. results/tables/04_cluster_summary_resolution_0.3.csv")
message("8. results/tables/04_cluster_summary_resolution_0.5.csv")
message("============================================================\n")