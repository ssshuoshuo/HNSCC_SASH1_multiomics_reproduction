# ============================================================
# 04b_standard_Seurat_PCA_UMAP_resolution_scan.R
#
# 目标：
# 1. 读取 03b 固定阈值 QC 后对象
# 2. 按论文公开方法采用标准 Seurat 主流程
# 3. LogNormalize → HVG → ScaleData → PCA → UMAP → clustering
# 4. 同时保存多个 resolution，便于后续选择最合适的版本
#
# 论文方法公开描述：
# QC → normalization → highly variable features → PCA → UMAP
# → unsupervised clustering → manual annotation
#
# 注意：
# 论文没有公开精确 QC cutoff、HVG 数、PC 数或 resolution；
#
# ============================================================

# ============================================================
# 用户配置说明
# ============================================================
# 运行前请检查以下设置：
# 1. project_dir：项目根目录。
# 2. raw_dir：原始数据目录。
# 3. object_dir：RDS对象输出目录。
# 4. table_dir：CSV和TXT结果输出目录。
# 5. figure_dir：PDF图输出目录。
# 6. 输入文件名：若本地文件名不同，请在对应input_file处修改。
# 7. 线程数、内存和运行位置：CopyKAT、Seurat聚类和Monocle3建议在服务器或高内存本地环境运行。
# ============================================================

# 可替换参数：
#
# variable_feature_number：
#   常用 2000–3000
#
# pca_dims：
#   常用 20–40
#
# clustering_resolutions：
#   肿瘤数据常见 0.2–0.8
#
# 使用建议：
# 当前 GSE215403 复现主线优先采用：
# 2000 HVG + 30 PCs + resolution 0.2 / 0.3 / 0.5
# ============================================================

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
# B. 项目路径
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
# C. 读取 03b QC 对象
# ============================================================

input_object_file <- file.path(
  object_dir,
  "03b_GSE215403_QC_reproduction_candidate.rds"
)

if (!file.exists(input_object_file)) {
  
  stop(
    paste0(
      "找不到 03b QC 对象：\n",
      input_object_file,
      "\n请先运行 03b_QC_reproduction_candidate.R"
    )
  )
}

sc <- readRDS(input_object_file)

DefaultAssay(sc) <- "RNA"

message("读取完成。")
message("当前细胞数：", ncol(sc))
message("当前基因数：", nrow(sc))

# ============================================================
# D. 固定参数
# ============================================================

# -----------------------------
# 复现候选参数
# -----------------------------

variable_feature_number <- 2000

n_pcs_total <- 50
n_dims_for_umap <- 30

clustering_resolutions <- c(
  0.2,
  0.3,
  0.5
)

set.seed(1234)

# ============================================================
# E. 合并 RNA layers
# ============================================================

# Seurat v5 合并样本后可能存在多个 counts layer。
# 标准化前先合并为统一 RNA layer。

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
# G. 高变基因
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
    "04b_variable_features_2000.csv"
  ),
  row.names = FALSE
)

message("前 20 个高变基因：")
print(head(variable_features, 20))

# ============================================================
# H. ScaleData：只缩放高变基因
# ============================================================

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
  ggtitle("04b PCA Elbow Plot")

ggsave(
  filename = file.path(
    figure_dir,
    "04b_PCA_elbow_plot.pdf"
  ),
  plot = p_elbow,
  width = 8,
  height = 5
)

ggsave(
  filename = file.path(
    figure_dir,
    "04b_PCA_elbow_plot.png"
  ),
  plot = p_elbow,
  width = 8,
  height = 5,
  dpi = 300
)

# ============================================================
# J. UMAP 与邻居图
# ============================================================

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

for (res in clustering_resolutions) {
  
  res_label <- format(
    res,
    nsmall = 1,
    trim = TRUE
  )
  
  message(
    "\n正在计算 resolution = ",
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
# L. 选择临时主分辨率
# ============================================================

# 这里只是为了生成统一诊断图。
# 最终使用哪个 resolution，
# 要等 05b 的 marker coherence 比较后决定。

primary_resolution <- 0.3

primary_cluster_column <- paste0(
  "cluster_res_",
  format(primary_resolution, nsmall = 1, trim = TRUE)
)

Idents(sc) <- sc[[primary_cluster_column, drop = TRUE]]

# ============================================================
# M. 输出 UMAP 图
# ============================================================

p_umap_sample <- DimPlot(
  object = sc,
  reduction = "umap_pca",
  group.by = "sample_id",
  pt.size = 0.15,
  raster = TRUE,
  shuffle = TRUE
) +
  ggtitle("04b PCA-UMAP colored by sample")

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
      "04b PCA-UMAP: resolution = ",
      primary_resolution
    )
  )

ggsave(
  filename = file.path(
    figure_dir,
    "04b_UMAP_sample_and_primary_cluster.pdf"
  ),
  plot = p_umap_sample + p_umap_primary_cluster,
  width = 16,
  height = 7
)

ggsave(
  filename = file.path(
    figure_dir,
    "04b_UMAP_sample_and_primary_cluster.png"
  ),
  plot = p_umap_sample + p_umap_primary_cluster,
  width = 16,
  height = 7,
  dpi = 300
)

# ============================================================
# N. 每个 resolution 单独保存 cluster UMAP 与统计表
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
        "04b PCA-UMAP: clustering resolution = ",
        res_label
      )
    )
  
  ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        "04b_UMAP_cluster_resolution_",
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
        "04b_cluster_summary_resolution_",
        res_label,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}

# ============================================================
# O. 保存对象
# ============================================================

sc$analysis_stage <- "standard_Seurat_PCA_UMAP_multi_resolution"

saveRDS(
  sc,
  file.path(
    object_dir,
    "04b_GSE215403_standard_Seurat_multi_resolution.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "04b_sessionInfo.txt"
  )
)

# ============================================================
# P. 完成提示
# ============================================================

message("\n============================================================")
message("04b_standard_Seurat_PCA_UMAP_resolution_scan.R 运行完成。")
message("")
message("已保存对象：")
message("results/objects/04b_GSE215403_standard_Seurat_multi_resolution.rds")
message("")
message("重点查看：")
message("1. results/figures/04b_PCA_elbow_plot.pdf")
message("2. results/figures/04b_UMAP_sample_and_primary_cluster.pdf")
message("3. results/figures/04b_UMAP_cluster_resolution_0.2.pdf")
message("4. results/figures/04b_UMAP_cluster_resolution_0.3.pdf")
message("5. results/figures/04b_UMAP_cluster_resolution_0.5.pdf")
message("6. results/tables/04b_cluster_summary_resolution_0.2.csv")
message("7. results/tables/04b_cluster_summary_resolution_0.3.csv")
message("8. results/tables/04b_cluster_summary_resolution_0.5.csv")
message("============================================================\n")