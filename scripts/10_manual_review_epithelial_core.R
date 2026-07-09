# 10_manual_review_epithelial_core.R

# 本脚本功能：
# 1. 读取malignant epithelial state characterization对象
# 2. 基于内部cluster marker结果进行manual-review epithelial-core refinement
# 3. 排除明确immune、macrophage、CAF-stromal和lineage-ambiguous内部cluster
# 4. 对保留的上皮核心细胞集重新构建轻量Seurat object
# 5. 对上皮核心细胞集重新执行NormalizeData、HVG、ScaleData、PCA、UMAP和聚类
# 6. 计算Cycling、Cancer-testis、Squamous differentiation和Basal epithelial相对program score
# 7. 比较SASH1、EMP1、MYH11、COL1A1在内部cluster和相对program中的表达
# 8. 保存manual-review epithelial-core对象和相关统计表

# 本项目专用数据：
# GSE215403
# 12个OSCC单细胞样本：
# OSCC, scB1, scB2, scB5, scB7, scB8,
# scB9, scB10, scB12, scB13, scB14, scB15
#
# 本脚本中的core subset表示：
# 经内部marker审查后保留的上皮核心细胞集。
#
# 该core subset用于辅助理解恶性上皮候选细胞的内部状态，
# 不直接替代08中的CopyKAT-supported final malignant call。
#
# 通用代码修改位置：
# 1. 换数据集时：
#    修改input_file和manual_cluster_review
#
# 2. 换内部聚类结果时：
#    修改manual_cluster_review中的internal_cluster和review_decision
#
# 3. 换program定义时：
#    修改cycling_genes、cancer_testis_genes、
#    squamous_differentiation_genes和basal_epithelial_genes
#
# 4. 换关注基因时：
#    修改core_genes对应的基因列表


# ============================================================
# A. 加载包
# ============================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
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
      paste(missing_packages, collapse = ", ")
    )
  )
}

library(Seurat)
library(SeuratObject)
library(Matrix)
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

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# C. 读取malignant epithelial state characterization对象
# ============================================================

input_file <- file.path(
  object_dir,
  "10_malignant_epithelial_state_characterization.rds"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "找不到输入对象：\n",
      input_file
    )
  )
}

sc <- readRDS(input_file)

DefaultAssay(sc) <- "RNA"

if (!"seurat_clusters" %in% colnames(sc@meta.data)) {
  stop(
    "输入对象中缺少seurat_clusters。"
  )
}

if (!"sample_id" %in% colnames(sc@meta.data)) {
  stop(
    "输入对象中缺少sample_id。"
  )
}

message(
  "10输入细胞数：",
  ncol(sc)
)

# ============================================================
# D. 基于内部marker审查的cluster决策表
# ============================================================

# manual_cluster_review记录每个内部cluster的保留或排除决定。
#
# Retain_epithelial_core：
# 保留为上皮核心细胞。
#
# Exclude_lineage_ambiguous：
# 排除lineage-ambiguous或混合信号细胞。
#
# Exclude_stromal_or_myeloid_like、Exclude_CAF_stromal_like、
# Exclude_T_cell_like、Exclude_macrophage_like：
# 排除非上皮谱系信号明显的内部cluster。

manual_cluster_review <- data.frame(
  internal_cluster = c(
    "0", "1", "2", "3", "4", "5", "6",
    "7", "8", "9", "10", "11"
  ),
  review_decision = c(
    "Retain_epithelial_core",
    "Retain_epithelial_core",
    "Retain_epithelial_core",
    "Retain_epithelial_core",
    "Retain_epithelial_core",
    "Retain_epithelial_core",
    "Exclude_lineage_ambiguous",
    "Retain_epithelial_core",
    "Exclude_stromal_or_myeloid_like",
    "Exclude_CAF_stromal_like",
    "Exclude_T_cell_like",
    "Exclude_macrophage_like"
  ),
  rationale = c(
    "Cancer-testis antigen-like epithelial program.",
    "Epithelial malignant state retained for downstream comparison.",
    "Interferon/chemokine-like epithelial state retained.",
    "Squamous differentiation/keratinization-like epithelial state retained.",
    "Epithelial malignant state retained for downstream comparison.",
    "Differentiated or secretory-like epithelial state retained.",
    "Mixed/lineage-ambiguous immune-associated signal.",
    "Differentiated or secretory-like epithelial state retained.",
    "Stromal or myeloid-associated marker signal.",
    "CAF/stromal-associated marker signal.",
    "T-cell-associated marker signal.",
    "Macrophage-associated marker signal."
  ),
  stringsAsFactors = FALSE
)

retained_clusters <- manual_cluster_review %>%
  filter(
    review_decision == "Retain_epithelial_core"
  ) %>%
  pull(internal_cluster)

excluded_clusters <- manual_cluster_review %>%
  filter(
    review_decision != "Retain_epithelial_core"
  ) %>%
  pull(internal_cluster)

write.csv(
  manual_cluster_review,
  file.path(
    table_dir,
    "10_manual_review_cluster_decision.csv"
  ),
  row.names = FALSE
)

message(
  "保留内部cluster：",
  paste(retained_clusters, collapse = ", ")
)

message(
  "排除内部cluster：",
  paste(excluded_clusters, collapse = ", ")
)

# ============================================================
# E. 写入审查状态并输出筛选概览
# ============================================================

sc$internal_cluster_10_input <- as.character(
  sc$seurat_clusters
)

decision_map <- setNames(
  manual_cluster_review$review_decision,
  manual_cluster_review$internal_cluster
)

sc$manual_review_status <- unname(
  decision_map[
    sc$internal_cluster_10_input
  ]
)

sc$manual_review_status <- factor(
  sc$manual_review_status,
  levels = c(
    "Retain_epithelial_core",
    "Exclude_lineage_ambiguous",
    "Exclude_stromal_or_myeloid_like",
    "Exclude_CAF_stromal_like",
    "Exclude_T_cell_like",
    "Exclude_macrophage_like"
  )
)

review_summary <- sc@meta.data %>%
  mutate(
    internal_cluster_10_input = as.character(
      internal_cluster_10_input
    )
  ) %>%
  count(
    internal_cluster_10_input,
    manual_review_status,
    name = "cell_number"
  ) %>%
  left_join(
    manual_cluster_review,
    by = c(
      "internal_cluster_10_input" = "internal_cluster"
    )
  ) %>%
  mutate(
    percent_of_all_input_cells = round(
      100 * cell_number / sum(cell_number),
      2
    )
  ) %>%
  arrange(
    suppressWarnings(
      as.numeric(internal_cluster_10_input)
    )
  )

write.csv(
  review_summary,
  file.path(
    table_dir,
    "10_manual_review_cell_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# F. 在现有UMAP上展示保留和排除决策
# ============================================================

p_review_status <- DimPlot(
  sc,
  reduction = "umap",
  group.by = "manual_review_status",
  pt.size = 0.3,
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  ggtitle(
    "Internal clusters: manual-review epithelial-core decision"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )

p_internal_cluster <- DimPlot(
  sc,
  reduction = "umap",
  group.by = "internal_cluster_10_input",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.3,
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  ggtitle(
    "Internal cluster identity"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(
    figure_dir,
    "10_manual_review_cluster_decision_UMAP.pdf"
  ),
  plot = p_review_status +
    p_internal_cluster,
  width = 16,
  height = 8
)

# ============================================================
# G. 提取上皮核心集并重新构建轻量对象
# ============================================================

retained_cells <- colnames(sc)[
  sc$internal_cluster_10_input %in% retained_clusters
]

if (length(retained_cells) < 100) {
  stop(
    "保留细胞数少于100，无法继续。"
  )
}

message(
  "manual-review epithelial-core cell number: ",
  length(retained_cells)
)

raw_counts <- LayerData(
  object = sc,
  assay = "RNA",
  layer = "counts"
)

core_counts <- raw_counts[
  ,
  retained_cells,
  drop = FALSE
]

core_metadata <- sc@meta.data[
  retained_cells,
  ,
  drop = FALSE
]

rm(raw_counts)
rm(sc)
gc()

sc_core <- CreateSeuratObject(
  counts = core_counts,
  meta.data = core_metadata,
  project = "manual_review_epithelial_core"
)

rm(core_counts)
gc()

DefaultAssay(sc_core) <- "RNA"

# ============================================================
# H. 上皮核心集重新降维与内部聚类
# ============================================================

set.seed(1234)

sc_core <- NormalizeData(
  sc_core,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

sc_core <- FindVariableFeatures(
  sc_core,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

sc_core <- ScaleData(
  sc_core,
  features = VariableFeatures(sc_core),
  verbose = FALSE
)

sc_core <- RunPCA(
  sc_core,
  features = VariableFeatures(sc_core),
  npcs = 30,
  verbose = FALSE
)

umap_dims <- 1:min(
  20,
  ncol(
    Embeddings(
      sc_core,
      reduction = "pca"
    )
  )
)

sc_core <- RunUMAP(
  sc_core,
  reduction = "pca",
  dims = umap_dims,
  seed.use = 1234,
  verbose = FALSE
)

sc_core <- FindNeighbors(
  sc_core,
  reduction = "pca",
  dims = umap_dims,
  verbose = FALSE
)

sc_core <- FindClusters(
  sc_core,
  resolution = 0.3,
  verbose = FALSE
)

sc_core$epithelial_core_cluster <- as.character(
  sc_core$seurat_clusters
)

# ============================================================
# I. 计算相对program score
# ============================================================

# 这里先在上皮核心集内计算每个program的平均表达，
# 再将每个program score转为z-score。
# 因此结果表示该program在上皮核心集内部的相对富集程度。

cycling_genes <- c(
  "MKI67", "TOP2A", "CDK1", "CDC20", "CDCA5",
  "BIRC5", "UBE2C", "NUSAP1", "KIF11", "KIF18A",
  "CCNB1", "CCNB2", "TYMS", "CENPF"
)

cancer_testis_genes <- c(
  "CTAG2", "MAGEA1", "MAGEA3", "MAGEA4",
  "MAGEA6", "MAGEB2", "PRAME", "SSX1",
  "SSX2", "TEX15", "SPANXC"
)

squamous_differentiation_genes <- c(
  "IVL", "KRT13", "KRT4", "KRT16",
  "SPRR1A", "SPRR1B", "SPRR2B", "SPRR2D",
  "KRTDAP", "FLG2", "LCE3D", "LCE3E"
)

basal_epithelial_genes <- c(
  "KRT5", "KRT14", "KRT17", "KRT15",
  "TP63", "KRT6A", "KRT6B", "KRT19",
  "EPCAM", "DSG3", "SFN", "CLDN4"
)

program_gene_sets <- list(
  Cycling = cycling_genes,
  Cancer_testis = cancer_testis_genes,
  Squamous_differentiation = squamous_differentiation_genes,
  Basal_epithelial = basal_epithelial_genes
)

program_gene_sets_present <- lapply(
  program_gene_sets,
  function(x) {
    intersect(
      x,
      rownames(sc_core)
    )
  }
)

program_gene_set_table <- data.frame(
  program = names(program_gene_sets_present),
  genes_available = vapply(
    program_gene_sets_present,
    length,
    integer(1)
  )
)

write.csv(
  program_gene_set_table,
  file.path(
    table_dir,
    "10_relative_program_gene_sets_available.csv"
  ),
  row.names = FALSE
)

expression_matrix <- LayerData(
  object = sc_core,
  assay = "RNA",
  layer = "data"
)

calculate_program_score <- function(
    matrix_object,
    genes
) {
  if (length(genes) == 0) {
    return(
      rep(
        NA_real_,
        ncol(matrix_object)
      )
    )
  }
  
  Matrix::colMeans(
    matrix_object[
      genes,
      ,
      drop = FALSE
    ]
  )
}

for (program_name in names(program_gene_sets_present)) {
  raw_score_name <- paste0(
    program_name,
    "_raw_score"
  )
  
  z_score_name <- paste0(
    program_name,
    "_relative_z"
  )
  
  raw_score <- calculate_program_score(
    expression_matrix,
    program_gene_sets_present[[program_name]]
  )
  
  raw_score <- raw_score[
    colnames(sc_core)
  ]
  
  sc_core[[raw_score_name]] <- raw_score
  
  if (
    all(is.na(raw_score)) ||
    sd(raw_score, na.rm = TRUE) == 0
  ) {
    sc_core[[z_score_name]] <- NA_real_
  } else {
    sc_core[[z_score_name]] <- as.numeric(
      scale(raw_score)
    )
  }
}

rm(expression_matrix)
gc()

relative_score_columns <- paste0(
  names(program_gene_sets_present),
  "_relative_z"
)

# ============================================================
# J. 以内部cluster为单位计算相对program富集
# ============================================================

core_cluster_program_scores <- sc_core@meta.data %>%
  mutate(
    epithelial_core_cluster = as.character(
      epithelial_core_cluster
    )
  ) %>%
  group_by(
    epithelial_core_cluster
  ) %>%
  summarise(
    cell_number = n(),
    across(
      all_of(relative_score_columns),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

program_long <- core_cluster_program_scores %>%
  pivot_longer(
    cols = all_of(relative_score_columns),
    names_to = "program_z",
    values_to = "mean_relative_z"
  ) %>%
  mutate(
    program = gsub(
      "_relative_z$",
      "",
      program_z
    )
  )

program_ranked <- program_long %>%
  group_by(
    epithelial_core_cluster
  ) %>%
  arrange(
    desc(mean_relative_z),
    .by_group = TRUE
  ) %>%
  mutate(
    within_cluster_rank = row_number()
  ) %>%
  ungroup()

program_top_two <- program_ranked %>%
  filter(
    within_cluster_rank <= 2
  ) %>%
  select(
    epithelial_core_cluster,
    program,
    mean_relative_z,
    within_cluster_rank
  ) %>%
  pivot_wider(
    names_from = within_cluster_rank,
    values_from = c(
      program,
      mean_relative_z
    ),
    names_glue = "{.value}_{within_cluster_rank}"
  ) %>%
  mutate(
    top_vs_second_margin =
      mean_relative_z_1 -
      mean_relative_z_2,
    relative_program_label = case_when(
      is.na(mean_relative_z_1) ~
        "Unassigned",
      mean_relative_z_1 < 0.20 ~
        "Mixed_epithelial_program",
      top_vs_second_margin < 0.15 ~
        "Mixed_epithelial_program",
      TRUE ~ program_1
    )
  ) %>%
  select(
    epithelial_core_cluster,
    relative_program_label,
    top_program = program_1,
    top_program_relative_z = mean_relative_z_1,
    second_program = program_2,
    second_program_relative_z = mean_relative_z_2,
    top_vs_second_margin
  )

core_cluster_program_summary <- core_cluster_program_scores %>%
  left_join(
    program_top_two,
    by = "epithelial_core_cluster"
  ) %>%
  arrange(
    suppressWarnings(
      as.numeric(epithelial_core_cluster)
    )
  )

write.csv(
  core_cluster_program_summary,
  file.path(
    table_dir,
    "10_epithelial_core_relative_program_summary.csv"
  ),
  row.names = FALSE
)

program_label_map <- setNames(
  core_cluster_program_summary$relative_program_label,
  core_cluster_program_summary$epithelial_core_cluster
)

sc_core$relative_program_label <- unname(
  program_label_map[
    sc_core$epithelial_core_cluster
  ]
)

sc_core$relative_program_label <- factor(
  sc_core$relative_program_label,
  levels = c(
    "Cycling",
    "Cancer_testis",
    "Squamous_differentiation",
    "Basal_epithelial",
    "Mixed_epithelial_program",
    "Unassigned"
  )
)

# ============================================================
# K. 输出上皮核心集组成表
# ============================================================

core_by_sample_cluster <- sc_core@meta.data %>%
  mutate(
    epithelial_core_cluster = as.character(
      epithelial_core_cluster
    )
  ) %>%
  count(
    sample_id,
    epithelial_core_cluster,
    relative_program_label,
    name = "cell_number"
  ) %>%
  group_by(
    sample_id
  ) %>%
  mutate(
    percent_within_sample = round(
      100 * cell_number / sum(cell_number),
      2
    )
  ) %>%
  ungroup() %>%
  arrange(
    sample_id,
    suppressWarnings(
      as.numeric(epithelial_core_cluster)
    )
  )

write.csv(
  core_by_sample_cluster,
  file.path(
    table_dir,
    "10_epithelial_core_by_sample_cluster_and_program.csv"
  ),
  row.names = FALSE
)

# ============================================================
# L. 图：上皮核心集UMAP
# ============================================================

p_core_sample <- DimPlot(
  sc_core,
  reduction = "umap",
  group.by = "sample_id",
  pt.size = 0.35,
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  ggtitle(
    "Manual-review epithelial-core subset: UMAP by sample"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )

p_core_cluster <- DimPlot(
  sc_core,
  reduction = "umap",
  group.by = "epithelial_core_cluster",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.35,
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  ggtitle(
    "Manual-review epithelial-core subset: internal clusters"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )

p_core_program <- DimPlot(
  sc_core,
  reduction = "umap",
  group.by = "relative_program_label",
  pt.size = 0.35,
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  ggtitle(
    "Manual-review epithelial-core subset: relative program labels"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(
    figure_dir,
    "10_epithelial_core_UMAP_sample_cluster_program.pdf"
  ),
  plot = p_core_sample +
    p_core_cluster +
    p_core_program +
    plot_layout(ncol = 3),
  width = 22,
  height = 7
)

# ============================================================
# M. 图：相对program score UMAP
# ============================================================

p_relative_program_scores <- FeaturePlot(
  sc_core,
  features = relative_score_columns,
  reduction = "umap",
  ncol = 2,
  order = TRUE,
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  plot_annotation(
    title = "Relative transcriptional program scores in epithelial-core subset"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "10_epithelial_core_relative_program_scores_UMAP.pdf"
  ),
  plot = p_relative_program_scores,
  width = 13,
  height = 10
)

# ============================================================
# N. 核心基因表达图
# ============================================================

core_genes <- intersect(
  c(
    "SASH1",
    "EMP1",
    "MYH11",
    "COL1A1"
  ),
  rownames(sc_core)
)

p_core_genes_by_cluster <- DotPlot(
  sc_core,
  features = core_genes,
  group.by = "epithelial_core_cluster",
  dot.scale = 8
) +
  RotatedAxis() +
  ggtitle(
    "Core genes across epithelial-core internal clusters"
  ) +
  theme_classic(base_size = 12)

p_core_genes_by_program <- DotPlot(
  sc_core,
  features = core_genes,
  group.by = "relative_program_label",
  dot.scale = 8
) +
  RotatedAxis() +
  ggtitle(
    "Core genes across relative epithelial programs"
  ) +
  theme_classic(base_size = 12)

ggsave(
  filename = file.path(
    figure_dir,
    "10_core_gene_expression_by_epithelial_core_state.pdf"
  ),
  plot = p_core_genes_by_cluster +
    p_core_genes_by_program,
  width = 16,
  height = 8
)

# ============================================================
# O. 输出sample×cluster×core gene汇总表
# ============================================================

core_expression <- FetchData(
  object = sc_core,
  vars = c(
    "sample_id",
    "epithelial_core_cluster",
    "relative_program_label",
    core_genes
  )
)

core_expression_summary <- core_expression %>%
  pivot_longer(
    cols = all_of(core_genes),
    names_to = "gene",
    values_to = "log_normalized_expression"
  ) %>%
  group_by(
    sample_id,
    epithelial_core_cluster,
    relative_program_label,
    gene
  ) %>%
  summarise(
    cell_number = n(),
    percent_expressed = round(
      100 * mean(log_normalized_expression > 0),
      2
    ),
    mean_log_normalized_expression = mean(
      log_normalized_expression
    ),
    .groups = "drop"
  ) %>%
  arrange(
    sample_id,
    epithelial_core_cluster,
    gene
  )

write.csv(
  core_expression_summary,
  file.path(
    table_dir,
    "10_core_gene_expression_by_sample_cluster_and_program.csv"
  ),
  row.names = FALSE
)

# ============================================================
# P. 保存对象和环境信息
# ============================================================

sc_core$analysis_stage <- "manual_review_epithelial_core_complete"

saveRDS(
  sc_core,
  file.path(
    object_dir,
    "10_manual_review_epithelial_core.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "10_sessionInfo.txt"
  )
)

# ============================================================
# Q. 最终提示
# ============================================================

message("\n============================================================")
message("10_manual_review_epithelial_core.R运行完成。")
message("")
message("已保存对象：")
message("results/objects/10_manual_review_epithelial_core.rds")
message("")
message("请重点查看：")
message("1. results/figures/10_manual_review_cluster_decision_UMAP.pdf")
message("2. results/figures/10_epithelial_core_UMAP_sample_cluster_program.pdf")
message("3. results/figures/10_epithelial_core_relative_program_scores_UMAP.pdf")
message("4. results/figures/10_core_gene_expression_by_epithelial_core_state.pdf")
message("5. results/tables/10_manual_review_cell_summary.csv")
message("6. results/tables/10_epithelial_core_relative_program_summary.csv")
message("7. results/tables/10_epithelial_core_by_sample_cluster_and_program.csv")
message("============================================================\n")