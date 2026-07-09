# ============================================================
# 07d_manual_review_epithelial_core_local.R
#
# 目的：
# 1. 基于07c内部cluster marker结果，对恶性候选细胞做透明的
#    manual-review epithelial-core refinement；
# 2. 排除明确immune / macrophage / CAF-stromal /
#    lineage-ambiguous internal clusters；
# 3. 对保留的上皮核心集重新UMAP；
# 4. 使用相对标准化program score，而非原始平均表达，
#    描述Cycling / Cancer-testis / Squamous / Basal程序；
# 5. 比较SASH1、EMP1、MYH11、COL1A1在不同内部cluster与
#    相对program中的表达。
#
# 注意：
# 这里的core subset是“经内部marker审查后保留的上皮核心细胞集”，
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


# ============================================================
# A. 包与路径
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
# B. 读取07c对象
# ============================================================

input_file <- file.path(
  object_dir,
  "07c_GSE215403_malignant_epithelial_state_characterization.rds"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "找不到07c对象：\n",
      input_file
    )
  )
}

sc <- readRDS(input_file)

DefaultAssay(sc) <- "RNA"

if (!"seurat_clusters" %in% colnames(sc@meta.data)) {
  stop(
    "07c对象中缺少seurat_clusters。"
  )
}

if (!"sample_id" %in% colnames(sc@meta.data)) {
  stop(
    "07c对象中缺少sample_id。"
  )
}

message(
  "07d输入细胞数：",
  ncol(sc)
)

# ============================================================
# C. 基于07c marker审查的cluster决策表
# ============================================================

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
    "07d_manual_review_cluster_decision.csv"
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
# D. 写入审查状态，并输出筛选概览
# ============================================================

sc$internal_cluster_07c <- as.character(
  sc$seurat_clusters
)

decision_map <- setNames(
  manual_cluster_review$review_decision,
  manual_cluster_review$internal_cluster
)

sc$manual_review_status <- unname(
  decision_map[
    sc$internal_cluster_07c
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
    internal_cluster_07c = as.character(
      internal_cluster_07c
    )
  ) %>%
  count(
    internal_cluster_07c,
    manual_review_status,
    name = "cell_number"
  ) %>%
  left_join(
    manual_cluster_review,
    by = c(
      "internal_cluster_07c" = "internal_cluster"
    )
  ) %>%
  mutate(
    percent_of_all_07c_cells = round(
      100 * cell_number / sum(cell_number),
      2
    )
  ) %>%
  arrange(
    suppressWarnings(
      as.numeric(internal_cluster_07c)
    )
  )

write.csv(
  review_summary,
  file.path(
    table_dir,
    "07d_manual_review_cell_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# E. 在07c现有UMAP上展示保留/排除决策
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
    "07c internal clusters: manual-review epithelial-core decision"
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
  group.by = "internal_cluster_07c",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.3,
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  ggtitle(
    "07c internal cluster identity"
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
    "07d_manual_review_cluster_decision_UMAP.pdf"
  ),
  plot = p_review_status +
    p_internal_cluster,
  width = 16,
  height = 8
)

# ============================================================
# F. 提取上皮核心集，并重新构建轻量对象
# ============================================================

retained_cells <- colnames(sc)[
  sc$internal_cluster_07c %in% retained_clusters
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
  project = "GSE215403_manual_review_epithelial_core"
)

rm(core_counts)
gc()

DefaultAssay(sc_core) <- "RNA"

# ============================================================
# G. 上皮核心集重新降维与内部聚类
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
# H. 相对program score
# ============================================================
#
# 与07c不同：
# 这里不再比较原始平均表达值。
# 而是先在上皮核心集内对每个program做z-score标准化，
# 再比较每个cluster的相对富集程度。
# ============================================================

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
    "07d_relative_program_gene_sets_available.csv"
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
# I. 以内部cluster为单位计算相对program富集
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
    "07d_epithelial_core_relative_program_summary.csv"
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
# J. 上皮核心集组成表
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
    "07d_epithelial_core_by_sample_cluster_and_program.csv"
  ),
  row.names = FALSE
)

# ============================================================
# K. 图：上皮核心集UMAP
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
    "07d_epithelial_core_UMAP_sample_cluster_program.pdf"
  ),
  plot = p_core_sample +
    p_core_cluster +
    p_core_program +
    plot_layout(ncol = 3),
  width = 22,
  height = 7
)

# ============================================================
# L. 图：相对program score UMAP
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
    "07d_epithelial_core_relative_program_scores_UMAP.pdf"
  ),
  plot = p_relative_program_scores,
  width = 13,
  height = 10
)

# ============================================================
# M. 核心基因表达
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
    "07d_core_gene_expression_by_epithelial_core_state.pdf"
  ),
  plot = p_core_genes_by_cluster +
    p_core_genes_by_program,
  width = 16,
  height = 8
)

# ============================================================
# N. Sample × cluster × core gene汇总
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
    "07d_core_gene_expression_by_sample_cluster_and_program.csv"
  ),
  row.names = FALSE
)

# ============================================================
# O. 保存对象与运行信息
# ============================================================

sc_core$analysis_stage <- "manual_review_epithelial_core_complete"

saveRDS(
  sc_core,
  file.path(
    object_dir,
    "07d_GSE215403_manual_review_epithelial_core.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "07d_sessionInfo.txt"
  )
)

message("\n============================================================")
message("07d_manual_review_epithelial_core_local.R运行完成。")
message("")
message("重点查看：")
message("1. results/figures/07d_manual_review_cluster_decision_UMAP.pdf")
message("2. results/figures/07d_epithelial_core_UMAP_sample_cluster_program.pdf")
message("3. results/figures/07d_epithelial_core_relative_program_scores_UMAP.pdf")
message("4. results/figures/07d_core_gene_expression_by_epithelial_core_state.pdf")
message("5. results/tables/07d_manual_review_cell_summary.csv")
message("6. results/tables/07d_epithelial_core_relative_program_summary.csv")
message("7. results/tables/07d_epithelial_core_by_sample_cluster_and_program.csv")
message("============================================================\n")