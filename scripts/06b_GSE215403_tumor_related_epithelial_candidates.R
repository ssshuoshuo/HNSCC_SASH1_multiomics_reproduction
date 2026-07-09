# ============================================================
#
# 功能：
# 1. 读取 05c 人工注释对象
# 2. 定义 tumor epithelial candidate、salivary epithelial reference
# 3. 不复制 Seurat 子对象，避免本地内存问题
# 4. 汇总 SASH1 / MYH11 / EMP1 / COL1A1 的表达
# 5. 输出后续 CNV / malignant call 所需的诊断表和图
#
# 注意：
# 最终恶性判定需要后续 CNV 或其他证据支持。
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
  "patchwork",
  "ggrastr"
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
# C. 读取 05c 对象
# ============================================================

plot_ready_file <- file.path(
  object_dir,
  "05c_GSE215403_manual_annotated_plot_ready.rds"
)

base_file <- file.path(
  object_dir,
  "05c_GSE215403_manual_annotated_before_malignant_call.rds"
)

input_candidates <- c(
  plot_ready_file,
  base_file
)

input_object_file <- input_candidates[
  file.exists(input_candidates)
][1]

if (is.na(input_object_file)) {
  
  stop(
    paste0(
      "找不到 05c 对象。\n",
      "请先运行 05c_manual_annotation_and_target_gene_summary.R"
    )
  )
}

sc <- readRDS(input_object_file)

DefaultAssay(sc) <- "RNA"

cluster_column <- "cluster_res_0.2"

required_metadata <- c(
  cluster_column,
  "celltype_manual",
  "sample_id"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(sc@meta.data)
)

if (length(missing_metadata) > 0) {
  
  stop(
    paste0(
      "缺少 metadata 列：",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

message("读取完成。")
message("当前细胞数：", ncol(sc))

# ============================================================
# D. 定义 malignant candidate 与 epithelial reference
# ============================================================
#
# 基于 05b / 05c 的人工 marker 注释：
#
# 2  = Differentiated epithelial tumor candidate
# 3  = Cycling epithelial tumor candidate
# 4  = Cancer-testis epithelial tumor candidate
# 6  = Epithelial tumor candidate
# 11 = Epithelial tumor candidate
#
# 15 = Salivary epithelial normal-like
#
# 通用项目：
# 换数据集时，这一段 cluster 编号必须重新判断，
# ============================================================

tumor_candidate_clusters <- c(
  "2",
  "3",
  "4",
  "6",
  "11"
)

salivary_reference_clusters <- c(
  "15"
)

cluster_vector <- as.character(
  sc[[cluster_column, drop = TRUE]]
)

diagnostic_status <- rep(
  "Other_cells",
  length(cluster_vector)
)

diagnostic_status[
  cluster_vector %in% tumor_candidate_clusters
] <- "Tumor_epithelial_candidate"

diagnostic_status[
  cluster_vector %in% salivary_reference_clusters
] <- "Salivary_epithelial_reference"

diagnostic_status <- factor(
  diagnostic_status,
  levels = c(
    "Tumor_epithelial_candidate",
    "Salivary_epithelial_reference",
    "Other_cells"
  )
)

names(diagnostic_status) <- colnames(sc)

sc[["malignant_status_diagnostic"]] <- diagnostic_status

# ============================================================
# E. 输出 candidate 群细胞数与样本组成
# ============================================================

candidate_summary <- sc@meta.data %>%
  mutate(
    cluster = as.character(
      .data[[cluster_column]]
    )
  ) %>%
  count(
    malignant_status_diagnostic,
    cluster,
    celltype_manual,
    name = "cell_number"
  ) %>%
  mutate(
    percent_within_status = round(
      100 * cell_number /
        sum(cell_number),
      2
    )
  ) %>%
  arrange(
    malignant_status_diagnostic,
    suppressWarnings(as.numeric(cluster))
  )

write.csv(
  candidate_summary,
  file.path(
    table_dir,
    "06b_malignant_candidate_cluster_summary.csv"
  ),
  row.names = FALSE
)

candidate_by_sample <- sc@meta.data %>%
  mutate(
    cluster = as.character(
      .data[[cluster_column]]
    )
  ) %>%
  filter(
    malignant_status_diagnostic != "Other_cells"
  ) %>%
  count(
    sample_id,
    malignant_status_diagnostic,
    cluster,
    name = "cell_number"
  ) %>%
  group_by(
    malignant_status_diagnostic,
    cluster
  ) %>%
  mutate(
    percent_within_cluster = round(
      100 * cell_number / sum(cell_number),
      2
    )
  ) %>%
  ungroup() %>%
  arrange(
    malignant_status_diagnostic,
    suppressWarnings(as.numeric(cluster)),
    desc(cell_number)
  )

write.csv(
  candidate_by_sample,
  file.path(
    table_dir,
    "06b_malignant_candidate_by_sample.csv"
  ),
  row.names = FALSE
)

print(candidate_summary)

# ============================================================
# F. 诊断状态 UMAP
# ============================================================

diagnostic_colors <- c(
  "Tumor_epithelial_candidate" = "#D73027",
  "Salivary_epithelial_reference" = "#7570B3",
  "Other_cells" = "#BDBDBD"
)

p_diagnostic_umap <- DimPlot(
  object = sc,
  reduction = "umap_pca",
  group.by = "malignant_status_diagnostic",
  cols = diagnostic_colors,
  pt.size = 0.22,
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  ggtitle("Tumor epithelial candidates and salivary epithelial reference") +
  theme_classic(base_size = 13) +
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
    "06b_tumor_candidate_and_salivary_reference_UMAP.pdf"
  ),
  plot = p_diagnostic_umap,
  width = 12,
  height = 8
)

ggsave(
  filename = file.path(
    figure_dir,
    "06b_tumor_candidate_and_salivary_reference_UMAP.png"
  ),
  plot = p_diagnostic_umap,
  width = 12,
  height = 8,
  dpi = 400
)

# ============================================================
# G. 定义诊断基因
# ============================================================
#
# target_genes：
# 论文中的四个核心基因。
#
# epithelial / salivary / tumor-state marker：
# 用于确认 candidate 与 normal-like reference 的身份。
#
# 这些 marker 用于诊断，不直接构成最终 CNV 判定。
# ============================================================

target_genes <- c(
  "SASH1",
  "MYH11",
  "EMP1",
  "COL1A1"
)

epithelial_markers <- c(
  "EPCAM",
  "KRT5",
  "KRT14",
  "KRT17",
  "KRT19",
  "TP63"
)

salivary_markers <- c(
  "MUC7",
  "PIP",
  "PRR4",
  "CRISP3",
  "STATH",
  "FDCSP"
)

tumor_state_markers <- c(
  "MKI67",
  "TOP2A",
  "CDC20",
  "CTAG2",
  "MAGEA4",
  "MAGEB2",
  "MAGEA1",
  "IVL",
  "SPRR1B"
)

diagnostic_genes <- unique(
  c(
    target_genes,
    epithelial_markers,
    salivary_markers,
    tumor_state_markers
  )
)

diagnostic_genes_present <- intersect(
  diagnostic_genes,
  rownames(sc)
)

target_genes_present <- intersect(
  target_genes,
  rownames(sc)
)

write.csv(
  data.frame(
    gene = diagnostic_genes,
    present_in_matrix = diagnostic_genes %in% rownames(sc),
    stringsAsFactors = FALSE
  ),
  file.path(
    table_dir,
    "06b_diagnostic_gene_check.csv"
  ),
  row.names = FALSE
)

# ============================================================
# H. 提取表达数据
# ============================================================
#
# 这里只提取少量 marker 和 metadata，
# 本地内存负担较低。
# ============================================================

expression_data <- FetchData(
  object = sc,
  vars = c(
    "sample_id",
    cluster_column,
    "celltype_manual",
    "malignant_status_diagnostic",
    diagnostic_genes_present
  )
)

expression_long <- expression_data %>%
  mutate(
    cluster = as.character(
      .data[[cluster_column]]
    )
  ) %>%
  select(
    -all_of(cluster_column)
  ) %>%
  pivot_longer(
    cols = all_of(diagnostic_genes_present),
    names_to = "gene",
    values_to = "expression"
  )

# ============================================================
# I. 四个核心基因：按诊断状态汇总
# ============================================================
#
#
# 真正统计比较时要按 sample 做 pseudobulk，
# 避免把每个细胞都当独立生物学重复。
# ============================================================

target_by_status <- expression_long %>%
  filter(
    gene %in% target_genes_present
  ) %>%
  group_by(
    malignant_status_diagnostic,
    gene
  ) %>%
  summarise(
    cell_number = n(),
    percent_expressed = round(
      100 * mean(expression > 0),
      2
    ),
    mean_expression = mean(expression),
    median_expression = median(expression),
    .groups = "drop"
  )

write.csv(
  target_by_status,
  file.path(
    table_dir,
    "06b_target_gene_expression_by_diagnostic_status.csv"
  ),
  row.names = FALSE
)

print(target_by_status)

# ============================================================
# J. 四个核心基因：sample-aware pseudobulk 汇总
# ============================================================
#
# 每个 sample × status × gene 得到一个平均表达值。
# 后续若做统计检验，应基于这张表，而非逐细胞检验。
# ============================================================

target_by_sample <- expression_long %>%
  filter(
    gene %in% target_genes_present,
    malignant_status_diagnostic != "Other_cells"
  ) %>%
  group_by(
    sample_id,
    malignant_status_diagnostic,
    gene
  ) %>%
  summarise(
    cell_number = n(),
    percent_expressed = round(
      100 * mean(expression > 0),
      2
    ),
    mean_expression = mean(expression),
    median_expression = median(expression),
    .groups = "drop"
  )

write.csv(
  target_by_sample,
  file.path(
    table_dir,
    "06b_target_gene_pseudobulk_by_sample.csv"
  ),
  row.names = FALSE
)

# ============================================================
# K. 目标基因 DotPlot 与 violin plot
# ============================================================

p_target_dotplot <- DotPlot(
  object = sc,
  features = target_genes_present,
  group.by = "malignant_status_diagnostic",
  assay = "RNA",
  dot.scale = 8
) +
  RotatedAxis() +
  ggtitle("Core genes by tumor-candidate diagnostic status") +
  theme_classic(base_size = 12)

ggsave(
  filename = file.path(
    figure_dir,
    "06b_target_genes_DotPlot_by_diagnostic_status.pdf"
  ),
  plot = p_target_dotplot,
  width = 12,
  height = 7
)

p_target_violin <- VlnPlot(
  object = sc,
  features = target_genes_present,
  group.by = "malignant_status_diagnostic",
  pt.size = 0,
  ncol = 2
) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

ggsave(
  filename = file.path(
    figure_dir,
    "06b_target_genes_VlnPlot_by_diagnostic_status.pdf"
  ),
  plot = p_target_violin,
  width = 14,
  height = 10
)

# ============================================================
# L. 关键诊断 marker DotPlot
# ============================================================

marker_dotplot_genes <- intersect(
  c(
    epithelial_markers,
    salivary_markers,
    tumor_state_markers
  ),
  rownames(sc)
)

p_marker_dotplot <- DotPlot(
  object = sc,
  features = marker_dotplot_genes,
  group.by = "malignant_status_diagnostic",
  assay = "RNA",
  dot.scale = 6
) +
  RotatedAxis() +
  ggtitle("Epithelial, salivary and tumor-state markers") +
  theme_classic(base_size = 11)

ggsave(
  filename = file.path(
    figure_dir,
    "06b_epithelial_salivary_tumor_marker_DotPlot.pdf"
  ),
  plot = p_marker_dotplot,
  width = 18,
  height = 7
)

# ============================================================
# M. 保存诊断对象
# ============================================================

sc$analysis_stage <- "malignant_candidate_diagnostic_before_CNV"

saveRDS(
  sc,
  file.path(
    object_dir,
    "06b_GSE215403_malignant_candidate_diagnostic.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "06b_sessionInfo.txt"
  )
)

# ============================================================
# N. 完成提示
# ============================================================

message("\n============================================================")
message("06b_malignant_candidate_diagnostic.R 运行完成。")
message("")
message("已保存对象：")
message("results/objects/06b_GSE215403_malignant_candidate_diagnostic.rds")
message("")
message("重点查看：")
message("1. results/figures/06b_tumor_candidate_and_salivary_reference_UMAP.pdf")
message("2. results/figures/06b_target_genes_DotPlot_by_diagnostic_status.pdf")
message("3. results/figures/06b_target_genes_VlnPlot_by_diagnostic_status.pdf")
message("4. results/figures/06b_epithelial_salivary_tumor_marker_DotPlot.pdf")
message("5. results/tables/06b_target_gene_expression_by_diagnostic_status.csv")
message("6. results/tables/06b_target_gene_pseudobulk_by_sample.csv")
message("============================================================\n")