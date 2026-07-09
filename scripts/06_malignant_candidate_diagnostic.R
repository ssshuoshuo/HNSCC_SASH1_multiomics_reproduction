# 06_malignant_candidate_diagnostic.R

# 本脚本功能：
# 1. 读取05人工注释后的Seurat object
# 2. 定义tumor epithelial candidate clusters
# 3. 定义salivary epithelial reference cluster
# 4. 写入malignant_status_diagnostic初步诊断标签
# 5. 汇总candidate和reference群的细胞数与样本组成
# 6. 汇总SASH1、MYH11、EMP1、COL1A1的表达
# 7. 输出target gene和diagnostic marker图表
# 8. 保存供后续CopyKAT和final malignant call使用的对象

# 本项目专用数据：
# GSE215403
# 12个OSCC单细胞样本：
# OSCC, scB1, scB2, scB5, scB7, scB8,
# scB9, scB10, scB12, scB13, scB14, scB15
#
# 本脚本定义的是tumor epithelial candidate和reference。
# 最终malignant/non-malignant判定不在本脚本完成。
# 后续步骤会结合CopyKAT inferred CNV和marker evidence进一步确认。
#
# 通用代码修改位置：
# 1. 换数据集时：
#    修改input_object_file和cluster_column
#
# 2. 换聚类分辨率时：
#    修改cluster_column
#
# 3. 换candidate/reference定义时：
#    修改tumor_candidate_clusters和salivary_reference_clusters
#
# 4. 换关注基因或诊断marker时：
#    修改target_genes、epithelial_markers、salivary_markers和tumor_state_markers


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
# C. 读取05对象
# ============================================================

input_object_file <- file.path(
  object_dir,
  "05_manual_annotated_before_malignant_call.rds"
)

if (!file.exists(input_object_file)) {
  
  stop(
    paste0(
      "找不到05对象：\n",
      input_object_file,
      "\n请先运行05_manual_annotation_and_target_gene_summary.R"
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
      "缺少metadata列：",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

message("读取完成。")
message("当前细胞数：", ncol(sc))

# ============================================================
# D. 定义tumor candidate和epithelial reference
# ============================================================

# 基于05的人工annotation和marker诊断结果：
#
# 2=Differentiated epithelial tumor candidate
# 3=Cycling epithelial tumor candidate
# 4=Cancer-testis epithelial tumor candidate
# 6=Epithelial tumor candidate
# 11=Epithelial tumor candidate
#
# 15=Salivary epithelial normal-like
#
# 这里的tumor candidate仍是候选标签。
# 后续需要CopyKAT inferred CNV结果进一步支持。

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
# E. 输出candidate细胞数与样本组成
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
    "06_malignant_candidate_cluster_summary.csv"
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
    "06_malignant_candidate_by_sample.csv"
  ),
  row.names = FALSE
)

print(candidate_summary)

# ============================================================
# F. 诊断状态UMAP
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
    "06_tumor_candidate_and_salivary_reference_UMAP.pdf"
  ),
  plot = p_diagnostic_umap,
  width = 12,
  height = 8
)

ggsave(
  filename = file.path(
    figure_dir,
    "06_tumor_candidate_and_salivary_reference_UMAP.png"
  ),
  plot = p_diagnostic_umap,
  width = 12,
  height = 8,
  dpi = 400
)

# ============================================================
# G. 定义诊断基因
# ============================================================

# target_genes：
# 论文主线关注的4个核心基因。
#
# epithelial_markers：
# 辅助确认上皮身份。
#
# salivary_markers：
# 辅助确认salivary/normal-like epithelial reference。
#
# tumor_state_markers：
# 辅助确认cycling、cancer-testis和分化上皮肿瘤状态。
#
# 这些marker用于诊断展示，不直接构成最终CNV判定。

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
    "06_diagnostic_gene_check.csv"
  ),
  row.names = FALSE
)

# ============================================================
# H. 提取表达数据
# ============================================================

# 这里只提取少量marker和metadata，避免复制完整表达矩阵。

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
# I. 核心基因按诊断状态汇总
# ============================================================

# 本表按diagnostic status汇总4个核心基因表达。
# 该表用于描述性比较。
# 更严格的统计比较应基于sample-aware pseudobulk表。

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
    "06_target_gene_expression_by_diagnostic_status.csv"
  ),
  row.names = FALSE
)

print(target_by_status)

# ============================================================
# J. 核心基因sample-aware pseudobulk汇总
# ============================================================

# 每个sample×status×gene得到一个平均表达值。
# 后续如果做统计检验，应优先基于这张表，
# 避免把每个细胞都当作独立生物学重复。

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
    "06_target_gene_pseudobulk_by_sample.csv"
  ),
  row.names = FALSE
)

# ============================================================
# K. 目标基因DotPlot与VlnPlot
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
    "06_target_genes_DotPlot_by_diagnostic_status.pdf"
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
    "06_target_genes_VlnPlot_by_diagnostic_status.pdf"
  ),
  plot = p_target_violin,
  width = 14,
  height = 10
)

# ============================================================
# L. 关键诊断marker DotPlot
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
    "06_epithelial_salivary_tumor_marker_DotPlot.pdf"
  ),
  plot = p_marker_dotplot,
  width = 18,
  height = 7
)

# ============================================================
# M. 保存诊断对象和环境信息
# ============================================================

sc$analysis_stage <- "malignant_candidate_diagnostic_before_CNV"

saveRDS(
  sc,
  file.path(
    object_dir,
    "06_malignant_candidate_diagnostic.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "06_sessionInfo.txt"
  )
)

# ============================================================
# N. 最终提示
# ============================================================

message("\n============================================================")
message("06_malignant_candidate_diagnostic.R 运行完成。")
message("")
message("已保存对象：")
message("results/objects/06_malignant_candidate_diagnostic.rds")
message("")
message("请重点查看：")
message("1. results/figures/06_tumor_candidate_and_salivary_reference_UMAP.pdf")
message("2. results/figures/06_target_genes_DotPlot_by_diagnostic_status.pdf")
message("3. results/figures/06_target_genes_VlnPlot_by_diagnostic_status.pdf")
message("4. results/figures/06_epithelial_salivary_tumor_marker_DotPlot.pdf")
message("5. results/tables/06_target_gene_expression_by_diagnostic_status.csv")
message("6. results/tables/06_target_gene_pseudobulk_by_sample.csv")
message("============================================================\n")