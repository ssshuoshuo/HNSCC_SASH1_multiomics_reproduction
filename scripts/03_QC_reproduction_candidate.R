# 03_QC_reproduction_candidate.R

# 本脚本功能：
# 1. 从02保存的原始Seurat object开始做正式QC过滤
# 2. 使用固定、可解释、可复现的QC阈值
# 3. 为每个细胞生成QC通过/未通过判定
# 4. 输出每个样本过滤前后的细胞数和保留率
# 5. 输出过滤前后的QC小提琴图
# 6. 保存供后续标准Seurat分析、细胞注释和恶性细胞分析使用的QC对象
# 7. 记录本步骤使用的QC参数和session信息

# 本项目专用数据：
# GSE215403
# 12个OSCC单细胞样本：
# OSCC, scB1, scB2, scB5, scB7, scB8,
# scB9, scB10, scB12, scB13, scB14, scB15
#
# 本脚本使用固定QC阈值：
# nFeature_RNA：200-8000
# nCount_RNA：500-100000
# percent.mt：≤20
#
# 通用代码修改位置：
# 1. 换数据集时：
#    修改input_object_file和sample_order
#
# 2. 换QC策略时：
#    修改analysis_preset对应的阈值
#
# 3. 换组织、平台或物种时：
#    重点重新检查min_nFeature、max_nFeature、
#    min_nCount、max_nCount和max_percent_mt


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
# C. 读取02原始对象
# ============================================================

input_object_file <- file.path(
  object_dir,
  "02_raw_before_QC_filtering.rds"
)

if (!file.exists(input_object_file)) {
  
  stop(
    paste0(
      "找不到02的原始对象：\n",
      input_object_file,
      "\n请先完整运行02_read_and_QC_scRNA.R"
    )
  )
}

sc_raw <- readRDS(input_object_file)

message("读取完成。")
message("原始细胞数：", ncol(sc_raw))
message("原始基因数：", nrow(sc_raw))

# 检查后续QC过滤所需的metadata列。
required_metadata <- c(
  "sample_id",
  "nFeature_RNA",
  "nCount_RNA",
  "percent.mt",
  "percent.ribo"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(sc_raw@meta.data)
)

if (length(missing_metadata) > 0) {
  
  stop(
    paste0(
      "缺少metadata列：",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

# ============================================================
# D. 固定样本顺序
# ============================================================

# 该顺序用于表格和图中的样本展示。
# 如果当前对象中出现额外样本，会自动追加到顺序末尾。

sample_order <- c(
  "OSCC",
  "scB1",
  "scB2",
  "scB5",
  "scB7",
  "scB8",
  "scB9",
  "scB10",
  "scB12",
  "scB13",
  "scB14",
  "scB15"
)

observed_samples <- unique(as.character(sc_raw$sample_id))

if (!all(observed_samples %in% sample_order)) {
  
  warning(
    "sample_order中可能缺少部分当前对象的样本名；",
    "将自动补到顺序末尾。"
  )
  
  sample_order <- c(
    sample_order,
    setdiff(observed_samples, sample_order)
  )
}

sc_raw$sample_id <- factor(
  sc_raw$sample_id,
  levels = sample_order
)

# ============================================================
# E. 设置QC参数
# ============================================================

analysis_preset <- "reproduction_candidate"

# reproduction_candidate：
# 当前复现流程使用的主QC参数。
# 目标是在保留主要OSCC细胞群的同时，
# 去除低复杂度、高线粒体比例和极端高UMI细胞。
#
# strict：
# 更严格的敏感性分析参数。
#
# custom：
# 预留给其他项目或手动调整参数使用。

if (analysis_preset == "reproduction_candidate") {
  
  min_nFeature <- 200
  max_nFeature <- 8000
  
  min_nCount <- 500
  max_nCount <- 100000
  
  max_percent_mt <- 20
}

if (analysis_preset == "strict") {
  
  min_nFeature <- 200
  max_nFeature <- 6000
  
  min_nCount <- 500
  max_nCount <- 80000
  
  max_percent_mt <- 20
}

if (analysis_preset == "custom") {
  
  min_nFeature <- 200
  max_nFeature <- 8000
  
  min_nCount <- 500
  max_nCount <- 100000
  
  max_percent_mt <- 20
}

message("\n============================================================")
message("当前QC preset：", analysis_preset)
message("nFeature_RNA：", min_nFeature, " - ", max_nFeature)
message("nCount_RNA：", min_nCount, " - ", max_nCount)
message("percent.mt：≤ ", max_percent_mt)
message("============================================================\n")

# ============================================================
# F. 生成每细胞QC判定
# ============================================================

# 每个细胞分别判断是否通过feature、UMI和线粒体比例阈值。
# qc_pass_reproduction为最终综合QC结果。

qc_cell_table <- sc_raw@meta.data %>%
  mutate(
    cell_id = rownames(sc_raw@meta.data),
    
    pass_min_feature = nFeature_RNA >= min_nFeature,
    pass_max_feature = nFeature_RNA <= max_nFeature,
    
    pass_min_count = nCount_RNA >= min_nCount,
    pass_max_count = nCount_RNA <= max_nCount,
    
    pass_mt = percent.mt <= max_percent_mt,
    
    qc_pass_reproduction = pass_min_feature &
      pass_max_feature &
      pass_min_count &
      pass_max_count &
      pass_mt
  )

# ============================================================
# G. 输出样本级QC统计表
# ============================================================

# 输出每个样本过滤前后的细胞数、各类未通过原因和保留率。
# fail_min_feature、fail_max_feature、fail_high_mt等列可用于判断
# 哪一种QC规则对不同样本影响最大。

qc_summary <- qc_cell_table %>%
  group_by(sample_id) %>%
  summarise(
    
    cells_before_QC = n(),
    
    fail_min_feature = sum(!pass_min_feature),
    fail_max_feature = sum(!pass_max_feature),
    
    fail_min_count = sum(!pass_min_count),
    fail_max_count = sum(!pass_max_count),
    
    fail_high_mt = sum(!pass_mt),
    
    cells_removed_total = sum(!qc_pass_reproduction),
    cells_after_QC = sum(qc_pass_reproduction),
    
    retention_rate_percent = round(
      100 * cells_after_QC / cells_before_QC,
      2
    ),
    
    median_nFeature_before = median(nFeature_RNA),
    median_nCount_before = median(nCount_RNA),
    median_percent_mt_before = median(percent.mt),
    
    .groups = "drop"
    
  ) %>%
  arrange(match(as.character(sample_id), sample_order))

print(qc_summary)

write.csv(
  qc_summary,
  file.path(
    table_dir,
    "03_QC_reproduction_candidate_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# H. 构建QC过滤后的Seurat object
# ============================================================

cells_to_keep <- qc_cell_table$cell_id[
  qc_cell_table$qc_pass_reproduction
]

sc_qc_reproduction <- subset(
  sc_raw,
  cells = cells_to_keep
)

sc_qc_reproduction$qc_preset <- analysis_preset
sc_qc_reproduction$qc_pass_reproduction <- TRUE
sc_qc_reproduction$analysis_stage <- "QC_reproduction_candidate"

message("\n============================================================")
message("03 QC完成")
message("原始总细胞数：", ncol(sc_raw))
message("03后总细胞数：", ncol(sc_qc_reproduction))
message(
  "总体保留率：",
  round(
    100 * ncol(sc_qc_reproduction) / ncol(sc_raw),
    2
  ),
  "%"
)
message("============================================================\n")

# ============================================================
# I. 过滤前后QC图
# ============================================================

p_feature_before <- VlnPlot(
  sc_raw,
  features = "nFeature_RNA",
  group.by = "sample_id",
  pt.size = 0
) +
  NoLegend() +
  ggtitle("Before QC: detected genes per cell")

p_feature_after <- VlnPlot(
  sc_qc_reproduction,
  features = "nFeature_RNA",
  group.by = "sample_id",
  pt.size = 0
) +
  NoLegend() +
  ggtitle("03 QC: detected genes per cell")

p_count_before <- VlnPlot(
  sc_raw,
  features = "nCount_RNA",
  group.by = "sample_id",
  pt.size = 0
) +
  NoLegend() +
  ggtitle("Before QC: total UMI counts per cell")

p_count_after <- VlnPlot(
  sc_qc_reproduction,
  features = "nCount_RNA",
  group.by = "sample_id",
  pt.size = 0
) +
  NoLegend() +
  ggtitle("03 QC: total UMI counts per cell")

p_mt_before <- VlnPlot(
  sc_raw,
  features = "percent.mt",
  group.by = "sample_id",
  pt.size = 0
) +
  NoLegend() +
  ggtitle("Before QC: mitochondrial percentage")

p_mt_after <- VlnPlot(
  sc_qc_reproduction,
  features = "percent.mt",
  group.by = "sample_id",
  pt.size = 0
) +
  NoLegend() +
  ggtitle("03 QC: mitochondrial percentage")

p_qc_before_after <- (
  p_feature_before + p_feature_after
) / (
  p_count_before + p_count_after
) / (
  p_mt_before + p_mt_after
)

ggsave(
  filename = file.path(
    figure_dir,
    "03_QC_before_after_reproduction_candidate.pdf"
  ),
  plot = p_qc_before_after,
  width = 20,
  height = 18
)

ggsave(
  filename = file.path(
    figure_dir,
    "03_QC_before_after_reproduction_candidate.png"
  ),
  plot = p_qc_before_after,
  width = 20,
  height = 18,
  dpi = 300
)

# ============================================================
# J. 保存对象、参数和环境信息
# ============================================================

qc_parameters <- data.frame(
  analysis_preset = analysis_preset,
  min_nFeature = min_nFeature,
  max_nFeature = max_nFeature,
  min_nCount = min_nCount,
  max_nCount = max_nCount,
  max_percent_mt = max_percent_mt,
  stringsAsFactors = FALSE
)

write.csv(
  qc_parameters,
  file.path(
    table_dir,
    "03_QC_reproduction_candidate_parameters.csv"
  ),
  row.names = FALSE
)

saveRDS(
  sc_qc_reproduction,
  file.path(
    object_dir,
    "03_QC_reproduction_candidate.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "03_sessionInfo.txt"
  )
)

# ============================================================
# K. 最终提示
# ============================================================

message("\n============================================================")
message("03_QC_reproduction_candidate.R 运行完成。")
message("")
message("已保存对象：")
message("results/objects/03_QC_reproduction_candidate.rds")
message("")
message("请重点查看：")
message("1. results/tables/03_QC_reproduction_candidate_summary.csv")
message("2. results/tables/03_QC_reproduction_candidate_parameters.csv")
message("3. results/figures/03_QC_before_after_reproduction_candidate.pdf")
message("============================================================\n")