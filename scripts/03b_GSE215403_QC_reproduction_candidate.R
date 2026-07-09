# ============================================================
# 03b_QC_reproduction_candidate.R
#
# 目的：
# 1. 从 02 原始对象重新做最终复现候选版 QC
# 3. 使用固定、可解释、可迁移的阈值
# 4. 保存供 04b / 05b / malignant-cell 分析使用的新对象
#
# 重要说明：
# 本脚本不覆盖 03 的结果。
#
#
# 本版 03b：
# 使用明确的固定阈值，方便复现、报告和跨项目迁移。
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

# 参数可替换区：
#
# analysis_preset:
#   "reproduction_candidate"
#   "strict"
#   "custom"
#
# 更换疾病 / 平台 / 组织时，最需要重新检查：
#   min_nFeature
#   max_nFeature
#   min_nCount
#   max_nCount
#   max_percent_mt
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
# C. 读取 02 原始对象
# ============================================================

input_object_file <- file.path(
  object_dir,
  "02_GSE215403_raw_before_QC_filtering.rds"
)

if (!file.exists(input_object_file)) {
  
  stop(
    paste0(
      "找不到 02 的原始对象：\n",
      input_object_file,
      "\n请先完整运行 02_read_and_QC_scRNA.R"
    )
  )
}

sc_raw <- readRDS(input_object_file)

message("读取完成。")
message("原始细胞数：", ncol(sc_raw))
message("原始基因数：", nrow(sc_raw))

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
      "缺少 metadata 列：",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

# ============================================================
# D. 固定样本顺序
# ============================================================

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
    "sample_order 中可能缺少部分当前对象的样本名；",
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
# E. QC 参数 preset
# ============================================================

analysis_preset <- "reproduction_candidate"

# ------------------------------------------------------------
# reproduction_candidate：
#
# 用于当前 GSE215403 复现候选版。
#
# 目的：
# 1. 保留主要 OSCC / epithelial 群；
# 2. 去除低复杂度和高 mt 细胞；
# 4. 为后续 doublet 诊断、Harmony、manual annotation 打基础。
# ------------------------------------------------------------

if (analysis_preset == "reproduction_candidate") {
  
  min_nFeature <- 200
  max_nFeature <- 8000
  
  min_nCount <- 500
  max_nCount <- 100000
  
  max_percent_mt <- 20
}

# ------------------------------------------------------------
# strict：
#
# 更严格的版本，适合后续发现高 UMI 尾部仍明显、
# 或 doublet / 异常小岛很多时进行敏感性分析。
# ------------------------------------------------------------

if (analysis_preset == "strict") {
  
  min_nFeature <- 200
  max_nFeature <- 6000
  
  min_nCount <- 500
  max_nCount <- 80000
  
  max_percent_mt <- 20
}

# ------------------------------------------------------------
# custom：
#
# 其他项目可直接在这里手动填自己的阈值。
# ------------------------------------------------------------

if (analysis_preset == "custom") {
  
  min_nFeature <- 200
  max_nFeature <- 8000
  
  min_nCount <- 500
  max_nCount <- 100000
  
  max_percent_mt <- 20
}

message("\n============================================================")
message("当前 QC preset：", analysis_preset)
message("nFeature_RNA：", min_nFeature, " - ", max_nFeature)
message("nCount_RNA：", min_nCount, " - ", max_nCount)
message("percent.mt：≤ ", max_percent_mt)
message("============================================================\n")

# ============================================================
# F. 生成每细胞 QC 判定
# ============================================================

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
# G. 样本级 QC 统计
# ============================================================

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
    "03b_QC_reproduction_candidate_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# ============================================================

old_qc_summary_file <- file.path(
  table_dir,
  "03_QC_filter_summary_by_sample.csv"
)

if (file.exists(old_qc_summary_file)) {
  
  old_qc_summary <- read.csv(
    old_qc_summary_file,
    stringsAsFactors = FALSE
  )
  
  qc_comparison <- qc_summary %>%
    mutate(
      sample_id = as.character(sample_id)
    ) %>%
    select(
      sample_id,
      cells_before_QC,
      cells_after_QC,
      retention_rate_percent
    ) %>%
    rename(
      cells_after_03b = cells_after_QC,
      retention_03b_percent = retention_rate_percent
    ) %>%
    left_join(
      old_qc_summary %>%
        select(
          sample_id,
          cells_after_QC,
          retention_rate_percent
        ) %>%
        rename(
          cells_after_old_03 = cells_after_QC,
          retention_old_03_percent = retention_rate_percent
        ),
      by = "sample_id"
    ) %>%
    mutate(
      difference_03b_minus_old03 =
        cells_after_03b - cells_after_old_03
    ) %>%
    arrange(match(sample_id, sample_order))
  
  write.csv(
    qc_comparison,
    file.path(
      table_dir,
      "03b_QC_comparison_old03_vs_reproduction_candidate.csv"
    ),
    row.names = FALSE
  )
  
  print(qc_comparison)
  
} else {
  
  warning(
    "没有找到旧版 03 的 summary 表，",
    "将跳过 03 vs 03b 对照。"
  )
  
  qc_comparison <- NULL
}

# ============================================================
# I. 构建新版 QC 对象
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
sc_qc_reproduction$analysis_stage <- "QC_reproduction_candidate_before_doublet"

message("\n============================================================")
message("03b QC 完成")
message("原始总细胞数：", ncol(sc_raw))
message("03b 后总细胞数：", ncol(sc_qc_reproduction))
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
# J. 过滤前后 QC 图
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
  ggtitle("03b QC: detected genes per cell")

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
  ggtitle("03b QC: total UMI counts per cell")

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
  ggtitle("03b QC: mitochondrial percentage")

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
    "03b_QC_before_after_reproduction_candidate.pdf"
  ),
  plot = p_qc_before_after,
  width = 20,
  height = 18
)

ggsave(
  filename = file.path(
    figure_dir,
    "03b_QC_before_after_reproduction_candidate.png"
  ),
  plot = p_qc_before_after,
  width = 20,
  height = 18,
  dpi = 300
)

# ============================================================
# ============================================================

if (!is.null(qc_comparison)) {
  
  comparison_long <- qc_comparison %>%
    select(
      sample_id,
      cells_before_QC,
      cells_after_old_03,
      cells_after_03b
    ) %>%
    pivot_longer(
      cols = c(
        cells_before_QC,
        cells_after_old_03,
        cells_after_03b
      ),
      names_to = "stage",
      values_to = "cell_number"
    ) %>%
    mutate(
      stage = recode(
        stage,
        "cells_before_QC" = "Before QC",
        "cells_after_old_03" = "Old 03 q99 QC",
        "cells_after_03b" = "03b reproduction candidate"
      ),
      sample_id = factor(
        sample_id,
        levels = sample_order
      )
    )
  
  p_qc_comparison <- ggplot(
    comparison_long,
    aes(
      x = sample_id,
      y = cell_number,
      fill = stage
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.8)
    ) +
    labs(
      title = "Cell number comparison: old 03 versus 03b QC",
      x = "Sample",
      y = "Number of cells",
      fill = NULL
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_text(
        hjust = 0.5
      )
    )
  
  ggsave(
    filename = file.path(
      figure_dir,
      "03b_cell_number_old03_vs_reproduction_candidate.pdf"
    ),
    plot = p_qc_comparison,
    width = 14,
    height = 7
  )
}

# ============================================================
# L. 保存对象与参数
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
    "03b_QC_reproduction_candidate_parameters.csv"
  ),
  row.names = FALSE
)

saveRDS(
  sc_qc_reproduction,
  file.path(
    object_dir,
    "03b_GSE215403_QC_reproduction_candidate.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "03b_sessionInfo.txt"
  )
)

# ============================================================
# M. 完成提示
# ============================================================

message("\n============================================================")
message("03b_QC_reproduction_candidate.R 运行完成。")
message("")
message("已保存对象：")
message("results/objects/03b_GSE215403_QC_reproduction_candidate.rds")
message("")
message("重点查看：")
message("1. results/tables/03b_QC_reproduction_candidate_summary.csv")
message("2. results/tables/03b_QC_comparison_old03_vs_reproduction_candidate.csv")
message("3. results/figures/03b_QC_before_after_reproduction_candidate.pdf")
message("4. results/figures/03b_cell_number_old03_vs_reproduction_candidate.pdf")
message("============================================================\n")