# 02_read_and_QC_scRNA.R

# 本脚本功能：
# 1. 读取01中整理好的12个标准10x单细胞文件夹
# 2. 为每个细胞添加sample_id、dataset、disease等metadata
# 3. 合并所有样本为一个原始Seurat object
# 4. 计算常用QC指标：
#    - nFeature_RNA：每个细胞检测到的基因数
#    - nCount_RNA：每个细胞的总UMI数
#    - percent.mt：线粒体基因比例
#    - percent.ribo：核糖体基因比例
# 5. 检查本项目关注基因是否存在于表达矩阵中
# 6. 输出每个样本的QC汇总统计表
# 7. 基于每个样本自身分布生成建议QC阈值
# 8. 输出QC小提琴图和散点图
# 9. 保存正式过滤前的原始Seurat object

# 本项目专用数据：
# GSE215403
# 12个OSCC单细胞样本：
# OSCC, scB1, scB2, scB5, scB7, scB8,
# scB9, scB10, scB12, scB13, scB14, scB15
#
# 本脚本只做读取、合并和QC描述。
# 正式低质量细胞过滤、双细胞检测、归一化、
# 整合、降维和聚类在后续脚本完成。
#
# 通用代码修改位置：
# 1. 换数据集时：
#    修改tenx_root、sample_ids、dataset_name、disease_name
#
# 2. 换物种时：
#    修改线粒体基因匹配规则：
#    人类通常为"^MT-"
#    小鼠通常为"^mt-"
#
# 3. 换关注基因时：
#    修改target_genes
#
# 4. 换10x features.tsv结构时：
#    检查Read10X()中的gene.column


# ============================================================
# A. 安装并加载本步骤需要的包
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 3600)

# 本脚本需要的R包：
#
# Seurat/SeuratObject：
# 读取10x数据、建立Seurat object、绘制QC图
#
# dplyr：
# 整理metadata和QC统计表
#
# ggplot2：
# 绘制QC散点图
#
# patchwork：
# 拼接多个Seurat图
#
# fs：
# 文件夹和文件路径操作

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "ggplot2",
  "patchwork",
  "fs"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    message("正在安装：", pkg)
    install.packages(pkg)
    
  }
}

library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(patchwork)
library(fs)

# ============================================================
# B. 项目路径与文件夹
# ============================================================

project_dir <- getwd()

tenx_root <- file.path(
  project_dir,
  "data",
  "processed",
  "scRNA_GSE215403",
  "10x_by_sample"
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

object_dir <- file.path(
  project_dir,
  "results",
  "objects"
)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# C. 定义样本信息
# ============================================================

# 本课题的12个原始样本ID。
# 这些名字必须与01中整理出的10x_by_sample文件夹名称完全一致。
#
# 通用修改位置：
# 换数据集时，直接替换sample_ids。
#
# 也可以从文件夹自动读取：
# sample_ids <- basename(list.dirs(tenx_root, recursive = FALSE))
#
# 本流程为了保证输出顺序稳定，手动指定样本顺序。

sample_ids <- c(
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

dataset_name <- "GSE215403"
disease_name <- "OSCC"

# -----------------------------
# C1. 检查每个样本文件夹是否存在
# -----------------------------

missing_samples <- sample_ids[
  !dir.exists(file.path(tenx_root, sample_ids))
]

if (length(missing_samples) > 0) {
  
  stop(
    paste0(
      "以下样本文件夹不存在：\n",
      paste(missing_samples, collapse = ", "),
      "\n请先运行01_download_and_prepare_scRNA.R"
    )
  )
}

# ============================================================
# D. 逐个读取10x数据并建立Seurat object
# ============================================================

# CreateSeuratObject初始门槛设置得较宽松：
#
# min.cells=3：
# 一个基因至少出现在3个细胞中才保留。
# 主要用于去除极端罕见、信息量很低的基因。
#
# min.features=100：
# 一个细胞至少检测到100个基因才进入原始对象。
#
# 本步骤不进行正式QC过滤。
# 正式筛选会在下一步根据样本分布和QC图完成。

seurat_list <- list()

for (sid in sample_ids) {
  
  message("\n========================================")
  message("正在读取样本：", sid)
  message("========================================")
  
  sample_dir <- file.path(
    tenx_root,
    sid
  )
  
  # Read10X()会读取标准10x文件：
  # matrix.mtx.gz
  # barcodes.tsv.gz
  # features.tsv.gz
  #
  # gene.column=2：
  # 10x features.tsv通常第二列是gene symbol。
  # 例如SASH1、COL1A1、EMP1等。
  #
  # 通用修改位置：
  # 如果第二列不是gene symbol，或数据使用Ensembl ID，
  # 需要修改gene.column，或在读取后进行ID转换。
  
  counts <- Read10X(
    data.dir = sample_dir,
    gene.column = 2
  )
  
  # 创建单个样本的Seurat object。
  obj <- CreateSeuratObject(
    counts = counts,
    project = dataset_name,
    min.cells = 3,
    min.features = 100
  )
  
  # 添加样本和数据集metadata。
  # 后续QC、Harmony整合、分样本作图和组成统计都依赖这些字段。
  
  obj$sample_id <- sid
  obj$dataset <- dataset_name
  obj$disease <- disease_name
  
  # 暂存到list，后续统一合并。
  seurat_list[[sid]] <- obj
  
  message(
    sid,
    "：保留 ",
    ncol(obj),
    " 个细胞；",
    nrow(obj),
    " 个基因"
  )
}

# ============================================================
# E. 合并所有样本
# ============================================================

# add.cell.ids会给细胞barcode加上样本前缀。
#
# 例如原始barcode：
# AAACCTGAGAGTAATC-1
#
# 合并后会变为：
# OSCC_AAACCTGAGAGTAATC-1
#
# 这一步很重要。
# 不同10x样本可能出现相同barcode。
# 加入样本前缀可以避免合并后细胞名称冲突。

sc_raw <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = sample_ids,
  project = dataset_name,
  merge.data = FALSE
)

message("\n========================================")
message("全部样本合并完成")
message("总细胞数：", ncol(sc_raw))
message("总基因数：", nrow(sc_raw))
message("========================================")

# ============================================================
# F. 计算QC指标
# ============================================================

# -----------------------------
# F1. 线粒体基因比例
# -----------------------------
#
# 人类线粒体基因通常以MT-开头：
# MT-CO1、MT-ND1、MT-ATP6等。
#
# 通用修改位置：
# 人类：pattern="^MT-"
# 小鼠：pattern="^mt-"
#
# 线粒体比例较高通常提示：
# - 细胞受损
# - RNA降解
# - 濒死细胞
#
# 具体阈值需要结合每个样本的真实分布决定。

sc_raw[["percent.mt"]] <- PercentageFeatureSet(
  sc_raw,
  pattern = "^MT-"
)

# -----------------------------
# F2. 核糖体基因比例
# -----------------------------
#
# 人类核糖体蛋白基因通常以RPL或RPS开头。
# 该指标可辅助判断细胞状态和样本差异。

sc_raw[["percent.ribo"]] <- PercentageFeatureSet(
  sc_raw,
  pattern = "^RP[SL]"
)

# -----------------------------
# F3. 检查目标基因是否存在
# -----------------------------
#
# 这一步确认本项目关注的基因符号是否能在表达矩阵中找到。
# 如果found_in_matrix为FALSE，后续FeaturePlot或DotPlot会找不到该基因。
#
# 通用修改位置：
# 换课题时，修改target_genes。

target_genes <- c(
  "SASH1",
  "MYH11",
  "EMP1",
  "COL1A1"
)

target_gene_check <- data.frame(
  gene = target_genes,
  found_in_matrix = target_genes %in% rownames(sc_raw),
  stringsAsFactors = FALSE
)

print(target_gene_check)

write.csv(
  target_gene_check,
  file.path(
    table_dir,
    "02_target_gene_check.csv"
  ),
  row.names = FALSE
)

if (!all(target_gene_check$found_in_matrix)) {
  
  warning(
    paste0(
      "有目标基因未在表达矩阵中找到。\n",
      "可能原因：gene symbol不一致、Ensembl ID、旧基因名，",
      "或Read10X(gene.column=2)不适合该数据。"
    )
  )
}

# ============================================================
# G. 输出样本级QC统计表
# ============================================================

# 输出每个样本的QC概览：
#
# raw_cells：
# 初步读取后保留的细胞数
#
# median/mean：
# 中位数和均值
#
# q01/q05/q95/q99：
# 分位数，用于观察低质量细胞和高UMI尾部
#
# 这些统计有助于判断：
# - 是否有测序深度明显偏低的样本
# - 是否有线粒体比例整体偏高的样本
# - 是否存在潜在doublet对应的高UMI/高feature尾部

meta_qc <- sc_raw@meta.data

sample_qc_summary <- meta_qc %>%
  mutate(
    cell_id = rownames(meta_qc)
  ) %>%
  group_by(sample_id) %>%
  summarise(
    raw_cells = n(),
    
    median_nFeature_RNA = median(nFeature_RNA),
    mean_nFeature_RNA = mean(nFeature_RNA),
    q01_nFeature_RNA = quantile(nFeature_RNA, 0.01),
    q05_nFeature_RNA = quantile(nFeature_RNA, 0.05),
    q95_nFeature_RNA = quantile(nFeature_RNA, 0.95),
    q99_nFeature_RNA = quantile(nFeature_RNA, 0.99),
    
    median_nCount_RNA = median(nCount_RNA),
    mean_nCount_RNA = mean(nCount_RNA),
    q01_nCount_RNA = quantile(nCount_RNA, 0.01),
    q05_nCount_RNA = quantile(nCount_RNA, 0.05),
    q95_nCount_RNA = quantile(nCount_RNA, 0.95),
    q99_nCount_RNA = quantile(nCount_RNA, 0.99),
    
    median_percent_mt = median(percent.mt),
    mean_percent_mt = mean(percent.mt),
    q01_percent_mt = quantile(percent.mt, 0.01),
    q05_percent_mt = quantile(percent.mt, 0.05),
    q95_percent_mt = quantile(percent.mt, 0.95),
    q99_percent_mt = quantile(percent.mt, 0.99),
    
    median_percent_ribo = median(percent.ribo),
    .groups = "drop"
  ) %>%
  arrange(match(sample_id, sample_ids))

print(sample_qc_summary)

write.csv(
  sample_qc_summary,
  file.path(
    table_dir,
    "02_sample_QC_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# H. 生成每个样本的建议QC阈值
# ============================================================

# 本表只是基于每个样本自身分布生成的建议范围。
# 下一步会结合QC图和该表确定正式筛选阈值。
#
# 建议规则：
#
# nFeature_RNA：
# 下限取max(200,1%分位数)
# 上限取99%分位数
#
# nCount_RNA：
# 下限取1%分位数
# 上限取99%分位数
#
# percent.mt：
# 上限取max(20,99%分位数)
#
# 不同肿瘤样本的测序深度、细胞状态和组织组成可能差异明显。
# 因此正式QC过滤不在本脚本中直接执行。

qc_threshold_suggestion <- meta_qc %>%
  group_by(sample_id) %>%
  summarise(
    suggested_min_nFeature = max(
      200,
      floor(quantile(nFeature_RNA, 0.01))
    ),
    
    suggested_max_nFeature = ceiling(
      quantile(nFeature_RNA, 0.99)
    ),
    
    suggested_min_nCount = floor(
      quantile(nCount_RNA, 0.01)
    ),
    
    suggested_max_nCount = ceiling(
      quantile(nCount_RNA, 0.99)
    ),
    
    suggested_max_percent_mt = max(
      20,
      ceiling(quantile(percent.mt, 0.99))
    ),
    
    stringsAsFactors = FALSE,
    .groups = "drop"
  ) %>%
  arrange(match(sample_id, sample_ids))

print(qc_threshold_suggestion)

write.csv(
  qc_threshold_suggestion,
  file.path(
    table_dir,
    "02_QC_threshold_suggestions.csv"
  ),
  row.names = FALSE
)

# ============================================================
# I. QC小提琴图
# ============================================================

# pt.size=0：
# 不显示单个细胞点。
# 细胞数量较多时，显示所有点会使图像过密且绘图较慢。
#
# group.by="sample_id"：
# 按样本分别展示QC分布。
# 不同病人的测序深度和组织质量可能不同，因此应先按样本查看。

p_vln_feature <- VlnPlot(
  sc_raw,
  features = "nFeature_RNA",
  group.by = "sample_id",
  pt.size = 0
) +
  NoLegend() +
  ggtitle("QC: detected genes per cell")

p_vln_count <- VlnPlot(
  sc_raw,
  features = "nCount_RNA",
  group.by = "sample_id",
  pt.size = 0
) +
  NoLegend() +
  ggtitle("QC: total UMI counts per cell")

p_vln_mt <- VlnPlot(
  sc_raw,
  features = "percent.mt",
  group.by = "sample_id",
  pt.size = 0
) +
  NoLegend() +
  ggtitle("QC: mitochondrial gene percentage")

p_vln_ribo <- VlnPlot(
  sc_raw,
  features = "percent.ribo",
  group.by = "sample_id",
  pt.size = 0
) +
  NoLegend() +
  ggtitle("QC: ribosomal gene percentage")

# patchwork的/表示上下排列。

p_qc_violin <- p_vln_feature /
  p_vln_count /
  p_vln_mt /
  p_vln_ribo

ggsave(
  filename = file.path(
    figure_dir,
    "02_QC_violin_by_sample.pdf"
  ),
  plot = p_qc_violin,
  width = 16,
  height = 18
)

ggsave(
  filename = file.path(
    figure_dir,
    "02_QC_violin_by_sample.png"
  ),
  plot = p_qc_violin,
  width = 16,
  height = 18,
  dpi = 300
)

# ============================================================
# J. QC散点图
# ============================================================

# 每个样本最多随机抽取5000个细胞用于QC散点图。
# 如果某个样本少于5000个细胞，则保留全部细胞。

set.seed(1234)

plot_meta <- meta_qc %>%
  mutate(
    cell_id = rownames(meta_qc)
  ) %>%
  group_by(sample_id) %>%
  group_modify(
    ~ slice_sample(
      .x,
      n = min(5000, nrow(.x))
    )
  ) %>%
  ungroup()

# nCount_RNA vs percent.mt：
# 高UMI且线粒体比例异常高的细胞，可能对应受损或濒死细胞。
#
# nCount_RNA vs nFeature_RNA：
# 极端高UMI且极端高feature的尾部，
# 可能提示潜在doublet或其他异常细胞。

# nCount_RNA vs percent.mt
# 点的颜色表示每个细胞检测到的基因数nFeature_RNA。

p_scatter_mt <- ggplot(
  plot_meta,
  aes(
    x = nCount_RNA,
    y = percent.mt,
    color = nFeature_RNA
  )
) +
  geom_point(
    alpha = 0.35,
    size = 0.35
  ) +
  facet_wrap(
    ~ sample_id,
    scales = "free"
  ) +
  scale_color_viridis_c(
    option = "D",
    name = "nFeature_RNA"
  ) +
  labs(
    title = "QC: UMI counts versus mitochondrial percentage",
    x = "nCount_RNA",
    y = "percent.mt"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 9)
  )

# nCount_RNA vs nFeature_RNA
# 点的颜色表示线粒体比例percent.mt。

p_scatter_feature <- ggplot(
  plot_meta,
  aes(
    x = nCount_RNA,
    y = nFeature_RNA,
    color = percent.mt
  )
) +
  geom_point(
    alpha = 0.35,
    size = 0.35
  ) +
  facet_wrap(
    ~ sample_id,
    scales = "free"
  ) +
  scale_color_viridis_c(
    option = "C",
    name = "percent.mt"
  ) +
  labs(
    title = "QC: UMI counts versus detected genes",
    x = "nCount_RNA",
    y = "nFeature_RNA"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 9)
  )

ggsave(
  filename = file.path(
    figure_dir,
    "02_QC_scatter_nCount_vs_mt.pdf"
  ),
  plot = p_scatter_mt,
  width = 16,
  height = 10
)

ggsave(
  filename = file.path(
    figure_dir,
    "02_QC_scatter_nCount_vs_feature.pdf"
  ),
  plot = p_scatter_feature,
  width = 16,
  height = 10
)

# ============================================================
# K. 保存原始Seurat object
# ============================================================

# 这个对象已经完成：
# 1. 读取每个样本的10x表达矩阵
# 2. 合并12个样本
# 3. 添加sample_id、dataset、disease metadata
# 4. 计算percent.mt和percent.ribo
#
# 这个对象还没有正式过滤低质量细胞。
# 后续如果QC阈值需要调整，可以从该对象重新开始。

saveRDS(
  sc_raw,
  file.path(
    object_dir,
    "02_raw_before_QC_filtering.rds"
  )
)

# ============================================================
# L. 输出session信息
# ============================================================

# 记录R版本和包版本。
# 后续写方法、复现分析和排查环境差异时会用到。

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "02_sessionInfo.txt"
  )
)

# ============================================================
# M. 最终提示
# ============================================================

message("\n============================================================")
message("02_read_and_QC_scRNA.R 运行完成。")
message("已保存原始Seurat object：")
message("results/objects/02_raw_before_QC_filtering.rds")
message("")
message("请重点查看：")
message("1. results/tables/02_sample_QC_summary.csv")
message("2. results/tables/02_QC_threshold_suggestions.csv")
message("3. results/tables/02_target_gene_check.csv")
message("4. results/figures/02_QC_violin_by_sample.pdf")
message("5. results/figures/02_QC_scatter_nCount_vs_mt.pdf")
message("6. results/figures/02_QC_scatter_nCount_vs_feature.pdf")
message("============================================================\n")