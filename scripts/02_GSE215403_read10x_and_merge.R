# ============================================================
# 02_read_and_QC_scRNA.R
#
# 项目：
# HNSCC / OSCC SASH1 单细胞 + 空间转录组复现与通用流程整理
#
# 本脚本功能：
# 1. 读取 01 中整理好的 12 个标准 10x 单细胞文件夹
# 2. 为每个细胞添加 sample_id 等 metadata
# 3. 合并为一个原始 Seurat object
# 4. 计算常用 QC 指标：
#    - nFeature_RNA：每个细胞检测到的基因数
#    - nCount_RNA：每个细胞的总 UMI 数
#    - percent.mt：线粒体基因比例
#    - percent.ribo：核糖体基因比例
# 5. 输出每个样本的 QC 汇总统计表
# 6. 画 QC 小提琴图和散点图
# 7. 基于每个样本自身分布，生成“建议 QC 阈值”
#
# - 正式过滤低质量细胞
# - 双细胞检测
# - NormalizeData / SCTransform
# - Harmony 整合
# - PCA / UMAP / clustering
#
# 原因：
# 先检查每个样本的真实 QC 分布，
# 再决定第三步如何正式过滤。
#
# ------------------------------------------------------------
# 通用代码修改位置：
#
# 1. tenx_root：
#    换成你自己的 10x 标准文件夹总目录
#
# 2. sample_ids：
#    换成你的样本名，或者直接从文件夹自动读取
#
# 3. mt_pattern：
#    人类通常为 "^MT-"
#    小鼠通常为 "^mt-"
#
# 4. disease / dataset：
#    换成自己的疾病和数据集名称
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
# A. 安装并加载本步骤需要的包
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 3600)

# 这一步真正需要的包：
# Seurat：读取 10x、建立对象、QC 图
# dplyr：整理统计表
# ggplot2：画散点图
# patchwork：拼多个 Seurat 图
# fs：文件夹操作

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
# B. 项目路径
# ============================================================

# getwd() 应该是 R Project 根目录。
# 例如：
# /Users/yaoshuo/Desktop/HNSCC_SASH1_reproduction

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

# 本课题专用：
# 这 12 个名字与 01 中整理出的文件夹名称必须完全一致。
#
# 通用项目：
# 可以直接替换为自己的 sample_id。
#
# 也可以写成：
# sample_ids <- basename(list.dirs(tenx_root, recursive = FALSE))
#
# 但这里为了输出顺序稳定，先手动指定。

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
# 检查每个样本文件夹是否存在
# -----------------------------

missing_samples <- sample_ids[
  !dir.exists(file.path(tenx_root, sample_ids))
]

if (length(missing_samples) > 0) {
  
  stop(
    paste0(
      "以下样本文件夹不存在：\n",
      paste(missing_samples, collapse = ", "),
      "\n请先运行 01_download_and_prepare_scRNA.R"
    )
  )
}

# ============================================================
# D. 逐个读取 10x 数据并建立 Seurat object
# ============================================================

# CreateSeuratObject 初始门槛故意设得宽松：
#
# min.cells = 3：
# 一个基因至少出现在 3 个细胞中才保留。
# 主要去掉极端罕见、无信息量的基因。
#
# min.features = 100：
# 一个细胞至少检测到 100 个基因才进入原始对象。
#
# 正式过滤会在下一步完成。

seurat_list <- list()

for (sid in sample_ids) {
  
  message("\n========================================")
  message("正在读取样本：", sid)
  message("========================================")
  
  sample_dir <- file.path(
    tenx_root,
    sid
  )
  
  # Read10X() 读取：
  # matrix.mtx.gz
  # barcodes.tsv.gz
  # features.tsv.gz
  #
  # gene.column = 2：
  # 10x features.tsv 通常第二列是 gene symbol。
  # 比如 SASH1、COL1A1、EMP1 等。
  #
  # 而是 Ensembl ID，就需要改 gene.column，
  # 或者后续做 ID 转换。
  
  counts <- Read10X(
    data.dir = sample_dir,
    gene.column = 2
  )
  
  # 创建单个样本的 Seurat object
  obj <- CreateSeuratObject(
    counts = counts,
    project = dataset_name,
    min.cells = 3,
    min.features = 100
  )
  
  # 给每个细胞增加 metadata。
  # 后续 Harmony、QC、分样本作图都依赖 sample_id。
  
  obj$sample_id <- sid
  obj$dataset <- dataset_name
  obj$disease <- disease_name
  
  # 暂存至 list
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

# add.cell.ids 会给细胞 barcode 加上样本前缀。
#
# 例如原始 barcode：
# AAACCTGAGAGTAATC-1
#
# 合并后会变为：
# OSCC_AAACCTGAGAGTAATC-1
#
# 这一步非常重要。
# 不同 10x 样本可能出现相同 barcode，
# 若没有前缀，合并后可能造成细胞名称冲突。

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
# F. 计算 QC 指标
# ============================================================

# -----------------------------
# F1. 线粒体基因比例
# -----------------------------
#
# 人类基因符号常用 MT- 开头：
# MT-CO1、MT-ND1、MT-ATP6 等。
#
# 通用修改：
# 人类：pattern = "^MT-"
# 小鼠：pattern = "^mt-"
#
# 线粒体比例过高通常意味着：
# - 细胞受损
# - RNA 降解
# - 濒死细胞
#
# 但阈值不应一刀切。
# 后面会按样本真实分布判断。

sc_raw[["percent.mt"]] <- PercentageFeatureSet(
  sc_raw,
  pattern = "^MT-"
)

# -----------------------------
# F2. 核糖体基因比例
# -----------------------------
#
# 人类核糖体蛋白基因通常以：
# RPL、RPS 开头。
#
# 但可辅助判断细胞状态和某些异常样本。

sc_raw[["percent.ribo"]] <- PercentageFeatureSet(
  sc_raw,
  pattern = "^RP[SL]"
)

# -----------------------------
# F3. 检查四个论文目标基因
# -----------------------------
#
# 这一步确认 gene symbol 是否正确读取。
# 若这里显示 FALSE，后续 FeaturePlot 会找不到基因。
#
# 通用项目：
# 修改 target_genes 即可。

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
      "可能原因：gene symbol 不一致、Ensembl ID、旧基因名，",
      "或 Read10X(gene.column = 2) 不适合该数据。"
    )
  )
}

# ============================================================
# G. 输出样本级 QC 统计表
# ============================================================

# 这里会输出每个样本：
#
# raw_cells：
# 初步读取后保留的细胞数
#
# median / mean：
# 中位数、均值
#
# q01 / q99：
# 1% 和 99% 分位数
#
# 这些可以帮助判断：
# - 是否有测序深度特别低的样本
# - 是否有 mt 比例整体异常的样本
# - 是否存在潜在 doublet 高 UMI / 高 feature 尾部

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
# H. 生成每个样本的建议 QC 阈值
# ============================================================

# 只是基于每个样本自身的分布，生成一个“建议范围”。
#
# 规则：
# - nFeature_RNA：
#   下限取 max(200, 1%分位数)
#   上限取 99%分位数
#
# - nCount_RNA：
#   下限取 1%分位数
#   上限取 99%分位数
#
# - percent.mt：
#   上限取 max(20, 99%分位数)
#
# nFeature_RNA > 200 & nFeature_RNA < 5000？
#
# 因为肿瘤组织、不同平台、不同样本的测序深度差异很大。
# 某些高质量肿瘤细胞可能天然有 6000-10000 个基因，
#
# 下一步我们会结合图和这张建议表，
# 再决定最终正式筛选规则。

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
# I. QC 小提琴图
# ============================================================

# pt.size = 0：
# 不画单个细胞的小点。
# 因为细胞太多，画出来会很密、也会很慢。
#
# group.by = "sample_id"：
# 每个样本单独看。
# 不建议把所有样本混在一起看 QC，
# 因为不同病人的测序深度可能明显不同。

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

# 拼成一张总图。
# patchwork 的 / 表示上下排列。

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
# J. QC 散点图
# ============================================================
# 每个 sample 最多随机抽取 5000 个细胞用于 QC 散点图。
# 若该样本少于 5000 个细胞，则保留全部细胞。

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
# 高 UMI 但 mt 异常高，可能是受损细胞或异常细胞。
#
# nCount_RNA vs nFeature_RNA：
# 极端高 UMI + 极端高 feature 的尾部，
# 有时提示潜在 doublet。

# nCount_RNA vs percent.mt
# 点的颜色表示每个细胞检测到的基因数 nFeature_RNA

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
# 点的颜色表示线粒体比例 percent.mt

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
# K. 保存原始 Seurat object
# ============================================================

# 这个对象是：
# 已读取 + 已合并 + 已计算 QC 指标
#
# 但还没有真正过滤细胞。
#
# 后面即使 QC 参数改了，也可以从这个对象重新开始，
# 不必重新读取所有 10x 文件。

saveRDS(
  sc_raw,
  file.path(
    object_dir,
    "02_GSE215403_raw_before_QC_filtering.rds"
  )
)

# ============================================================
# L. 输出 session 信息
# ============================================================

# 记录 R 版本和包版本。
# 后面写方法、复现、排错时很有用。

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
message("已保存原始 Seurat object：")
message("results/objects/02_GSE215403_raw_before_QC_filtering.rds")
message("")
message("请重点查看：")
message("1. results/tables/02_sample_QC_summary.csv")
message("2. results/tables/02_QC_threshold_suggestions.csv")
message("3. results/tables/02_target_gene_check.csv")
message("4. results/figures/02_QC_violin_by_sample.pdf")
message("5. results/figures/02_QC_scatter_nCount_vs_mt.pdf")
message("6. results/figures/02_QC_scatter_nCount_vs_feature.pdf")
message("============================================================\n")