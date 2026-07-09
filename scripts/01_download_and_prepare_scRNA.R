# 01_download_and_prepare_scRNA.R

# 本脚本功能：
# 1. 建立项目文件夹
# 2. 安装后续分析需要的核心R包
# 3. 下载原始单细胞数据
# 4. 解压GEO RAW.tar文件
# 5. 将散装的GEO 10x文件整理为标准Read10X()文件夹
# 6. 输出每个样本的文件完整性检查表

# 本项目专用数据：
# GSE215403
# 12个OSCC单细胞样本：
# OSCC, scB1, scB2, scB5, scB7, scB8,
# scB9, scB10, scB12, scB13, scB14, scB15
#
# 通用代码修改位置：
# 1. 换数据集时：
#    修改gse_id、raw_tar_url、expected_samples
#
# 2. 换GEO文件命名规则时：
#    修改filename_pattern


# ============================================================
# A. 项目路径与文件夹
# ============================================================

# getwd()应当是你的R Project根目录，可以换成本项目文件夹。
# 我现在的项目路径应为：
# /Users/yaoshuo/Desktop/HNSCC_SASH1_reproduction

project_dir <- getwd()

data_dir <- file.path(project_dir, "data")

raw_dir <- file.path(
  data_dir,
  "raw",
  "scRNA_GSE215403"
)

processed_dir <- file.path(
  data_dir,
  "processed",
  "scRNA_GSE215403"
)

table_dir <- file.path(
  project_dir,
  "results",
  "tables"
)

object_dir <- file.path(
  project_dir,
  "results",
  "objects"
)

figure_dir <- file.path(
  project_dir,
  "results",
  "figures"
)

config_dir <- file.path(
  project_dir,
  "config"
)

script_dir <- file.path(
  project_dir,
  "scripts"
)

# 自动建立需要的文件夹。
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# B. 安装与检查R包
# ============================================================

# 设置CRAN下载镜像。
# 只影响install.packages()从哪里下载，不影响分析结果。
options(repos = c(CRAN = "https://cloud.r-project.org"))

# 设置下载超时。
# GEO大文件或安装大型R包时，避免默认超时过短。
options(timeout = 3600)

# -----------------------------
# B1. 安装BiocManager
# -----------------------------

# BiocManager用于后续安装Bioconductor包。
# 例如：
# scDblFinder、clusterProfiler、org.Hs.eg.db、
# SingleCellExperiment、monocle3等。
#
# 目前这一步不强制安装所有Bioconductor包，
# 因为它们依赖多、体积大，等后续真正用到时再装。

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# -----------------------------
# B2. 安装核心CRAN包
# -----------------------------

# 这些包后续一定会用到，因此一次装好。
#
# Seurat：
# 单细胞、空间转录组主要分析包
#
# harmony：
# 多样本整合与批次校正
#
# dplyr/tidyr/stringr：
# 数据整理
#
# fs：
# 文件夹与文件操作
#
# ggplot2/patchwork/cowplot：
# 绘图与拼图
#
# Matrix：
# 单细胞稀疏矩阵
#
# data.table：
# 大型文本表读取
#
# future：
# 后续可用于并行计算
#
# R.utils：
# 处理压缩文件时备用

cran_packages <- c(
  "Seurat",
  "SeuratObject",
  "harmony",
  "Matrix",
  "data.table",
  "dplyr",
  "tidyr",
  "stringr",
  "fs",
  "ggplot2",
  "patchwork",
  "cowplot",
  "future",
  "R.utils",
  "curl"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("正在安装 R 包：", pkg)
    install.packages(pkg)
  } else {
    message("已安装，跳过：", pkg)
  }
}


# -----------------------------
# B3. 检查关键包
# -----------------------------

package_check <- data.frame(
  package = cran_packages,
  installed = sapply(
    cran_packages,
    requireNamespace,
    quietly = TRUE
  ),
  stringsAsFactors = FALSE
)

print(package_check)

write.csv(
  package_check,
  file.path(
    table_dir,
    "01_core_package_check.csv"
  ),
  row.names = FALSE
)

# 加载本脚本后面实际会使用的包。
library(dplyr)
library(tidyr)
library(stringr)

# ============================================================
# C. 下载GSE215403原始数据
# ============================================================

# 本课题专用GEO编号。
# 换数据集时改这里。

gse_id <- "GSE215403"

# GEO补充文件下载地址。
# GSE215403_RAW.tar大小约268.5MB。
#
# 通用替换提示：
# 不同GSE的路径层级通常是：
# GSE215nnn/GSE215403
#
# 例如GSE123456通常可能位于：
# GSE123nnn/GSE123456

raw_tar_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE215nnn/",
  gse_id,
  "/suppl/",
  gse_id,
  "_RAW.tar"
)

raw_tar_file <- file.path(
  raw_dir,
  paste0(gse_id, "_RAW.tar")
)

# -----------------------------
# C1. 下载逻辑
# -----------------------------

# 注意：
# NCBI有时网络较慢。
# 如果R下载反复中断，优先用浏览器或Terminal curl下载，
# 然后放进raw_dir即可。

if (!file.exists(raw_tar_file)) {
  
  message("开始下载：", gse_id)
  message("下载地址：", raw_tar_url)
  
  download.file(
    url = raw_tar_url,
    destfile = raw_tar_file,
    method = "libcurl",
    mode = "wb",
    quiet = FALSE
  )
  
} else {
  
  message("检测到文件已存在，跳过下载：")
  message(raw_tar_file)
}

# -----------------------------
# C2. 检查下载文件大小
# -----------------------------

# 正常约为268.5MB。
# 若只有几KB、几MB，通常说明下载不完整。

tar_size_mb <- file.info(raw_tar_file)$size / 1024^2

message("RAW.tar 文件大小：", round(tar_size_mb, 2), " MB")

if (is.na(tar_size_mb) || tar_size_mb < 200) {
  
  warning(
    paste0(
      "RAW.tar 文件可能未完整下载。\n",
      "当前大小：", round(tar_size_mb, 2), " MB\n",
      "预期大小约：268.5 MB"
    )
  )
}

# ============================================================
# D. 解压GSE215403_RAW.tar
# ============================================================

# GEO解压后得到36个文件：
# 12个样本×3个文件
#
# 每个样本有：
# - barcodes.tsv.gz
# - features.tsv.gz
# - matrix.mtx.gz

raw_expression_files <- list.files(
  processed_dir,
  pattern = "\\.(mtx|tsv)\\.gz$",
  recursive = FALSE,
  full.names = TRUE
)

if (length(raw_expression_files) == 0) {
  message("开始解压 RAW.tar：")
  utils::untar(
    tarfile = raw_tar_file,
    exdir = processed_dir
  )
  message("解压完成。")
} else {
  message("检测到已解压的表达矩阵文件，跳过解压。")
}

# -----------------------------
# D1. 输出原始文件清单
# -----------------------------

all_raw_files <- list.files(
  processed_dir,
  recursive = FALSE,
  full.names = TRUE
)

file_table <- data.frame(
  full_path = all_raw_files,
  file_name = basename(all_raw_files),
  stringsAsFactors = FALSE
)

write.csv(
  file_table,
  file.path(
    table_dir,
    "01_raw_file_list.csv"
  ),
  row.names = FALSE
)

message("解压后的原始文件数：", nrow(file_table))

# ============================================================
# E. 识别样本与建立标准10x文件夹
# ============================================================

# 本课题的12个原始样本ID。
#
# 直接保留GEO原始sampleID，
# 后面作为Seurat metadata中的sample_id使用。

expected_samples <- c(
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

# -----------------------------
# E1. 解析GEO原始文件名
# -----------------------------

# 本数据集文件名规则：
#
# GSM6634869_OSCC_barcodes.tsv.gz
# GSM6634869_OSCC_features.tsv.gz
# GSM6634869_OSCC_matrix.mtx.gz
#
# 正则表达式含义：
#
# ^GSM\\d+_
# 文件以GSM编号_开头
#
# (.+?)
# 提取中间的样本名，例如OSCC或scB1
#
# (barcodes...|features...|matrix...)
# 提取文件类型
#
# 通用修改位置：
# 换GEO数据集时，最有可能需要改这里的filename_pattern。

filename_pattern <- paste0(
  "^GSM\\d+_(.+?)_",
  "(barcodes\\.tsv\\.gz|features\\.tsv\\.gz|matrix\\.mtx\\.gz)$"
)

parsed <- stringr::str_match(
  file_table$file_name,
  filename_pattern
)

manifest <- file_table %>%
  dplyr::mutate(
    sample_id = parsed[, 2],
    file_type = parsed[, 3]
  ) %>%
  dplyr::filter(
    !is.na(sample_id),
    !is.na(file_type)
  ) %>%
  dplyr::arrange(sample_id, file_type)

# -----------------------------
# E2. 检查每个样本是否都有三类文件
# -----------------------------

sample_check <- manifest %>%
  dplyr::count(sample_id, file_type) %>%
  tidyr::pivot_wider(
    names_from = file_type,
    values_from = n,
    values_fill = 0
  )

print(sample_check)

write.csv(
  sample_check,
  file.path(
    table_dir,
    "01_raw_sample_file_check.csv"
  ),
  row.names = FALSE
)

# -----------------------------
# E3. 检查检测到的样本是否与预期一致
# -----------------------------

detected_samples <- sort(unique(manifest$sample_id))

message("检测到的样本：")
print(detected_samples)

if (!setequal(expected_samples, detected_samples)) {
  
  stop(
    paste0(
      "检测到的样本ID与预期不一致。\n",
      "请查看：results/tables/01_raw_sample_file_check.csv"
    )
  )
}

# -----------------------------
# E4. 建立标准10x文件夹
# -----------------------------

# 最终结构：
#
# data/processed/scRNA_GSE215403/10x_by_sample/
# ├── OSCC/
# │   ├── matrix.mtx.gz
# │   ├── barcodes.tsv.gz
# │   └── features.tsv.gz
# ├── scB1/
# │   ├── matrix.mtx.gz
# │   ├── barcodes.tsv.gz
# │   └── features.tsv.gz
# ...
#
# 后续可以直接：
# Seurat::Read10X(data.dir = ".../10x_by_sample/OSCC")

tenx_root <- file.path(
  processed_dir,
  "10x_by_sample"
)

# 为保证可以重复运行，每次重新建立标准10x文件夹。

if (dir.exists(tenx_root)) {
  unlink(
    tenx_root,
    recursive = TRUE,
    force = TRUE
  )
}

dir.create(
  tenx_root,
  recursive = TRUE,
  showWarnings = FALSE
)

# -----------------------------
# E5. 复制每个样本的3个文件
# -----------------------------

for (sid in expected_samples) {
  
  sample_files <- manifest %>%
    dplyr::filter(.data$sample_id == sid)
  
  matrix_file <- sample_files %>%
    dplyr::filter(.data$file_type == "matrix.mtx.gz") %>%
    dplyr::pull(full_path)
  
  barcode_file <- sample_files %>%
    dplyr::filter(.data$file_type == "barcodes.tsv.gz") %>%
    dplyr::pull(full_path)
  
  feature_file <- sample_files %>%
    dplyr::filter(.data$file_type == "features.tsv.gz") %>%
    dplyr::pull(full_path)
  
  # 每种文件必须恰好找到一个。
  
  if (
    length(matrix_file) != 1 ||
    length(barcode_file) != 1 ||
    length(feature_file) != 1
  ) {
    
    stop(
      paste0(
        "样本 ", sid, " 文件不完整或匹配重复。\n",
        "matrix = ", length(matrix_file), "\n",
        "barcodes = ", length(barcode_file), "\n",
        "features = ", length(feature_file)
      )
    )
  }
  
  sample_outdir <- file.path(
    tenx_root,
    sid
  )
  
  dir.create(
    sample_outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # 统一标准10x文件名。
  # Read10X()才能直接识别。
  
  file.copy(
    matrix_file,
    file.path(sample_outdir, "matrix.mtx.gz"),
    overwrite = TRUE
  )
  
  file.copy(
    barcode_file,
    file.path(sample_outdir, "barcodes.tsv.gz"),
    overwrite = TRUE
  )
  
  file.copy(
    feature_file,
    file.path(sample_outdir, "features.tsv.gz"),
    overwrite = TRUE
  )
  
  message("已整理样本：", sid)
}

# -----------------------------
# E6. 最终检查12个样本文件夹
# -----------------------------

final_check <- lapply(expected_samples, function(sid) {
  
  sample_dir <- file.path(
    tenx_root,
    sid
  )
  
  files <- list.files(
    sample_dir,
    full.names = FALSE
  )
  
  data.frame(
    sample_id = sid,
    has_matrix = "matrix.mtx.gz" %in% files,
    has_barcodes = "barcodes.tsv.gz" %in% files,
    has_features = "features.tsv.gz" %in% files,
    stringsAsFactors = FALSE
  )
}) %>%
  dplyr::bind_rows()

print(final_check)

write.csv(
  final_check,
  file.path(
    table_dir,
    "01_10x_folder_check.csv"
  ),
  row.names = FALSE
)

# -----------------------------
# E7. 保存基础样本metadata
# -----------------------------

# 目前只记录：
# sample_id、数据集、疾病
#
# 后面拿到patient clinical metadata后，
# 可以在这个csv里补充：
# stage、sex、age、site、treatment等列。

sample_metadata <- data.frame(
  sample_id = expected_samples,
  dataset = "GSE215403",
  disease = "OSCC",
  stringsAsFactors = FALSE
)

write.csv(
  sample_metadata,
  file.path(
    config_dir,
    "GSE215403_sample_metadata.csv"
  ),
  row.names = FALSE
)

# ============================================================
# F. 最终提示
# ============================================================

message("\n============================================================")
message("01_download_and_prepare_scRNA.R 运行完成。")
message("下一步可以开始读取 12 个样本并做 QC。")
message("请重点确认以下文件：")
message("1. results/tables/01_10x_folder_check.csv")
message("2. results/tables/01_raw_sample_file_check.csv")
message("3. config/GSE215403_sample_metadata.csv")
message("============================================================\n")