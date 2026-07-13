# 17_scTenifold_virtual_KO_OE.R

# 本脚本功能：
# 1. 读取已经完成注释的OSCC单细胞Seurat对象
# 2. 优先筛选核心/扩展恶性focus细胞
# 3. 使用scTenifoldKnk进行目标基因虚拟敲除分析
# 4. 使用scTenifoldNet构建WT网络并模拟目标基因虚拟过表达
# 5. 输出差异调控基因、显著结果表、火山图和Top基因条形图
# 6. 保存实际参与分析的细胞对象、运行参数和sessionInfo
#
# 本项目默认目标基因：
# SASH1、COL1A1、EMP1、MYH11
#
# 本项目默认分析对象：
# Core malignant focus
# Extended malignant focus
#
# 重要说明：
# 1. 本脚本属于基因调控网络虚拟扰动分析，不是拟时序分析
# 2. 虚拟敲除不是CRISPR真实敲除实验
# 3. 虚拟过表达是网络边权扰动模拟，不等同于真实湿实验过表达
# 4. COL1A1、MYH11等基因如果在恶性细胞中表达不足，脚本会自动跳过并记录原因
#
# 通用代码修改位置：
# 1. 换项目路径时：
#    修改project_dir
#
# 2. 换输入Seurat对象时：
#    修改input_seurat_candidates
#
# 3. 换目标基因时：
#    修改target_genes
#
# 4. 换分析细胞群时：
#    修改analysis_mode、malignant_status_keep、
#    malignant_status_column_candidates和tumor_celltype_patterns
#
# 5. 调整最低细胞数时：
#    修改params$minimum_cells
#
# 6. 调整目标基因最低表达细胞数时：
#    修改params$minimum_target_expressing_cells
#
# 7. 调整网络参数时：
#    修改params中的nc_nNet、nc_nCells、n_hvg、
#    overexpression_factor、seed和n_cores
#
# 8. 缺少scTenifold包时：
#    将install_missing_packages改为TRUE
#
# ============================================================
# A. 加载R包
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(Matrix)
})

# ============================================================
# B. 项目路径、输入文件和输出文件夹
# ============================================================

# getwd()应当是R Project根目录。
# 为避免从错误路径运行，本项目默认使用绝对路径。
# 换电脑或移动项目文件夹时，只需要修改这里。

project_dir <- "/Users/yaoshuo/Desktop/HNSCC_SASH1_reproduction"

# -----------------------------
# B1. 输入Seurat对象候选路径
# -----------------------------

# 脚本会按顺序选择第一个真实存在的文件。
#
# 推荐优先使用包含恶性focus标签的11阶段对象。
# 如果11阶段完整RDS不存在，则退回05阶段人工注释对象。
#
# 通用修改位置：
# 换项目或换对象时修改下面的候选路径。

input_seurat_candidates <- c(
  file.path(project_dir, "results/objects/11_global_trajectory_Seurat.rds"),
  file.path(project_dir, "results/objects/05_manual_annotated_plot_ready.rds")
)

# -----------------------------
# B2. 目标基因和分析细胞群
# -----------------------------

# 待进行虚拟敲除和虚拟过表达的基因。
# 可以只保留SASH1，也可以加入其他候选基因。

target_genes <- c("SASH1", "COL1A1", "EMP1", "MYH11")
analysis_mode <- "malignant_focus"

malignant_status_column_candidates <- c(
  "malignant_focus_status",
  "malignant_focus_label",
  "Figure5_forced_status",
  "trajectory_status",
  "focus_label",
  "malignant_state"
)

malignant_status_keep <- c(
  "Core malignant focus",
  "Extended malignant focus",
  "Core_malignant_focus",
  "Extended_malignant_focus"
)

celltype_column_candidates <- c(
  "celltype_plot",
  "celltype_manual",
  "cell_type",
  "celltype",
  "annotation",
  "Cell_Type"
)

tumor_celltype_patterns <- c(
  "tumor", "malignant", "epithelial", "ct-antigen", "cycling"
)

cluster_column_candidates <- c(
  "seurat_clusters",
  "cluster",
  "Cluster",
  "cluster_id",
  "Cluster_ID"
)

malignant_focus_clusters_keep <- c(
  "4",
  "6",
  "11"
)

# ============================================================
# C. 分析参数
# ============================================================

# n_hvg：
# 用于构建网络的高变基因数量
#
# gene_min_cells：
# 一个基因至少需要在多少细胞中检测到
#
# minimum_cells：
# 筛选后的分析细胞总数下限
#
# minimum_target_expressing_cells：
# 目标基因至少需要在多少细胞中表达
#
# nc_nNet：
# scTenifold构建的网络重复数
#
# nc_nCells：
# 每个网络抽样使用的细胞数
#
# overexpression_factor：
# 虚拟过表达时目标基因输出边权的倍数
#
# 通用修改位置：
# 数据量变化时优先调整minimum_cells、nc_nCells和n_hvg。

params <- list(
  n_hvg = 2500L,
  gene_min_cells = 5L,
  minimum_cells = 120L,
  minimum_target_expressing_cells = 10L,
  nc_nNet = 5L,
  nc_nCells = 300L,
  qc_mtThreshold = 0.20,
  qc_minLSize = 500L,
  ma_nDim = 2L,
  pval_threshold = 0.05,
  overexpression_factor = 2,
  seed = 123L,
  n_cores = 1L,
  top_n_plot = 20L
)

# ============================================================
# D. 安装与检查scTenifold相关R包
# ============================================================

# scTenifoldNet包含C++和Fortran相关链接依赖。
# 在Apple Silicon Mac上从GitHub源码安装前，需要先安装：
#
# 1. Apple Command Line Tools
#    Terminal运行：
#    xcode-select --install
#
# 2. 与R 4.5匹配的GNU Fortran
#    推荐安装R官方提供的：
#    gfortran-14.2-universal.pkg
#
# R官方工具页面：
# https://mac.r-project.org/tools/
#
# 第一次安装依赖时设为TRUE。
# 安装成功后改回FALSE，避免每次运行都重复检查和安装。

install_missing_packages <- TRUE

# -----------------------------
# D1. 检查macOS源码编译工具
# -----------------------------

check_macos_build_tools <- function() {
  
  if (Sys.info()[["sysname"]] != "Darwin") {
    return(invisible(TRUE))
  }
  
  xcode_path <- tryCatch(
    system2(
      command = "xcode-select",
      args = "-p",
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) character(0)
  )
  
  if (
    length(xcode_path) == 0 ||
    any(grepl("error", xcode_path, ignore.case = TRUE))
  ) {
    stop(
      paste0(
        "未检测到Apple Command Line Tools。\n",
        "请先退出R，在Terminal运行：\n",
        "xcode-select --install\n",
        "安装完成后重新打开RStudio。"
      )
    )
  }
  
  gfortran_path <- Sys.which("gfortran")
  
  additional_gfortran_candidates <- c(
    "/opt/gfortran/bin/gfortran",
    "/usr/local/gfortran/bin/gfortran",
    "/opt/R/arm64/gfortran/bin/gfortran"
  )
  
  if (!nzchar(gfortran_path)) {
    
    existing_candidates <- additional_gfortran_candidates[
      file.exists(additional_gfortran_candidates)
    ]
    
    if (length(existing_candidates) > 0) {
      gfortran_path <- existing_candidates[1]
      
      Sys.setenv(
        PATH = paste(
          dirname(gfortran_path),
          Sys.getenv("PATH"),
          sep = .Platform$path.sep
        )
      )
    }
  }
  
  if (!nzchar(gfortran_path) || !file.exists(gfortran_path)) {
    stop(
      paste0(
        "未检测到GNU Fortran。\n",
        "scTenifoldNet 1.4需要从源码编译。\n",
        "请安装与R 4.5匹配的gfortran-14.2-universal.pkg：\n",
        "https://mac.r-project.org/tools/gfortran-14.2-universal.pkg\n",
        "安装后彻底关闭并重新打开RStudio。"
      )
    )
  }
  
  gfortran_version <- tryCatch(
    system2(
      command = gfortran_path,
      args = "--version",
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) character(0)
  )
  
  message("检测到Command Line Tools：", xcode_path[1])
  message("检测到gfortran：", gfortran_path)
  
  if (length(gfortran_version) > 0) {
    message("gfortran版本：", gfortran_version[1])
  }
  
  invisible(TRUE)
}

# -----------------------------
# D2. 安装基础安装工具
# -----------------------------

install_scTenifold_dependencies <- function() {
  
  if (!isTRUE(install_missing_packages)) {
    return(invisible(NULL))
  }
  
  options(
    repos = c(CRAN = "https://cloud.r-project.org"),
    timeout = 3600
  )
  
  cran_packages <- c(
    "remotes",
    "Rcpp",
    "RcppArmadillo",
    "RSpectra",
    "RhpcBLASctl",
    "furrr",
    "future",
    "cli",
    "pbapply",
    "enrichR",
    "igraph",
    "reshape2"
  )
  
  missing_cran <- cran_packages[
    !vapply(
      cran_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]
  
  if (length(missing_cran) > 0) {
    message(
      "安装缺失CRAN依赖：",
      paste(missing_cran, collapse = ", ")
    )
    
    install.packages(
      missing_cran,
      dependencies = TRUE,
      type = "binary"
    )
  }
  
  check_macos_build_tools()
  
  # 如果当前安装的是旧版1.3，需要移除后从GitHub安装1.4。
  if (requireNamespace("scTenifoldNet", quietly = TRUE)) {
    
    current_net_version <- packageVersion("scTenifoldNet")
    
    message(
      "当前scTenifoldNet版本：",
      as.character(current_net_version)
    )
    
    if (current_net_version < "1.4") {
      message("移除旧版scTenifoldNet。")
      remove.packages("scTenifoldNet")
    }
  }
  
  if (
    !requireNamespace("scTenifoldNet", quietly = TRUE) ||
    packageVersion("scTenifoldNet") < "1.4"
  ) {
    
    message("从GitHub安装scTenifoldNet 1.4。")
    
    remotes::install_github(
      repo = "cailab-tamu/scTenifoldNet",
      dependencies = FALSE,
      upgrade = "never",
      force = TRUE,
      build_vignettes = FALSE
    )
  }
  
  if (!requireNamespace("scTenifoldNet", quietly = TRUE)) {
    stop("scTenifoldNet安装失败。")
  }
  
  if (packageVersion("scTenifoldNet") < "1.4") {
    stop(
      "scTenifoldNet版本仍低于1.4：",
      as.character(packageVersion("scTenifoldNet"))
    )
  }
  
  if (requireNamespace("scTenifoldKnk", quietly = TRUE)) {
    message(
      "当前scTenifoldKnk版本：",
      as.character(packageVersion("scTenifoldKnk"))
    )
  }
  
  if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
    
    message("从GitHub安装scTenifoldKnk。")
    
    remotes::install_github(
      repo = "cailab-tamu/scTenifoldKnk",
      dependencies = FALSE,
      upgrade = "never",
      force = TRUE,
      build_vignettes = FALSE
    )
  }
  
  invisible(TRUE)
}

install_scTenifold_dependencies()

# -----------------------------
# D3. 最终检查包和版本
# -----------------------------

required_namespaces <- c(
  "scTenifoldKnk",
  "scTenifoldNet",
  "RSpectra"
)

missing_namespaces <- required_namespaces[
  !vapply(
    required_namespaces,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_namespaces) > 0) {
  stop(
    paste0(
      "缺少以下R包：",
      paste(missing_namespaces, collapse = ", "),
      "\n请先完成Command Line Tools和GNU Fortran安装，",
      "然后重新运行本脚本D部分。"
    )
  )
}

if (packageVersion("scTenifoldNet") < "1.4") {
  stop(
    paste0(
      "scTenifoldNet版本过低：",
      as.character(packageVersion("scTenifoldNet")),
      "\n当前scTenifoldKnk需要scTenifoldNet >= 1.4。"
    )
  )
}

message(
  "scTenifoldNet版本：",
  as.character(packageVersion("scTenifoldNet"))
)

message(
  "scTenifoldKnk版本：",
  as.character(packageVersion("scTenifoldKnk"))
)

message(
  "RSpectra版本：",
  as.character(packageVersion("RSpectra"))
)

# ============================================================
# D4. 输出文件夹
# ============================================================

# 图统一保存到results/figures。
# CSV和TXT统一保存到results/tables。
# RDS对象统一保存到results/objects。

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

ko_figure_dir <- file.path(
  figure_dir,
  "17_scTenifold_KO"
)

oe_figure_dir <- file.path(
  figure_dir,
  "17_scTenifold_OE"
)

ko_table_dir <- file.path(
  table_dir,
  "17_scTenifold_KO"
)

oe_table_dir <- file.path(
  table_dir,
  "17_scTenifold_OE"
)

oe_object_dir <- file.path(
  object_dir,
  "17_scTenifold_OE_networks"
)

required_output_dirs <- c(
  table_dir,
  figure_dir,
  object_dir,
  ko_figure_dir,
  oe_figure_dir,
  ko_table_dir,
  oe_table_dir,
  oe_object_dir
)

for (path_use in required_output_dirs) {
  dir.create(
    path_use,
    recursive = TRUE,
    showWarnings = FALSE
  )
}
# ============================================================
# E. 通用辅助函数
# ============================================================

# -----------------------------
# E1. 输入文件与metadata列识别
# -----------------------------

first_existing_file <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    stop("未找到输入Seurat对象：\n", paste(paths, collapse = "\n"))
  }
  existing[1]
}

first_existing_column <- function(metadata, candidates) {
  found <- candidates[candidates %in% colnames(metadata)]
  if (length(found) == 0) NA_character_ else found[1]
}

get_count_matrix <- function(seurat_obj, assay_name) {
  tryCatch(
    GetAssayData(seurat_obj, assay = assay_name, layer = "counts"),
    error = function(e) {
      GetAssayData(seurat_obj, assay = assay_name, slot = "counts")
    }
  )
}

write_utf8_csv <- function(x, file_name) {
  write.csv(
    x,
    file = file_name,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

write_utf8_text <- function(lines, file_name) {
  con <- file(file_name, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con)
}

# -----------------------------
# E2. 筛选实际用于分析的细胞
# -----------------------------

# 优先根据恶性focus状态列筛选。
# 若无法找到恶性focus列或细胞数不足，
# 则根据cell type名称中的tumor、malignant、
# epithelial等关键词进行回退筛选。

select_analysis_object <- function(seurat_obj) {
  
  metadata <- seurat_obj@meta.data
  
  status_column <- first_existing_column(
    metadata,
    malignant_status_column_candidates
  )
  
  cluster_column <- first_existing_column(
    metadata,
    cluster_column_candidates
  )
  
  celltype_column <- first_existing_column(
    metadata,
    celltype_column_candidates
  )
  
  selected_cells <- character(0)
  source_used <- NA_character_
  
  if (
    identical(analysis_mode, "malignant_focus") &&
    !is.na(status_column)
  ) {
    
    status_values <- as.character(
      metadata[[status_column]]
    )
    
    selected_cells <- rownames(metadata)[
      status_values %in% malignant_status_keep
    ]
    
    source_used <- paste0(
      "metadata_status:",
      status_column
    )
  }
  
  if (
    identical(analysis_mode, "malignant_focus") &&
    length(selected_cells) < params$minimum_cells &&
    !is.na(cluster_column)
  ) {
    
    cluster_values <- as.character(
      metadata[[cluster_column]]
    )
    
    selected_cells <- rownames(metadata)[
      cluster_values %in% malignant_focus_clusters_keep
    ]
    
    source_used <- paste0(
      "cluster_fallback:",
      cluster_column,
      ":clusters_",
      paste(
        malignant_focus_clusters_keep,
        collapse = "_"
      )
    )
  }
  
  if (
    identical(analysis_mode, "malignant_focus") &&
    length(selected_cells) < params$minimum_cells &&
    !is.na(celltype_column)
  ) {
    
    celltype_values <- as.character(
      metadata[[celltype_column]]
    )
    
    keep_by_pattern <- Reduce(
      `|`,
      lapply(
        tumor_celltype_patterns,
        function(pattern) {
          grepl(
            pattern,
            celltype_values,
            ignore.case = TRUE
          )
        }
      )
    )
    
    selected_cells <- rownames(metadata)[
      keep_by_pattern
    ]
    
    source_used <- paste0(
      "tumor_celltype_fallback:",
      celltype_column
    )
  }
  
  if (length(selected_cells) < params$minimum_cells) {
    stop(
      "筛选后仅剩",
      length(selected_cells),
      "个细胞，低于minimum_cells=",
      params$minimum_cells
    )
  }
  
  list(
    object = subset(
      seurat_obj,
      cells = selected_cells
    ),
    source = source_used,
    status_column = status_column,
    cluster_column = cluster_column,
    celltype_column = celltype_column,
    selected_clusters = paste(
      malignant_focus_clusters_keep,
      collapse = ","
    )
  )
}

# -----------------------------
# E3. 构建每个目标基因的网络输入矩阵
# -----------------------------

# 保留目标基因和高变基因。
# 自动过滤低检出基因和全零细胞。
# 如果目标基因表达覆盖不足，则跳过该基因。

prepare_network_matrix <- function(count_matrix, hvg_genes, target_gene) {
  selected_genes <- unique(c(target_gene, hvg_genes))
  selected_genes <- selected_genes[selected_genes %in% rownames(count_matrix)]

  if (!target_gene %in% selected_genes) {
    stop("目标基因不在表达矩阵中：", target_gene)
  }

  matrix_use <- count_matrix[selected_genes, , drop = FALSE]
  detected_cells <- Matrix::rowSums(matrix_use > 0)
  keep_genes <- detected_cells >= params$gene_min_cells
  keep_genes[rownames(matrix_use) == target_gene] <- TRUE
  matrix_use <- matrix_use[keep_genes, , drop = FALSE]

  keep_cells <- Matrix::colSums(matrix_use) > 0
  matrix_use <- matrix_use[, keep_cells, drop = FALSE]

  target_expressing_cells <- sum(matrix_use[target_gene, ] > 0)

  if (ncol(matrix_use) < params$minimum_cells) {
    stop("过滤后细胞数不足：", ncol(matrix_use))
  }

  if (target_expressing_cells < params$minimum_target_expressing_cells) {
    stop(
      "目标基因", target_gene, "仅在", target_expressing_cells,
      "个细胞表达，低于阈值", params$minimum_target_expressing_cells
    )
  }

  dense_matrix <- as.matrix(matrix_use)
  storage.mode(dense_matrix) <- "double"

  list(
    matrix = dense_matrix,
    target_expressing_cells = target_expressing_cells,
    target_expression_pct = target_expressing_cells / ncol(dense_matrix) * 100
  )
}

# -----------------------------
# E4. WT网络与扰动网络的manifold alignment
# -----------------------------

# X为WT网络，Y为KO/OE扰动网络。
# 输出两个网络在共同低维空间中的嵌入结果。

safe_manifold_alignment <- function(X, Y, d = 2L) {
  shared_genes <- Reduce(
    intersect,
    list(rownames(X), colnames(X), rownames(Y), colnames(Y))
  )

  if (length(shared_genes) < 3L) {
    stop("共享基因数量不足，无法进行manifold alignment。")
  }

  X <- X[shared_genes, shared_genes, drop = FALSE]
  Y <- Y[shared_genes, shared_genes, drop = FALSE]

  identity_matrix <- diag(length(shared_genes))
  weighted_x <- X + 1
  weighted_y <- Y + 1
  weighted_xy <- 0.9 *
    (sum(weighted_x) + sum(weighted_y)) /
    (2 * sum(identity_matrix)) *
    identity_matrix

  graph_matrix <- rbind(
    cbind(weighted_x, weighted_xy),
    cbind(t(weighted_xy), weighted_y)
  )

  graph_matrix <- -graph_matrix
  diag(graph_matrix) <- 0
  diag(graph_matrix) <- -apply(graph_matrix, 2, sum)

  requested_d <- as.integer(d)
  eig_k <- min(max(requested_d * 4L, 6L), nrow(graph_matrix) - 1L)

  eig_result <- suppressWarnings(
    RSpectra::eigs(graph_matrix, eig_k, which = "SR")
  )

  eig_values <- as.numeric(eig_result$values)
  eig_vectors <- eig_result$vectors
  ordering <- order(eig_values)
  eig_values <- eig_values[ordering]
  eig_vectors <- eig_vectors[, ordering, drop = FALSE]

  keep_dimensions <- which(is.finite(eig_values) & eig_values > 1e-08)
  if (length(keep_dimensions) < requested_d) {
    keep_dimensions <- which(is.finite(eig_values) & eig_values > -1e-08)
  }
  if (length(keep_dimensions) == 0) {
    stop("manifold alignment未获得有效特征向量。")
  }

  aligned_d <- min(requested_d, length(keep_dimensions))
  aligned_network <- eig_vectors[
    ,
    keep_dimensions[seq_len(aligned_d)],
    drop = FALSE
  ]

  colnames(aligned_network) <- paste0("NLMA_", seq_len(aligned_d))
  rownames(aligned_network) <- c(
    paste0("WT_", shared_genes),
    paste0("Perturbed_", shared_genes)
  )

  aligned_network
}

# ============================================================
# F. 虚拟敲除与虚拟过表达核心函数
# ============================================================

# -----------------------------
# F1. 虚拟敲除
# -----------------------------

# 使用scTenifoldKnk直接模拟目标基因KO，
# 返回差异调控结果表。

run_virtual_knockout <- function(
    network_matrix,
    target_gene
) {
  
  cells_for_network <- min(
    as.integer(params$nc_nCells),
    ncol(network_matrix)
  )
  
  set.seed(params$seed)
  
  knockout_function <- getExportedValue(
    "scTenifoldKnk",
    "scTenifoldKnk"
  )
  
  supported_arguments <- names(
    formals(knockout_function)
  )
  
  knockout_arguments <- list(
    countMatrix = as.data.frame(network_matrix),
    gKO = target_gene
  )
  
  optional_arguments <- list(
    qc = TRUE,
    qc_minLibSize = params$qc_minLSize,
    qc_maxMTratio = params$qc_mtThreshold,
    nc_nNet = params$nc_nNet,
    nc_nCells = cells_for_network,
    ma_nDim = params$ma_nDim,
    nCores = params$n_cores
  )
  
  for (argument_name in names(optional_arguments)) {
    
    if (argument_name %in% supported_arguments) {
      knockout_arguments[[argument_name]] <-
        optional_arguments[[argument_name]]
    }
  }
  
  message(
    "scTenifoldKnk实际使用参数：",
    paste(
      names(knockout_arguments),
      collapse = ", "
    )
  )
  
  result <- do.call(
    knockout_function,
    knockout_arguments
  )
  
  if (
    !is.list(result) ||
    is.null(result$diffRegulation)
  ) {
    stop(
      "scTenifoldKnk未返回diffRegulation结果。"
    )
  }
  
  result$diffRegulation
}

# -----------------------------
# F2. 虚拟过表达
# -----------------------------

# 先通过scTenifoldNet构建WT网络，
# 再将目标基因的正向输出边乘以overexpression_factor，
# 最后比较WT和OE网络的差异调控。
#
# 为避免引入强人为假设：
# 如果目标基因在WT网络中没有正向输出边，
# 本脚本不会使用相关性强行补边，而是记录失败原因。

run_dregulation_compatible <- function(
    aligned_network,
    target_gene
) {
  
  dregulation_function <- getFromNamespace(
    "dRegulation",
    "scTenifoldKnk"
  )
  
  supported_arguments <- names(
    formals(dregulation_function)
  )
  
  if (length(supported_arguments) == 0) {
    stop(
      "无法识别dRegulation函数参数。"
    )
  }
  
  dregulation_arguments <- list()
  
  alignment_argument_candidates <- c(
    "MA",
    "manifoldAlignment",
    "alignedNet",
    "X",
    "x"
  )
  
  alignment_argument <- alignment_argument_candidates[
    alignment_argument_candidates %in%
      supported_arguments
  ]
  
  if (length(alignment_argument) > 0) {
    
    dregulation_arguments[[alignment_argument[1]]
    ] <- aligned_network
    
  } else {
    
    dregulation_arguments[[supported_arguments[1]]
    ] <- aligned_network
  }
  
  target_argument_candidates <- c(
    "gKO",
    "gOE",
    "gene",
    "targetGene",
    "target_gene"
  )
  
  target_argument <- target_argument_candidates[
    target_argument_candidates %in%
      supported_arguments
  ]
  
  if (length(target_argument) > 0) {
    
    dregulation_arguments[[target_argument[1]]
    ] <- target_gene
  }
  
  message(
    "dRegulation实际使用参数：",
    paste(
      names(dregulation_arguments),
      collapse = ", "
    )
  )
  
  do.call(
    dregulation_function,
    dregulation_arguments
  )
}

run_dregulation_compatible <- function(
    aligned_network
) {
  
  dregulation_function <- getFromNamespace(
    "dRegulation",
    "scTenifoldKnk"
  )
  
  supported_arguments <- names(
    formals(dregulation_function)
  )
  
  message(
    "dRegulation函数参数：",
    paste(
      supported_arguments,
      collapse = ", "
    )
  )
  
  if (
    length(supported_arguments) == 1L &&
    identical(
      supported_arguments[1],
      "manifoldOutput"
    )
  ) {
    return(
      dregulation_function(
        manifoldOutput = aligned_network
      )
    )
  }
  
  do.call(
    dregulation_function,
    setNames(
      list(aligned_network),
      supported_arguments[1]
    )
  )
}

run_virtual_overexpression <- function(network_matrix, target_gene) {
  strict_direction <- getFromNamespace("strictDirection", "scTenifoldKnk")

  cells_for_network <- min(params$nc_nCells, ncol(network_matrix))
  set.seed(params$seed)

  wt_networks <- scTenifoldNet::makeNetworks(
    X = network_matrix,
    q = 0.9,
    nNet = params$nc_nNet,
    nCells = cells_for_network,
    scaleScores = TRUE,
    symmetric = FALSE,
    nComp = 3L,
    nCores = params$n_cores
  )

  wt_tensor <- scTenifoldNet::tensorDecomposition(
    xList = wt_networks,
    K = 3L,
    maxError = 1e-05,
    maxIter = 1000L,
    nDecimal = 3L
  )

  wt <- strict_direction(wt_tensor$X, lambda = 0)
  wt <- as.matrix(wt)
  diag(wt) <- 0
  wt <- t(wt)

  oe <- wt
  outgoing_edges <- oe[target_gene, , drop = TRUE]
  positive_edges <- is.finite(outgoing_edges) & outgoing_edges > 0

  if (!any(positive_edges)) {
    stop(
      "目标基因", target_gene,
      "在WT网络中没有正向输出边；不进行人为相关性补边，已跳过OE。"
    )
  }

  oe[target_gene, positive_edges] <- outgoing_edges[positive_edges] *
    params$overexpression_factor
  diag(oe) <- 0

  aligned <- safe_manifold_alignment(wt, oe, d = params$ma_nDim)
  differential_regulation <-
    run_dregulation_compatible(
      aligned_network = aligned,
      target_gene = target_gene
    )

  list(
    differential_regulation = differential_regulation,
    wild_type_network = wt,
    overexpression_network = oe
  )
}

# -----------------------------
# F3. 统一KO和OE结果格式
# -----------------------------

standardize_result <- function(result_table, target_gene, perturbation) {
  result_table <- as.data.frame(result_table, stringsAsFactors = FALSE)

  if (!"gene" %in% colnames(result_table)) {
    stop("结果缺少gene列。")
  }

  if (!"p.adj" %in% colnames(result_table)) {
    if ("p.value" %in% colnames(result_table)) {
      result_table$p.adj <- p.adjust(result_table$p.value, method = "BH")
    } else {
      stop("结果缺少p.adj或p.value列。")
    }
  }

  if (!"FC" %in% colnames(result_table)) result_table$FC <- NA_real_
  if (!"Z" %in% colnames(result_table)) result_table$Z <- NA_real_

  result_table %>%
    filter(gene != target_gene) %>%
    mutate(
      target_gene = target_gene,
      perturbation = perturbation,
      significant = p.adj < params$pval_threshold,
      direction = case_when(
        is.finite(FC) & FC > 0 ~ "Activated",
        is.finite(FC) & FC < 0 ~ "Suppressed",
        TRUE ~ "Unresolved"
      ),
      minus_log10_fdr = -log10(pmax(p.adj, .Machine$double.xmin))
    ) %>%
    arrange(p.adj, desc(abs(FC)))
}

# -----------------------------
# F4. 绘制火山图和Top基因条形图
# -----------------------------

save_result_plots <- function(result_table, target_gene, perturbation, prefix) {
  y_cap <- suppressWarnings(
    quantile(
      result_table$minus_log10_fdr[
        is.finite(result_table$minus_log10_fdr)
      ],
      0.995,
      na.rm = TRUE
    )
  )

  if (!is.finite(y_cap)) y_cap <- 5
  y_cap <- min(max(y_cap, -log10(params$pval_threshold) * 1.2, 5), 50)

  plot_table <- result_table %>%
    mutate(minus_log10_fdr_plot = pmin(minus_log10_fdr, y_cap))

  label_table <- plot_table %>%
    filter(significant, is.finite(Z)) %>%
    slice_head(n = 25)

  p_scatter <- ggplot(
    plot_table,
    aes(x = Z, y = minus_log10_fdr_plot, color = direction)
  ) +
    geom_point(alpha = 0.7, size = 1.3) +
    geom_hline(
      yintercept = -log10(params$pval_threshold),
      linetype = "dashed"
    ) +
    scale_color_manual(
      values = c(
        "Activated" = "#D95F02",
        "Suppressed" = "#1B9E77",
        "Unresolved" = "grey70"
      )
    ) +
    coord_cartesian(ylim = c(0, y_cap)) +
    labs(
      title = paste0(target_gene, " ", perturbation, " Differential Regulation"),
      x = "Z-score",
      y = "-log10(FDR)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank()
    )

  if (nrow(label_table) > 0) {
    p_scatter <- p_scatter +
      geom_text_repel(
        data = label_table,
        aes(label = gene),
        size = 3,
        max.overlaps = 50,
        color = "black"
      )
  }

  ggsave(
    paste0(prefix, "_volcano.pdf"),
    p_scatter,
    width = 7,
    height = 5.5
  )

  top_table <- result_table %>%
    filter(is.finite(FC)) %>%
    slice_head(n = params$top_n_plot)

  if (nrow(top_table) > 0) {
    p_bar <- ggplot(
      top_table,
      aes(x = reorder(gene, FC), y = FC, fill = direction)
    ) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(
        values = c(
          "Activated" = "#D95F02",
          "Suppressed" = "#1B9E77",
          "Unresolved" = "grey70"
        )
      ) +
      labs(
        title = paste0(target_gene, " ", perturbation, " Top Genes"),
        x = "Gene",
        y = "FC"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_blank()
      )

    ggsave(
      paste0(prefix, "_top_genes_barplot.pdf"),
      p_bar,
      width = 7,
      height = 6
    )
  }
}

# ============================================================
# G. 读取Seurat对象并准备分析数据
# ============================================================

# -----------------------------
# G1. 读取输入对象
# -----------------------------

input_file <- first_existing_file(input_seurat_candidates)
message("读取Seurat对象：", input_file)
seurat_obj <- readRDS(input_file)

selection <- select_analysis_object(seurat_obj)
analysis_obj <- selection$object
message("细胞筛选来源：", selection$source)
message("用于分析的细胞数：", ncol(analysis_obj))

assay_to_use <- if ("RNA" %in% names(analysis_obj@assays)) {
  "RNA"
} else {
  DefaultAssay(analysis_obj)
}

DefaultAssay(analysis_obj) <- assay_to_use
count_matrix <- get_count_matrix(analysis_obj, assay_to_use)

analysis_obj <- FindVariableFeatures(
  analysis_obj,
  selection.method = "vst",
  nfeatures = params$n_hvg,
  verbose = FALSE
)

hvg_genes <- VariableFeatures(analysis_obj)

runtime_lines <- c(
  paste0("Input_Seurat: ", input_file),
  paste0("Analysis_Mode: ", analysis_mode),
  paste0("Selection_Source: ", selection$source),
  paste0("Cells_Used: ", ncol(analysis_obj)),
  paste0("Assay_Used: ", assay_to_use),
  paste0("Target_Genes: ", paste(target_genes, collapse = ", ")),
  paste0("n_hvg: ", params$n_hvg),
  paste0("nc_nNet: ", params$nc_nNet),
  paste0("nc_nCells: ", params$nc_nCells),
  paste0("overexpression_factor: ", params$overexpression_factor),
  paste0("seed: ", params$seed),
  paste0(
    "scTenifoldKnk_version: ",
    as.character(packageVersion("scTenifoldKnk"))
  ),
  paste0(
    "scTenifoldNet_version: ",
    as.character(packageVersion("scTenifoldNet"))
  )
)

write_utf8_text(
  runtime_lines,
  file.path(table_dir, "17_runtime_config.txt")
)

summary_rows <- list()
combined_results <- list()

# ============================================================
# H. 逐个目标基因运行KO和OE
# ============================================================

# 每个基因分别进行：
# 1. 输入矩阵检查
# 2. 虚拟敲除
# 3. 虚拟过表达
# 4. 保存结果表和图
# 5. 汇总运行状态

for (target_gene in target_genes) {
  message("开始分析：", target_gene)

  if (!target_gene %in% rownames(count_matrix)) {
    summary_rows[[paste0(target_gene, "_missing")]] <- data.frame(
      target_gene = target_gene,
      perturbation = "Not_run",
      status = "Skipped",
      reason = "Gene_not_found",
      stringsAsFactors = FALSE
    )
    next
  }

  prepared <- tryCatch(
    prepare_network_matrix(count_matrix, hvg_genes, target_gene),
    error = function(e) {
      summary_rows[[paste0(target_gene, "_prepare")]] <<- data.frame(
        target_gene = target_gene,
        perturbation = "Not_run",
        status = "Skipped",
        reason = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      NULL
    }
  )

  if (is.null(prepared)) next

  for (perturbation in c("KO", "OE")) {
    
    target_figure_dir <- if (perturbation == "KO") {
      file.path(
        ko_figure_dir,
        paste0("Gene_", target_gene)
      )
    } else {
      file.path(
        oe_figure_dir,
        paste0("Gene_", target_gene)
      )
    }
    
    target_table_dir <- if (perturbation == "KO") {
      file.path(
        ko_table_dir,
        paste0("Gene_", target_gene)
      )
    } else {
      file.path(
        oe_table_dir,
        paste0("Gene_", target_gene)
      )
    }
    
    target_object_dir <- if (perturbation == "OE") {
      file.path(
        oe_object_dir,
        paste0("Gene_", target_gene)
      )
    } else {
      NULL
    }
    
    dir.create(
      target_figure_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    dir.create(
      target_table_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    if (!is.null(target_object_dir)) {
      dir.create(
        target_object_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )
    }

    result <- tryCatch(
      {
        raw_result <- if (perturbation == "KO") {
          run_virtual_knockout(prepared$matrix, target_gene)
        } else {
          run_virtual_overexpression(prepared$matrix, target_gene)
        }

        diff_table <- if (perturbation == "KO") {
          raw_result
        } else {
          raw_result$differential_regulation
        }

        standardized <- standardize_result(
          diff_table,
          target_gene,
          perturbation
        )

        write_utf8_csv(
          standardized,
          file.path(
            target_table_dir,
            paste0(
              target_gene,
              "_",
              perturbation,
              "_all_results.csv"
            )
          )
        )

        write_utf8_csv(
          standardized %>%
            filter(significant),
          file.path(
            target_table_dir,
            paste0(
              target_gene,
              "_",
              perturbation,
              "_significant_results.csv"
            )
          )
        )

        save_result_plots(
          result_table = standardized,
          target_gene = target_gene,
          perturbation = perturbation,
          prefix = file.path(
            target_figure_dir,
            paste0(
              target_gene,
              "_",
              perturbation
            )
          )
        )

        if (perturbation == "OE") {
          saveRDS(
            list(
              wild_type_network =
                raw_result$wild_type_network,
              overexpression_network =
                raw_result$overexpression_network
            ),
            file.path(
              target_object_dir,
              paste0(
                target_gene,
                "_OE_networks.rds"
              )
            )
          )
        }

        summary_rows[[paste0(target_gene, "_", perturbation)]] <- data.frame(
          target_gene = target_gene,
          perturbation = perturbation,
          status = "Success",
          reason = NA_character_,
          cells_used = ncol(prepared$matrix),
          genes_used = nrow(prepared$matrix),
          target_expressing_cells = prepared$target_expressing_cells,
          target_expression_pct = prepared$target_expression_pct,
          significant_genes = sum(standardized$significant, na.rm = TRUE),
          stringsAsFactors = FALSE
        )

        combined_results[[paste0(target_gene, "_", perturbation)]] <- standardized
        standardized
      },
      error = function(e) {
        write_utf8_text(
          conditionMessage(e),
          file.path(
            target_table_dir,
            paste0(
              target_gene,
              "_",
              perturbation,
              "_error.txt"
            )
          )
        )

        summary_rows[[paste0(target_gene, "_", perturbation)]] <<- data.frame(
          target_gene = target_gene,
          perturbation = perturbation,
          status = "Failed",
          reason = conditionMessage(e),
          cells_used = ncol(prepared$matrix),
          genes_used = nrow(prepared$matrix),
          target_expressing_cells = prepared$target_expressing_cells,
          target_expression_pct = prepared$target_expression_pct,
          significant_genes = NA_integer_,
          stringsAsFactors = FALSE
        )
        NULL
      }
    )
  }
}

# ============================================================
# I. 汇总输出
# ============================================================

summary_table <- bind_rows(summary_rows)
write_utf8_csv(
  summary_table,
  file.path(table_dir, "17_scTenifold_KO_OE_summary.csv")
)

if (length(combined_results) > 0) {
  write_utf8_csv(
    bind_rows(combined_results, .id = "analysis_id"),
    file.path(table_dir, "17_scTenifold_KO_OE_all_results.csv")
  )
}

saveRDS(
  analysis_obj,
  file.path(object_dir, "17_scTenifold_analysis_input_Seurat.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(table_dir, "17_sessionInfo.txt")
)

output_check <- data.frame(
  output = c(
    "17_scTenifold_KO_OE_summary.csv",
    "17_scTenifold_KO_OE_all_results.csv",
    "17_runtime_config.txt",
    "17_sessionInfo.txt",
    "17_scTenifold_analysis_input_Seurat.rds",
    "17_scTenifold_KO_figure_directory",
    "17_scTenifold_OE_figure_directory",
    "17_scTenifold_KO_table_directory",
    "17_scTenifold_OE_table_directory",
    "17_scTenifold_OE_network_directory"
  ),
  exists = c(
    file.exists(
      file.path(
        table_dir,
        "17_scTenifold_KO_OE_summary.csv"
      )
    ),
    file.exists(
      file.path(
        table_dir,
        "17_scTenifold_KO_OE_all_results.csv"
      )
    ),
    file.exists(
      file.path(
        table_dir,
        "17_runtime_config.txt"
      )
    ),
    file.exists(
      file.path(
        table_dir,
        "17_sessionInfo.txt"
      )
    ),
    file.exists(
      file.path(
        object_dir,
        "17_scTenifold_analysis_input_Seurat.rds"
      )
    ),
    dir.exists(ko_figure_dir),
    dir.exists(oe_figure_dir),
    dir.exists(ko_table_dir),
    dir.exists(oe_table_dir),
    dir.exists(oe_object_dir)
  ),
  stringsAsFactors = FALSE
)

write_utf8_csv(
  output_check,
  file.path(
    table_dir,
    "17_output_file_check.csv"
  )
)

# ============================================================
# J. 最终提示
# ============================================================

message("\n============================================================")
message("17_scTenifold_virtual_KO_OE.R运行完成。")
message("请重点检查以下文件和目录：")
message("1. results/tables/17_scTenifold_KO_OE_summary.csv")
message("2. results/tables/17_scTenifold_KO_OE_all_results.csv")
message("3. results/tables/17_runtime_config.txt")
message("4. results/tables/17_output_file_check.csv")
message("5. results/figures/17_scTenifold_KO/")
message("6. results/figures/17_scTenifold_OE/")
message("7. results/tables/17_scTenifold_KO/")
message("8. results/tables/17_scTenifold_OE/")
message("9. results/objects/17_scTenifold_OE_networks/")
message("10. results/objects/17_scTenifold_analysis_input_Seurat.rds")
message("============================================================\n")