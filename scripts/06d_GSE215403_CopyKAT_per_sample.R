# ============================================================
# 06d_CopyKAT_malignant_call.R
#
# 目标：
# 1. 按 sample_id 单独运行 CopyKAT
# 2. 使用同一样本内的非 tumor epithelial candidate 细胞
#    作为已知 normal reference
# 3. 获得 CopyKAT 的 diploid / aneuploid prediction
# 4. 回填到原始 Seurat metadata
# 5. 重点汇总 cluster 2 / 3 / 4 / 6 / 11 的预测结果
#
# 说明：
# CopyKAT prediction 是基于 scRNA-seq 推断的 CNV/aneuploidy 结果，
#
# 后续 malignant cell 定义将结合：
# - tumor epithelial marker
# - CopyKAT aneuploid prediction
# - sample-level consistency
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
# A. R library 与包
# ============================================================

options(timeout = 3600)

user_r_library <- Sys.getenv("R_LIBS_USER")

if (nzchar(user_r_library)) {
  
  dir.create(
    user_r_library,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  .libPaths(
    c(
      user_r_library,
      .libPaths()
    )
  )
}

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "copykat",
  "dplyr",
  "tidyr",
  "ggplot2",
  "Matrix"
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
      "缺少 R 包：",
      paste(missing_packages, collapse = ", "),
      "\n请先完成 CopyKAT 安装。"
    )
  )
}

library(Seurat)
library(SeuratObject)
library(copykat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)

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

copykat_dir <- file.path(
  project_dir,
  "results",
  "copykat"
)

dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(copykat_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# C. 读取 06b 对象
# ============================================================

input_object_file <- file.path(
  object_dir,
  "06b_GSE215403_malignant_candidate_diagnostic.rds"
)

if (!file.exists(input_object_file)) {
  
  stop(
    paste0(
      "找不到对象：\n",
      input_object_file
    )
  )
}

sc <- readRDS(input_object_file)

DefaultAssay(sc) <- "RNA"

cluster_column <- "cluster_res_0.2"
sample_column <- "sample_id"
status_column <- "malignant_status_diagnostic"

required_metadata <- c(
  cluster_column,
  sample_column,
  status_column,
  "celltype_manual"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(sc@meta.data)
)

if (length(missing_metadata) > 0) {
  
  stop(
    paste0(
      "缺少 metadata：",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

# ============================================================
# D. 参数
# ============================================================

tumor_candidate_clusters <- c(
  "2",
  "3",
  "4",
  "6",
  "11"
)

# 每个 sample 至少需要：
# - 50 个 tumor epithelial candidate
# - 50 个 known normal reference cells
#
# 不满足时跳过，避免 CopyKAT 在极小样本中产生不稳定预测。

min_tumor_cells <- 50
min_normal_cells <- 50

# Eddie 当前申请一个 CPU slot，因此设为 1。
# 以后改成多核 qlogin 后，可相应提高。

copykat_cores <- 1

set.seed(1234)

# ============================================================
# E. 获取 raw count matrix
# ============================================================

raw_counts <- LayerData(
  object = sc,
  assay = "RNA",
  layer = "counts"
)

message("Raw count matrix：")
message("Genes: ", nrow(raw_counts))
message("Cells: ", ncol(raw_counts))

# ============================================================
# F. 标记每个细胞的候选状态
# ============================================================

meta <- sc@meta.data %>%
  mutate(
    cell_barcode = rownames(sc@meta.data),
    cluster = as.character(
      .data[[cluster_column]]
    ),
    is_tumor_candidate = cluster %in% tumor_candidate_clusters,
    is_known_normal = !is_tumor_candidate
  )

sample_ids <- sort(
  unique(
    as.character(meta[[sample_column]])
  )
)

sample_summary <- meta %>%
  group_by(
    .data[[sample_column]]
  ) %>%
  summarise(
    total_cells = n(),
    tumor_candidate_cells = sum(is_tumor_candidate),
    known_normal_cells = sum(is_known_normal),
    .groups = "drop"
  ) %>%
  mutate(
    will_run_copykat =
      tumor_candidate_cells >= min_tumor_cells &
      known_normal_cells >= min_normal_cells
  )

write.csv(
  sample_summary,
  file.path(
    table_dir,
    "06d_CopyKAT_sample_input_summary.csv"
  ),
  row.names = FALSE
)

print(sample_summary)

# ============================================================
# G. 按 sample 独立运行 CopyKAT
# ============================================================
#
# 关键点：
# 1. 每个 sample 独立运行，避免 patient-to-patient 差异混入。
# 3. non-tumor candidate 细胞作为已知 normal reference。
# ============================================================

prediction_list <- list()
run_summary_list <- list()

for (current_sample in sample_ids) {
  
  message("\n============================================================")
  message("开始 CopyKAT：", current_sample)
  message("============================================================")
  
  current_meta <- meta %>%
    filter(
      .data[[sample_column]] == current_sample
    )
  
  current_cells <- current_meta$cell_barcode
  
  current_tumor_cells <- current_meta$cell_barcode[
    current_meta$is_tumor_candidate
  ]
  
  current_normal_cells <- current_meta$cell_barcode[
    current_meta$is_known_normal
  ]
  
  n_tumor <- length(current_tumor_cells)
  n_normal <- length(current_normal_cells)
  n_total <- length(current_cells)
  # ------------------------------------------------------------
# 断点续跑：
# 若该 sample 已存在完整且细胞数匹配的 prediction 文件，
# 则直接读取并跳过 CopyKAT 重算。
# ------------------------------------------------------------

existing_prediction_file <- file.path(
  copykat_dir,
  paste0(
    "06d_CopyKAT_",
    current_sample
  ),
  paste0(
    "06d_CopyKAT_",
    current_sample,
    "_cell_prediction.csv"
  )
)

if (file.exists(existing_prediction_file)) {
  
  existing_prediction <- read.csv(
    existing_prediction_file,
    stringsAsFactors = FALSE
  )
  
  required_prediction_columns <- c(
    "cell_barcode",
    "copykat_prediction",
    "copykat_sample"
  )
  
  prediction_is_complete <- all(
    required_prediction_columns %in%
      colnames(existing_prediction)
  ) &&
    nrow(existing_prediction) == n_total &&
    all(current_cells %in% existing_prediction$cell_barcode)
  
  if (prediction_is_complete) {
    
    prediction_list[[current_sample]] <- existing_prediction
    
    run_summary_list[[current_sample]] <- data.frame(
      sample_id = current_sample,
      total_cells = n_total,
      tumor_candidate_cells = n_tumor,
      known_normal_cells = n_normal,
      CopyKAT_status = "reused_existing_result",
      note = "existing complete prediction reused",
      stringsAsFactors = FALSE
    )
    
    message(
      "检测到已完成 CopyKAT 输出，跳过重算：",
      current_sample
    )
    
    next
  }
}
  current_status <- "not_run"
  current_note <- ""
  
  if (n_tumor < min_tumor_cells) {
    
    current_note <- paste0(
      "tumor candidate cells < ",
      min_tumor_cells
    )
    
  } else if (n_normal < min_normal_cells) {
    
    current_note <- paste0(
      "known normal cells < ",
      min_normal_cells
    )
    
  } else {
    
    current_output_dir <- file.path(
      copykat_dir,
      paste0(
        "06d_CopyKAT_",
        current_sample
      )
    )
    
    dir.create(
      current_output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    old_working_dir <- getwd()
    setwd(current_output_dir)
    
    current_raw_counts <- raw_counts[
      ,
      current_cells,
      drop = FALSE
    ]
    
    # CopyKAT 官方示例使用 raw UMI matrix。
    # 输入需要 dense matrix，因此仅在单个 sample 内转换。
    
    current_raw_matrix <- as.matrix(
      current_raw_counts
    )
    
    rm(current_raw_counts)
    gc()
    
    copykat_result <- tryCatch(
      {
        copykat(
          rawmat = current_raw_matrix,
          id.type = "S",
          ngene.chr = 5,
          win.size = 25,
          KS.cut = 0.1,
          sam.name = current_sample,
          distance = "euclidean",
          norm.cell.names = current_normal_cells,
          output.seg = "FALSE",
          plot.genes = TRUE,
          genome = "hg20",
          n.cores = copykat_cores
        )
      },
      error = function(e) {
        
        message(
          "CopyKAT 失败：",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    rm(current_raw_matrix)
    gc()
    
    setwd(old_working_dir)
    
    if (is.null(copykat_result)) {
      
      current_status <- "failed"
      current_note <- "CopyKAT returned NULL or error"
      
    } else {
      
      current_status <- "completed"
      current_note <- "success"
      
      saveRDS(
        copykat_result,
        file.path(
          current_output_dir,
          paste0(
            "06d_CopyKAT_",
            current_sample,
            "_result.rds"
          )
        )
      )
      
      current_prediction <- as.data.frame(
        copykat_result$prediction,
        stringsAsFactors = FALSE
      )
      
      # CopyKAT 通常输出 cell.names 和 copykat.pred。
      # 此处检查列名，避免静默写错。
      
      if (!all(
        c(
          "cell.names",
          "copykat.pred"
        ) %in% colnames(current_prediction)
      )) {
        
        current_status <- "failed"
        current_note <- paste0(
          "unexpected prediction columns: ",
          paste(
            colnames(current_prediction),
            collapse = ", "
          )
        )
        
      } else {
        
        current_prediction <- current_prediction %>%
          transmute(
            cell_barcode = cell.names,
            copykat_prediction = copykat.pred,
            copykat_sample = current_sample
          )
        
        prediction_list[[current_sample]] <- current_prediction
        
        write.csv(
          current_prediction,
          file.path(
            current_output_dir,
            paste0(
              "06d_CopyKAT_",
              current_sample,
              "_cell_prediction.csv"
            )
          ),
          row.names = FALSE
        )
      }
    }
  }
  
  run_summary_list[[current_sample]] <- data.frame(
    sample_id = current_sample,
    total_cells = n_total,
    tumor_candidate_cells = n_tumor,
    known_normal_cells = n_normal,
    CopyKAT_status = current_status,
    note = current_note,
    stringsAsFactors = FALSE
  )
  
  gc()
}

run_summary <- bind_rows(
  run_summary_list
)

write.csv(
  run_summary,
  file.path(
    table_dir,
    "06d_CopyKAT_run_summary.csv"
  ),
  row.names = FALSE
)

print(run_summary)

# ============================================================
# H. 合并 CopyKAT prediction
# ============================================================

if (length(prediction_list) == 0) {
  
  stop(
    "没有任何 sample 成功完成 CopyKAT。请检查 06d_CopyKAT_run_summary.csv。"
  )
}

copykat_prediction_all <- bind_rows(
  prediction_list
)

write.csv(
  copykat_prediction_all,
  file.path(
    table_dir,
    "06d_CopyKAT_all_cell_predictions.csv"
  ),
  row.names = FALSE
)

prediction_vector <- copykat_prediction_all$copykat_prediction
names(prediction_vector) <- copykat_prediction_all$cell_barcode

sc$copykat_prediction <- prediction_vector[
  colnames(sc)
]

sc$copykat_prediction <- factor(
  sc$copykat_prediction,
  levels = c(
    "aneuploid",
    "diploid",
    "not.defined"
  )
)

sc$copykat_prediction[
  is.na(sc$copykat_prediction)
] <- "not.run"

sc$copykat_prediction <- factor(
  sc$copykat_prediction,
  levels = c(
    "aneuploid",
    "diploid",
    "not.defined",
    "not.run"
  )
)

# ============================================================
# I. 汇总 tumor epithelial candidate 的 CopyKAT 结果
# ============================================================

candidate_prediction_summary <- sc@meta.data %>%
  mutate(
    cell_barcode = rownames(sc@meta.data),
    cluster = as.character(
      .data[[cluster_column]]
    )
  ) %>%
  filter(
    cluster %in% tumor_candidate_clusters
  ) %>%
  count(
    sample_id,
    cluster,
    celltype_manual,
    copykat_prediction,
    name = "cell_number"
  ) %>%
  group_by(
    sample_id,
    cluster
  ) %>%
  mutate(
    percent_within_sample_cluster = round(
      100 * cell_number / sum(cell_number),
      2
    )
  ) %>%
  ungroup() %>%
  arrange(
    sample_id,
    suppressWarnings(as.numeric(cluster)),
    desc(cell_number)
  )

write.csv(
  candidate_prediction_summary,
  file.path(
    table_dir,
    "06d_CopyKAT_tumor_candidate_prediction_by_sample_cluster.csv"
  ),
  row.names = FALSE
)

candidate_prediction_overall <- sc@meta.data %>%
  mutate(
    cluster = as.character(
      .data[[cluster_column]]
    )
  ) %>%
  filter(
    cluster %in% tumor_candidate_clusters
  ) %>%
  count(
    cluster,
    celltype_manual,
    copykat_prediction,
    name = "cell_number"
  ) %>%
  group_by(
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
    suppressWarnings(as.numeric(cluster)),
    desc(cell_number)
  )

write.csv(
  candidate_prediction_overall,
  file.path(
    table_dir,
    "06d_CopyKAT_tumor_candidate_prediction_overall.csv"
  ),
  row.names = FALSE
)

print(candidate_prediction_overall)

# ============================================================
# J. UMAP：CopyKAT prediction
# ============================================================

copykat_colors <- c(
  "aneuploid" = "#D73027",
  "diploid" = "#4575B4",
  "not.defined" = "#BDBDBD",
  "not.run" = "#FEE08B"
)

p_copykat_umap <- DimPlot(
  object = sc,
  reduction = "umap_pca",
  group.by = "copykat_prediction",
  cols = copykat_colors,
  pt.size = 0.22,
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  ggtitle("CopyKAT prediction across all analyzed samples") +
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
    "06d_CopyKAT_prediction_UMAP.pdf"
  ),
  plot = p_copykat_umap,
  width = 12,
  height = 8
)

# ============================================================
# K. UMAP：只突出 tumor epithelial candidate 中的 aneuploid
# ============================================================

candidate_copykat_status <- rep(
  "Other_cells",
  ncol(sc)
)

names(candidate_copykat_status) <- colnames(sc)

candidate_cells <- colnames(sc)[
  as.character(
    sc[[cluster_column, drop = TRUE]]
  ) %in% tumor_candidate_clusters
]

candidate_copykat_status[
  candidate_cells
] <- "Tumor_candidate_non_aneuploid"

candidate_copykat_status[
  candidate_cells[
    sc$copykat_prediction[
      candidate_cells
    ] == "aneuploid"
  ]
] <- "Tumor_candidate_aneuploid"

sc$candidate_copykat_status <- factor(
  candidate_copykat_status,
  levels = c(
    "Tumor_candidate_aneuploid",
    "Tumor_candidate_non_aneuploid",
    "Other_cells"
  )
)

candidate_copykat_colors <- c(
  "Tumor_candidate_aneuploid" = "#D73027",
  "Tumor_candidate_non_aneuploid" = "#FDAE61",
  "Other_cells" = "#D9D9D9"
)

p_candidate_copykat_umap <- DimPlot(
  object = sc,
  reduction = "umap_pca",
  group.by = "candidate_copykat_status",
  cols = candidate_copykat_colors,
  pt.size = 0.22,
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  ggtitle("Aneuploid tumor epithelial candidates predicted by CopyKAT") +
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
    "06d_CopyKAT_aneuploid_tumor_candidate_UMAP.pdf"
  ),
  plot = p_candidate_copykat_umap,
  width = 12,
  height = 8
)

# ============================================================
# L. 定义严格 malignant call
# ============================================================
#
# 严格定义：
# 1. 位于 tumor epithelial candidate cluster
# 2. CopyKAT predicted aneuploid
#
# diploid / not.defined candidate 不直接删除，
# 仅标为 non-confirmed，保留用于敏感性分析。
# ============================================================

sc$malignant_call_copykat <- "Other_cells"

sc$malignant_call_copykat[
  candidate_cells
] <- "Tumor_candidate_not_confirmed"

sc$malignant_call_copykat[
  candidate_cells[
    sc$copykat_prediction[
      candidate_cells
    ] == "aneuploid"
  ]
] <- "Malignant_copykat_aneuploid"

sc$malignant_call_copykat <- factor(
  sc$malignant_call_copykat,
  levels = c(
    "Malignant_copykat_aneuploid",
    "Tumor_candidate_not_confirmed",
    "Other_cells"
  )
)

# ============================================================
# M. 保存对象与 session 信息
# ============================================================

sc$analysis_stage <- "CopyKAT_samplewise_malignant_call_completed"

saveRDS(
  sc,
  file.path(
    object_dir,
    "06d_GSE215403_CopyKAT_malignant_call.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "06d_CopyKAT_sessionInfo.txt"
  )
)

# ============================================================
# N. 完成提示
# ============================================================

message("\n============================================================")
message("06d_CopyKAT_malignant_call.R 运行完成。")
message("")
message("重点查看：")
message("1. results/tables/06d_CopyKAT_run_summary.csv")
message("2. results/tables/06d_CopyKAT_tumor_candidate_prediction_overall.csv")
message("3. results/tables/06d_CopyKAT_tumor_candidate_prediction_by_sample_cluster.csv")
message("4. results/figures/06d_CopyKAT_prediction_UMAP.pdf")
message("5. results/figures/06d_CopyKAT_aneuploid_tumor_candidate_UMAP.pdf")
message("============================================================\n")
