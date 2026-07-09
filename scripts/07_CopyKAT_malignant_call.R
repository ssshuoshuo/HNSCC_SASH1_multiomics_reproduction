# 07_CopyKAT_malignant_call.R

# 本脚本功能：
# 1. 读取06定义tumor epithelial candidate后的Seurat object
# 2. 按sample_id单独运行CopyKAT
# 3. 使用同一样本内的非tumor epithelial candidate细胞作为known normal reference
# 4. 获得CopyKAT的diploid/aneuploid prediction
# 5. 将CopyKAT预测结果回填到Seurat metadata
# 6. 重点汇总cluster 2、3、4、6、11的CopyKAT预测结果
# 7. 输出CopyKAT诊断UMAP和candidate aneuploid UMAP
# 8. 保存供后续final malignant call使用的对象

# 本项目专用数据：
# GSE215403
# 12个OSCC单细胞样本：
# OSCC, scB1, scB2, scB5, scB7, scB8,
# scB9, scB10, scB12, scB13, scB14, scB15
#
# CopyKAT prediction是基于scRNA-seq推断的CNV/aneuploidy结果。
# 后续final malignant cell定义会结合：
# - tumor epithelial candidate identity
# - CopyKAT aneuploid prediction
# - sample-level consistency
#
# 通用代码修改位置：
# 1. 换数据集时：
#    修改input_object_file、cluster_column和sample_column
#
# 2. 换candidate定义时：
#    修改tumor_candidate_clusters
#
# 3. 调整CopyKAT运行门槛时：
#    修改min_tumor_cells和min_normal_cells
#
# 4. 调整CopyKAT参数时：
#    修改copykat()中的ngene.chr、win.size、KS.cut、genome和n.cores


# ============================================================
# A. 加载包
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
      "缺少R包：",
      paste(missing_packages, collapse = ", "),
      "\n请先完成CopyKAT安装。"
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
# C. 读取06对象
# ============================================================

input_object_file <- file.path(
  object_dir,
  "06_malignant_candidate_diagnostic.rds"
)

if (!file.exists(input_object_file)) {
  
  stop(
    paste0(
      "找不到06对象：\n",
      input_object_file,
      "\n请先运行06_malignant_candidate_diagnostic.R"
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
      "缺少metadata：",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

# ============================================================
# D. 设置CopyKAT输入参数
# ============================================================

tumor_candidate_clusters <- c(
  "2",
  "3",
  "4",
  "6",
  "11"
)

# 每个sample至少需要：
# - 50个tumor epithelial candidate
# - 50个known normal reference cells
#
# 不满足时跳过，避免CopyKAT在极小样本中产生不稳定预测。

min_tumor_cells <- 50
min_normal_cells <- 50

# 当前设置为单核运行。
# 如果服务器或本地环境允许多核，可以相应提高。

copykat_cores <- 1

set.seed(1234)

# ============================================================
# E. 获取raw count matrix
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
    "07_CopyKAT_sample_input_summary.csv"
  ),
  row.names = FALSE
)

print(sample_summary)

# ============================================================
# G. 按sample独立运行CopyKAT
# ============================================================

# 关键点：
# 1. 每个sample独立运行，避免patient-to-patient差异混入。
# 2. non-tumor candidate细胞作为known normal reference。
# 3. 如果某个sample已经有完整prediction文件，则直接复用。

prediction_list <- list()
run_summary_list <- list()

for (current_sample in sample_ids) {
  
  message("\n============================================================")
  message("开始CopyKAT：", current_sample)
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
  # 如果该sample已存在完整且细胞数匹配的prediction文件，
  # 则直接读取并跳过CopyKAT重算。
  # ------------------------------------------------------------
  
  existing_prediction_file <- file.path(
    copykat_dir,
    paste0(
      "07_CopyKAT_",
      current_sample
    ),
    paste0(
      "07_CopyKAT_",
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
        "检测到已完成CopyKAT输出，跳过重算：",
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
        "07_CopyKAT_",
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
    
    # CopyKAT官方示例使用raw UMI matrix。
    # 输入需要dense matrix，因此仅在单个sample内转换。
    
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
          "CopyKAT失败：",
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
            "07_CopyKAT_",
            current_sample,
            "_result.rds"
          )
        )
      )
      
      current_prediction <- as.data.frame(
        copykat_result$prediction,
        stringsAsFactors = FALSE
      )
      
      # CopyKAT通常输出cell.names和copykat.pred。
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
              "07_CopyKAT_",
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
    "07_CopyKAT_run_summary.csv"
  ),
  row.names = FALSE
)

print(run_summary)

# ============================================================
# H. 合并CopyKAT prediction
# ============================================================

if (length(prediction_list) == 0) {
  
  stop(
    "没有任何sample成功完成CopyKAT。请检查07_CopyKAT_run_summary.csv。"
  )
}

copykat_prediction_all <- bind_rows(
  prediction_list
)

write.csv(
  copykat_prediction_all,
  file.path(
    table_dir,
    "07_CopyKAT_all_cell_predictions.csv"
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
# I. 汇总tumor epithelial candidate的CopyKAT结果
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
    "07_CopyKAT_tumor_candidate_prediction_by_sample_cluster.csv"
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
    "07_CopyKAT_tumor_candidate_prediction_overall.csv"
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
    "07_CopyKAT_prediction_UMAP.pdf"
  ),
  plot = p_copykat_umap,
  width = 12,
  height = 8
)

# ============================================================
# K. UMAP：突出tumor epithelial candidate中的aneuploid细胞
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
    "07_CopyKAT_aneuploid_tumor_candidate_UMAP.pdf"
  ),
  plot = p_candidate_copykat_umap,
  width = 12,
  height = 8
)

# ============================================================
# L. 定义严格malignant call
# ============================================================

# 严格定义：
# 1. 位于tumor epithelial candidate cluster
# 2. CopyKAT predicted aneuploid
#
# diploid/not.defined candidate不直接删除。
# 它们标记为non-confirmed，保留用于敏感性分析。

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
# M. 保存对象和环境信息
# ============================================================

sc$analysis_stage <- "CopyKAT_samplewise_malignant_call_completed"

saveRDS(
  sc,
  file.path(
    object_dir,
    "07_CopyKAT_malignant_call.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    table_dir,
    "07_CopyKAT_sessionInfo.txt"
  )
)

# ============================================================
# N. 最终提示
# ============================================================

message("\n============================================================")
message("07_CopyKAT_malignant_call.R 运行完成。")
message("")
message("已保存对象：")
message("results/objects/07_CopyKAT_malignant_call.rds")
message("")
message("请重点查看：")
message("1. results/tables/07_CopyKAT_run_summary.csv")
message("2. results/tables/07_CopyKAT_tumor_candidate_prediction_overall.csv")
message("3. results/tables/07_CopyKAT_tumor_candidate_prediction_by_sample_cluster.csv")
message("4. results/figures/07_CopyKAT_prediction_UMAP.pdf")
message("5. results/figures/07_CopyKAT_aneuploid_tumor_candidate_UMAP.pdf")
message("============================================================\n")