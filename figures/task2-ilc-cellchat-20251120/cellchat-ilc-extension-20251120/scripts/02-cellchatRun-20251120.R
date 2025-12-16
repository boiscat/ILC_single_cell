#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(future)
})

parse_args <- function(arg_vec) {
  res <- list()
  i <- 1
  while (i <= length(arg_vec)) {
    token <- arg_vec[[i]]
    if (startsWith(token, "--")) {
      token_clean <- substring(token, 3)
      if (grepl("=", token_clean, fixed = TRUE)) {
        parts <- strsplit(token_clean, "=", fixed = TRUE)[[1]]
        key <- parts[1]
        value <- if (length(parts) > 1 && nzchar(parts[2])) parts[2] else TRUE
      } else {
        key <- token_clean
        value <- TRUE
        if (i < length(arg_vec)) {
          nxt <- arg_vec[[i + 1]]
          if (!startsWith(nxt, "--")) {
            value <- nxt
            i <- i + 1
          }
        }
      }
      res[[key]] <- value
    }
    i <- i + 1
  }
  res
}
`%||%` <- function(a, b) if (!is.null(a)) a else b
to_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  if (is.logical(x)) return(x)
  val <- tolower(as.character(x))
  val %in% c("1", "true", "yes")
}
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", cmd_args[grep(file_arg, cmd_args)])
  if (length(script_path) > 0) {
    return(normalizePath(dirname(script_path)))
  }
  normalizePath(getwd())
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
opt$tag <- opt$tag %||% format(Sys.Date(), "%Y%m%d")
opt$conditions <- opt$conditions %||% "lesion,non_lesion,Healthy"
opt$condition_col <- opt$condition_col %||% "Site_Status"
opt$max_cells_tsub <- as.numeric(opt$max_cells_tsub %||% 2000)
opt$max_cells_all <- as.numeric(opt$max_cells_all %||% 1500)
opt$min_cells_condition <- as.numeric(opt$min_cells_condition %||% 200)
opt$skip_global <- to_bool(opt$skip_global, FALSE)
opt$reuse_existing <- to_bool(opt$reuse_existing, FALSE)

script_dir <- get_script_dir()
base_dir <- normalizePath(file.path(script_dir, ".."))
input_dir <- file.path(base_dir, "inputs")
output_dir <- file.path(base_dir, "outputs", "objects")
log_dir <- file.path(base_dir, "logs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

full_path <- file.path(input_dir, sprintf("full_cellchat_ready_%s.rds", opt$tag))
ilc_subset_path <- file.path(input_dir, sprintf("ilc_tcell_subset_%s.rds", opt$tag))
stopifnot(file.exists(full_path), file.exists(ilc_subset_path))

plan("sequential")
options(future.globals.maxSize = 6e9)
data("CellChatDB.human")
cellchat_db <- CellChatDB.human

meta_value <- function(meta, col) {
  if (!col %in% colnames(meta)) return(NULL)
  val <- meta[[col]]
  if (is.data.frame(val)) {
    if (col %in% colnames(val)) {
      val <- val[[col]]
    } else if (ncol(val) == 1) {
      val <- val[[1]]
    } else {
      stop(sprintf("列 %s 是 data.frame，无法自动解析。", col))
    }
  }
  as.character(val)
}

downsample_by_group <- function(obj, group_col, max_cells = 2000) {
  meta <- obj@meta.data
  stopifnot(group_col %in% colnames(meta))
  groups <- as.character(meta[[group_col]])
  meta$group_tmp <- groups
  meta <- meta[!is.na(meta$group_tmp), , drop = FALSE]
  keep <- unlist(lapply(split(rownames(meta), meta$group_tmp), function(ids) {
    if (length(ids) <= max_cells) return(ids)
    sample(ids, max_cells)
  }))
  subset(obj, cells = keep)
}

build_cellchat <- function(obj, group_col, min_cells = 1, population_size = FALSE) {
  meta <- obj@meta.data
  stopifnot(group_col %in% colnames(meta))
  meta$group_tmp <- as.character(meta[[group_col]])
  keep <- !is.na(meta$group_tmp)
  obj <- subset(obj, cells = rownames(meta)[keep])
  data.use <- GetAssayData(obj, assay = "RNA", slot = "data")
  meta.use <- obj@meta.data

  cc <- createCellChat(object = data.use, meta = meta.use, group.by = group_col)
  cc@DB <- cellchat_db
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc@options$population.size <- population_size
  cc <- computeCommunProb(cc, population.size = population_size)
  cc <- filterCommunication(cc, min.cells = min_cells, nonFilter.keep = TRUE)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  invisible(rankNet(cc, mode = "single"))
  cc
}

subset_condition <- function(obj, condition, cond_col) {
  meta <- obj@meta.data
  cond_vec <- meta_value(meta, cond_col)
  if (is.null(cond_vec)) {
    stop("对象中缺少列：", cond_col)
  }
  keep <- !is.na(cond_vec) & cond_vec == condition
  if (!any(keep)) return(NULL)
  subset(obj, cells = rownames(meta)[keep])
}

message(">>> 载入对象")
full_obj <- readRDS(full_path)
ilc_subset_obj <- readRDS(ilc_subset_path)

cond_vec_full <- meta_value(full_obj@meta.data, opt$condition_col)
if (is.null(cond_vec_full)) {
  stop("full 对象缺少列：", opt$condition_col)
}
cond_vec_ilc <- meta_value(ilc_subset_obj@meta.data, opt$condition_col)
if (is.null(cond_vec_ilc)) {
  stop("ILC/T cell 对象缺少列：", opt$condition_col)
}

available_conditions <- sort(unique(cond_vec_full))
target_conditions <- strsplit(opt$conditions, ",")[[1]]
target_conditions <- trimws(target_conditions)
target_conditions <- target_conditions[target_conditions != ""]
target_conditions <- unique(target_conditions)
target_conditions <- intersect(target_conditions, available_conditions)
if (length(target_conditions) == 0) {
  stop("没有匹配的 Site_Status，可用值：", paste(available_conditions, collapse = ", "))
}

summary_df <- data.frame(
  Site_Status = target_conditions,
  n_cells_full = as.numeric(table(factor(cond_vec_full, levels = target_conditions))),
  n_cells_ilc = as.numeric(table(factor(cond_vec_ilc, levels = target_conditions))),
  stringsAsFactors = FALSE
)
summary_path <- file.path(log_dir, sprintf("cellchat_condition_summary_%s.csv", opt$tag))
write.csv(summary_df, summary_path, row.names = FALSE)
message("条件细胞数量已记录：", summary_path)

extract_and_save <- function(cellchat_obj, prefix, condition) {
  suffix <- if (condition == "all" || is.null(condition)) "" else paste0("_", condition)
  out_csv <- file.path(output_dir, sprintf("%s%s_communications_%s.csv", prefix, suffix, opt$tag))
  comm <- subsetCommunication(cellchat_obj)
  write.csv(comm, out_csv, row.names = FALSE)
  message("保存通讯表：", out_csv)
}

cellchat_tsub_all <- NULL
cellchat_all_all <- NULL

message(">>> 构建 ILC + T cell CellChat（全部）")
if (!opt$skip_global) {
  ilc_ds <- downsample_by_group(ilc_subset_obj, "cellchat_celltype", max_cells = opt$max_cells_tsub)
  cellchat_tsub_all <- build_cellchat(ilc_ds, "cellchat_celltype", min_cells = 1, population_size = FALSE)
  tcell_groups <- sort(unique(ilc_ds@meta.data$cellchat_celltype[grepl("^Tcell_", ilc_ds@meta.data$cellchat_celltype)]))
  weight_mat <- cellchat_tsub_all@net$weight
  if ("ILC" %in% rownames(weight_mat) && length(tcell_groups) > 0) {
    available_targets <- intersect(tcell_groups, colnames(weight_mat))
    ilc_summary <- data.frame(
      target = available_targets,
      weight = as.numeric(weight_mat["ILC", available_targets]),
      stringsAsFactors = FALSE
    )
    if (nrow(ilc_summary) > 0) {
      max_weight <- max(ilc_summary$weight, na.rm = TRUE)
      ilc_summary$normalized_weight <- if (max_weight > 0) ilc_summary$weight / max_weight else 0
      ilc_summary$is_cd4_naive <- ilc_summary$target == "Tcell_CD4_NAIVE"
      summary_path <- file.path(output_dir, sprintf("ilc_to_tcell_weights_%s.csv", opt$tag))
      write.csv(ilc_summary, summary_path, row.names = FALSE)
      message("保存 ILC->T cell 权重表（包含 CD4-naive）：", summary_path)
      cd4_row <- ilc_summary[ilc_summary$is_cd4_naive, , drop = FALSE]
      if (nrow(cd4_row) == 1) {
        message(sprintf("ILC->CD4-naive 权重 %.4f（归一化 %.3f）", cd4_row$weight, cd4_row$normalized_weight))
      } else {
        message("未找到 Tcell_CD4_NAIVE 作为目标，无法记录 CD4-naive 关系。")
      }
    }
  }
  tcell_out <- file.path(output_dir, sprintf("cellchat_ILC_Tsub_%s.rds", opt$tag))
  saveRDS(cellchat_tsub_all, tcell_out)

  full_ds <- downsample_by_group(full_obj, "cellchat_celltype", max_cells = opt$max_cells_all)
  cellchat_all_all <- build_cellchat(full_ds, "cellchat_celltype", min_cells = 5, population_size = FALSE)
  all_out <- file.path(output_dir, sprintf("cellchat_ILC_all_%s.rds", opt$tag))
  saveRDS(cellchat_all_all, all_out)

  group_file <- file.path(output_dir, sprintf("cellchat_tcell_groups_%s.txt", opt$tag))
  writeLines(tcell_groups, con = group_file)
  message("T cell 组列表：", group_file)

  extract_and_save(cellchat_tsub_all, "cellchat_ILC_Tsub", "all")
  extract_and_save(cellchat_all_all, "cellchat_ILC_all", "all")
} else {
  message(">>> 跳过全量对象构建，使用现有 RDS。")
  tsub_file <- file.path(output_dir, sprintf("cellchat_ILC_Tsub_%s.rds", opt$tag))
  all_file <- file.path(output_dir, sprintf("cellchat_ILC_all_%s.rds", opt$tag))
  if (file.exists(tsub_file)) {
    cellchat_tsub_all <- readRDS(tsub_file)
    extract_and_save(cellchat_tsub_all, "cellchat_ILC_Tsub", "all")
  } else {
    warning("未找到 ", tsub_file, "，无法生成全局 Tsub 通讯表。")
  }
  if (file.exists(all_file)) {
    cellchat_all_all <- readRDS(all_file)
    extract_and_save(cellchat_all_all, "cellchat_ILC_all", "all")
  } else {
    warning("未找到 ", all_file, "，无法生成全局全细胞通讯表。")
  }
}

for (cond in target_conditions) {
  message(">>> 处理条件：", cond)
  cond_out <- file.path(output_dir, sprintf("cellchat_ILC_Tsub_%s_%s.rds", cond, opt$tag))
  build_tsub <- TRUE
  if (opt$reuse_existing && file.exists(cond_out)) {
    message("  - 发现既有 ILC/T cell 对象，跳过重新构建：", cond_out)
    build_tsub <- FALSE
  }
  if (build_tsub) {
    ilc_cond <- subset_condition(ilc_subset_obj, cond, opt$condition_col)
    if (is.null(ilc_cond) || ncol(ilc_cond) < opt$min_cells_condition) {
      message("  - ILC/T cell 子集细胞数不足，跳过：", cond)
    } else {
      ilc_cond <- downsample_by_group(ilc_cond, "cellchat_celltype", max_cells = opt$max_cells_tsub)
      cellchat_cond <- build_cellchat(ilc_cond, "cellchat_celltype", min_cells = 1, population_size = FALSE)
      saveRDS(cellchat_cond, cond_out)
      message("  - ILC/T cell CellChat 已保存：", cond_out)
    }
  }
  if (file.exists(cond_out)) {
    cc_obj <- readRDS(cond_out)
    extract_and_save(cc_obj, "cellchat_ILC_Tsub", cond)
  }

  cond_all_out <- file.path(output_dir, sprintf("cellchat_ILC_all_%s_%s.rds", cond, opt$tag))
  build_all <- TRUE
  if (opt$reuse_existing && file.exists(cond_all_out)) {
    message("  - 发现既有全细胞对象，跳过重新构建：", cond_all_out)
    build_all <- FALSE
  }
  if (build_all) {
    full_cond <- subset_condition(full_obj, cond, opt$condition_col)
    if (is.null(full_cond) || ncol(full_cond) < opt$min_cells_condition) {
      message("  - 全细胞子集细胞数不足，跳过：", cond)
    } else {
      full_cond <- downsample_by_group(full_cond, "cellchat_celltype", max_cells = opt$max_cells_all)
      cellchat_full_cond <- build_cellchat(full_cond, "cellchat_celltype", min_cells = 5, population_size = FALSE)
      saveRDS(cellchat_full_cond, cond_all_out)
      message("  - 全细胞 CellChat 已保存：", cond_all_out)
    }
  }
  if (file.exists(cond_all_out)) {
    cc_obj <- readRDS(cond_all_out)
    extract_and_save(cc_obj, "cellchat_ILC_all", cond)
  }
}

session_path <- file.path(log_dir, sprintf("sessionInfo_cellchat_%s.txt", opt$tag))
writeLines(capture.output(sessionInfo()), con = session_path)
message("sessionInfo 保存：", session_path)
message("完成 CellChat 运行。")
