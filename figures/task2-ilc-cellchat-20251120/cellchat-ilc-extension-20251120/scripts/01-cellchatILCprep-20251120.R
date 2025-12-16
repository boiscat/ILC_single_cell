#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

parse_args <- function(arg_vec) {
  res <- list()
  i <- 1
  while (i <= length(arg_vec)) {
    token <- arg_vec[[i]]
    if (startsWith(token, "--")) {
      key <- substring(token, 3)
      value <- TRUE
      if (i < length(arg_vec)) {
        nxt <- arg_vec[[i + 1]]
        if (!startsWith(nxt, "--")) {
          value <- nxt
          i <- i + 1
        }
      }
      res[[key]] <- value
    }
    i <- i + 1
  }
  res
}
`%||%` <- function(a, b) if (!is.null(a)) a else b

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
script_dir <- get_script_dir()
base_dir <- normalizePath(file.path(script_dir, ".."))
input_dir <- file.path(base_dir, "inputs")
log_dir <- file.path(base_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

default_full <- file.path(input_dir, "full_obj_saved_1117.rds")
tcell_candidates <- c(
  "subset_obj_tcell_subtyped_20251204.rds",
  "subset_obj_tcell_subtyped_20251203.rds",
  "subset_obj_Tcell_clustered_annoed_1117.rds"
)
tcell_default <- NULL
for (cand in tcell_candidates) {
  cand_path <- file.path(input_dir, cand)
  if (file.exists(cand_path)) {
    tcell_default <- cand_path
    break
  }
}

opt$full_rds <- opt$full_rds %||% default_full
opt$tcell_rds <- opt$tcell_rds %||% tcell_default
if (is.null(opt$tcell_rds)) {
  stop("inputs/ 下没有找到新的 T cell 注释文件，请通过 --tcell-rds 指定。")
}
stopifnot(file.exists(opt$full_rds), file.exists(opt$tcell_rds))

message(">>> 使用 full 对象：", opt$full_rds)
message(">>> 使用 T cell 对象：", opt$tcell_rds)
full_obj <- readRDS(opt$full_rds)
tcell_obj <- readRDS(opt$tcell_rds)

required_cols <- c("cell_anno", "Site_Status")
missing_cols <- setdiff(required_cols, colnames(full_obj@meta.data))
if (length(missing_cols) > 0) {
  stop("full_obj 缺少必要列：", paste(missing_cols, collapse = ", "))
}

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

major_vec <- meta_value(full_obj@meta.data, "major_annotation")
if (is.null(major_vec)) {
  warning("full_obj 缺少 major_annotation 列，改用 cell_anno 作为粗分类。")
  major_vec <- meta_value(full_obj@meta.data, "cell_anno")
}
anno_vec <- meta_value(full_obj@meta.data, "cell_anno")

subtype_col <- if ("Tcell_subtype" %in% colnames(tcell_obj@meta.data)) {
  "Tcell_subtype"
} else if ("tcell_subtype" %in% colnames(tcell_obj@meta.data)) {
  "tcell_subtype"
} else {
  stop("T cell 对象缺少 tcell_subtype / Tcell_subtype 列")
}

message(">>> 匹配 T cell 亚型信息")
shared_cells <- intersect(colnames(full_obj), colnames(tcell_obj))
if (length(shared_cells) == 0) {
  stop("full 与 T cell 对象没有共享细胞名称，无法映射亚型")
}

map_subtype <- function(x) {
  clean <- iconv(x, to = "ASCII//TRANSLIT")
  clean <- toupper(clean)
  clean <- gsub("[^A-Z0-9]+", "_", clean)
  clean <- gsub("_+", "_", clean)
  clean <- gsub("^_|_$", "", clean)
  ifelse(is.na(clean) | clean == "", NA_character_, paste0("Tcell_", clean))
}

tcell_labels <- map_subtype(tcell_obj@meta.data[shared_cells, subtype_col, drop = TRUE])
full_meta <- full_obj@meta.data
full_meta$Tcell_subtype_label <- NA_character_
full_meta[shared_cells, "Tcell_subtype_label"] <- tcell_labels
site_vec <- meta_value(full_meta, "Site_Status")
if (is.null(site_vec)) {
  stop("meta data 中没有 Site_Status 列，无法继续。")
}
full_meta$Site_Status <- site_vec

is_tcell <- major_vec == "T_cells" | anno_vec %in% c("T", "T_cells")
full_meta$cellchat_celltype <- as.character(major_vec)
full_meta$cellchat_celltype[is.na(full_meta$cellchat_celltype)] <- anno_vec[is.na(full_meta$cellchat_celltype)]

full_meta$cellchat_celltype[is_tcell] <- ifelse(
  is.na(full_meta$Tcell_subtype_label[is_tcell]),
  "Tcell_unspecified",
  full_meta$Tcell_subtype_label[is_tcell]
)

bad_label <- grepl("unknown|unassigned|unspecified", full_meta$cellchat_celltype, ignore.case = TRUE)
bad_label[is.na(bad_label)] <- FALSE
message(sprintf("在 cellchat_celltype 中识别到 %d 个未知/未分配/unspecified 记录，后续子集构建中将排除。", sum(bad_label)))

full_clustering <- meta_value(full_meta, "full_clustering")
if (!is.null(full_clustering)) {
  ilc_idx <- grepl("ILC", full_clustering, ignore.case = TRUE)
  full_meta$cellchat_celltype[ilc_idx] <- "ILC"
}

full_obj@meta.data <- full_meta

message(">>> 保存全量对象（附加 CellChat 元信息）")
full_ready_path <- file.path(input_dir, sprintf("full_cellchat_ready_%s.rds", opt$tag))
saveRDS(full_obj, full_ready_path)

message(">>> 提取 ILC + T cell 亚型子集")
is_ilc <- full_meta$cellchat_celltype == "ILC"
is_tcell <- grepl("^Tcell_", full_meta$cellchat_celltype)
is_ilc[is.na(is_ilc)] <- FALSE
is_tcell[is.na(is_tcell)] <- FALSE
keep_cells <- (is_ilc | is_tcell) & !bad_label
keep_cells[is.na(keep_cells)] <- FALSE
ilc_cells <- rownames(full_meta)[keep_cells]
ilc_tcell_obj <- subset(full_obj, cells = ilc_cells)
ilc_subset_path <- file.path(input_dir, sprintf("ilc_tcell_subset_%s.rds", opt$tag))
saveRDS(ilc_tcell_obj, ilc_subset_path)
cd4_total <- sum(full_meta$cellchat_celltype == "Tcell_CD4_NAIVE", na.rm = TRUE)
cd4_subset <- sum(ilc_tcell_obj@meta.data$cellchat_celltype == "Tcell_CD4_NAIVE", na.rm = TRUE)
message(sprintf("CD4-naive 细胞共 %d 个，其中 %d 个被保留进 CellChat 子集，以增强 ILC ↔ CD4-naive 关系的可追踪度。", cd4_total, cd4_subset))

message(">>> 汇总细胞数量")
cell_counts <- as.data.frame(table(full_meta$cellchat_celltype), stringsAsFactors = FALSE)
colnames(cell_counts) <- c("cellchat_celltype", "n_cells")
count_path <- file.path(log_dir, sprintf("cellchat_prep_counts_%s.csv", opt$tag))
write.csv(cell_counts, count_path, row.names = FALSE)

condition_counts <- as.data.frame(table(full_meta$Site_Status, full_meta$cellchat_celltype))
colnames(condition_counts) <- c("Site_Status", "cellchat_celltype", "n_cells")
condition_path <- file.path(log_dir, sprintf("cellchat_prep_condition_counts_%s.csv", opt$tag))
write.csv(condition_counts, condition_path, row.names = FALSE)

message("完成：", full_ready_path)
message("完成：", ilc_subset_path)
message("计数表：", count_path)
message("分组计数：", condition_path)
