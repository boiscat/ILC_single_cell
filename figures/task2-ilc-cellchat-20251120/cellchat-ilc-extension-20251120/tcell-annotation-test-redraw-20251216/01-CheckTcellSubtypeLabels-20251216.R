#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

pick_latest_input <- function(input_dir) {
  files <- list.files(
    input_dir,
    pattern = "^subset_obj_tcell_subtyped_[0-9]{8}\\.rds$",
    full.names = TRUE
  )
  if (length(files) == 0) return(NULL)
  files[order(files, decreasing = TRUE)][[1]]
}

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) args[[1]] else NULL

base_dir <- normalizePath(file.path(getwd()))
input_dir <- file.path(base_dir, "inputs")
if (is.null(input_path) || !nzchar(input_path)) {
  input_path <- pick_latest_input(input_dir)
}
if (is.null(input_path) || !file.exists(input_path)) {
  stop("找不到输入文件：请传入 inputs/subset_obj_tcell_subtyped_*.rds 作为第一个参数")
}

obj <- readRDS(input_path)
group_col <- if ("tcell_subtype" %in% colnames(obj@meta.data)) {
  "tcell_subtype"
} else if ("Tcell_subtype" %in% colnames(obj@meta.data)) {
  "Tcell_subtype"
} else {
  stop("对象 meta.data 缺少 tcell_subtype / Tcell_subtype 列")
}

label_map <- c(
  "CD4-naive" = "Naive CD4+T cell",
  "Cytotoxic CD8 T" = "Cytotoxic CD8+T cell",
  "exhausted CD8 T" = "Exhuasted CD8+T cell",
  "iNKT" = "iNKT cell",
  "Th17" = "Th17",
  "Treg" = "Treg",
  "TSCM" = "Tscm"
)

old_vals <- as.character(obj@meta.data[[group_col]])
new_vals <- old_vals
hit <- new_vals %in% names(label_map)
new_vals[hit] <- unname(label_map[new_vals[hit]])

expected <- unname(label_map)
missing <- setdiff(expected, unique(new_vals))
message("输入：", input_path)
message("分组列：", group_col)
message("旧标签出现：", paste(sort(unique(intersect(names(label_map), old_vals))), collapse = ", "))
message("新标签出现：", paste(sort(unique(intersect(expected, new_vals))), collapse = ", "))
if (length(missing) > 0) {
  warning("以下目标标签在该输入对象中未出现（通常表示该亚群在当前数据里不存在，可忽略）：", paste(missing, collapse = ", "))
}
message("OK：标签映射逻辑可运行。")

