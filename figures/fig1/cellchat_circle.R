#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
})

# 解析简化后的参数
parse_args <- function(arg_vec) {
  tag <- format(Sys.Date(), "%Y%m%d")
  input_tag <- NULL
  bubble_width <- 8.5
  bubble_height <- 9

  get_one <- function(prefix) {
    hit <- grep(paste0("^", prefix, "="), arg_vec, value = TRUE)
    if (length(hit) == 1) sub(paste0("^", prefix, "="), "", hit) else NULL
  }

  # 获取基本参数
  input_tag <- get_one("--input-tag") %||% tag
  bubble_width_arg <- get_one("--bubble-width") %||% bubble_width
  bubble_height_arg <- get_one("--bubble-height") %||% bubble_height

  list(input_tag = input_tag, bubble_width = as.numeric(bubble_width_arg), bubble_height = as.numeric(bubble_height_arg))
}

# 获取脚本路径
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
script_dir <- get_script_dir()
base_dir <- normalizePath(file.path(script_dir, ".."))
object_dir <- file.path(base_dir, "outputs", "objects")
plot_dir <- file.path(base_dir, "outputs", "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# 读取 CellChat 对象
tcell_obj_path <- file.path(object_dir, sprintf("cellchat_ILC_Tsub_%s.rds", opt$input_tag))
all_obj_path <- file.path(object_dir, sprintf("cellchat_ILC_all_%s.rds", opt$input_tag))
stopifnot(file.exists(tcell_obj_path), file.exists(all_obj_path))

cellchat_tsub <- readRDS(tcell_obj_path)
cellchat_all <- readRDS(all_obj_path)

# 创建基础颜色映射
group_cols <- setNames(RColorBrewer::brewer.pal(8, "Set2"), levels(cellchat_all@idents))

# 选择要展示的图形类型（仅保留 bubble 图）
run_bubble <- TRUE

# 绘制 Bubble 图
if (run_bubble) {
  bubble_args <- list(
    object = cellchat_tsub,
    sources.use = "ILC",
    targets.use = levels(cellchat_tsub@idents),
    angle.x = 45,
    remove.isolate = FALSE
  )
  
  # 创建 Bubble 图
  bubble_plot <- do.call(netVisual_bubble, bubble_args) +
    scale_x_discrete(labels = function(x) gsub("Tcell_", "", x)) +  # 简化标签
    scale_color_gradientn(colors = c("#FCB8B8", "#E4151D"), name = "Commun. Prob.") +
    scale_fill_gradientn(colors = c("#FCB8B8", "#E4151D"), name = "Commun. Prob.") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "ILC to T cell pathways",
      subtitle = "Data shown as communication probability"
    )

  # 保存为 PDF 和 PNG
  bubble_pdf <- file.path(plot_dir, sprintf("Fig_ILC_to_Tcell_bubble_%s.pdf", opt$input_tag))
  bubble_png <- file.path(plot_dir, sprintf("Fig_ILC_to_Tcell_bubble_%s.png", opt$input_tag))
  ggsave(bubble_pdf, plot = bubble_plot, width = opt$bubble_width, height = opt$bubble_height, dpi = 600, units = "in")
  ggsave(bubble_png, plot = bubble_plot, width = opt$bubble_width, height = opt$bubble_height, dpi = 600, units = "in")
  
  message("Bubble 图已保存：", bubble_pdf, " / ", bubble_png)
}

