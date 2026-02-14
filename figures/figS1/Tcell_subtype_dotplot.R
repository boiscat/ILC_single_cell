#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# 简化后的参数解析函数
parse_args <- function(arg_vec) {
  res <- list(
    tag = format(Sys.Date(), "%Y%m%d"),
    input = NULL,
    outdir = NULL,
    group_col = NULL,
    width_dot = 11,
    height_dot = 5.5
  )
  
  i <- 1
  while (i <= length(arg_vec)) {
    token <- arg_vec[[i]]
    if (startsWith(token, "--")) {
      token_clean <- substring(token, 3)
      key <- NA_character_
      value <- NULL
      if (grepl("=", token_clean, fixed = TRUE)) {
        parts <- strsplit(token_clean, "=", fixed = TRUE)[[1]]
        key <- parts[1]
        if (length(parts) > 1 && nzchar(parts[2])) {
          value <- parts[2]
        }
      } else {
        key <- token_clean
        if (i < length(arg_vec)) {
          next_token <- arg_vec[[i + 1]]
          if (!startsWith(next_token, "--")) {
            value <- next_token
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

# 获取脚本所在目录
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", cmd_args[grep(file_arg, cmd_args)])
  if (length(script_path) > 0) {
    return(normalizePath(dirname(script_path)))
  }
  normalizePath(getwd())
}

# 读取输入对象路径
opt <- parse_args(commandArgs(trailingOnly = TRUE))
script_dir <- get_script_dir()
base_dir <- normalizePath(file.path(script_dir, ".."))
input_dir <- file.path(base_dir, "inputs")
plot_dir <- opt$outdir %||% file.path(base_dir, "outputs", "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# 获取输入文件路径
input_path <- opt$input
if (is.null(input_path) || identical(input_path, TRUE) || !nzchar(input_path)) {
  input_path <- file.path(input_dir, sprintf("subset_obj_tcell_subtyped_%s.rds", opt$tag))
  if (!file.exists(input_path)) {
    stop("找不到输入文件。请确保输入路径存在。")
  }
}

obj <- readRDS(input_path)

# 获取分组列并重命名
group_col <- opt$group_col
if (is.null(group_col) || identical(group_col, TRUE) || !nzchar(group_col)) {
  stop("请提供分组列名 --group_col")
}

# Marker genes
marker_genes <- c(
  "CD3D", "CD3E", "CD3G", "CD4", "LEF1", "TCF7", "CD8A", "CD8B", "NKG7", "LAG3", 
  "TIGIT", "ZBTB16", "RIPK2", "IL17A", "IL17F", "FOXP3", "SELL", "CD28", "FAS", 
  "TRDV1", "TRGV3", "TRGV4", "TRGV8"
)

# 创建 DotPlot
p_dot <- DotPlot(obj, features = marker_genes) +
  scale_color_gradientn(
    values = rescale(c(-2, -0.1, 0, 0.1, 2)), 
    colors = c("#5E92CD", "#ADC9E6", "white", "#F39894", "#E71F19"),
    limits = c(-2, 2), oob = squish
  ) +
  RotatedAxis() +
  ggtitle(sprintf("DotPlot of marker genes across %s", group_col)) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)
  )

# 保存 DotPlot 图
dot_file_png <- file.path(plot_dir, "DotPlot_tcell_subtype_markers.png")
dot_file_pdf <- file.path(plot_dir, "DotPlot_tcell_subtype_markers.pdf")

ggsave(dot_file_pdf, plot = p_dot, width = as.numeric(opt$width_dot), height = as.numeric(opt$height_dot), units = "in", dpi = 600, bg = "white")
ggsave(dot_file_png, plot = p_dot, width = as.numeric(opt$width_dot), height = as.numeric(opt$height_dot), units = "in", dpi = 600, bg = "white")

# 输出路径
message("DotPlot 输出：", dot_file_png)
