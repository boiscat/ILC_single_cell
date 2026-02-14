#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0 && !identical(a, "") && !is.na(a)) a else b
}

parse_args <- function(arg_vec) {
  res <- list(
    tag = format(Sys.Date(), "%Y%m%d"),
    input = NULL,
    outdir = NULL,
    group_col = NULL,
    width_dim = 7.5,
    height_dim = 5.5,
    width_dot = 11,
    height_dot = 5
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
      if (!is.na(key)) {
        if (is.null(value)) value <- TRUE
        res[[key]] <- value
      }
    }
    i <- i + 1
  }
  res
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

pick_latest_input <- function(input_dir) {
  files <- list.files(
    input_dir,
    pattern = "^subset_obj_tcell_subtyped_[0-9]{8}\\.rds$",
    full.names = TRUE
  )
  if (length(files) == 0) return(NULL)
  files[order(files, decreasing = TRUE)][[1]]
}

rename_tcell_subtype <- function(obj, group_col) {
  if (!group_col %in% colnames(obj@meta.data)) {
    stop("对象 meta.data 中找不到列：", group_col)
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
  desired_order <- c(
    "Naive CD4+T cell",
    "Cytotoxic CD8+T cell",
    "Exhuasted CD8+T cell",
    "iNKT cell",
    "Th17",
    "Treg",
    "Tscm"
  )

  old_vals <- as.character(obj@meta.data[[group_col]])
  new_vals <- old_vals
  hit <- new_vals %in% names(label_map)
  new_vals[hit] <- unname(label_map[new_vals[hit]])

  present_new <- unique(new_vals)
  levels_final <- c(desired_order, setdiff(sort(present_new), desired_order))
  obj@meta.data[[group_col]] <- factor(new_vals, levels = levels_final)
  Idents(obj) <- group_col

  missing_old <- setdiff(names(label_map), unique(old_vals))
  if (length(missing_old) > 0) {
    warning("输入对象中未出现的旧标签（可忽略）：", paste(missing_old, collapse = ", "))
  }
  missing_new <- setdiff(desired_order, unique(new_vals))
  if (length(missing_new) > 0) {
    warning("输入对象中未出现的目标标签（可忽略）：", paste(missing_new, collapse = ", "))
  }
  obj
}
#3K颜色
make_original_dimplot_palette <- function(levels_vec) {
  # 颜色需与原始 `subset_obj_Tcell_DimPlot_by_cellanno_nature_palette.png` 对齐
 old_palette <- c(
    "CD4-naive" = "#EC7576",
    "Cytotoxic CD8 T" = "#fcb3b3",
    "exhausted CD8 T" = "#C6A3CE",
    "iNKT" = "#6DA443",
    "Th17" = "#FC8D62",
    "Treg" = "#14A0BD",
    "TSCM" = "#5E92CD",
    "γδ T" = "#e5e55aff"
  )
  rename_map <- c(
    "CD4-naive" = "Naive CD4+T cell",
    "Cytotoxic CD8 T" = "Cytotoxic CD8+T cell",
    "exhausted CD8 T" = "Exhuasted CD8+T cell",
    "iNKT" = "iNKT cell",
    "Th17" = "Th17",
    "Treg" = "Treg",
    "TSCM" = "Tscm",
    "γδ T" = "γδ T"
  )
  new_palette <- old_palette
  names(new_palette) <- unname(rename_map[names(old_palette)])

  pal <- new_palette[intersect(levels_vec, names(new_palette))]
  extras <- setdiff(levels_vec, names(new_palette))
  if (length(extras) > 0) {
    extra_cols <- grDevices::hcl.colors(length(extras), palette = "Dark 3")
    pal <- c(pal, setNames(extra_cols, extras))
  }
  pal
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
script_dir <- get_script_dir()
base_dir <- normalizePath(file.path(script_dir, ".."))
input_dir <- file.path(base_dir, "inputs")
plot_dir <- opt$outdir %||% file.path(base_dir, "outputs", "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- opt$input
if (is.null(input_path) || identical(input_path, TRUE) || !nzchar(input_path)) {
  input_path <- file.path(input_dir, sprintf("subset_obj_tcell_subtyped_%s.rds", opt$tag))
  if (!file.exists(input_path)) {
    input_path <- pick_latest_input(input_dir)
  }
}
if (is.null(input_path) || !file.exists(input_path)) {
  stop("找不到输入文件。请用 --input 指定，或确保 inputs/ 下存在 subset_obj_tcell_subtyped_*.rds")
}

obj <- readRDS(input_path)

group_col <- opt$group_col
if (is.null(group_col) || identical(group_col, TRUE) || !nzchar(group_col)) {
  if ("tcell_subtype" %in% colnames(obj@meta.data)) {
    group_col <- "tcell_subtype"
  } else if ("Tcell_subtype" %in% colnames(obj@meta.data)) {
    group_col <- "Tcell_subtype"
  } else {
    stop("对象 meta.data 缺少 tcell_subtype / Tcell_subtype 列，请用 --group_col 指定。")
  }
}

obj <- rename_tcell_subtype(obj, group_col)
level_vec <- levels(obj@meta.data[[group_col]])
pal <- make_original_dimplot_palette(level_vec)

dim_file_png <- file.path(plot_dir, "subset_obj_Tcell_DimPlot_by_cellanno_nature_palette.png")
dim_file_pdf <- file.path(plot_dir, "subset_obj_Tcell_DimPlot_by_cellanno_nature_palette_6.pdf")
dot_file_png <- file.path(plot_dir, "subset_obj_DotPlot_tcell_subtype_markers.png")
dot_file_pdf <- file.path(plot_dir, "subset_obj_DotPlot_tcell_subtype_markers7.pdf")

p_dim <- DimPlot(
  obj,
  reduction = "umap",
  group.by = group_col,
  cols = pal,
  label = TRUE,
  repel = TRUE
) +
  ggtitle(group_col) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    text = element_text(family = "Helvetica")
  )

ggsave(
  dim_file_pdf,
  plot = p_dim,
  device = cairo_pdf, 
  width = as.numeric(opt$width_dim),
  height = as.numeric(opt$height_dim),
  units = "in",
  dpi = 600,
  bg = "white"
)

marker_genes <- c(
  "CD3D", "CD3E", "CD3G",
  "CD4", "LEF1", "TCF7",
  "CD8A", "CD8B", "NKG7",
  "LAG3", "TIGIT",
  "ZBTB16", "RIPK2",
  "IL17A", "IL17F",
  "FOXP3", "SELL", "CD28", "FAS",
  "TRDV1", "TRGV3", "TRGV4", "TRGV8"
)

# p_dot <- DotPlot(obj, features = marker_genes, cols = c("#5E92CD","#E71F19")) +
#   RotatedAxis() +
#   ggtitle(sprintf("DotPlot of marker genes across %s", group_col)) +
#   theme(
#     plot.title = element_text(face = "bold", hjust = 0.5),
#     text = element_text(family = "Helvetica"),
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)
#   )
library(scales)
p_dot <- DotPlot(obj, features = marker_genes) +
  scale_color_gradientn(
    # values 定义颜色在 scale 中的相对位置（0到1之间）
    # 我们把 -2, 0, 2 映射到 0, 0.5, 1
    # 如果想让白色更窄，可以将转折点设置得离 0.5 更近
    values = rescale(c(-2, -0.1, 0, 0.1, 2)), 
    # 对应的颜色也要补齐：蓝 -> 浅蓝 -> 白 -> 浅红 -> 红
    colors = c("#5E92CD", "#ADC9E6", "white", "#F39894", "#E71F19"),
    limits = c(-2, 2),
    oob = squish
  ) +
  # scale_color_gradient2(
  #   low = "#5E92CD",
  #   mid = "#e5c5f2ff",
  #   high = "#E71F19",
  #   midpoint =0 
  # ) +
  RotatedAxis() +
  ggtitle(sprintf("DotPlot of marker genes across %s", group_col)) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)
  )



ggsave(
  dot_file_pdf,
  plot = p_dot,
  width = as.numeric(opt$width_dot),
  height = as.numeric(opt$height_dot),
  units = "in",
  dpi = 600,
  bg = "white"
)

message("输入对象：", input_path)
message("使用分组列：", group_col)
message("DimPlot 输出：", dim_file_png)
message("DotPlot 输出：", dot_file_png)
