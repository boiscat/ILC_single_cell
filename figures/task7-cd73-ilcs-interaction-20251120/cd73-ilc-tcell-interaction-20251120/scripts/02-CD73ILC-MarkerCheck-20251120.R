#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(patchwork)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "CD73阴性ILC对象RDS"),
  make_option(c("-o", "--output_dir"), type = "character", default = "outputs", help = "输出目录根路径"),
  make_option(c("--tag"), type = "character", default = format(Sys.Date(), "%Y%m%d"), help = "输出tag"),
  make_option(c("--feature_list"), type = "character",
              default = "NT5E,IL7R,IL1RL1,KLRB1,GATA3,IL7"),
  make_option(c("--group_by"), type = "character", default = "Site_Status",
              help = "绘图分组列（默认Site_Status）")
)

opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$input), file.exists(opt$input))

feature_vec <- unique(trimws(strsplit(opt$feature_list, ",")[[1]]))
feature_vec <- feature_vec[nchar(feature_vec) > 0]
stopifnot(length(feature_vec) > 0)

plots_dir <- file.path(opt$output_dir, "plots")
tables_dir <- file.path(opt$output_dir, "tables")
objects_dir <- file.path(opt$output_dir, "objects")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)

message("读取对象: ", opt$input)
obj <- readRDS(opt$input)
meta <- obj@meta.data
stopifnot(opt$group_by %in% colnames(meta))

ntu <- FetchData(obj, vars = unique(c("NT5E", opt$group_by)))
meta$NT5E_expr <- ntu$NT5E

summary_tbl <- meta %>%
  group_by(.data[[opt$group_by]]) %>%
  summarise(
    n_cells = dplyr::n(),
    mean_NT5E = mean(NT5E_expr),
    median_NT5E = median(NT5E_expr),
    pct_zero = mean(NT5E_expr == 0) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

summary_path <- file.path(tables_dir, sprintf("cd73neg_nt5e_stats_%s.csv", opt$tag))
write_csv(summary_tbl, summary_path)
message("写出 NT5E 统计: ", summary_path)

vln <- VlnPlot(
  obj,
  features = "NT5E",
  group.by = opt$group_by,
  pt.size = 0.05
) + ggtitle(sprintf("NT5E expression (%s)", opt$group_by)) +
  theme_minimal(base_size = 12)

marker_features <- intersect(feature_vec, rownames(obj))
if (length(marker_features) == 0) {
  warning("输入的feature_list在对象中未找到，跳过DotPlot")
  marker_plot <- ggplot() + ggtitle("No marker features available")
} else {
  marker_plot <- DotPlot(
    obj,
    features = marker_features,
    group.by = opt$group_by
  ) + RotatedAxis() +
    ggtitle(sprintf("CD73- ILC marker panel (%s)", opt$group_by)) +
    theme_minimal(base_size = 11)
}

combined <- vln / marker_plot + plot_layout(heights = c(1, 1.2))

plot_basename <- sprintf("cd73neg_qc_%s", opt$tag)
plot_pdf <- file.path(plots_dir, paste0(plot_basename, ".pdf"))
plot_png <- file.path(plots_dir, paste0(plot_basename, ".png"))
ggsave(plot_pdf, combined, width = 10, height = 10, units = "in")
ggsave(plot_png, combined, width = 10, height = 10, units = "in", dpi = 600)
message("图形输出: ", plot_pdf)

session_path <- file.path(objects_dir, sprintf("cd73neg_markercheck_sessionInfo_%s.txt", opt$tag))
writeLines(capture.output(sessionInfo()), con = session_path)
message("SessionInfo: ", session_path)
