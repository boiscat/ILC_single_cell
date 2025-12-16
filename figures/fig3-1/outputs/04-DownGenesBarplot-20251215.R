#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default) {
  hit <- which(args == flag)
  if (length(hit) == 1 && hit < length(args)) {
    return(args[[hit + 1]])
  }
  default
}

de_file <- get_arg("--de", "CD73_ILC_DE_results.csv")
down_genes_file <- get_arg("--down_genes", "CD73_ILC_DE_down_genes.txt")
out_png <- get_arg("--out", "CD73_ILC_DE_barplot_down_p0.05.png")

if (!file.exists(de_file)) stop("找不到差异结果文件: ", de_file)
if (!file.exists(down_genes_file)) stop("找不到下调基因文件: ", down_genes_file)

de_tbl <- read.csv(de_file, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("gene", "avg_log2FC")
missing_cols <- setdiff(required_cols, colnames(de_tbl))
if (length(missing_cols) > 0) {
  stop("差异结果缺少列: ", paste(missing_cols, collapse = ", "), "；文件=", de_file)
}

down_genes <- readLines(down_genes_file, warn = FALSE) |>
  trimws() |>
  (\(x) x[nzchar(x)])() |>
  unique()

if (length(down_genes) == 0) stop("下调基因文件为空: ", down_genes_file)

plot_data <- de_tbl |>
  filter(gene %in% down_genes) |>
  filter(avg_log2FC < 0) |>
  arrange(avg_log2FC) |>
  mutate(
    Group = "Down",
    gene = factor(gene, levels = unique(gene[order(avg_log2FC, decreasing = FALSE)]))
  )

if (nrow(plot_data) == 0) {
  stop("在差异表中没有匹配到任何下调基因；请检查gene命名是否一致。")
}

height_in <- max(6, 1 + 0.18 * nrow(plot_data))

p <- ggplot(plot_data, aes(x = gene, y = avg_log2FC, fill = Group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  coord_flip() +
  labs(
    title = sprintf("Human (down genes p<0.05, n=%d)", nrow(plot_data)),
    x = "gene",
    y = "Avg_log2FC",
    fill = "Group"
  ) +
  scale_fill_manual(values = c(Down = "#4DBBD5")) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
  )

ggsave(out_png, p, width = 4, height = height_in, dpi = 300)
message("绘图保存: ", out_png)
