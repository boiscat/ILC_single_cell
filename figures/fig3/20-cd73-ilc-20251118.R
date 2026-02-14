#!/usr/bin/env Rscript

# 10-cd73-ilc-20251118.R
# 目的：统计 full_obj_saved_1117.rds 中 CD73+ ILC 的数量与比例，并输出图表

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

obj_path <- "figures/fig3/full_obj_saved_1117.rds"
output_dir <- "outputs/cd73_ilc_stats_20251118"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("读取对象：", obj_path)
obj <- readRDS(obj_path)

meta_df <- obj@meta.data %>%
  rownames_to_column(var = "cell_barcode")

ilc_meta <- meta_df %>%
  filter(cell_anno == "ILC")

if (nrow(ilc_meta) == 0) {
  stop("对象中不存在 ILC 细胞，无法继续")
}

rna_data <- SeuratObject::GetAssayData(obj[["RNA"]], layer = "data")

if (!"NT5E" %in% rownames(rna_data)) {
  stop("RNA 数据中不存在 NT5E (CD73) 基因")
}

cd73_expr <- as.numeric(rna_data["NT5E", ilc_meta$cell_barcode, drop = TRUE])
ilc_meta <- ilc_meta %>%
  mutate(
    NT5E_expr = cd73_expr,
    CD73_positive = NT5E_expr > 0
  )

cd73_ilc <- ilc_meta %>% filter(CD73_positive)

if (nrow(cd73_ilc) == 0) {
  stop("没有检测到 CD73+ ILC 细胞，请检查阈值或数据")
}

status_summary <- cd73_ilc %>%
  group_by(Status) %>%
  summarise(cd73_ilc_n = n(), .groups = "drop") %>%
  mutate(fraction = cd73_ilc_n / sum(cd73_ilc_n)) %>%
  arrange(desc(cd73_ilc_n))

site_summary <- cd73_ilc %>%
  group_by(Site_Status) %>%
  summarise(cd73_ilc_n = n(), .groups = "drop") %>%
  mutate(fraction = cd73_ilc_n / sum(cd73_ilc_n)) %>%
  arrange(desc(cd73_ilc_n))

write.table(status_summary,
            file = file.path(output_dir, "cd73_ilc_by_status.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(site_summary,
            file = file.path(output_dir, "cd73_ilc_by_site_status.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_fraction <- function(df, x_col, x_label, title_text, y_label) {
  x_sym <- rlang::sym(x_col)
  ggplot(df, aes(x = !!x_sym, y = fraction, fill = !!x_sym)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = paste0(round(fraction * 100, 1), "% (n=", cd73_ilc_n, ")")),
              vjust = -0.3, size = 3) +
    scale_y_continuous(labels = function(x) paste0(x * 100, "%"),
                       expand = expansion(mult = c(0, 0.15))) +
    labs(x = x_label, y = y_label, title = title_text, fill = x_label) +
    theme_minimal(base_size = 12)
}

status_plot <- plot_fraction(
  status_summary,
  x_col = "Status",
  x_label = "Status",
  title_text = "CD73+ ILC 在 Healthy vs Psoriasis 的比例",
  y_label = "百分比"
)

site_plot <- plot_fraction(
  site_summary,
  x_col = "Site_Status",
  x_label = "Site_Status",
  title_text = "CD73+ ILC 在皮损 vs 非皮损的比例",
  y_label = "百分比"
)

ggsave(filename = file.path(output_dir, "cd73_ilc_status_ratio.png"),
       plot = status_plot, width = 5, height = 4, dpi = 300)

ggsave(filename = file.path(output_dir, "cd73_ilc_site_ratio.png"),
       plot = site_plot, width = 5, height = 4, dpi = 300)

message("CD73+ ILC 统计完成，总数量：", nrow(cd73_ilc))
message("结果输出目录：", output_dir)
