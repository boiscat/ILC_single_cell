#!/usr/bin/env Rscript

setwd("/lustre/home/yangyb/youmz/data")

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(RColorBrewer)
})

set.seed(20251110)
options(stringsAsFactors = FALSE)

input_rds <- "outputs/healthy_pso_merge_20251110/healthy_pso_merged_samples_20251110.rds"
output_dir <- "outputs/healthy_pso_merge_20251110"
figure_dir <- "figures"
filtered_rds <- file.path(output_dir, "healthy_pso_harmony_integrated_20251110.rds")
qc_stats_tsv <- file.path(output_dir, "healthy_pso_qc_summary_20251110.tsv")
sample_stats_tsv <- file.path(output_dir, "healthy_pso_major_annotations_20251110.tsv")
figure_pdf <- file.path(figure_dir, "nature_healthy_pso_harmony_20251110.pdf")
figure_png <- file.path(figure_dir, "nature_healthy_pso_harmony_20251110.png")

if (!file.exists(input_rds)) {
  stop(sprintf("输入文件不存在：%s", input_rds))
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

message("读取合并对象 ...")
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"


message("执行质控筛选 ...")
qc_filters <- with(
  obj@meta.data,
  nFeature_RNA >= 500 &
    nFeature_RNA <= 6000 &
    nCount_RNA >= 1000 &
    nCount_RNA <= 60000 &
    percent.mito <= 0.2
)
cells_keep <- names(which(qc_filters))
if (length(cells_keep) == 0) {
  stop("质控条件过严，未保留任何细胞。")
}
filtered_obj <- subset(obj, cells = cells_keep)

qc_summary <- tibble::tibble(
  metric = c(
    "cells_before", "cells_after", "samples_before", "samples_after",
    "median_nFeature_after", "median_nCount_after", "median_pct_mito_after"
  ),
  value = c(
    ncol(obj), ncol(filtered_obj),
    length(unique(obj$sample_id)), length(unique(filtered_obj$sample_id)),
    median(filtered_obj$nFeature_RNA), median(filtered_obj$nCount_RNA),
    median(filtered_obj$percent.mito)
  )
)
readr::write_tsv(qc_summary, qc_stats_tsv)

message("标准Seurat预处理 ...")
filtered_obj <- NormalizeData(filtered_obj, verbose = FALSE)
filtered_obj <- FindVariableFeatures(filtered_obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
filtered_obj <- ScaleData(filtered_obj, vars.to.regress = c("percent.mito", "nCount_RNA"), verbose = FALSE)
filtered_obj <- RunPCA(filtered_obj, npcs = 50, verbose = FALSE)

message("运行Harmony批次校正（按 sample_id） ...")
filtered_obj <- RunHarmony(
  object = filtered_obj,
  group.by.vars = "sample_id",
  reduction = "pca",
  dims.use = 1:50,
  assay.use = "RNA",
  verbose = FALSE
)

message("降维与聚类 ...")
filtered_obj <- RunUMAP(filtered_obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
filtered_obj <- FindNeighbors(filtered_obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
filtered_obj <- FindClusters(filtered_obj, resolution = 0.4, verbose = FALSE)

message("计算大类注释分数 ...")
marker_panels <- list(
  T_cells = c("CD3D", "CD3E", "CD2"),
  NK_cells = c("NKG7", "GNLY", "PRF1"),
  B_cells = c("MS4A1", "CD79A", "CD74"),
  Plasma = c("SDC1", "MZB1", "XBP1"),
  Myeloid = c("LYZ", "LST1", "S100A8"),
  Dendritic = c("ITGAX", "IRF8", "CLEC10A"),
  Keratinocyte = c("KRT1", "KRT14", "KRT10"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN"),
  Endothelial = c("PECAM1", "VWF", "CLDN5"),
  Melanocyte = c("PMEL", "TYR", "MLANA")
)

marker_panels <- lapply(marker_panels, function(genes) {
  intersect(genes, rownames(filtered_obj))
})
non_empty_sets <- lengths(marker_panels) > 0
marker_panels <- marker_panels[non_empty_sets]

if (length(marker_panels) == 0) {
  stop("所有 marker set 均为空，请检查基因命名。")
}

filtered_obj <- AddModuleScore(filtered_obj, features = marker_panels, name = "typeScore", search = FALSE, assay = "RNA", seed = 20251110)
score_cols <- grep("^typeScore", colnames(filtered_obj@meta.data), value = TRUE)
panel_names <- names(marker_panels)
if (length(score_cols) != length(panel_names)) {
  stop("模块分数列与marker集合数量不匹配。")
}
colnames(filtered_obj@meta.data)[match(score_cols, colnames(filtered_obj@meta.data))] <- paste0("typeScore_", panel_names)
score_cols <- paste0("typeScore_", panel_names)

filtered_obj$major_annotation <- apply(filtered_obj@meta.data[, score_cols, drop = FALSE], 1, function(x) {
  panel_names[which.max(x)]
})

message("生成Nature风格图 ...")
umap_df <- as.data.frame(Embeddings(filtered_obj, "umap"))
umap_df$major_annotation <- filtered_obj$major_annotation
umap_df$sample_id <- filtered_obj$sample_id
palette_nature <- c(
  "T_cells" = "#1f77b4",
  "NK_cells" = "#17becf",
  "B_cells" = "#9467bd",
  "Plasma" = "#c5b0d5",
  "Myeloid" = "#ff7f0e",
  "Dendritic" = "#ffbb78",
  "Keratinocyte" = "#2ca02c",
  "Fibroblast" = "#8c564b",
  "Endothelial" = "#e377c2",
  "Melanocyte" = "#7f7f7f"
)
missing_colors <- setdiff(unique(umap_df$major_annotation), names(palette_nature))
if (length(missing_colors) > 0) {
  extra_cols <- RColorBrewer::brewer.pal(max(length(missing_colors), 3), "Set2")
  palette_nature <- c(palette_nature, setNames(extra_cols[seq_along(missing_colors)], missing_colors))
}

nature_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3)
  )

umap_plot <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = major_annotation)) +
  geom_point(size = 0.3, alpha = 0.6, stroke = 0) +
  scale_color_manual(values = palette_nature) +
  nature_theme +
  ggtitle("Harmony UMAP | Major Cell Classes")

composition_df <- filtered_obj@meta.data %>%
  count(sample_id, major_annotation, name = "cell_count") %>%
  group_by(sample_id) %>%
  mutate(fraction = cell_count / sum(cell_count)) %>%
  ungroup()

composition_plot <- ggplot(composition_df, aes(x = sample_id, y = fraction, fill = major_annotation)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = palette_nature) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.title.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    legend.position = "none"
  ) +
  ylab("Fraction per sample")

final_plot <- umap_plot + composition_plot + plot_layout(widths = c(2, 1))
ggsave(figure_pdf, final_plot, width = 10, height = 4.2, device = cairo_pdf)
ggsave(figure_png, final_plot, width = 10, height = 4.2, dpi = 400)

message("写出对象与统计 ...")
readr::write_tsv(
  composition_df %>% arrange(desc(cell_count)),
  sample_stats_tsv
)

filtered_obj@misc$analysis_info <- list(
  script = basename("05-harmony-healthy-pso-20251110.R"),
  created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  qc_thresholds = list(
    min_genes = 500,
    max_genes = 6000,
    min_counts = 1000,
    max_counts = 60000,
    max_pct_mito = 0.2
  )
)

saveRDS(filtered_obj, filtered_rds)

message("完成 Harmony 批次校正与可视化。")



