#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
})

set.seed(20251110)
options(stringsAsFactors = FALSE)

write_tsv <- function(df, path) {
  utils::write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
}

input_rds <- "data/outputs/healthy_pso_merge_20251110/Dhealthy_pso_merged_samples_20251110.rds"
output_dir <- "outputs/healthy_pso_merge_20251110"
figure_dir <- "figures"
de_dir <- file.path(output_dir, "de_results")
filtered_rds <- file.path(output_dir, "healthy_pso_harmony_integrated_de_20251110.rds")
qc_stats_tsv <- file.path(output_dir, "healthy_pso_qc_summary_20251110.tsv")
sample_stats_tsv <- file.path(output_dir, "healthy_pso_major_annotations_20251110.tsv")
cluster_marker_tsv <- file.path(de_dir, "healthy_pso_cluster_markers_20251110.tsv")
major_marker_tsv <- file.path(de_dir, "healthy_pso_major_markers_20251110.tsv")
figure_pdf <- file.path(figure_dir, "nature_healthy_pso_harmony_20251110.pdf")
figure_png <- file.path(figure_dir, "nature_healthy_pso_harmony_20251110.png")

qc_params <- list(
  min_genes = 400,
  max_genes = 8000,
  min_counts = 800,
  max_counts = 80000,
  max_pct_mito = 0.25
)

pca_dims <- 50
neighbors_dims <- 1:30
harmony_dims <- 1:50
cluster_resolution <- 0.4

if (!file.exists(input_rds)) {
  stop(sprintf("输入文件不存在：%s", input_rds))
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(de_dir)) {
  dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
}
getwd()
message("读取合并对象 ...")
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"

required_cols <- c("sample_id", "percent.mito")
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
  stop(sprintf("对象缺少必要元数据列：%s", paste(missing_cols, collapse = ", ")))
}

message("执行宽阈值质控 ...")
qc_filters <- with(
  obj@meta.data,
  nFeature_RNA >= qc_params$min_genes &
    nFeature_RNA <= qc_params$max_genes &
    nCount_RNA >= qc_params$min_counts &
    nCount_RNA <= qc_params$max_counts &
    percent.mito <= qc_params$max_pct_mito
)
names(qc_filters) <- colnames(obj)
qc_filters[is.na(qc_filters)] <- FALSE
cells_keep <- names(qc_filters)[qc_filters]
if (length(cells_keep) == 0) {
  stop("质控条件过宽或元数据异常，未保留任何细胞。")
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
write_tsv(qc_summary, qc_stats_tsv)

message("标准 Seurat 预处理 + PCA ...")
filtered_obj <- NormalizeData(filtered_obj, verbose = FALSE)
filtered_obj <- FindVariableFeatures(filtered_obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
filtered_obj <- ScaleData(filtered_obj, vars.to.regress = c("percent.mito", "nCount_RNA"), verbose = FALSE)
filtered_obj <- RunPCA(filtered_obj, npcs = pca_dims, verbose = FALSE)

message("运行 Harmony 批次校正 (sample_id) ...")
filtered_obj <- RunHarmony(
  object = filtered_obj,
  group.by.vars = "sample_id",
  reduction.use = "pca",
  dims.use = harmony_dims,
  assay.use = "RNA",
  verbose = FALSE
)

message("降维、建图与聚类 ...")
filtered_obj <- RunUMAP(filtered_obj, reduction = "harmony", dims = neighbors_dims, verbose = FALSE)
filtered_obj <- FindNeighbors(filtered_obj, reduction = "harmony", dims = neighbors_dims, verbose = FALSE)
filtered_obj <- FindClusters(filtered_obj, resolution = cluster_resolution, verbose = FALSE)

message("大类 marker 注释 ...")
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

marker_panels <- lapply(marker_panels, function(genes) intersect(genes, rownames(filtered_obj)))
marker_panels <- marker_panels[lengths(marker_panels) > 0]
if (length(marker_panels) == 0) {
  stop("所有 marker set 均为空，请检查基因命名。")
}

filtered_obj <- AddModuleScore(
  filtered_obj,
  features = marker_panels,
  name = "typeScore",
  search = FALSE,
  assay = "RNA",
  seed = 20251110
)
score_cols <- grep("^typeScore", colnames(filtered_obj@meta.data), value = TRUE)
stopifnot(length(score_cols) == length(marker_panels))
panel_names <- names(marker_panels)
colnames(filtered_obj@meta.data)[match(score_cols, colnames(filtered_obj@meta.data))] <- paste0("typeScore_", panel_names)
score_cols <- paste0("typeScore_", panel_names)

filtered_obj$major_annotation <- apply(filtered_obj@meta.data[, score_cols, drop = FALSE], 1, function(x) {
  panel_names[which.max(x)]
})

message("构建 Nature 风格 UMAP & 组成图 ...")
umap_df <- as.data.frame(Embeddings(filtered_obj, "umap"))
if (ncol(umap_df) < 2) {
  stop("UMAP 维度少于2，无法绘制散点图。")
}
colnames(umap_df)[seq_len(2)] <- paste0("UMAP_", seq_len(2))
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
  extra_cols <- RColorBrewer::brewer.pal(max(3, length(missing_colors)), "Set2")
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

message("差异表达分析 (cluster & major) ...")
cluster_markers <- FindAllMarkers(
  filtered_obj,
  assay = "RNA",
  features = VariableFeatures(filtered_obj),
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox",
  only.pos = FALSE,
  return.thresh = 0.01,
  max.cells.per.ident = 5000
)
cluster_markers <- cluster_markers %>% arrange(cluster, desc(avg_log2FC))
write_tsv(cluster_markers, cluster_marker_tsv)

orig_idents <- Idents(filtered_obj)
Idents(filtered_obj) <- filtered_obj$major_annotation
major_markers <- FindAllMarkers(
  filtered_obj,
  assay = "RNA",
  features = VariableFeatures(filtered_obj),
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox",
  only.pos = FALSE,
  return.thresh = 0.01,
  max.cells.per.ident = 5000
)
major_markers <- major_markers %>% arrange(cluster, desc(avg_log2FC))
write_tsv(major_markers, major_marker_tsv)
Idents(filtered_obj) <- orig_idents

message("写出对象与统计 ...")
write_tsv(
  composition_df %>% arrange(sample_id, desc(cell_count)),
  sample_stats_tsv
)

filtered_obj@misc$analysis_info <- list(
  script = basename("06-harmony-de-20251110.R"),
  created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  environment = list(r_module = "r-pasta-20250929"),
  qc_thresholds = qc_params,
  dims = list(pca = pca_dims, harmony = harmony_dims, clustering_dims = neighbors_dims),
  cluster_resolution = cluster_resolution,
  differential_expression = list(
    test = "wilcox",
    logfc_threshold = 0.25,
    min_pct = 0.1,
    return_thresh = 0.01,
    max_cells_per_ident = 5000
  )
)

saveRDS(filtered_obj, filtered_rds)

message("完成 Harmony 批次校正 + 差异表达 + Nature 图。")


