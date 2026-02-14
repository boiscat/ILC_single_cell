suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(FNN)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop("参数缺少取值：", flag)
  args[[hit + 1]]
}

in_rds <- get_arg(
  "--in_rds",
  "/lustre/home/yangyb/youms/data/outputs/healthy_pso_merge_20251110/healthy_pso_harmony_integrated_de_20251110.rds"
)
out_dir <- get_arg("--out_dir", "/lustre/home/yangyb/youms/data/figures/fig1")
k <- as.integer(get_arg("--k", "30"))
mono_thresh <- as.numeric(get_arg("--mono_thresh", "0.6"))

if (!file.exists(in_rds)) stop("输入RDS不存在：", in_rds)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_info <- file.path(out_dir, "figure1_clean_run_info_20260106.txt")
session_info <- file.path(out_dir, "figure1_clean_sessionInfo_20260106.txt")
writeLines(
  c(
    paste0("in_rds\t", in_rds),
    paste0("in_rds_md5\t", unname(tools::md5sum(in_rds))),
    paste0("out_dir\t", out_dir),
    paste0("k\t", k),
    paste0("mono_thresh\t", mono_thresh),
    paste0("run_time\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  ),
  run_info
)

set.seed(20260106)
obj <- readRDS(in_rds)

required_cols <- c("cell_anno", "seurat_clusters", "nCount_RNA", "nFeature_RNA", "percent.mito", "percent.HB")
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) stop("meta.data 缺少必要列：", paste(missing_cols, collapse = ", "))
if (!("umap" %in% Reductions(obj))) stop("对象里找不到 umap reduction（Reductions(obj) 不包含 'umap'）")

# 复现 figure1_DimPlot_cell_anno_7.pdf 的注释修正逻辑（figures/fig1/04-Fig1DimPlotCellAnnoRename-20251216.R）
ann <- obj@meta.data$cell_anno
clust <- obj@meta.data$seurat_clusters
if (is.factor(ann)) ann <- as.character(ann)
if (is.factor(clust)) clust <- as.character(clust)
ann[clust == "5" | clust == 5] <- "ILG"
ann[clust == "6" | clust == 6] <- "NK"
ann[ann %in% c("T/NK", "TNK")] <- "T"
ann[ann == "ILG"] <- "ILC"
ann[ann == "Myeloid"] <- "Monocyte/Macrophage"
desired_order <- c("FB", "ILC", "KC", "Mast cells", "Monocyte/Macrophage", "NK", "T")
other_levels <- setdiff(sort(unique(ann)), desired_order)
obj@meta.data$cell_anno_fixed <- factor(ann, levels = c(desired_order, other_levels))

# 保持与Fig1一致的细胞类型过滤
obj <- subset(obj, subset = !(cell_anno_fixed %in% c("VE", "LE", "Melanocytes", "Pericyte")))

meta <- obj@meta.data
umap <- Embeddings(obj, "umap")
labels <- as.character(meta$cell_anno_fixed)
mono_label <- "Monocyte/Macrophage"

knn <- FNN::get.knn(umap, k = k)
idx <- knn$nn.index
neighbor_labels <- matrix(labels[idx], nrow = nrow(idx), ncol = ncol(idx))
neighbor_mono_frac <- rowMeans(neighbor_labels == mono_label)
neighbor_same_frac <- rowMeans(neighbor_labels == matrix(labels, nrow = nrow(idx), ncol = ncol(idx)))

meta$umap_neighbor_mono_frac <- neighbor_mono_frac
meta$umap_neighbor_same_frac <- neighbor_same_frac
meta$embedded_in_mono <- (labels != mono_label) & (neighbor_mono_frac >= mono_thresh)
meta$mono_on_boundary <- (labels == mono_label) & (neighbor_same_frac < mono_thresh)
meta$low_quality_like <- (meta$percent.mito >= 0.10) | (meta$nFeature_RNA <= 500)
obj@meta.data <- meta

total_cells <- ncol(obj)
rm_mode1 <- meta$embedded_in_mono
rm_mode2 <- meta$embedded_in_mono | meta$mono_on_boundary

summarise_rm <- function(rm_flag, tag) {
  rm_n <- sum(rm_flag, na.rm = TRUE)
  keep_n <- total_cells - rm_n
  rm_pct <- rm_n / total_cells
  out_tsv <- file.path(out_dir, paste0("figure1_removed_cells_", tag, "_20260106.tsv"))
  rm_df <- meta %>%
    mutate(cell = rownames(meta)) %>%
    filter(rm_flag) %>%
    select(
      cell,
      cell_anno_fixed,
      embedded_in_mono,
      mono_on_boundary,
      low_quality_like,
      nCount_RNA,
      nFeature_RNA,
      percent.mito,
      percent.HB,
      umap_neighbor_mono_frac,
      umap_neighbor_same_frac
    ) %>%
    arrange(desc(embedded_in_mono), desc(mono_on_boundary), desc(percent.mito))
  write.table(rm_df, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  list(rm_n = rm_n, keep_n = keep_n, rm_pct = rm_pct, out_tsv = out_tsv)
}

res1 <- summarise_rm(rm_mode1, "embedded_in_mono")
res2 <- summarise_rm(rm_mode2, "embedded_plus_boundary")

writeLines(
  c(
    paste0("total_cells_after_Fig1_subset\t", total_cells),
    paste0("remove_embedded_in_mono_n\t", res1$rm_n),
    paste0("remove_embedded_in_mono_pct\t", sprintf("%.4f", res1$rm_pct)),
    paste0("remove_embedded_plus_boundary_n\t", res2$rm_n),
    paste0("remove_embedded_plus_boundary_pct\t", sprintf("%.4f", res2$rm_pct))
  ),
  file.path(out_dir, "figure1_clean_summary_20260106.tsv")
)

anno_colors <- c(
  "FB" = "#F0904D",
  "ILC" = "#EC7576",
  "KC" = "#14A0BD",
  "Mast cells" = "#6DA443",
  "Monocyte/Macrophage" = "#C6A3CE",
  "NK" = "#F7C8C8",
  "T" = "#5E92CD"
)

pdf_device <- if (capabilities("cairo")) cairo_pdf else pdf

plot_clean <- function(obj_in, pdf_path, title) {
  pdf_device(pdf_path, width = 6.5, height = 5.0)
  print(
    DimPlot(
      obj_in,
      group.by = "cell_anno_fixed",
      cols = anno_colors,
      label = TRUE,
      pt.size = 0.05,
      label.size = 5
    ) +
      ggtitle(title) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 10)
      )
  )
  dev.off()
}

obj_clean1 <- subset(obj, cells = rownames(meta)[!rm_mode1])
obj_clean2 <- subset(obj, cells = rownames(meta)[!rm_mode2])

out_pdf1 <- file.path(out_dir, "figure1_DimPlot_cell_anno_7.clean_embedded_20260106.pdf")
out_pdf2 <- file.path(out_dir, "figure1_DimPlot_cell_anno_7.clean_embedded_boundary_20260106.pdf")

plot_clean(
  obj_clean1,
  out_pdf1,
  sprintf("Cell annotation (clean: remove embedded-in-mono, %.2f%% removed)", res1$rm_pct * 100)
)
plot_clean(
  obj_clean2,
  out_pdf2,
  sprintf("Cell annotation (clean: remove embedded+boundary, %.2f%% removed)", res2$rm_pct * 100)
)

writeLines(capture.output(sessionInfo()), session_info)

message("完成：Fig1 清理版 DimPlot 已输出")
message("- 总细胞数(按Fig1过滤后)：", total_cells)
message("- 方案1：去掉 embedded_in_mono：", res1$rm_n, " (", sprintf("%.2f%%", res1$rm_pct * 100), ")")
message("  PDF：", out_pdf1)
message("- 方案2：去掉 embedded_in_mono + mono_on_boundary：", res2$rm_n, " (", sprintf("%.2f%%", res2$rm_pct * 100), ")")
message("  PDF：", out_pdf2)
message("- 细胞名单：", res1$out_tsv, " / ", res2$out_tsv)
