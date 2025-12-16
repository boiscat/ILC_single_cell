#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Seurat))

# 加载 PSO 选定样本对象
load("outputs/lesion_nonlesion_20251108/pso_seurat_selected_samples_20251108.RData")
if (!exists("pso_seurat_selected_samples")) {
  stop("对象 pso_seurat_selected_samples 未在加载的文件中找到。")
}

# 加载 Healthy 6 个样本的子集对象
load("healthy_subset_6samples_20251108.RData")
if (!exists("subset_obj")) {
  stop("对象 subset_obj 未在加载的文件中找到。")
}

# 合并两个 Seurat 对象，并为不同来源加上标识前缀
merged_obj <- merge(
  x = pso_seurat_selected_samples,
  y = subset_obj,
  add.cell.ids = c("PSO", "HEALTHY"),
  project = "pso_healthy_subset"
)

cat("PSO 细胞数:", ncol(pso_seurat_selected_samples), "\n")
cat("Healthy 细胞数:", ncol(subset_obj), "\n")
cat("合并后细胞数:", ncol(merged_obj), "\n")

save(merged_obj, file = "pso_healthy_selected_20251108.RData")
