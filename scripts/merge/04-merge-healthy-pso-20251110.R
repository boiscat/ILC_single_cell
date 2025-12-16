#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

options(stringsAsFactors = FALSE)

healthy_rdata <- "outputs/healthy_subset_20251110/healthy_seurat_selected_samples_20251110.RData"
healthy_object_name <- "healthy_seurat_selected_samples"
pso_rdata <- "outputs/lesion_nonlesion_20251108/pso_seurat_selected_samples_20251108.RData"
pso_object_name <- "pso_seurat_selected_samples"
output_dir <- "outputs/healthy_pso_merge_20251110"
merged_object_name <- "healthy_pso_merged_samples"
merged_rdata <- file.path(output_dir, "healthy_pso_merged_samples_20251110.RData")
merged_rds <- file.path(output_dir, "healthy_pso_merged_samples_20251110.rds")
sample_summary_tsv <- file.path(output_dir, "healthy_pso_sample_counts_20251110.tsv")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

message("加载 Healthy 子集...")
load(healthy_rdata)
if (!exists(healthy_object_name, inherits = FALSE)) {
  stop(sprintf("对象 %s 在 %s 中不存在", healthy_object_name, healthy_rdata))
}
healthy_obj <- get(healthy_object_name)

message("加载 PSO 子集...")
load(pso_rdata)
if (!exists(pso_object_name, inherits = FALSE)) {
  stop(sprintf("对象 %s 在 %s 中不存在", pso_object_name, pso_rdata))
}
pso_obj <- get(pso_object_name)

if (!"sample_id" %in% colnames(healthy_obj@meta.data)) {
  stop("healthy 对象缺少 sample_id 字段")
}
if (!"sample_id" %in% colnames(pso_obj@meta.data)) {
  stop("pso 对象缺少 sample_id 字段")
}

message("合并对象...")
merged_obj <- merge(
  x = pso_obj,
  y = healthy_obj,
  add.cell.ids = c("PSO", "HEALTHY"),
  project = "healthy_pso_merge_20251110"
)

sample_counts <- sort(table(merged_obj@meta.data$sample_id), decreasing = TRUE)
sample_summary <- data.frame(
  sample_id = names(sample_counts),
  cell_count = as.integer(sample_counts),
  stringsAsFactors = FALSE
)

merged_obj@misc$merge_info <- list(
  healthy_source = healthy_rdata,
  pso_source = pso_rdata,
  healthy_cells = ncol(healthy_obj),
  pso_cells = ncol(pso_obj),
  total_cells = ncol(merged_obj),
  total_samples = nrow(sample_summary),
  created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  sample_counts = sample_summary
)

message(sprintf("PSO 细胞数: %d", ncol(pso_obj)))
message(sprintf("Healthy 细胞数: %d", ncol(healthy_obj)))
message(sprintf("合并后细胞数: %d", ncol(merged_obj)))
message(sprintf("合并后样本数: %d", nrow(sample_summary)))

message("写出样本统计 TSV ...")
write.table(sample_summary, file = sample_summary_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

assign(merged_object_name, merged_obj)

message("保存 RData ...")
save(list = merged_object_name, file = merged_rdata)

message("保存 RDS ...")
saveRDS(merged_obj, file = merged_rds)

message("处理完成。")
