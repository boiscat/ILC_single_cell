#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

options(stringsAsFactors = FALSE)

samples_to_keep <- c(
  "SKN8104896","SKN8104902","SKN8105194","SKN8105200",
  "SKN8090576","SKN8090578","SKN8090580","SKN8090582",
  "SKN8090588","SKN8090590","SKN8090592","SKN8090594",
  "SKN8090600","SKN8090602","SKN8090604","SKN8090606"
)

input_rdata <- "pso_seurat_obj.RData"
object_name <- "pso_seurat_obj"
output_dir <- "outputs/lesion_nonlesion_20251108"
subset_rdata <- file.path(output_dir, "pso_seurat_selected_samples_20251108.RData")
subset_rds <- file.path(output_dir, "pso_seurat_selected_samples_20251108.rds")
subset_object_name <- "pso_seurat_selected_samples"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

message("加载Seurat对象...")
load(input_rdata)
if (!exists(object_name, inherits = FALSE)) {
  stop(sprintf("对象 %s 在 %s 中不存在", object_name, input_rdata))
}
seurat_obj <- get(object_name)

if (!"sample_id" %in% colnames(seurat_obj@meta.data)) {
  stop("meta.data 中缺少 sample_id 字段，无法按样本筛选")
}

cells_to_keep <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data$sample_id %in% samples_to_keep]
if (length(cells_to_keep) == 0) {
  stop("指定的样本在对象中均无细胞，无法创建子集")
}

missing_samples <- setdiff(samples_to_keep, unique(seurat_obj@meta.data$sample_id))
message(sprintf("本次筛选共找到 %d 个细胞。", length(cells_to_keep)))
if (length(missing_samples) > 0) {
  warning(sprintf("以下样本在对象中不存在：%s", paste(missing_samples, collapse = ", ")))
}

subset_obj <- subset(seurat_obj, cells = cells_to_keep)
subset_obj@misc$subset_info <- list(
  samples_requested = samples_to_keep,
  samples_missing = missing_samples,
  total_cells = ncol(subset_obj),
  created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

assign(subset_object_name, subset_obj)

message("保存 RData ...")
save(list = subset_object_name, file = subset_rdata)
message("保存 RDS ...")
saveRDS(subset_obj, file = subset_rds)

message(sprintf("完成：%s 和 %s", subset_rdata, subset_rds))
