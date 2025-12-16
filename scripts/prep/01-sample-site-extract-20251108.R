#!/usr/bin/env Rscript

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
metadata_output <- file.path(output_dir, "selected_samples_site_metadata.tsv")
summary_output <- file.path(output_dir, "selected_samples_site_summary.tsv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load(input_rdata)
if (!exists(object_name)) {
  stop(sprintf("对象 %s 在 %s 中不存在", object_name, input_rdata))
}

seurat_obj <- get(object_name)
meta_df <- seurat_obj@meta.data
meta_df$cell_barcode <- colnames(seurat_obj)

selected_meta <- meta_df[meta_df$sample_id %in% samples_to_keep, , drop = FALSE]
if (nrow(selected_meta) == 0) {
  stop("指定样本在Seurat对象中不存在任何细胞")
}

cols_of_interest <- c(
  "cell_barcode","sample_id","Site","Tissue","Status",
  "Location","Sex","Age","stage","Enrichment","full_clustering"
)
cols_available <- intersect(cols_of_interest, colnames(selected_meta))
selected_meta_out <- selected_meta[, cols_available, drop = FALSE]

write.table(
  selected_meta_out,
  file = metadata_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 汇总每个样本在lesion/non_lesion的细胞数量
counts_df <- aggregate(
  cell_barcode ~ sample_id + Site,
  data = selected_meta,
  FUN = length
)
names(counts_df)[names(counts_df) == "cell_barcode"] <- "cell_count"

site_levels <- c("lesion", "non_lesion")
full_grid <- expand.grid(
  sample_id = samples_to_keep,
  Site = site_levels,
  stringsAsFactors = FALSE
)
summary_df <- merge(full_grid, counts_df, by = c("sample_id", "Site"), all.x = TRUE)
summary_df$cell_count[is.na(summary_df$cell_count)] <- 0
summary_df$has_cells <- summary_df$cell_count > 0

write.table(
  summary_df[order(summary_df$sample_id, summary_df$Site), ],
  file = summary_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("完成：", metadata_output)
message("完成：", summary_output)
