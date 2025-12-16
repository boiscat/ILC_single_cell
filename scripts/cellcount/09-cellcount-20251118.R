#!/usr/bin/env Rscript

# 09-cellcount-20251118.R
# 功能：统计 Seurat 对象中每个样本的细胞量，并映射 Healthy/Psoriasis 及 lesion/nonlesion 信息

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})

obj_path <- "figures/fig3/full_obj_saved_1117.rds"
output_dir <- "outputs/sample_condition_counts_20251118"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("读取对象：", obj_path)
obj <- readRDS(obj_path)
meta_df <- obj@meta.data %>%
  rownames_to_column(var = "cell_barcode")

summarise_field <- function(values) {
  uniq_vals <- unique(values)
  paste(uniq_vals, collapse = ";")
}

sample_summary <- meta_df %>%
  group_by(sample_id) %>%
  summarise(
    Status = summarise_field(Status),
    Site_Status = summarise_field(Site_Status),
    Site = summarise_field(Site),
    Tissue = summarise_field(Tissue),
    cell_count = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(cell_count))

condition_summary <- sample_summary %>%
  group_by(Status, Site_Status) %>%
  summarise(
    sample_n = n(),
    cell_n = sum(cell_count),
    .groups = "drop"
  ) %>%
  arrange(desc(cell_n))

status_summary <- sample_summary %>%
  group_by(Status) %>%
  summarise(
    sample_n = n(),
    cell_n = sum(cell_count),
    .groups = "drop"
  ) %>%
  arrange(desc(cell_n))

write.table(sample_summary, file.path(output_dir, "sample_cell_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(condition_summary, file.path(output_dir, "status_sitestats_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(status_summary, file.path(output_dir, "status_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("结果保存到目录：", output_dir)
message("样本数量：", nrow(sample_summary), ", 总细胞数：", sum(sample_summary$cell_count))
