#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("用法: cd73-interaction-test1-20251120.R <input_rds> <output_rds>")
}

input_rds <- args[[1]]
output_rds <- args[[2]]
stopifnot(file.exists(input_rds))

obj <- readRDS(input_rds)
set.seed(2025)
subset_cells <- sample(colnames(obj), min(500, ncol(obj)))
mini_obj <- subset(obj, cells = subset_cells)
saveRDS(mini_obj, output_rds)
cat("测试对象已保存：", output_rds, "\n")
