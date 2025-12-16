#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
})
args_full <- commandArgs(trailingOnly = FALSE)
file_flag <- "--file="
script_path <- NULL
for (arg in args_full) {
  if (startsWith(arg, file_flag)) {
    script_path <- normalizePath(sub(file_flag, "", arg))
  }
}
if (is.null(script_path)) {
  script_path <- normalizePath(getwd())
}
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
pipeline_root <- file.path(project_root, "task6-cd73pos-ilc-interaction-20251120", "cd73pos-ilc-interaction-pipeline-20251120")
obj_path <- file.path(pipeline_root, "inputs", "full_obj_lesionFixed_20251120.rds")
if (!file.exists(obj_path)) {
  obj_path <- file.path(project_root, "task1-lesion-fix-20251120", "lesion-fix-rerun-20251120", "outputs", "objects", "full_obj_lesionFixed_20251120.rds")
}
stopifnot(file.exists(obj_path))
obj <- readRDS(obj_path)
stopifnot("NT5E" %in% rownames(obj))
expr_tbl <- FetchData(obj, vars = c("NT5E", "cell_anno", "Site_Status"))
ilc_tbl <- expr_tbl[expr_tbl$cell_anno == "ILC", , drop = FALSE]
stopifnot(nrow(ilc_tbl) > 0)
stopifnot(all(ilc_tbl$NT5E >= 0))
set.seed(20251120)
sample_n <- min(200, nrow(ilc_tbl))
sample_rows <- sample(seq_len(nrow(ilc_tbl)), size = sample_n)
sample_tbl <- ilc_tbl[sample_rows, ]
candidate_thresholds <- c(0, 0.1, 0.25)
cat("# CD73ILC-test1-20251120\n")
cat("sample_size=", sample_n, " of ", nrow(ilc_tbl), "\n", sep = "")
cat("Site_Status breakdown (sample):\n")
print(table(sample_tbl$Site_Status))
for (th in candidate_thresholds) {
  frac <- mean(sample_tbl$NT5E > th)
  cat(sprintf("threshold > %.2f => %.2f%% positive\n", th, frac * 100))
}
stopifnot(mean(sample_tbl$NT5E > 0) <= 0.05)
cat("Test assertions passed.\n")
