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
args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) > 0) args[[1]] else "20251120"
message("[INFO] 开始读取 full_obj 进行 CD73 (NT5E) 闸门扫描，tag=", tag)
obj_path <- file.path(project_root, "task1-lesion-fix-20251120", "lesion-fix-rerun-20251120", "outputs", "objects", "full_obj_lesionFixed_20251120.rds")
stopifnot(file.exists(obj_path))
obj <- readRDS(obj_path)
stopifnot("NT5E" %in% rownames(obj), "cell_anno" %in% colnames(obj@meta.data))
expr_tbl <- FetchData(obj, vars = c("NT5E", "cell_anno", "Site_Status"))
ilc_tbl <- expr_tbl[expr_tbl$cell_anno == "ILC", , drop = FALSE]
stopifnot(nrow(ilc_tbl) > 0)
site_levels <- c("Healthy", "non_lesion", "lesion")
ilc_tbl$Site_Status <- factor(ilc_tbl$Site_Status, levels = site_levels)
site_totals <- as.data.frame(table(ilc_tbl$Site_Status), stringsAsFactors = FALSE)
colnames(site_totals) <- c("Site_Status", "total_ilc")
thresholds <- data.frame(
  label = c(">0", ">=0.1", ">=0.25", ">=0.5"),
  value = c(0, 0.1, 0.25, 0.5),
  strict = c(TRUE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)
records <- list()
for (i in seq_len(nrow(thresholds))) {
  th <- thresholds[i, ]
  mask <- if (th$strict) ilc_tbl$NT5E > th$value else ilc_tbl$NT5E >= th$value
  ilc_tbl$cd73_positive <- mask
  counts <- aggregate(cd73_positive ~ Site_Status, data = ilc_tbl, FUN = sum)
  merged <- merge(site_totals, counts, by = "Site_Status", all.x = TRUE)
  merged$cd73_positive[is.na(merged$cd73_positive)] <- 0
  merged$threshold <- th$label
  merged$fraction <- ifelse(merged$total_ilc > 0, merged$cd73_positive / merged$total_ilc, NA)
  records[[i]] <- merged
}
res <- do.call(rbind, records)
res <- res[order(res$threshold, res$Site_Status), ]
all_totals <- data.frame(Site_Status = "All", total_ilc = sum(site_totals$total_ilc))
for (i in seq_len(nrow(thresholds))) {
  th <- thresholds[i, ]
  mask <- if (th$strict) ilc_tbl$NT5E > th$value else ilc_tbl$NT5E >= th$value
  pos_total <- sum(mask)
  frac <- ifelse(all_totals$total_ilc > 0, pos_total / all_totals$total_ilc, NA)
  res <- rbind(res, data.frame(
    Site_Status = "All",
    total_ilc = all_totals$total_ilc,
    cd73_positive = pos_total,
    threshold = thresholds$label[i],
    fraction = frac
  ))
}
out_dir <- file.path(project_root, "task6-cd73pos-ilc-interaction-20251120", "cd73pos-ilc-interaction-pipeline-20251120", "outputs", "tables")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
out_file <- file.path(out_dir, sprintf("CD73_ILC_threshold_counts_%s.tsv", tag))
write.table(res, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("[INFO] 输出写入: ", out_file)
summary_file <- file.path(project_root, "task6-cd73pos-ilc-interaction-20251120", "cd73pos-ilc-interaction-pipeline-20251120", "logs", sprintf("stage1_threshold_summary_%s.log", tag))
summary_lines <- c(
  paste0("total_ILC=", sum(site_totals$total_ilc)),
  paste0("thresholds=", paste(thresholds$label, collapse = ","))
)
writeLines(summary_lines, con = summary_file)
