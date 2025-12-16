#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
})

# resolve paths
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

# parse arguments
args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) > 0) args[[1]] else "20251120"
threshold <- if (length(args) > 1) as.numeric(args[[2]]) else 0
if (!is.finite(threshold)) {
  stop("threshold must be numeric")
}
threshold_mode <- if (length(args) > 2) args[[3]] else "gt"
if (!threshold_mode %in% c("gt", "gte")) {
  stop("threshold mode must be 'gt' or 'gte'")
}
threshold_label <- if (threshold_mode == "gt") paste0(">", threshold) else paste0(">=", threshold)
message("[INFO] tag=", tag, " threshold=", threshold_label)

# input paths
full_path <- file.path(project_root, "task1-lesion-fix-20251120", "lesion-fix-rerun-20251120", "outputs", "objects", "full_obj_lesionFixed_20251120.rds")
tcell_path <- file.path(project_root, "task1-lesion-fix-20251120", "lesion-fix-rerun-20251120", "outputs", "objects", "tcell_obj_lesionFixed_20251120.rds")
stopifnot(file.exists(full_path), file.exists(tcell_path))

# load objects
full_obj <- readRDS(full_path)
tcell_obj <- readRDS(tcell_path)
if (!"NT5E" %in% rownames(full_obj)) {
  stop("NT5E not found in full_obj")
}
if (!"cell_anno" %in% colnames(full_obj@meta.data)) {
  stop("cell_anno missing in full_obj meta")
}
if (!"Site_Status" %in% colnames(full_obj@meta.data)) {
  stop("Site_Status missing in full_obj meta")
}
if (!"Site_Status" %in% colnames(tcell_obj@meta.data)) {
  stop("Site_Status missing in tcell_obj meta")
}

# fetch NT5E information
expr_tbl <- FetchData(full_obj, vars = c("NT5E", "cell_anno", "Site_Status"))
expr_tbl$cell_id <- rownames(expr_tbl)
ilc_tbl <- expr_tbl[expr_tbl$cell_anno == "ILC", , drop = FALSE]
stopifnot(nrow(ilc_tbl) > 0)
if (threshold_mode == "gt") {
  ilc_tbl$cd73_positive <- ilc_tbl$NT5E > threshold
} else {
  ilc_tbl$cd73_positive <- ilc_tbl$NT5E >= threshold
}
cd73_cells <- ilc_tbl$cell_id[ilc_tbl$cd73_positive]
if (length(cd73_cells) == 0) {
  stop("No CD73+ ILC identified with threshold ", threshold_label)
}

# prepare directories
outputs_objects <- file.path(pipeline_root, "outputs", "objects")
outputs_tables <- file.path(pipeline_root, "outputs", "tables")
if (!dir.exists(outputs_objects)) dir.create(outputs_objects, recursive = TRUE)
if (!dir.exists(outputs_tables)) dir.create(outputs_tables, recursive = TRUE)

# prepare tcell subtype mapping
tcell_meta <- tcell_obj@meta.data
tcell_meta$cell_id <- rownames(tcell_meta)
if (!"tcell_subtype" %in% colnames(tcell_meta)) {
  stop("tcell_subtype column not found in tcell_obj meta")
}
present_in_full <- intersect(tcell_meta$cell_id, colnames(full_obj))
if (length(present_in_full) == 0) {
  stop("No overlapping cells between tcell_obj and full_obj")
}

conditions <- c("Healthy", "non_lesion", "lesion")
summary_records <- list()
tcell_subtype_counts <- as.data.frame(table(tcell_meta$Site_Status, tcell_meta$tcell_subtype), stringsAsFactors = FALSE)
colnames(tcell_subtype_counts) <- c("Site_Status", "tcell_subtype", "cells")

# total ILC counts per condition for reference
site_totals <- as.data.frame(table(ilc_tbl$Site_Status), stringsAsFactors = FALSE)
colnames(site_totals) <- c("Site_Status", "total_ilc")

for (cond in conditions) {
  cond_ilc <- cd73_cells[ilc_tbl$Site_Status[match(cd73_cells, ilc_tbl$cell_id)] == cond]
  cond_ilc <- cond_ilc[!is.na(cond_ilc)]
  cond_tcells <- tcell_meta$cell_id[tcell_meta$Site_Status == cond]
  cond_tcells <- cond_tcells[cond_tcells %in% colnames(full_obj)]
  combined <- unique(c(cond_ilc, cond_tcells))
  combined <- combined[combined %in% colnames(full_obj)]
  if (length(combined) == 0) {
    warning("No cells found for condition ", cond, " - skipping")
    next
  }
  subset_obj <- subset(full_obj, cells = combined)
  subset_obj$interaction_condition <- cond
  subset_obj$cd73_threshold <- threshold_label
  subset_obj$interaction_group <- NA_character_
  subset_obj$tcell_subtype <- NA_character_
  matched <- intersect(cond_tcells, colnames(subset_obj))
  if (length(matched) > 0) {
    subset_obj@meta.data[matched, "tcell_subtype"] <- as.character(tcell_meta[matched, "tcell_subtype"])
  }
  is_ilc <- subset_obj@meta.data$cell_anno == "ILC"
  subset_obj@meta.data$interaction_group[is_ilc] <- "CD73pos_ILC"
  subset_obj@meta.data$interaction_group[!is_ilc] <- ifelse(
    !is.na(subset_obj@meta.data$tcell_subtype[!is_ilc]) & nzchar(subset_obj@meta.data$tcell_subtype[!is_ilc]),
    subset_obj@meta.data$tcell_subtype[!is_ilc],
    subset_obj@meta.data$cell_anno[!is_ilc]
  )
  out_file <- file.path(outputs_objects, sprintf("CD73_ILC_Tcell_%s_%s.rds", cond, tag))
  saveRDS(subset_obj, out_file)
  total_ilc_cond <- site_totals$total_ilc[site_totals$Site_Status == cond]
  total_ilc_cond <- ifelse(length(total_ilc_cond) == 0, NA_integer_, total_ilc_cond)
  summary_records[[cond]] <- data.frame(
    Site_Status = cond,
    total_cd73_ilc = length(cond_ilc),
    total_ilc = total_ilc_cond,
    total_tcells = length(cond_tcells),
    total_combined = length(combined),
    threshold = threshold_label,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_records)
summary_df$fraction_cd73_ilc <- ifelse(!is.na(summary_df$total_ilc) & summary_df$total_ilc > 0,
                                       summary_df$total_cd73_ilc / summary_df$total_ilc, NA_real_)
summary_path <- file.path(outputs_tables, sprintf("CD73_ILC_Tcell_conditionSummary_%s.tsv", tag))
write.table(summary_df, file = summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

tcell_counts_path <- file.path(outputs_tables, sprintf("Tcell_subtype_condition_counts_%s.tsv", tag))
write.table(tcell_subtype_counts, file = tcell_counts_path, sep = "\t", quote = FALSE, row.names = FALSE)

message("[INFO] 写出对象: ", paste(basename(list.files(outputs_objects, pattern = sprintf("%s.rds$", tag))), collapse = ", "))
message("[INFO] 条件摘要: ", summary_path)
message("[INFO] T cell subtype 统计: ", tcell_counts_path)
