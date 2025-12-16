#!/usr/bin/env Rscript

# 目标：为 fig1 的 Tissue 图统计 sample 级别与 patient 级别的数量（可复现）。
# 输入：
# - fig1_ILC_proportion_by_Sample_Tissue.csv（sample级别，一行一个sample）
# - patient_ILC_proportion_by_Tissue.csv（patient级别，一行一个 patient×Tissue）
# - meta_id_patient.tsv（sample_id -> patient 的映射，用于交叉核对）
#
# 兼容两种工作目录：
# 1) 直接在 fig1 目录运行：Rscript 03-CountSamplePatient-20251215.R
# 2) 在上一级 figures 目录运行：Rscript fig1/03-CountSamplePatient-20251215.R

pick_path <- function(candidates) {
  for (p in candidates) {
    if (file.exists(p)) return(p)
  }
  stop("找不到输入文件，尝试过：\n- ", paste(candidates, collapse = "\n- "))
}

sample_path <- pick_path(c("fig1_ILC_proportion_by_Sample_Tissue.csv",
                           file.path("fig1", "fig1_ILC_proportion_by_Sample_Tissue.csv")))
patient_path <- pick_path(c("patient_ILC_proportion_by_Tissue.csv",
                            file.path("fig1", "patient_ILC_proportion_by_Tissue.csv")))
meta_path <- pick_path(c("meta_id_patient.tsv",
                         file.path("fig1", "meta_id_patient.tsv")))

samples <- read.csv(sample_path, stringsAsFactors = FALSE, check.names = FALSE)
patients <- read.csv(patient_path, stringsAsFactors = FALSE, check.names = FALSE)
meta <- read.delim(meta_path, stringsAsFactors = FALSE, check.names = FALSE)

stopifnot("sample_id" %in% names(samples), "Tissue" %in% names(samples))
stopifnot("patient" %in% names(patients), "Tissue" %in% names(patients))
stopifnot("sample_id" %in% names(meta), "patient" %in% names(meta))

samples2 <- samples[!is.na(samples$sample_id) & !is.na(samples$Tissue), , drop = FALSE]
patients2 <- patients[!is.na(patients$patient) & !is.na(patients$Tissue), , drop = FALSE]
meta2 <- meta[!is.na(meta$sample_id) & !is.na(meta$patient), , drop = FALSE]

sample_n_total <- length(unique(samples2$sample_id))
patient_n_total <- length(unique(patients2$patient))

sample_n_by_tissue <- tapply(samples2$sample_id, samples2$Tissue, function(x) length(unique(x)))
patient_n_by_tissue <- tapply(patients2$patient, patients2$Tissue, function(x) length(unique(x)))

all_tissues <- sort(unique(c(names(sample_n_by_tissue), names(patient_n_by_tissue))))
out <- data.frame(
  Tissue = all_tissues,
  sample_n = as.integer(sample_n_by_tissue[all_tissues]),
  patient_n = as.integer(patient_n_by_tissue[all_tissues]),
  stringsAsFactors = FALSE
)

out_dir <- dirname(sample_path)
out_path <- file.path(out_dir, "fig1_ILC_Tissue_sample_patient_counts.tsv")
write.table(out, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)

sample_list_path <- file.path(out_dir, "fig1_ILC_Tissue_sample_list.tsv")
patient_list_path <- file.path(out_dir, "fig1_ILC_Tissue_patient_list.tsv")
map_path <- file.path(out_dir, "fig1_ILC_Tissue_sample_patient_map.tsv")

sample_list <- samples2[, c("Tissue", "sample_id"), drop = FALSE]
sample_list <- sample_list[order(sample_list$Tissue, sample_list$sample_id), , drop = FALSE]
sample_list <- unique(sample_list)
write.table(sample_list, file = sample_list_path, sep = "\t", quote = FALSE, row.names = FALSE)

patient_list <- patients2[, c("Tissue", "patient"), drop = FALSE]
patient_list <- patient_list[order(patient_list$Tissue, patient_list$patient), , drop = FALSE]
patient_list <- unique(patient_list)
write.table(patient_list, file = patient_list_path, sep = "\t", quote = FALSE, row.names = FALSE)

map <- merge(sample_list, meta2[, c("sample_id", "patient"), drop = FALSE], by = "sample_id", all.x = TRUE)
map <- map[, c("Tissue", "sample_id", "patient"), drop = FALSE]
map <- map[order(map$Tissue, map$patient, map$sample_id), , drop = FALSE]
write.table(map, file = map_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("=== Tissue 维度样本数/病人数统计 ===\n")
cat("输入（sample级别）：", sample_path, "\n", sep = "")
cat("输入（patient级别）：", patient_path, "\n", sep = "")
cat("输出：", out_path, "\n\n", sep = "")

cat("Sample 级别（unique sample_id）：总数 =", sample_n_total, "\n", sep = "")
for (t in all_tissues) {
  cat("  - ", t, ": ", out$sample_n[out$Tissue == t], "\n", sep = "")
}

cat("\nPatient 级别（unique patient）：总数 =", patient_n_total, "\n", sep = "")
for (t in all_tissues) {
  cat("  - ", t, ": ", out$patient_n[out$Tissue == t], "\n", sep = "")
}

# 交叉核对：映射表里的病人数与样本数是否一致（只作提示，不中断）
cat("\n[核对] meta_id_patient.tsv：unique sample_id =", length(unique(meta2$sample_id)),
    ", unique patient =", length(unique(meta2$patient)), "\n", sep = "")

cat("\n=== 详细清单输出 ===\n")
cat("Sample 清单：", sample_list_path, "\n", sep = "")
cat("Patient 清单：", patient_list_path, "\n", sep = "")
cat("Sample->Patient 映射：", map_path, "\n", sep = "")
