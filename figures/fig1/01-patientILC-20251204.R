#!/usr/bin/env Rscript

# 读取样本层面 ILC 比例与个体映射，以便在个体级别重新汇总组织/部位下的比例。
library(readr)
library(dplyr)
library(tidyr)

samples <- read_csv("fig1/fig1_ILC_proportion_by_Sample_Tissue.csv", show_col_types = FALSE)
meta <- read_tsv("fig1/meta_id_patient.tsv", show_col_types = FALSE)

# 确保样本和映射表包含 sample_id。
stopifnot("sample_id" %in% names(samples), "sample_id" %in% names(meta))

data <- left_join(samples, meta, by = "sample_id")
if (any(is.na(data$patient))) {
  warning("存在无法映射的 sample_id：", paste(unique(data$sample_id[is.na(data$patient)]), collapse = ", "))
}

# 按个体 + Tissue 汇总 counts，并重新计算比例，保证权重跟随细胞数。
patient_tissue <- data %>%
  filter(!is.na(Tissue)) %>%
  group_by(patient, Tissue) %>%
  summarise(
    total_count = sum(Total_Count, na.rm = TRUE),
    ilc_count = sum(ILC_Count, na.rm = TRUE),
    sample_n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    ilc_proportion = if_else(total_count > 0, ilc_count / total_count, NA_real_)
  )

# 按个体 + Site_Status 汇总，兼容 lesion / non_lesion / Healthy 等标签。
patient_site_status <- data %>%
  filter(!is.na(Site_Status)) %>%
  group_by(patient, Site_Status) %>%
  summarise(
    total_count = sum(Total_Count, na.rm = TRUE),
    ilc_count = sum(ILC_Count, na.rm = TRUE),
    sample_n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    ilc_proportion = if_else(total_count > 0, ilc_count / total_count, NA_real_)
  )

patient_status <- data %>%
  filter(!is.na(Status)) %>%
  group_by(patient, Status) %>%
  summarise(
    total_count = sum(Total_Count, na.rm = TRUE),
    ilc_count = sum(ILC_Count, na.rm = TRUE),
    sample_n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    ilc_proportion = if_else(total_count > 0, ilc_count / total_count, NA_real_)
  )

# 计算每个患者的 Dermis / Epidermis 比例并保留各自的比例。
patient_ratio <- patient_tissue %>%
  select(patient, Tissue, ilc_proportion) %>%
  pivot_wider(names_from = Tissue, values_from = ilc_proportion) %>%
  mutate(
    dermis_prop = Dermis,
    epidermis_prop = Epidermis,
    dermis_over_epidermis = {
      valid <- !is.na(Dermis) & !is.na(Epidermis) & Epidermis > 0
      ratio <- rep(NA_real_, n())
      ratio[valid] <- Dermis[valid] / Epidermis[valid]
      ratio
    }
  )

# 保存结果，供后续人工验证或绘图使用。
write_csv(patient_tissue, "fig1/patient_ILC_proportion_by_Tissue.csv")
write_csv(patient_site_status, "fig1/patient_ILC_proportion_by_Site_Status.csv")
write_csv(patient_ratio, "fig1/patient_ILC_ratio_Dermis_over_Epidermis_by_Patient.csv")
write_csv(patient_status, "fig1/patient_ILC_proportion_by_Status.csv")
