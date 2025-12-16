#!/usr/bin/env Rscript

# 依托 01-patientILC-20251204.R 的汇总，补全 fig1 原始样式的患者层面图。
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

theme_classic2 <- function(base_size = 14, ...) {
  theme_classic(base_size = base_size, ...) +
    theme(panel.background = element_blank())
}

patient_tissue <- read_csv("fig1/patient_ILC_proportion_by_Tissue.csv", show_col_types = FALSE)
patient_site <- read_csv("fig1/patient_ILC_proportion_by_Site_Status.csv", show_col_types = FALSE)
patient_status <- read_csv("fig1/patient_ILC_proportion_by_Status.csv", show_col_types = FALSE)
patient_ratio <- read_csv("fig1/patient_ILC_ratio_Dermis_over_Epidermis_by_Patient.csv", show_col_types = FALSE)
sample_tissue <- read_csv("fig1/fig1_ILC_proportion_by_Sample_Tissue.csv", show_col_types = FALSE)

out_dir <- "fig1/ILC"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

theme_fig1 <- theme_classic2(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# 组织图：Dermis vs Epidermis
tissue_colors <- c(Dermis = "#C4D9F1", Epidermis = "#F8D6B7")
tissue_summary <- patient_tissue %>%
  group_by(Tissue) %>%
  summarise(
    mean_prop = mean(ilc_proportion, na.rm = TRUE),
    se_prop = sd(ilc_proportion, na.rm = TRUE) / sqrt(sum(!is.na(ilc_proportion))),
    .groups = "drop"
  ) %>%
  mutate(
    lower = pmax((mean_prop - se_prop) * 100, 0),
    upper = (mean_prop + se_prop) * 100
  )

tissue_plot <- ggplot(tissue_summary, aes(x = Tissue, y = mean_prop * 100, fill = Tissue)) +
  geom_col(width = 0.65, color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.18) +
  geom_jitter(
    data = patient_tissue,
    aes(x = Tissue, y = ilc_proportion * 100),
    width = 0.18, height = 0, size = 2, color = "black"
  ) +
  scale_fill_manual(values = tissue_colors) +
  scale_y_continuous(
    labels = percent_format(scale = 1),
    limits = c(0, 25),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = "Tissue",
    y = "ILC proportion",
    title = "ILC proportion by Tissue (patient-level)"
  ) +
  theme_fig1 +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "fig1_ILC_proportion_Tissue_mean_SE.png"), tissue_plot, width = 5, height = 4.5, dpi = 150)

# Status 图：Healthy vs Psoriasis
status_colors <- c(Healthy = "#C4D9F1", Psoriasis = "#F6BDC1")
status_summary <- patient_status %>%
  group_by(Status) %>%
  summarise(
    mean_prop = mean(ilc_proportion, na.rm = TRUE),
    se_prop = sd(ilc_proportion, na.rm = TRUE) / sqrt(sum(!is.na(ilc_proportion))),
    .groups = "drop"
  ) %>%
  mutate(
    lower = pmax((mean_prop - se_prop) * 100, 0),
    upper = (mean_prop + se_prop) * 100
  )

status_plot <- ggplot(status_summary, aes(x = Status, y = mean_prop * 100, fill = Status)) +
  geom_col(width = 0.65, color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.18) +
  geom_jitter(
    data = patient_status,
    aes(x = Status, y = ilc_proportion * 100),
    width = 0.18, height = 0, size = 2, color = "black"
  ) +
  scale_fill_manual(values = status_colors) +
  scale_y_continuous(
    labels = percent_format(scale = 1),
    limits = c(0, 25),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = "Status",
    y = "ILC proportion",
    title = "ILC proportion by Status (patient-level)"
  ) +
  theme_fig1 +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "fig1_ILC_proportion_Status_mean_SE.png"), status_plot, width = 5, height = 4.5, dpi = 150)

# Lesion / non_lesion 图：保持原图风格
lesion_colors <- c(lesion = "#D8C8F2", non_lesion = "#C4D9F1")
lesion_summary <- patient_site %>%
  filter(Site_Status %in% c("lesion", "non_lesion")) %>%
  group_by(Site_Status) %>%
  summarise(
    mean_prop = mean(ilc_proportion, na.rm = TRUE),
    se_prop = sd(ilc_proportion, na.rm = TRUE) / sqrt(sum(!is.na(ilc_proportion))),
    .groups = "drop"
  ) %>%
  mutate(
    lower = pmax((mean_prop - se_prop) * 100, 0),
    upper = (mean_prop + se_prop) * 100
  )

lesion_plot <- ggplot(lesion_summary, aes(x = Site_Status, y = mean_prop * 100, fill = Site_Status)) +
  geom_col(width = 0.65, color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.18) +
  geom_jitter(
    data = patient_site %>% filter(Site_Status %in% c("lesion", "non_lesion")),
    aes(x = Site_Status, y = ilc_proportion * 100),
    width = 0.18, height = 0, size = 2, color = "black"
  ) +
  scale_fill_manual(values = lesion_colors) +
  scale_y_continuous(
    labels = percent_format(scale = 1),
    limits = c(0, 25),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = "Site Status",
    y = "ILC proportion",
    title = "ILC proportion by Site Status (patient-level)"
  ) +
  theme_fig1 +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0)
  )

ggsave(file.path(out_dir, "fig1_ILC_proportion_Site_Status_mean_SE.png"), lesion_plot, width = 5, height = 4.5, dpi = 150)

# 样本层面：按 sample_id 展示真皮/表皮样本的 ILC 比例（ILC_Count / Total_Count）
sample_tissue_plot_data <- sample_tissue %>%
  filter(!is.na(sample_id), !is.na(Tissue), !is.na(ILC_Proportion)) %>%
  mutate(
    Tissue = factor(Tissue, levels = c("Dermis", "Epidermis"))
  ) %>%
  arrange(Tissue, desc(ILC_Proportion)) %>%
  mutate(sample_id = factor(sample_id, levels = unique(sample_id)))

sample_ymax <- max(sample_tissue_plot_data$ILC_Proportion, na.rm = TRUE)

sample_plot <- ggplot(sample_tissue_plot_data, aes(x = sample_id, y = ILC_Proportion, fill = Tissue)) +
  geom_col(width = 0.75, color = "black") +
  scale_fill_manual(values = tissue_colors) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    limits = c(0, min(0.25, sample_ymax * 1.15)),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = "Sample",
    y = "ILC proportion",
    title = "ILC proportion by Sample",
    fill = "Tissue"
  ) +
  theme_fig1 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

ggsave(file.path(out_dir, "fig1_ILC_proportion_by_Sample_Tissue_sample.png"), sample_plot, width = 10.5, height = 4.8, dpi = 150)

# 样本层面：按 Tissue 汇总（均值±SE + 散点）
sample_tissue_summary <- sample_tissue %>%
  filter(!is.na(Tissue), !is.na(ILC_Proportion)) %>%
  group_by(Tissue) %>%
  summarise(
    mean_prop = mean(ILC_Proportion, na.rm = TRUE),
    se_prop = sd(ILC_Proportion, na.rm = TRUE) / sqrt(sum(!is.na(ILC_Proportion))),
    .groups = "drop"
  ) %>%
  mutate(
    lower = pmax((mean_prop - se_prop) * 100, 0),
    upper = pmin((mean_prop + se_prop) * 100, 100)
  )

sample_tissue_mean_plot <- ggplot(sample_tissue_summary, aes(x = Tissue, y = mean_prop * 100, fill = Tissue)) +
  geom_col(width = 0.65, color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.18) +
  geom_jitter(
    data = sample_tissue %>% filter(!is.na(Tissue), !is.na(ILC_Proportion)),
    aes(x = Tissue, y = ILC_Proportion * 100),
    width = 0.18, height = 0, size = 2, color = "black"
  ) +
  scale_fill_manual(values = tissue_colors) +
  scale_y_continuous(
    labels = percent_format(scale = 1),
    limits = c(0, 25),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = "Tissue",
    y = "ILC proportion",
    title = "ILC proportion by Tissue (sample-level)"
  ) +
  theme_fig1 +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "fig1_ILC_proportion_Tissue_sample_mean_SE.png"), sample_tissue_mean_plot, width = 5, height = 4.5, dpi = 150)

# Dermis / Epidermis 比值
ratio_plot_data <- patient_ratio %>%
  filter(!is.na(dermis_over_epidermis)) %>%
  mutate(patient = factor(patient, levels = patient))

ratio_ymax <- max(ratio_plot_data$dermis_over_epidermis, na.rm = TRUE) * 1.15

ratio_plot <- ggplot(ratio_plot_data, aes(x = patient, y = dermis_over_epidermis)) +
  geom_col(fill = "#A8C7E7", color = "black") +
  geom_point(color = "black", size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)), limits = c(0, ratio_ymax)) +
  labs(
    x = "Patient",
    y = "Dermis / Epidermis ILC proportion ratio",
    title = "ILC proportion ratio Dermis / Epidermis by Patient"
  ) +
  theme_fig1 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "fig1_ILC_ratio_Dermis_over_Epidermis_by_Patient.png"), ratio_plot, width = 6, height = 4.5, dpi = 150)

# ----------- 修改：堆叠图改为并列柱状图 ----------- #

plot_ilc_by_condition_unstacked <- function(data, group_col, colors, order = NULL, title, legend_title, filename) {
  aggregated <- data %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(group_value = .data[[group_col]]) %>%
    summarise(ilc_cells = sum(ilc_count, na.rm = TRUE), .groups = "drop") %>%
    mutate(proportion = ilc_cells / sum(ilc_cells, na.rm = TRUE))

  group_levels <- if (!is.null(order)) {
    intersect(order, unique(aggregated$group_value))
  } else {
    unique(aggregated$group_value)
  }

  aggregated <- aggregated %>%
    mutate(group_value = factor(group_value, levels = group_levels)) %>%
    filter(!is.na(group_value)) %>%
    arrange(group_value)

  if (nrow(aggregated) == 0) {
    warning("没有可用的 ILC 计数用于 ", title)
    return(invisible(NULL))
  }

  plot <- ggplot(aggregated, aes(x = group_value, y = proportion, fill = group_value)) +
    geom_col(width = 0.55, color = "black") +
    geom_text(
      aes(label = scales::percent(proportion, accuracy = 0.1)),
      vjust = -0.5, 
      size = 3.6
    ) +
    {
      color_map <- colors[group_levels]
      color_map <- color_map[!is.na(color_map)]
      scale_fill_manual(values = color_map, breaks = names(color_map))
    } +
    scale_y_continuous(
      labels = label_percent(accuracy = 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.12))
    ) +
    labs(
      x = legend_title,
      y = "Proportion of total ILC cells",
      title = title,
      fill = legend_title
    ) +
    theme_fig1 +
    theme(
      axis.text.x = element_text(),
      axis.ticks.x = element_line(),
      legend.position = "top"
    )

  ggsave(file.path(out_dir, filename), plot, width = 4.2, height = 4.2, dpi = 150)
}

plot_ilc_by_condition_unstacked(
  patient_tissue,
  group_col = "Tissue",
  colors = tissue_colors,
  order = c("Dermis", "Epidermis"),
  title = "Total ILC cells by Tissue",
  legend_title = "Tissue",
  filename = "fig1_ILC_total_ILC_Tissue_unstacked.png"
)

plot_ilc_by_condition_unstacked(
  patient_status,
  group_col = "Status",
  colors = status_colors,
  order = c("Healthy", "Psoriasis"),
  title = "Total ILC cells by Status",
  legend_title = "Status",
  filename = "fig1_ILC_total_ILC_Status_unstacked.png"
)

plot_ilc_by_condition_unstacked(
  patient_site,
  group_col = "Site_Status",
  colors = c(lesion_colors, Healthy = "#C4D9F1"),
  order = c("lesion", "non_lesion", "Healthy"),
  title = "Total ILC cells by Site Status",
  legend_title = "Site Status",
  filename = "fig1_ILC_total_ILC_Site_Status_unstacked.png"
)

plot_ilc_by_condition_unstacked(
  patient_site %>% filter(Site_Status %in% c("lesion", "non_lesion")),
  group_col = "Site_Status",
  colors = lesion_colors[c("lesion", "non_lesion")],
  order = c("non_lesion", "lesion"),
  title = "Total ILC cells: Lesion vs non-lesion",
  legend_title = "Site Status",
  filename = "fig1_ILC_total_ILC_lesion_vs_nonlesion_unstacked.png"
)

# ---------------------------------------------- #

# ----------- 新增：ILC vs Other cells 堆叠图改为并列柱状图 ----------- #

plot_ilc_vs_other_unstacked <- function(data, group_col, order = NULL, title, x_label, filename) {
  aggregated <- data %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(group_value = .data[[group_col]]) %>%
    summarise(
      total_cells = sum(total_count, na.rm = TRUE),
      ilc_cells = sum(ilc_count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(other_cells = pmax(total_cells - ilc_cells, 0)) %>%
    select(group_value, ilc_cells, other_cells)

  group_levels <- if (!is.null(order)) {
    intersect(order, unique(aggregated$group_value))
  } else {
    unique(aggregated$group_value)
  }

  plot_data <- aggregated %>%
    mutate(group_value = factor(group_value, levels = group_levels)) %>%
    filter(!is.na(group_value)) %>%
    arrange(group_value) %>%
    tidyr::pivot_longer(
      cols = c(ilc_cells, other_cells),
      names_to = "cell_type",
      values_to = "count"
    ) %>%
    mutate(
      cell_type = recode(
        cell_type,
        ilc_cells = "ILC",
        other_cells = "Other cells"
      ),
      cell_type = factor(cell_type, levels = c("Other cells", "ILC"))
    ) %>%
    group_by(group_value) %>%
    mutate(
      total = sum(count, na.rm = TRUE),
      proportion = if_else(total > 0, count / total, NA_real_)
    ) %>%
    ungroup()

  if (nrow(plot_data) == 0) {
    warning("没有可用的 total_count/ilc_count 用于 ", title)
    return(invisible(NULL))
  }

  dodge <- position_dodge(width = 0.7)
  plot <- ggplot(plot_data, aes(x = group_value, y = proportion, fill = cell_type)) +
    geom_col(width = 0.6, color = "black", position = dodge) +
    geom_text(
      aes(label = scales::percent(proportion, accuracy = 0.1)),
      position = dodge,
      vjust = -0.5,
      size = 3.4
    ) +
    scale_fill_manual(values = c(`Other cells` = "#E8E8E8", ILC = "#78BCE6")) +
    scale_y_continuous(
      labels = label_percent(accuracy = 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.12))
    ) +
    labs(
      x = x_label,
      y = "细胞占比",
      title = title,
      fill = "细胞类别"
    ) +
    theme_fig1 +
    theme(
      axis.text.x = element_text(),
      axis.ticks.x = element_line(),
      legend.position = "top"
    )

  ggsave(file.path(out_dir, filename), plot, width = 5.6, height = 4.6, dpi = 150)
}

plot_ilc_vs_other_unstacked(
  patient_tissue,
  group_col = "Tissue",
  order = c("Dermis", "Epidermis"),
  title = "按 Tissue 汇总的 ILC 与其他细胞数量",
  x_label = "Tissue",
  filename = "fig1_ILC_total_counts_Tissue_unstacked.png"
)

plot_ilc_vs_other_unstacked(
  patient_site,
  group_col = "Site_Status",
  order = c("lesion", "non_lesion", "Healthy"),
  title = "按 Site Status 汇总的 ILC 与其他细胞数量",
  x_label = "Site Status",
  filename = "fig1_ILC_total_counts_Site_Status_unstacked.png"
)

# -------------------------------------------------------- #


library(tidyverse)
library(scales) # 确保加载scales用于百分比格式化

cd73_condition_summary_path <- "todolist拆分-20251120/任务6-CD73阳性ILC互作-20251120/CD73阳性ILC互作-全流程-20251120/outputs/tables/CD73_ILC_Tcell_conditionSummary_20251120.tsv"

if (file.exists(cd73_condition_summary_path)) {
  cd73_condition <- read_tsv(cd73_condition_summary_path, show_col_types = FALSE)
  
  # 数据处理：计算占比
  condition_plot_data <- cd73_condition %>%
    filter(Site_Status %in% c("Healthy", "lesion", "non_lesion")) %>%
    mutate(condition = if_else(Site_Status == "Healthy", "Healthy", "Psoriasis")) %>%
    group_by(condition) %>%
    summarise(cd73_cells = sum(total_cd73_ilc, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      condition = factor(condition, levels = c("Healthy", "Psoriasis")),
      # 计算相对比例
      proportion = cd73_cells / sum(cd73_cells, na.rm = TRUE), 
      # 生成标签文本 (例如: "45%")
      label = percent(proportion, accuracy = 1) 
    )

  if (nrow(condition_plot_data) > 0) {
    condition_plot <- ggplot(condition_plot_data, aes(x = condition, y = proportion, fill = condition)) +
      geom_col(width = 0.55, color = "black") +
      # --- 关键部分：添加顶部比例标签 ---
      geom_text(
        aes(label = label), 
        vjust = -0.5,    # 负值将文字向上推，使其位于柱子上方
        size = 4, 
        color = "black",
        fontface = "bold" # 可选：加粗字体
      ) +
      # ----------------------------------
      scale_fill_manual(values = status_colors) +
      scale_y_continuous(
        labels = label_percent(accuracy = 1),
        limits = c(0, 1), # 纵坐标固定0-1
        # 重要：增加顶部留白(0.1)，防止文字被图片边缘切掉
        expand = expansion(mult = c(0, 0.1)) 
      ) +
      labs(
        x = "",
        y = "Proportion of CD73+ ILC", # 注意：这里你的逻辑是CD73在两组间的分布占比
        title = "CD73+ ILC distribution between Healthy and Psoriasis",
        fill = "Condition"
      ) +
      theme_fig1 +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5) # 标题居中
      )
    
    # 打印图片
    print(condition_plot)
  }
}




















# CD73+ vs CD73- 比例图
cd73_summary_path <- "fig3-1/outputs/CD73_ILC_DE_summary.csv"
if (file.exists(cd73_summary_path)) {
  cd73_summary <- read_csv(cd73_summary_path, show_col_types = FALSE)
  if (nrow(cd73_summary) > 0) {
    cd73_ratio <- tibble(
      status = c("CD73_neg", "CD73_pos"),
      count = c(cd73_summary$cd73_neg[1], cd73_summary$cd73_pos[1])
    ) %>%
      mutate(
        status = factor(status, levels = c("CD73_neg", "CD73_pos")),
        proportion = count / sum(count, na.rm = TRUE)
      )

    cd73_plot <- ggplot(cd73_ratio, aes(x = status, y = proportion, fill = status)) +
      geom_col(color = "black", width = 0.55) +
      geom_text(aes(label = scales::percent(proportion, accuracy = 1)), vjust = -0.5, size = 4, color = "black") +
      scale_fill_manual(values = c(CD73_neg = "#C4D9F1", CD73_pos = "#F6BDC1")) +
      scale_y_continuous(
        labels = label_percent(accuracy = 1),
        limits = c(0, 1),
        expand = expansion(mult = c(0, 0.08))
      ) +
      labs(
        x = "",
        y = "Proportion of ILC",
        title = "CD73+ vs CD73- ILCs",
        fill = "CD73 status"
      ) +
      theme_fig1 +
      theme(
        axis.text.x = element_text(size = 12),
        axis.ticks.x = element_line(),
        legend.position = "top"
      )

    ggsave(file.path(out_dir, "fig1_ILC_CD73_pos_neg_ratio_unstacked.png"), cd73_plot, width = 4.5, height = 4.5, dpi = 150)
  } else {
    warning("CD73 summary file exists but is empty: ", cd73_summary_path)
  }
} else {
  warning("CD73 summary file not found: ", cd73_summary_path)
}

cd73_condition_summary_path <- "todolist拆分-20251120/任务6-CD73阳性ILC互作-20251120/CD73阳性ILC互作-全流程-20251120/outputs/tables/CD73_ILC_Tcell_conditionSummary_20251120.tsv"
if (file.exists(cd73_condition_summary_path)) {
  cd73_condition <- read_tsv(cd73_condition_summary_path, show_col_types = FALSE)
    condition_plot_data <- cd73_condition %>%
      filter(Site_Status %in% c("Healthy", "lesion", "non_lesion")) %>%
      mutate(condition = if_else(Site_Status == "Healthy", "Healthy", "Psoriasis")) %>%
      group_by(condition) %>%
      summarise(cd73_cells = sum(total_cd73_ilc, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        condition = factor(condition, levels = c("Healthy", "Psoriasis")),
        proportion = cd73_cells / sum(cd73_cells, na.rm = TRUE),
        label = percent(proportion, accuracy = 1)
      )

  if (nrow(condition_plot_data) > 0) {
    condition_plot <- ggplot(condition_plot_data, aes(x = condition, y = proportion, fill = condition)) +
      geom_col(width = 0.55, color = "black") +
      geom_text(aes(label = label), vjust = -0.6, size = 4, color = "black") +
      scale_fill_manual(values = status_colors) +
      scale_y_continuous(
        labels = label_percent(accuracy = 1),
        limits = c(0, 1),
        expand = expansion(mult = c(0, 0.06))
      ) +
      labs(
        x = "",
        y = "Proportion of CD73+ ILC",
        title = "CD73+ ILC distribution between Healthy and Psoriasis",
        fill = "Condition"
      ) +
      theme_fig1 +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 12)
      )

    ggsave(file.path(out_dir, "fig1_ILC_CD73_positive_by_condition_unstacked.png"), condition_plot, width = 4.5, height = 4.5, dpi = 150)
  } else {
    warning("CD73 condition summary has no usable rows: ", cd73_condition_summary_path)
  }
} else {
  warning("CD73 condition summary file not found: ", cd73_condition_summary_path)
}
