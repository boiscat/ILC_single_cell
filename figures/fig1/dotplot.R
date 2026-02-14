suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

in_rds <- "/lustre/home/yangyb/youms/data/outputs/healthy_pso_merge_20251110/healthy_pso_harmony_integrated_de_20251110.rds"
out_dir <- "/lustre/home/yangyb/youms/data/figures/fig1"
out_png <- file.path(out_dir, "figure1_DimPlot_cell_anno.png")
out_pngpdf <- file.path(out_dir, "figure1_DimPlot_cell_anno_7.pdf")
out_dotplot <- file.path(out_dir, "fig1_marker_dotplot_by_cell_anno.png")
out_site_status <- file.path(out_dir, "fig1_cell_type_proportion_by_Site_Status.png")
out_tissue <- file.path(out_dir, "fig1_cell_type_proportion_by_Tissue.png")
out_dotplotpdf <- file.path(out_dir, "fig1_marker_dotplot_by_cell_anno_8.pdf")
out_site_statuspdf <- file.path(out_dir, "fig1_cell_type_proportion_by_Site_Status_7.pdf")
out_tissuepdf <- file.path(out_dir, "fig1_cell_type_proportion_by_Tissue_7.pdf")

if (!file.exists(in_rds)) stop("输入RDS不存在：", in_rds)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(in_rds)

if (!all(c("cell_anno", "seurat_clusters") %in% colnames(obj@meta.data))) {
  stop("meta.data 缺少必要列：cell_anno 或 seurat_clusters")
}

ann <- obj@meta.data$cell_anno
clust <- obj@meta.data$seurat_clusters
if (is.factor(ann)) ann <- as.character(ann)
if (is.factor(clust)) clust <- as.character(clust)
ann[clust == "5" | clust == 5] <- "ILG"
ann[clust == "6" | clust == 6] <- "NK"
ann[ann %in% c("T/NK", "TNK")] <- "T"
ann[ann == "ILG"] <- "ILC"

# 仅改显示名称：Myeloid -> Monocyte/Macrophage
ann[ann == "Myeloid"] <- "Monocyte/Macrophage"
desired_order <- c("FB", "ILC", "KC", "Mast cells", "Monocyte/Macrophage", "NK", "T")
other_levels <- setdiff(sort(unique(ann)), desired_order)
obj@meta.data$cell_anno <- factor(ann, levels = c(desired_order, other_levels))

obj <- subset(obj, subset = !(cell_anno %in% c("VE", "LE", "Melanocytes", "Pericyte")))
#2A批量图颜色
anno_colors <- c(
  "FB" = "#F0904D",
  "ILC" ="#EC7576",
  "KC" = "#14A0BD",
  "Mast cells" = "#6DA443",
  "Monocyte/Macrophage" = "#C6A3CE",
  "NK" = "#F7C8C8",
  "T" =  "#5E92CD"
)

# 如果出现未覆盖的新类别，自动补色（保证脚本不因新水平报错）
cell_lvls <- sort(unique(as.character(obj@meta.data$cell_anno)))
missing_lvls <- setdiff(cell_lvls, names(anno_colors))
if (length(missing_lvls) > 0) {
  extra_cols <- scales::hue_pal()(length(missing_lvls))
  names(extra_cols) <- missing_lvls
  anno_colors <- c(anno_colors, extra_cols)
}

backup_png <- function(path) {
  if (!file.exists(path)) return(invisible(FALSE))
  base <- tools::file_path_sans_ext(basename(path))
  ext <- tools::file_ext(path)
  ext_suffix <- if (nzchar(ext)) paste0(".", ext) else ""
  bak <- file.path(dirname(path), sprintf("%s.bak-%s%s", base, format(Sys.Date(), "%Y%m%d"), ext_suffix))
  if (!file.exists(bak)) file.copy(path, bak)
  invisible(TRUE)
}

backup_png(out_png)
#png(out_png, width = 2600, height = 2000, res = 400)
pdf_device <- if (capabilities("cairo")) cairo_pdf else pdf
# 开启 PDF 设备
pdf_device(
  out_pngpdf, 
  width = 6.5, 
  height = 5.0
)
print(
  DimPlot(
    obj,
    group.by = "cell_anno",
    cols = anno_colors,
    label = TRUE,
    pt.size = 0.05,
    label.size = 5,
  ) +
    ggtitle("Cell annotation") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      plot.margin = margin(10, 10, 10, 10)
    )
)
dev.off()

marker_genes <- c(
  "COL1A1", "COL1A2", "DCN", "KLRB1", "IL7R",
  "KRT14", "KRT5", "KRT1", "KRT10",
  "TPSAB1", "TPSB2", "KIT", "LYZ", "CST3", "LST1",
  "AREG", "NKG7",
  "CD3E", "CD3G"
)
library(scales)
backup_png(out_dotplot)
#png(out_dotplot, width = 2000, height = 900, res = 200)
pdf_device(
  out_dotplotpdf, 
  width = 10.0, 
  height = 4.5
)
print(
  DotPlot(
    obj,
    features = marker_genes,
    group.by = "cell_anno", 
    #cols = c("#5E92CD", "#E71F19"),
    dot.scale = 7
  ) + 
    # 使用 gradient2 实现三色渐变
  scale_color_gradientn(
    # values 定义颜色在 scale 中的相对位置（0到1之间）
    # 我们把 -2, 0, 2 映射到 0, 0.5, 1
    # 如果想让白色更窄，可以将转折点设置得离 0.5 更近
    values = rescale(c(-2, -0.1, 0, 0.1, 2)), 
    # 对应的颜色也要补齐：蓝 -> 浅蓝 -> 白 -> 浅红 -> 红
    colors = c("#5E92CD", "#ADC9E6", "white", "#F39894", "#E71F19"),
    limits = c(-2, 2),
    oob = squish
  ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 14),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
    ) +
    ggtitle("DotPlot of key marker genes by cell_anno")
)
dev.off()

if (!("Site_Status" %in% colnames(obj@meta.data))) {
  obj@meta.data <- obj@meta.data %>%
    mutate(
      Site_Status = case_when(
        Status == "Healthy" ~ "Healthy",
        Status == "Psoriasis" & Site == "lesion" ~ "lesion",
        Status == "Psoriasis" & Site == "non_lesion" ~ "non_lesion",
        TRUE ~ "Unknown"
      )
    )
}
if ("Site_Status" %in% colnames(obj@meta.data)) {
  obj@meta.data$Site_Status <- factor(obj@meta.data$Site_Status, levels = c("Healthy", "lesion", "non_lesion", "Unknown"))
}
if ("Tissue" %in% colnames(obj@meta.data)) {
  obj@meta.data$Tissue <- factor(obj@meta.data$Tissue, levels = c("Dermis", "Epidermis"))
}

fill_colors <- anno_colors

cell_prop_df <- obj@meta.data %>%
  group_by(Site_Status, cell_anno) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(Site_Status) %>%
  mutate(prop = cell_count / sum(cell_count)) %>%
  ungroup()

backup_png(out_site_status)
#png(out_site_status, width = 1800, height = 1500, res = 300)
pdf_device(
  out_site_statuspdf, 
  width = 6.0, 
  height = 5.0
)
print(
  ggplot(cell_prop_df, aes(x = Site_Status, y = prop, fill = cell_anno)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    scale_fill_manual(values = fill_colors) +
    ylab("Proportion of cells") +
    xlab("Site_Status") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    ) +
    ggtitle("Cell type proportion by Site_Status") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1))
)
dev.off()

tissue_prop_df <- obj@meta.data %>%
  group_by(Tissue, cell_anno) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(Tissue) %>%
  mutate(prop = cell_count / sum(cell_count)) %>%
  ungroup()

backup_png(out_tissue)
#png(out_tissue, width = 1800, height = 1500, res = 300)
pdf_device(
  out_tissuepdf, 
  width = 6.0, 
  height = 5.0
)
print(
  ggplot(tissue_prop_df, aes(x = Tissue, y = prop, fill = cell_anno)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    scale_fill_manual(values = fill_colors) +
    ylab("Proportion of cells") +
    xlab("Tissue") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    ) +
    ggtitle("Cell type proportion by Tissue") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1))
)
dev.off()

message("已输出：")
message("- ", out_png)
message("- ", out_dotplot)
message("- ", out_site_status)
message("- ", out_tissue)
