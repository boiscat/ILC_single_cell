
library(Seurat)

# 可视化obj对象中CD4基因的表达特征

obj <- readRDS("/lustre/home/yangyb/youmz/data/outputs/healthy_pso_merge_20251110/healthy_pso_harmony_integrated_de_20251110.rds")

# 设定想要展示的关键marker

marker_genes <- c(
  "COL1A1", "COL1A2", "DCN", "KLRB1", "IL7R", 
  #"GNLY","CD8A","KLRD1", 
  "KRT14", "KRT5", "KRT1", "KRT10",

  "TPSAB1", "TPSB2", "KIT", "LYZ", "CST3", "LST1",
  "AREG", "NKG7",
  
 #  "IL23R", 
  #"IL17A", "IL17F",
  # 下面是新加上的marker
  'CD3E','CD3G'
  #"PECAM1", "VWF", "CLDN5", "ITGAX", "CLEC10A", "CD207"
)

# 提取 T cell（设定 cell_anno/T/NK/T/NK/为T，常见label: "T", "T cell", "T/NK"等，按你的实际注释名称）
# 这里假设"T"为T细胞的注释标签，或者你应根据实际注释名称调整

# 首先确认 cell_anno 字段在 meta.data 中，并转换为字符向量
tcell_cells <- colnames(obj)[as.character(obj$cell_anno) == "T"]
tcell_obj <- subset(obj, cells = tcell_cells)


# 过滤去除核糖体基因影响的T细胞对象
# 核糖体基因通常以 "RPS" 和 "RPL" 开头

# 获取T细胞表达矩阵行名（基因名）
genes <- rownames(tcell_obj)
# 找到核糖体基因
ribosomal_genes <- grep("^RPS|^RPL", genes, value = TRUE)
# 去除这些基因的数据
tcell_obj <- subset(tcell_obj, features = setdiff(genes, ribosomal_genes))

# 过滤掉线粒体基因的影响
# 常见人类线粒体基因以 "MT-" 开头
mt_genes <- grep("^MT-", genes, value = TRUE)
genes_no_ribo_mt <- setdiff(genes, c(ribosomal_genes, mt_genes))
tcell_obj <- subset(tcell_obj, features = genes_no_ribo_mt)


# 可以保存T细胞亚群子对象
#saveRDS(tcell_obj, file = file.path(getwd(), "tcell_obj.rds"))

# 可选：输出T细胞数量
cat("Total T cells:", ncol(tcell_obj), "\n")

# 标准流程：归一化、筛选高变基因、标准化、PCA、邻接计算与聚类、UMAP
subset_obj <- NormalizeData(tcell_obj)
subset_obj <- FindVariableFeatures(subset_obj, nfeatures = 3000)  # 增加高变基因数
subset_obj <- ScaleData(subset_obj)
subset_obj <- RunPCA(subset_obj, npcs = 30)

# 增加紧密度可以通过减小邻居数k.param（更小的k意味着更细致的邻域）
subset_obj <- FindNeighbors(subset_obj, dims = 1:23, k.param = 50)

subset_obj <- FindClusters(subset_obj, resolution = 0.3)
subset_obj <- RunUMAP(subset_obj, dims = 1:23)


# 画UMAP降维分群图（DimPlot），按seurat_clusters分组并保存图片
png(file.path(out_dir, "subset_obj_umap_cluster.png"), width = 1800, height = 1400, res = 200)
print(
  DimPlot(subset_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
    ggtitle("UMAP Cluster Plot (subset_obj)")
)
dev.off()


# 绘制FeaturePlot，显示T细胞亚群对象（subset_obj）中CD4的表达情况
png(file.path(out_dir, "_featureplot_CD4.png"), width = 1800, height = 1400, res = 200)
print(
  FeaturePlot(obj, features = "CD4", reduction = "umap") +
    ggtitle("CD4 expression in T cell subset (subset_obj)")
)
dev.off()





subset_obj <- readRDS("/lustre/home/yangyb/youmz/data/subset_obj_tcell_subtyped_20251203.rds")


saveRDS(subset_obj, file = file.path(out_dir, "subset_obj_tcell_subtyped_20251204.rds"))


out_dir <- "/lustre/home/yangyb/youmz/data/figures/fig2"




# 任务：找出第16簇（seurat_clusters == 16）对应的细胞，对其标记为“Th17”亚群

# 初始化tcell_subtype列为NA
subset_obj$tcell_subtype <- NA


# 确保seurat_clusters是字符型（取决于Seurat版本，有时是数字，有时是字符）
cluster_ids <- as.character(subset_obj$seurat_clusters)
cells_th17 <- colnames(subset_obj)[cluster_ids == "16"]

# 定义0为 "T memory stem cell"
cells_tmstem <- colnames(subset_obj)[cluster_ids == "0"]
subset_obj$tcell_subtype[cells_tmstem] <- "TSCM"

# 标记1簇细胞为"iNKT"
cells_inkt <- colnames(subset_obj)[cluster_ids == "1"]
subset_obj$tcell_subtype[cells_inkt] <- "iNKT"


# 标记1簇细胞为"γδ T"
cells_gd_t <- colnames(subset_obj)[cluster_ids == "10"]
subset_obj$tcell_subtype[cells_gd_t] <- "γδ T"


# 标记16簇细胞为"Th17"
subset_obj$tcell_subtype[cells_th17] <- "Th17"

cells_eff_cd8 <- colnames(subset_obj)[cluster_ids %in% c("6", "5", "4")]
subset_obj$tcell_subtype[cells_eff_cd8] <- "Cytotoxic CD8 T"

# 标记4簇细胞为"PRF1high Cytotoxic CD8 T"
cells_prf1high_cd8 <- colnames(subset_obj)[cluster_ids == "4"]
subset_obj$tcell_subtype[cells_prf1high_cd8] <- "exhausted CD8 T"



cells_treg <- colnames(subset_obj)[cluster_ids %in% c("2", "9", "11")]
subset_obj$tcell_subtype[cells_treg] <- "Treg"


# DotPlot for NKT markers: TRAV10 (Vα24-Jα18), TRBV25-1, ZBTB16 (PLZF), KLRB1 (CD161), NKG7, GNLY, GZMB, PRF1, CXCR6, SLC4A10, TYROBP, KLRD1 (CD94), KLRF1

nkt_dotplot_genes <- c(
  "CD8A",'CD8B',
  "TRAV10",      # invariant NKT alpha chain
  "TRBV25-1",    # canonical iNKT beta chain
  "ZBTB16",      # PLZF key TF
  "KLRB1",       # CD161
  "NKG7",
  "GNLY", "GZMB", "PRF1",  # cytotoxic markers
  "CXCR6",
  "SLC4A10",     # common scRNA NKT marker
  "TYROBP",
  "KLRD1",       # CD94
  "KLRF1"
)

# Plot and save DotPlot for NKT-related markers grouped by tcell_subtype
png(file.path(out_dir, "subset_obj_DotPlot_NKT_markers.png"), width = 2600, height = 1400, res = 200)
print(
  DotPlot(
    subset_obj,
    features = nkt_dotplot_genes,
    group.by = "seurat_clusters",
    cols = c("lightgrey", "#0072B2")
  ) + RotatedAxis()
)
dev.off()





# 将剩余没有标注的细胞（tcell_subtype仍为NA的细胞）赋值为""
subset_obj$tcell_subtype[is.na(subset_obj$tcell_subtype)] <- ""

cat("Cell counts per tcell_subtype:\n")
print(tcell_subtype_counts)
subset_obj$tcell_subtype <- as.character(subset_obj$tcell_subtype)

# 绘制CD25、CD127、CXCR3在t细胞子集对象(subset_obj)上的FeaturePlot并保存图片
# 剔除3,8,12
exclude_clusters <- c("3", "12")
cells_to_exclude <- colnames(subset_obj)[cluster_ids %in% exclude_clusters]
subset_obj <- subset(subset_obj, cells = setdiff(colnames(subset_obj), cells_to_exclude))
cluster_ids <- as.character(subset_obj$seurat_clusters)  # 重新获取去除后的cluster_ids

# 剔除 tcell_subtype 为 "Unknown" 的细胞
cells_to_remove_unknown <- colnames(subset_obj)[subset_obj$tcell_subtype == "Unknown"]
if (length(cells_to_remove_unknown) > 0) {
  subset_obj <- subset(subset_obj, cells = setdiff(colnames(subset_obj), cells_to_remove_unknown))
  cluster_ids <- as.character(subset_obj$seurat_clusters)  # 重新更新cluster_ids
}


# 绘制带tcell_subtype的DimPlot并保存图片
png(file.path(out_dir, "subset_obj_DimPlot_tcell_subtype.png"), width = 1800, height = 1400, res = 200)
print(
  DimPlot(subset_obj, group.by = "tcell_subtype", reduction = "umap", label = TRUE)
)
dev.off()



# 将tcell_subtype为"Th"的细胞重新注释为"CD4-naive"
subset_obj$tcell_subtype[subset_obj$tcell_subtype == "Th"] <- "CD4-naive"


library(tidyverse)


# 确认marker_genes都在对象里
genes_to_plot <- c("CD3D","CD3E","CD3G","CD4",'LEF1',"TCF7","CD8A","CD8B",'NKG7',"LAG3","TIGIT","ZBTB16",

 "RIPK2","IL17A","IL17F","FOXP3","SELL",'CD28','FAS',

 
 "TRDV1","TRGV3", "TRGV4", "TRGV8"
 
 ) %>% unique()


png(file.path(out_dir, "subset_obj_DotPlot_tcell_subtype_markers.png"), width = 2800, height = 1300, res = 300)
print(
  DotPlot(subset_obj, features = genes_to_plot, group.by = "tcell_subtype", cols = c("lightgrey", "darkred")) +
    RotatedAxis() +
    ggtitle("DotPlot of marker genes across tcell_subtype")
)
dev.off()


png(file.path(out_dir, "subset_obj_DotPlot_Tcell_markers.png"), width = 3500, height = 1600, res = 200)
print(
  DotPlot(subset_obj, features = genes_to_plot) +
    ggtitle("DotPlot of T cell subset markers (subset_obj)")
)
dev.off()




# Nature级别的配色：使用Nature经典配色（近似于其杂志图常用的色条），可用RColorBrewer的"Dark2"、"Set1"或自定义色板
# 下面创建适合细胞亚群的色板并可视化legend以便检查

# 更加浅、明亮的配色方案（可进一步调整明度/色彩饱和度）
palette_use <- c(
  "#6FCF97", # 深一点绿色
  "#8FAADC", # 深一点紫色
  "#FFE066", # 深一点明黄
  "#B7E089", # 深一点黄绿
  "#FFAD8D", # 深一点橙色
  "#E1A2D7", # 深一点粉紫色
  "#85CFCB", # 深一点青绿
  "#FFE9A3", # 深一点米黄
  "#B7B5D2"  # 深一点灰紫
)

# 指定每个T细胞亚型的颜色（顺序需与你的tcell_subtype分组一致，否则可根据levels自动匹配）
if ("tcell_subtype" %in% colnames(subset_obj@meta.data)) {
  tcell_types <- levels(factor(subset_obj@meta.data$tcell_subtype))
  
  # 保证色板长度与分组数相符，超出则重复色板
  palette_use <- rep(palette_use, length.out = length(tcell_types))
  names(palette_use) <- tcell_types
  
  # 示例：绘制UMAP，应用Nature配色
  png(file.path(out_dir, "subset_obj_Tcell_DimPlot_by_cellanno_nature_palette.png"), width = 2300, height = 1600, res = 300)
  print(
    DimPlot(
      subset_obj, 
      reduction = "umap", 
      group.by = "tcell_subtype", 
      label = TRUE, 
      label.size = 3, 
      cols = palette_use
    ) 
  )
  dev.off()
}

library(ggplot2)
library(dplyr)

# 保存 subset_obj 为 RDS 文件，便于后续加载使用
saveRDS(subset_obj, file = file.path(out_dir, "subset_obj_tcell_subtyped_20251203.rds"))



# 假定 sample_id, tcell_subtype, Site_Status 都在 meta.data 里
# 列名矫正，你的数据列可能叫 Site_Status 或 condition，请确保一致
meta <- subset_obj@meta.data

# 统计subset_obj的各TcellSubtype在每个sample的数量
cell_count_tcell_sample <- meta %>%
  group_by(sample_id, TcellSubtype = tcell_subtype, Condition = Site_Status) %>%
  summarise(Freq = n(), .groups = "drop") %>%
  mutate(Source = "TcellSubset")


# obj中的ILC计数（从原始obj而不是subset_obj中拿ILC）
meta_obj <- obj@meta.data
if ("ILC" %in% meta_obj$cell_anno) {
  ilc_sample <- meta_obj %>%
    filter(cell_anno == "ILC") %>%
    group_by(sample_id, Condition = Site_Status) %>%
    summarise(Freq = n(), .groups = "drop") %>%
    mutate(TcellSubtype = "ILC", Source = "ILCfromObj")
} else {
  ilc_sample <- NULL
}

# 合并所有统计结果
cell_count_long <- bind_rows(
  cell_count_tcell_sample, ilc_sample
) %>% 
  mutate(
    TcellSubtype = factor(TcellSubtype, levels = unique(c(names(palette_use), "ILC")))
  )

# 检查数据
write.table(cell_count_long, file = file.path(out_dir, "debug_TcellSubtype_sample_counts.tsv"), sep = "\t", quote = F, row.names = F)

# 以样本为单位汇总: 计算每个Condition、T细胞亚型在各样本的数量分布，画均值±标准差
summary_per_type <- cell_count_long %>%
  group_by(TcellSubtype, Condition) %>%
  summarise(
    mean_count = mean(Freq),
    sd_count = sd(Freq),
    n_sample = n(),
    .groups = "drop"
  )

# 保证亚型配色一致
subtype_levels <- levels(cell_count_long$TcellSubtype)
condition_levels <- unique(subset_obj@meta.data$Site_Status)
# 按条形图示例配色（黑边黑点极差散点）重新设定配色，浅色为条形填充，黑色点和边：
bar_fill_colors <- c(
  "lesion" = "#cee2fb",
  "non_lesion" = "#dbd5f6",
  "Healthy" = "#f6e8d2" # 如有其他分组，补充
)
bar_fill_colors <- bar_fill_colors[condition_levels]
names(bar_fill_colors) <- condition_levels
point_shape_map <- c(
  "Healthy" = 16,      # circle
  "lesion" = 15,       # square
  "non_lesion" = 17,   # triangle
  "other" = 18         # diamond
)

# Condition 顺序统一排序为: Healthy, non_lesion, lesion
condition_levels <- c("Healthy", "non_lesion", "lesion")
cell_count_long$Condition <- factor(cell_count_long$Condition, levels = condition_levels)
summary_per_type$Condition <- factor(summary_per_type$Condition, levels = condition_levels)




# 绘图, 注意 y轴为每样本平均数，加误差线
p <- ggplot(summary_per_type, aes(x = TcellSubtype, y = mean_count, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width=0.8), width = 0.7, color = "black") +
geom_errorbar(aes(
  ymin = pmax(mean_count - sd_count, 0), 
  ymax = mean_count + sd_count
), position = position_dodge(width=0.8), width = 0.3, linewidth = 0.7)+
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = bar_fill_colors) +
  labs(
    x = "",
    y = "Cell count per sample (mean ± SD)",
    fill = "Condition"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = "right"
  )

ggsave(
  filename = file.path(out_dir, "subset_obj_Tcell_CountBar_BySubtype_sampleMeanSD.png"),
  plot = p,
  width = 9,
  height = 4,
  dpi = 300
)


#saveRDS(subset_obj, file = file.path(out_dir, "subset_obj_Tcell_clustered_annoed_1117.rds"))

saveRDS(obj, file = file.path(out_dir, "full_obj_saved_1117.rds"))
