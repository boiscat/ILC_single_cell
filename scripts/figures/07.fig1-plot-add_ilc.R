


obj <- readRDS("/lustre/home/yangyb/youmz/data/figures/todolist拆分-20251120/任务1-lesion修正-20251120/lesion修正-功能复跑-20251120/outputs/objects/full_obj_lesionFixed_20251120.rds")

#===============================
# 区分真皮层和表皮层，按不同样本统计 ILC 细胞比例
#===============================

library(dplyr)
library(ggplot2)
library(tidyr)

# 这里假定：
# - 真皮/表皮信息保存在 meta.data 的 `Tissue` 列（例如 "Dermis" / "Epidermis" 或中文“真皮”/“表皮”）
# - 每个样本的名称在 `sample_id` 列
# - ILC 细胞由 `cell_anno == "ILC"` 标记

meta_all <- obj@meta.data

# 1）先统计每个「样本 × 组织层（真皮/表皮）」的总细胞数
total_by_sample_layer <- meta_all %>%
  filter(!is.na(Tissue), !is.na(sample_id)) %>%
  group_by(sample_id, Tissue) %>%
  summarise(Total_Count = n(), .groups = "drop")

# 2）再统计每个「样本 × 组织层」中的 ILC 细胞数
ilc_by_sample_layer <- meta_all %>%
  filter(cell_anno == "ILC",
         !is.na(Tissue), !is.na(sample_id)) %>%
  group_by(sample_id, Tissue,Status,Site_Status) %>%
  summarise(ILC_Count = n(), .groups = "drop")



# 3）合并并计算比例：ILC 数 / 该样本该层的总细胞数
# 修正 left_join() 出错：只用两个共有的 key 字段进行合并
ilc_stat_df <- total_by_sample_layer %>%
  left_join(ilc_by_sample_layer,
            by = c("sample_id", "Tissue")) %>%
  mutate(
    ILC_Count = ifelse(is.na(ILC_Count), 0, ILC_Count),
    ILC_Proportion = ILC_Count / Total_Count
  )

cat("按样本、真皮/表皮分层的 ILC 比例统计:\n")
print(ilc_stat_df)

# 如果后面要复用，也可以保存一个表格
out_dir_ilc <- "/lustre/home/yangyb/youmz/data/figures/fig1"
if (!dir.exists(out_dir_ilc)) {
  dir.create(out_dir_ilc, recursive = TRUE, showWarnings = FALSE)
}
write.csv(ilc_stat_df,
          file = file.path(out_dir_ilc, "fig1_ILC_proportion_by_Sample_Tissue.csv"),
          row.names = FALSE)

# 4）画图：每个样本，真皮/表皮分开的 ILC 比例柱状图
png(file.path(out_dir_ilc, "fig1_ILC_proportion_by_Sample_Tissue.png"),
    width = 1800, height = 900, res = 200)
print(
  ggplot(ilc_stat_df,
         aes(x = sample_id, y = ILC_Proportion, fill = Tissue)) +
    geom_bar(stat = "identity", position = "dodge") +
    ylab("ILC proportion") +
    xlab("Sample") +
    ggtitle("Proportion of ILC in Dermis/Epidermis by Sample") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.title = element_blank()
    )
)
dev.off()

# 5）去掉 ILC 比例过低的样本，做“真皮 vs 表皮”带误差线的图（类似示例配色）

# 可以根据需要调整这个阈值，比如 0.005 / 0.01
prop_cutoff <- 0.005

ilc_filtered <- ilc_stat_df %>%
  dplyr::filter(ILC_Proportion >= prop_cutoff)

# 为每个组织层计算均值和标准误
ilc_tissue_summary <- ilc_filtered %>%
  dplyr::group_by(Tissue) %>%
  dplyr::summarise(
    mean_prop = mean(ILC_Proportion),
    sd_prop   = sd(ILC_Proportion),
    n         = dplyr::n(),
    se_prop   = sd_prop / sqrt(n),
    .groups   = "drop"
  )

# 自定义真皮 / 表皮配色（可以根据需要微调）
tissue_colors <- c(
  "Dermis"    = "#c9e3ff",  # 浅蓝
  "Epidermis" = "#f8e2c0"   # 浅橘
)

png(file.path(out_dir_ilc, "fig1_ILC_proportion_Tissue_mean_SE.png"),
    width = 1200, height = 1000, res = 200)
print(
  ggplot(ilc_tissue_summary,
         aes(x = Tissue, y = mean_prop, fill = Tissue)) +
    # 柱子（显示均值）
    geom_col(width = 0.6, color = "black") +
    # 误差线（均值 ± 标准误）
    geom_errorbar(
      aes(ymin = mean_prop - se_prop,
          ymax = mean_prop + se_prop),
      width = 0.13, size = 0.7
    ) +
    # 叠加每个样本的散点
    geom_jitter(
      data = ilc_filtered,
      aes(x = Tissue, y = ILC_Proportion),
      width = 0.08, height = 0,
      shape = 21, size = 2.5,
      fill = "black", color = "black",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = tissue_colors) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    xlab("") +
    ylab("ILC proportion") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      axis.title  = element_text(size = 14),
      plot.title  = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = "none"
    )
)
dev.off()


# 分析Status，统计不同Status（如Healthy、Psoriasis等）下ILC比例的均值和标准误，并画图
# 1. 汇总数据
ilc_status_summary <- ilc_filtered %>%
  dplyr::group_by(Status) %>%
  dplyr::summarise(
    mean_prop = mean(ILC_Proportion),
    sd_prop   = sd(ILC_Proportion),
    n         = dplyr::n(),
    se_prop   = sd_prop / sqrt(n),
    .groups   = "drop"
  )

cat("不同Status分组下ILC比例均值与标准误：\n")
print(ilc_status_summary)

# 2. 可视化: 各Status下ILC比例均值及标准误，并对每个样本添加散点
status_colors <- c(
  "Healthy"   = "#c9e3ff",
  "Psoriasis" = "#FFE0E0"
)
# 如果需要适配其它Status，在这里补充色板

png(file.path(out_dir_ilc, "fig1_ILC_proportion_Status_mean_SE.png"),
    width = 1200, height = 1000, res = 200)
print(
  ggplot(ilc_status_summary, aes(x = Status, y = mean_prop, fill = Status)) +
    geom_col(width = 0.6, color = "black") +
    geom_errorbar(
      aes(ymin = mean_prop - se_prop,
          ymax = mean_prop + se_prop),
      width = 0.13, size = 0.7
    ) +
    # 添加每个样本的散点
    geom_jitter(
      data = ilc_filtered,
      aes(x = Status, y = ILC_Proportion),
      width = 0.08, height = 0,
      shape = 21, size = 2.5,
      fill = "black", color = "black",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = status_colors) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    xlab("") +
    ylab("ILC proportion") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      axis.title  = element_text(size = 14),
      plot.title  = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = "none"
    )
)
dev.off()


# 分析Site_Status，统计不同Site_Status（如lesion、non_lesion）下ILC比例的均值和标准误，并画图
# 1. 汇总数据
ilc_site_status_summary <- ilc_filtered %>%
  dplyr::group_by(Site_Status) %>%
  dplyr::summarise(
    mean_prop = mean(ILC_Proportion),
    sd_prop   = sd(ILC_Proportion),
    n         = dplyr::n(),
    se_prop   = sd_prop / sqrt(n),
    .groups   = "drop"
  )

cat("不同Site_Status分组下ILC比例均值与标准误：\n")
print(ilc_site_status_summary)

# 可视化包含Healthy的样本

cat("不同Site_Status分组下ILC比例均值与标准误（含Healthy）：\n")
print(ilc_site_status_summary)

site_status_colors <- c(
  "lesion"      = "#c9e3ff",
  "non_lesion"  = "#FFE0E0",
  "Healthy"     = "#bcee8a"
)

# 只可视化lesion和non_lesion，不包含Healthy
ilc_site_status_summary_plot <- ilc_site_status_summary %>%
  filter(Site_Status != "Healthy")

ilc_filtered_plot <- ilc_filtered %>%
  filter(Site_Status != "Healthy")

# 绘图
png(file.path(out_dir_ilc, "fig1_ILC_proportion_Site_Status_mean_SE.png"),
    width = 1200, height = 1000, res = 200)
print(
  ggplot(ilc_site_status_summary_plot, aes(x = Site_Status, y = mean_prop, fill = Site_Status)) +
    geom_col(width = 0.6, color = "black") +
    geom_errorbar(
      aes(ymin = mean_prop - se_prop,
          ymax = mean_prop + se_prop),
      width = 0.13, size = 0.7
    ) +
    # 添加每个样本的散点，只包含lesion和non_lesion
    geom_jitter(
      data = ilc_filtered_plot,
      aes(x = Site_Status, y = ILC_Proportion),
      width = 0.08, height = 0,
      shape = 21, size = 2.5,
      fill = "black", color = "black",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = site_status_colors) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 0.07),
      expand = expansion(mult = c(0, 0.01))
    ) +
    xlab("") +
    ylab("ILC proportion") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      axis.title  = element_text(size = 14),
      plot.title  = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = "none"
    )
)
dev.off()

# 根据cluster绘制UMAP降维图，按seurat_clusters着色，并保存图片
png(file.path(out_dir, "fig1_DimPlot_by_cluster.png"), width = 1800, height = 1400, res = 200)
print(
  DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4) +
    ggtitle("UMAP by cluster")
)
dev.off()




anno_colors <- c(
  "#6197AD",
    "#E6674C",
  "#8887B4",
  "#B6D890",
  "#C176B0",
  "#8E5476",
  "#8C9FD1",
  "#F8A895"
)

names(anno_colors) <- sort(unique(obj$cell_anno))


png(
  file.path(out_dir, "figure1_DimPlot_cell_anno.png"),width = 2600, height = 2000, res = 400)
DimPlot(
  obj,
  group.by = "cell_anno",
  cols = anno_colors,
  label = TRUE,
   pt.size = 0.05,
  label.size = 5
) +
  ggtitle("Cell annotation") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

dev.off()

# 检查对象中是否有差异表达分析的结果（通常保存在obj@misc或obj@markers、obj@misc$de等）
# 下面代码会打印obj对象中和差异表达有关的元素名



#=========================

# 画各cell_anno在不同Site_Status（如Healthy/lesion/non_lesion）中的堆积柱状图

library(dplyr)
library(ggplot2)


# 确保obj有Site_Status这个变量（如果没有则根据Status和Site生成）
if (!("Site_Status" %in% colnames(obj@meta.data))) {
  obj@meta.data <- obj@meta.data %>%
    dplyr::mutate(
      Site_Status = dplyr::case_when(
        Status == "Healthy" ~ "Healthy",
        Status == "Psoriasis" & Site == "lesion" ~ "lesion",
        Status == "Psoriasis" & Site == "non_lesion" ~ "non_lesion",
        TRUE ~ "Unknown"
      )
    )
}


# 先统计每个Site_Status分组下不同cell_anno的细胞数量和比例
cell_prop_df <- obj@meta.data %>%
  group_by(Site_Status, cell_anno) %>%
  summarise(cell_count = n()) %>%
  ungroup() %>%
  group_by(Site_Status) %>%
  mutate(
    prop = cell_count / sum(cell_count)
  )

# 如有自定义cell_anno色板，可以按需设定；默认使用DistinctColors等
if (exists("anno_colors")) {
  fill_colors <- anno_colors
} else {
  cell_anno_lvls <- sort(unique(obj@meta.data$cell_anno))
  fill_colors <- scales::hue_pal()(length(cell_anno_lvls))
  names(fill_colors) <- cell_anno_lvls
}

png(file.path(out_dir, "fig1_cell_type_proportion_by_Site_Status.png"), width = 1400, height = 1500, res = 300)
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


# 绘制每个Tissue下不同cell_anno的比例堆积柱状图

tissue_prop_df <- obj@meta.data %>%
  group_by(Tissue, cell_anno) %>%
  summarise(cell_count = n()) %>%
  ungroup() %>%
  group_by(Tissue) %>%
  mutate(
    prop = cell_count / sum(cell_count)
  )

png(file.path(out_dir, "fig1_cell_type_proportion_by_Tissue.png"), width = 1400, height = 1500, res = 300)
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




table(obj$cell_anno )

ilc_cells <- colnames(obj)[obj$cell_anno == "ILC"]

# 提取ILC细胞的NT5E表达量
cd73_exp <- FetchData(obj, vars = "NT5E")[ilc_cells, , drop = FALSE]

# 按NT5E(即CD73)表达是否>0分组，统计ILC中CD73+与CD73-数量
cd73_pos_ilc <- ilc_cells[cd73_exp$NT5E > 0]
cd73_neg_ilc <- ilc_cells[cd73_exp$NT5E <= 0]

cat("CD73+ ILC 数量:", length(cd73_pos_ilc), "\n")
cat("CD73- ILC 数量:", length(cd73_neg_ilc), "\n")


# 统计 cell_anno 每个类型的细胞数
cellanno_counts <- table(obj$cell_anno)
print(cellanno_counts)

# 如需以数据框形式保存和查看
cellanno_counts_df <- as.data.frame(cellanno_counts)
colnames(cellanno_counts_df) <- c("cell_anno", "cell_count")
print(cellanno_counts_df)


colnames(obj@meta.data)
# 统计不同细胞类型在不同Tissue中的细胞数
celltype_tissue_table <- table(obj$cell_anno, obj$Tissue)
print(celltype_tissue_table)

# 如果需要更易读的数据框
celltype_tissue_df <- as.data.frame(celltype_tissue_table)
colnames(celltype_tissue_df) <- c("cell_anno", "Tissue", "cell_count")
print(celltype_tissue_df)
# 或者，如果Tissue是meta.data中的一列可以这样查看其所有取值
unique(obj@meta.data$Tissue)


# 绘制CD4、LCK、IL2RA、CD28的FeaturePlot并保存图片
genes_to_plot <- c("TRAC,CCR7,SELL,TCF7,IL7R,GPR183,LEF1, LTB")
for (gene in genes_to_plot) { 
  png(file.path(out_dir, paste0("fig1_FeaturePlot_", gene, ".png")), width = 1800, height = 1400, res = 200)
  print(FeaturePlot(obj, features = gene, max.cutoff = "q20") + ggtitle(paste0(gene, " FeaturePlot")))
  dev.off()
}


# 对不同cluster做差异表达分析

# 假定 Seurat object 的 cluster 信息 在 "seurat_clusters"
# 这里对每个cluster vs 其他所有cluster做差异分析
# 默认使用FindMarkers, 可以调整min.pct, logfc.threshold等参数

# 存结果的list
cluster_markers_list <- list()
clusters <- sort(unique(obj$seurat_clusters))
for (clust in clusters) {
  cat("Doing DE: cluster", clust, "vs rest\n")
  res <- FindMarkers(
    object = obj,
    ident.1 = clust,
    group.by = "seurat_clusters",
    logfc.threshold = 0.25,  # 常用默认阈值，根据需求可调
    min.pct = 0.1,
    only.pos = FALSE
  )
  cluster_markers_list[[as.character(clust)]] <- res
  # 可以保存到csv
  write.csv(res, file = file.path(out_dir, paste0("cluster", clust, "_markers.csv")))
}

# 如果需要把所有cluster的marker gene信息合并浏览，可以参考如下方法：
# 添加cluster列并合并所有，利于后续集中化筛选
markers_df <- do.call(rbind, lapply(names(cluster_markers_list), function(cl) {
  df <- cluster_markers_list[[cl]]
  df$cluster <- cl
  df$gene <- rownames(df)
  df
}))
write.csv(markers_df, file = file.path(out_dir, "all_clusters_markers_combined.csv"), row.names = FALSE)




# 绘制聚类（cluster）结果的UMAP图，并保存图片

png(file.path(out_dir, "fig1_umap_seurat_clusters.png"), width = 2000, height = 1600, res = 200)
print(
  DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 5) +
    ggtitle("UMAP by Seurat Clusters")
)
dev.off()

png(file.path(out_dir, "fig1_VlnPlot_nFeature_nCount_percentmt_by_clusters.png"), width = 2200, height = 1400, res = 200)
print(
  VlnPlot(
    obj, 
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
    group.by = "seurat_clusters", 
    pt.size = 0.05, 
    ncol = 3
  ) + 
    ggtitle("Quality Metrics by Seurat Cluster") +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
)
dev.off()


# 提取seurat_clusters为4、11、14的细胞，重新聚类分析




selected_clusters <- c(4, 11, 14)
selected_cells <- WhichCells(obj, idents = selected_clusters)  # 默认idents用seurat_clusters
subset_obj <- subset(obj, cells = selected_cells)

# 标准流程：归一化、筛选高变基因、标准化、PCA、邻接计算与聚类、UMAP
subset_obj <- NormalizeData(subset_obj)
subset_obj <- FindVariableFeatures(subset_obj, nfeatures = 4000)  # 增加高变基因数
subset_obj <- ScaleData(subset_obj)
subset_obj <- RunPCA(subset_obj, npcs = 30)

subset_obj <- FindNeighbors(subset_obj, dims = 1:30)
subset_obj <- FindClusters(subset_obj, resolution = 0.5)
subset_obj <- RunUMAP(subset_obj, dims = 1:30)


# 绘制IL17A和RORC在subset_obj中的FeaturePlot并保存图片

genes_to_plot <- c("IL17A", "FOXP3")
for (gene in genes_to_plot) {
  png(file.path(out_dir, paste0("subset_obj_FeaturePlot_", gene, ".png")), width = 1800, height = 1400, res = 200)
  print(
    FeaturePlot(subset_obj, features = gene, max.cutoff = "q20") +
      ggtitle(paste0(gene, " FeaturePlot (subset_obj)"))
  )
  dev.off()
}


#saveRDS(subset_obj, file = file.path(out_dir, "subset.rds"))


