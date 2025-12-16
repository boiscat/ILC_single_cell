

obj <- readRDS("/lustre/home/yangyb/youmz/data/outputs/healthy_pso_merge_20251110/healthy_pso_harmony_integrated_de_20251110.rds")


 table(obj@meta.data$cell_anno )
# 对cell_anno中的特定簇重命名：簇5改为"ILG"，簇6改为"NK"
# 假定cell_anno为obj@meta.data$cell_anno
# 首先确认簇的命名方式（字符还是数字），如果是因子则先转为字符

# 将seurat_clusters为5和6的细胞的cell_anno重命名
if ("cell_anno" %in% colnames(obj@meta.data) && "seurat_clusters" %in% colnames(obj@meta.data)) {
  ann <- obj@meta.data$cell_anno
  clust <- obj@meta.data$seurat_clusters
  # 如果cell_anno是因子，则转为字符
  if (is.factor(ann)) ann <- as.character(ann)
  # 如果seurat_clusters是因子，则转为字符/数值
  if (is.factor(clust)) clust <- as.character(clust)
  # 处理数字型和字符型cluster一致性
  # 对于seurat_clusters为5和6分别赋予新名字
  ann[clust == "5" | clust == 5] <- "ILG"
  ann[clust == "6" | clust == 6] <- "NK"
  obj@meta.data$cell_anno <- factor(ann)
}

# cell_anno的T/NK改为T
if ("cell_anno" %in% colnames(obj@meta.data)) {
  # 提取cell_anno列，如果是factor先转为字符向量
  ann <- obj@meta.data$cell_anno
  if (is.factor(ann)) ann <- as.character(ann)
  # 把"TNK"或"T/NK"等（包含T和NK的标签）统一改为"T"
  ann[ann == "T/NK" | ann == "TNK"] <- "T"
    ann[ann == "ILG"] <- "ILC"
  obj@meta.data$cell_anno <- factor(ann)
}



# 去掉细胞注释为 "VE" 或 "LE" 的细胞
obj <- subset(obj, subset = !(cell_anno %in% c("VE", "LE","Melanocytes","Pericyte")))


out_dir <- "/lustre/home/yangyb/youmz/data/figures/fig1"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# 计算总细胞数
cat("Total cells:", ncol(obj), "\n")







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


# DotPlot展示，按 cell_anno 聚类分组（如需更换分组可修改 group.by 参数）
png(file.path(out_dir, "fig1_marker_dotplot_by_cell_anno.png"), width = 2000, height = 900, res = 200)
DotPlot(
  obj, 
  features = marker_genes, 
  group.by = "cell_anno", 
  cols = c("lightgrey", "#E41A1C"),
  dot.scale = 7
) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 14),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  ) +
  ggtitle("DotPlot of key marker genes by cell_anno")
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


