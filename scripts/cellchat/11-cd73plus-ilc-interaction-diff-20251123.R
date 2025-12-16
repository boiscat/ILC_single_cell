#!/usr/bin/env Rscript

# 11-cd73plus-ilc-interaction-diff-20251123.R
# 目的：统计/比较 CD73+ ILC 相关的 CellChat 互作差异，并输出图表

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(future)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

obj_path <- "/lustre/home/yangyb/youmz/data/figures/todolist拆分-20251120/任务1-lesion修正-20251120/lesion修正-功能复跑-20251120/outputs/objects/full_obj_lesionFixed_20251120.rds"

output_dir <- "/lustre/home/yangyb/youmz/data/figures/todolist拆分-20251120/任务5-ilc+-交互差异"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("读取对象：", obj_path)
obj <- readRDS(obj_path)




CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB # 默认使用全部数据库


# 提取ILC细胞
ilc_meta <- obj@meta.data %>%
  tibble::rownames_to_column(var = "cell_barcode") %>%
  dplyr::filter(cell_anno == "ILC")
if (nrow(ilc_meta) == 0) {
  stop("对象中不存在 ILC 细胞，无法继续")
}
ilc_obj <- subset(obj, cells = ilc_meta$cell_barcode)

# 提取 ILC 中 NT5E (CD73) 阳性细胞
rna_data <- SeuratObject::GetAssayData(ilc_obj[["RNA"]], layer = "data")
nt5e_expr <- as.numeric(rna_data["NT5E", colnames(rna_data), drop = TRUE])

ilc_meta <- ilc_meta %>%
  mutate(NT5E_expr = nt5e_expr,
         CD73_positive = NT5E_expr > 0)

pos_ilc_barcodes <- rownames(ilc_meta %>% filter(CD73_positive))
neg_ilc_barcodes <- rownames(ilc_meta %>% filter(!CD73_positive))


other_cells_barcodes <- rownames(obj@meta.data %>% filter(cell_anno != "ILC"))
# 列表 A：CD73+ ILCs + 所有其他细胞
cells_for_pos_group <- c(pos_ilc_barcodes, other_cells_barcodes)
# 列表 B：CD73- ILCs + 所有其他细胞
cells_for_neg_group <- c(neg_ilc_barcodes, other_cells_barcodes)

# === 运行 Group 1: CD73 Positive ===
# 提取子集
seurat_pos <- subset(obj, cells = cells_for_pos_group)

cc_pos@DB <- CellChatDB.use



# 创建 CellChat 对象
cc_pos <- createCellChat(object = seurat_pos, group.by = "cell_anno") 
cc_pos <- subsetData(cc_pos)
# 标准流程 (这里简化写，你需要跑完完整流程)
cc_pos <- cc_pos %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions() %>% 
  computeCommunProb() %>% 
  filterCommunication() %>% 
  computeCommunProbPathway() %>% 
  aggregateNet()



# === 运行 Group 2: CD73 Negative ===
# 提取子集
seurat_neg <- subset(obj, cells = cells_for_neg_group)

# 创建 CellChat 对象
cc_neg <- createCellChat(object = seurat_neg, group.by = "cell_anno")
# 赋值数据库
cc_neg@DB <- CellChatDB.use

# 跑基本处理流程
cc_neg <- subsetData(cc_neg)
cc_neg <- cc_neg %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions() %>% 
  computeCommunProb() %>% 
  filterCommunication() %>% 
  computeCommunProbPathway() %>% 
  aggregateNet()




# 1. 合并对象
object.list <- list(CD73_Pos = cc_pos, CD73_Neg = cc_neg)
cellchat <- mergeCellChat(object.list, 
                          add.names = names(object.list),
                          cell.prefix = TRUE) # <-- 关键修正



cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")

cellchat@net$diff$prob
cellchat@net$diff$weight


# 总体数量差异
png(filename = file.path(output_dir, "diffInteraction_count.png"), width = 1200, height = 1000, res = 150)
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
dev.off()

# 总体强度差异
png(filename = file.path(output_dir, "diffInteraction_weight.png"), width = 1200, height = 1000, res = 150)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
dev.off()

# 通路热图
png(filename = file.path(output_dir, "diffInteraction_heatmap.png"), width = 1400, height = 1200, res = 150)
netVisual_heatmap(cellchat)
dev.off()



library(CellChat)

# --- 1. 比较总体互作数量和强度的差异 ---
# 这部分您已运行，作为参考保留
# par(mfrow = c(1, 2), xpd=TRUE)
# netVisual_diffInteraction(cellchat, weight.scale = T)
# netVisual_diffInteraclibrary(CellChat)

# --- 2. 通路贡献差异比较 (Pathogenicity Score) ---
# 找出在 Pos 和 Neg 组中，哪些信号通路的总强度变化最大


netVisual_heatmap(cellchat, color.heatmap = "RdBu")

summary(cellchat@net$diff$prob)
