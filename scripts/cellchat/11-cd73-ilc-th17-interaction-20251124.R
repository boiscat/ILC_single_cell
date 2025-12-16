#!/usr/bin/env Rscript

# 11-cd73-ilc-th17-interaction-20251124.R
# 目的：分析 CD73+ ILC 和 Th17 细胞的 CellChat 交互

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(future)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

obj_path <- "/lustre/home/yangyb/youmz/data/figures/todolist拆分-20251120/任务1-lesion修正-20251120/lesion修正-功能复跑-20251120/outputs/objects/full_obj_lesionFixed_20251120.rds"

output_dir <- "/lustre/home/yangyb/youmz/data/figures/todolist拆分-20251120/任务4-ilc-th17交互"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("读取对象：", obj_path)
obj <- readRDS(obj_path)

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

# 仅选取 ILC 中 CD73+（NT5E高表达）细胞
ilc_meta <- obj@meta.data %>%
  tibble::rownames_to_column(var = "cell_barcode") %>%
  dplyr::filter(cell_anno == "ILC")
if (nrow(ilc_meta) == 0) {
  stop("对象中不存在 ILC 细胞，无法继续")
}
ilc_obj <- subset(obj, cells = ilc_meta$cell_barcode)
rna_data <- SeuratObject::GetAssayData(ilc_obj[["RNA"]], layer = "data")
nt5e_expr <- as.numeric(rna_data["NT5E", colnames(rna_data), drop = TRUE])
ilc_meta <- ilc_meta %>%
  mutate(NT5E_expr = nt5e_expr,
         CD73_positive = NT5E_expr > 0)


pos_ilc_barcodes <- (ilc_meta %>% filter(CD73_positive))$cell_barcode

if (length(pos_ilc_barcodes) == 0) {
  stop("未检测到 CD73+ ILC 细胞")
}

# 提取 FOXP3 阳性细胞
foXP3_expr <- SeuratObject::GetAssayData(obj[["RNA"]], layer = "data")["FOXP3", , drop = TRUE]
foXP3_positive_barcodes <- names(foXP3_expr)[foXP3_expr > 0]
if (length(foXP3_positive_barcodes) == 0) {
  stop("未检测到 FOXP3 阳性细胞")
}
th17_barcodes <- foXP3_positive_barcodes

# 只保留 CD73+ ILC 和 Th17 两类细胞，检查barcode数量差异
selected_barcodes <- c(pos_ilc_barcodes, th17_barcodes)
selected_barcodes <- unique(selected_barcodes)  # 移除重复barcodes

message("CD73+ ILC 细胞数: ", length(pos_ilc_barcodes))
message("Th17 细胞数: ", length(th17_barcodes))
message("合并后 unique barcode 数: ", length(selected_barcodes))
message("原始 barcode 总数: ", length(unique(c(pos_ilc_barcodes, th17_barcodes))))
# 检查总数为何小于单项之和
overlap_barcodes <- intersect(pos_ilc_barcodes, th17_barcodes)
message("CD73+ ILC 与 Th17 细胞重叠的 barcode 数量: ", length(overlap_barcodes))

obj2 <- subset(obj, cells = selected_barcodes)
message("obj2 subset 后细胞数: ", ncol(obj2))

# 重新定义分组
obj2$cellchat_group <- NA
obj2$cellchat_group[colnames(obj2) %in% th17_barcodes] <- "Th17"
obj2$cellchat_group[colnames(obj2) %in% pos_ilc_barcodes] <- "CD73pos_ILC"


obj2 <- subset(obj2, cells = colnames(obj2)[!is.na(obj2$cellchat_group)])

# 检查分组情况
# CellChat 分析
message("创建 CellChat 对象...")
cc_obj <- createCellChat(object = obj2, group.by = "cellchat_group")
CellChatDB <- CellChatDB.human
cc_obj@DB <- CellChatDB
cc_obj <- subsetData(cc_obj)

message("运行 CellChat 分析流程...")
future::plan("sequential")

cc_obj <- identifyOverExpressedGenes(cc_obj, 
                                     thresh.p = 0.05,        # P值阈值
                                     thresh.fc = 0.1)        # 表达倍数差异阈值


cc_obj <- identifyOverExpressedInteractions(cc_obj)
cc_obj <- computeCommunProb(cc_obj, raw.use = TRUE, population.size = FALSE)
cc_obj <- filterCommunication(cc_obj, min.cells = 3)
cc_obj <- computeCommunProbPathway(cc_obj)
cc_obj <- aggregateNet(cc_obj)

# 保存分析对象
save_path <- file.path(output_dir, "cellchat_cd73posILC_th17_20251122.rds")
saveRDS(cc_obj, file = save_path)
message("CellChat 对象已保存：", save_path)

# 提取通讯结果
comm_df <- subsetCommunication(cc_obj)
comm_path <- file.path(output_dir, "cellchat_cd73posILC_th17_communications_20251122.csv")
write.csv(comm_df, comm_path, row.names = FALSE)
message("通讯结果已保存：", comm_path)



# --- 1. 确认 CellChat 对象中的细胞分组 ---
# 确保分组名是正确的，并且存在于 cc_obj@DB 中
sources.id <- "CD73pos_ILC"
targets.id <- "Th17"

comm_df <- subsetCommunication(cc_obj)


message(paste("当前分析的 Source Group:", sources.id))
message(paste("当前分析的 Target Group:", targets.id))
# --- 2. 气泡图 (Bubble Plot): 展示具体的配体-受体对并保存 ---
# 统计 ILC(CD73+) -> Th17 所有配体-受体对（表格文本），直接输出前几行并存文件



# 气泡图是最常用且信息量最大的图，用于展示 ILC(CD73+) 发出的具体分子信号。
# 关键修复：移除 signaling.filter，并使用 thresh = 0 强制绘制所有互作。
p1 <- netVisual_bubble(cc_obj, 
                       sources.use = sources.id,       # 指定信号发出细胞
                       targets.use = targets.id,       # 指定信号接收细胞
                       thresh = 0,                     # <-- 修正：强制包含所有互作（非显著也包含）
                       remove.isolate = FALSE,         
                       angle.x = 90) +                 
  labs(title = paste0("ILC(CD73+) -> Th17(Foxp3+) Specific L-R Interactions (All Inferred)")) +
  theme(plot.title = element_text(hjust = 0.5))

print(p1)
bubble_plot_path <- file.path(output_dir, "cellchat_bubble_cd73posILC_th17_20251122.pdf")
ggsave(filename = bubble_plot_path, plot = p1, width = 9, height = 7)
message("气泡图已保存：", bubble_plot_path)


# --- 3. 和弦图 (Chord Diagram): 展示信号通路的总强度并保存 ---

p2 <- netVisual_chord_gene(cc_obj, 
                           sources.use = sources.id, 
                           targets.use = targets.id, 
                           signaling.filter = FALSE,   # <-- 强制绘制所有互作
                           lab.cex = 1.0,               
                           legend.pos.x = 20,
                           title.name = paste0("Pathway Contribution: ILC(CD73+) to Th17(Foxp3+)"))

print(p2)
chord_plot_path <- file.path(output_dir, "cellchat_chord_cd73posILC_th17_20251122.pdf")
ggsave(filename = chord_plot_path, plot = p2, width = 8, height = 8)
message("和弦图已保存：", chord_plot_path)


# --- 4. 通路气泡图 (Pathway Bubble Plot): 看哪些通路最重要并保存 ---

p3 <- plot_pathway_bubble(cc_obj, 
                          sources.use = sources.id, 
                          targets.use = targets.id,
                          title = "Key Signaling Pathways in ILC(CD73+) <-> Th17(Foxp3+) Communication")
print(p3)
pathway_bubble_path <- file.path(output_dir, "cellchat_pathway_bubble_cd73posILC_th17_20251122.pdf")
ggsave(filename = pathway_bubble_path, plot = p3, width = 10, height = 7)
message("通路气泡图已保存：", pathway_bubble_path)


# 假设我们想看 MIF 通路的详细互作图，并保存
pathway.name <- "MIF" # 替换为你感兴趣的通路名称
if (pathway.name %in% cc_obj@netP$pathways) {
  p4 <- netVisual_aggregate(cc_obj, 
                            signaling = pathway.name, 
                            sources.use = sources.id, 
                            targets.use = targets.id,
                            layout = "chord")
  print(p4)
  pathway_chord_path <- file.path(output_dir, paste0("cellchat_", pathway.name, "_chord_cd73posILC_th17_20251122.pdf"))
  ggsave(filename = pathway_chord_path, plot = p4, width = 7, height = 7)
  message(paste0("已绘制并保存特定通路 (", pathway.name, ") 的和弦图: ", pathway_chord_path))
} else {
  message(paste0("通路 '", pathway.name, "' 在当前通讯中不活跃，无法绘制特定通路图。"))
}


