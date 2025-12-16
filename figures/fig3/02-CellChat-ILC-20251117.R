#!/usr/bin/env Rscript
# 目标：用CellChat比较Healthy / non_lesion / lesion间的ILC相关通讯，突出腺苷节点（CD39/ENTPD1, CD73/NT5E, ADORA1/2A/2B/3）

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(future)
  library(ggplot2)
  library(patchwork)
})

message(">>> 载入Seurat对象")
seurat_path <- "full_obj_saved_1117.rds"
obj <- readRDS(seurat_path)

meta <- obj@meta.data
required_cols <- c("Site_Status", "major_annotation", "cell_anno", "full_clustering")
stopifnot(all(required_cols %in% colnames(meta)))

# 只保留有分组信息的细胞
obj <- subset(obj, subset = !is.na(Site_Status))

# 控制规模：每个状态、每个major_annotation最多保留2000细胞，加快计算
set.seed(123)
downsample_per_group <- function(seurat_obj, max_cells = 2000) {
  meta <- seurat_obj@meta.data
  keep <- unlist(lapply(split(rownames(meta), meta$cell_anno), function(ids) {
    if (length(ids) <= max_cells) return(ids)
    sample(ids, max_cells)
  }))
  subset(seurat_obj, cells = keep)
}

# CellChat数据库（人）
data("CellChatDB.human")
cellchat_db <- CellChatDB.human

plan("sequential")

# 允许控制人口规模调整
use_population_size <- TRUE

# 只保留ILC及与ILC互作的意义：在下游画图时过滤 edges 涉及 ILC（发送或接收），此处不过滤，以免丢失配对细胞类型
# 保持群体标签为 major_annotation

# 辅助函数：对单个状态运行CellChat全流程
run_cellchat <- function(seurat_obj, status_label) {
  message(">>> 构建 CellChat: ", status_label)
  data.use <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  meta.use <- seurat_obj@meta.data
  # 创建 CellChat 对象
  cc <- createCellChat(object = data.use, meta = meta.use, group.by = "cell_anno")
  cc@DB <- cellchat_db
  cc <- subsetData(cc)                       # 预过滤
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, raw.use = TRUE, population.size = use_population_size)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc@options$levels <- levels(factor(meta.use$major_annotation))
  return(cc)
}

status_levels <- c("Healthy", "non_lesion", "lesion")
cellchat_list <- list()
for (st in status_levels) {
  sub_obj <- subset(obj, subset = Site_Status == st)
  if (ncol(sub_obj) < 100) {
    warning("分组 ", st, " 细胞数过少：", ncol(sub_obj))
  }
  sub_obj <- downsample_per_group(sub_obj, max_cells = 2000)
  cellchat_list[[st]] <- run_cellchat(sub_obj, st)
}

saveRDS(cellchat_list, file = "outputs/cellchat_list_ILC_allstatus.rds")
message(">>> 保存 cellchat 列表：outputs/cellchat_list_ILC_allstatus.rds")

# 整理腺苷相关通讯
adenosine_genes <- c("ENTPD1", "NT5E", "ADORA1", "ADORA2A", "ADORA2B", "ADORA3")
pattern <- paste(adenosine_genes, collapse = "|")

extract_comm <- function(cc, status_label) {
  df <- subsetCommunication(cc)
  df <- df[
    grepl(pattern, df$interaction_name, ignore.case = TRUE) |
      grepl(pattern, df$ligand, ignore.case = TRUE) |
      grepl(pattern, df$receptor, ignore.case = TRUE),
    ]
  if (nrow(df) == 0) {
    return(df)
  }
  df$Site_Status <- status_label
  df$pair <- paste(df$source, "→", df$target)
  df
}

comm_list <- lapply(names(cellchat_list), function(nm) {
  extract_comm(cellchat_list[[nm]], nm)
})
comm_list <- comm_list[sapply(comm_list, nrow) > 0]

if (length(comm_list) == 0) {
  stop("未在CellChat结果中找到腺苷相关通讯")
}

comm_tbl <- do.call(rbind, comm_list)

# 全部状态合并的 CellChat（便于筛选 ILC 全局交互）
message(">>> 构建 CellChat: All")
obj_all <- downsample_per_group(obj, max_cells = 2000)
cellchat_all <- run_cellchat(obj_all, "All")
saveRDS(cellchat_all, file = "outputs/cellchat_total_ILC_all.rds")
message(">>> 保存总体 cellchat 对象：outputs/cellchat_total_ILC_all.rds")

all_comm <- subsetCommunication(cellchat_all)
all_comm$pair <- paste(all_comm$source, "→", all_comm$target)
write.csv(all_comm, file = "outputs/CellChat_total_all_interactions.csv", row.names = FALSE)

ilc_comm <- all_comm[all_comm$source == "ILC" | all_comm$target == "ILC", ]
if (nrow(ilc_comm) > 0) {
  write.csv(ilc_comm, file = "outputs/CellChat_total_ILC_interactions.csv", row.names = FALSE)
  ilc_comm$neglogP <- -log10(pmax(ilc_comm$pval, 1e-6))
  ilc_comm$interaction_name <- factor(ilc_comm$interaction_name, levels = sort(unique(ilc_comm$interaction_name)))
  ilc_comm$pair <- factor(ilc_comm$pair, levels = sort(unique(ilc_comm$pair)))
  pdf("outputs/CellChat_total_ILC_interactions_bubble.pdf", width = 9, height = 5)
  print(
    ggplot(ilc_comm, aes(x = pair, y = interaction_name)) +
      geom_point(aes(size = neglogP, fill = prob), shape = 21, color = "grey30", alpha = 0.85) +
      scale_fill_gradientn(colors = c("#2E86AB", "#F3DF6C", "#D7263D"), name = "Commun. prob") +
      scale_size(range = c(2, 7), name = "-log10(p)") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Cell pair (含ILC)", y = "Ligand-Receptor", title = "全部样本：ILC相关通讯")
  )
  dev.off()
} else {
  message("All状态下未检测到涉及ILC的通讯")
}

adenosine_all <- extract_comm(cellchat_all, "All")
if (nrow(adenosine_all) > 0) {
  write.csv(adenosine_all, "outputs/CellChat_total_ILC_adenosine.csv", row.names = FALSE)
}

write.csv(comm_tbl, file = "outputs/CellChat_ILC_adenosine_communications.csv", row.names = FALSE)

healthy_tbl <- comm_tbl[comm_tbl$Site_Status == "Healthy", ]
comp_tbl <- comm_tbl[comm_tbl$Site_Status != "Healthy", ]

diff_tbl <- merge(
  comp_tbl,
  healthy_tbl[, c("interaction_name", "pair", "prob", "pval")],
  by = c("interaction_name", "pair"),
  all.x = TRUE,
  suffixes = c("", "_Healthy")
)

diff_tbl$prob_Healthy[is.na(diff_tbl$prob_Healthy)] <- 0
diff_tbl$pval_Healthy[is.na(diff_tbl$pval_Healthy)] <- 1
diff_tbl$prob_diff <- diff_tbl$prob - diff_tbl$prob_Healthy
diff_tbl$direction <- ifelse(diff_tbl$prob_diff >= 0, "Increased signaling", "Decreased signaling")
diff_tbl$direction <- factor(diff_tbl$direction, levels = c("Increased signaling", "Decreased signaling"))
diff_tbl$status_vs <- paste0(diff_tbl$Site_Status, " vs Healthy")
diff_tbl$status_vs <- factor(diff_tbl$status_vs, levels = c("non_lesion vs Healthy", "lesion vs Healthy"))
diff_tbl$neglogP <- -log10(pmax(diff_tbl$pval, 1e-6))
diff_tbl$abs_prob_diff <- abs(diff_tbl$prob_diff)
diff_tbl$signif_cat <- ifelse(diff_tbl$pval < 0.01, "p < 0.01",
                              ifelse(diff_tbl$pval < 0.05, "0.01 ≤ p < 0.05", NA))
diff_tbl$signif_cat <- factor(diff_tbl$signif_cat, levels = c("0.01 ≤ p < 0.05", "p < 0.01"))
diff_tbl$pair_status <- paste0(diff_tbl$pair, "\n", diff_tbl$status_vs)
diff_tbl$rank_within <- ave(
  diff_tbl$abs_prob_diff,
  interaction(diff_tbl$status_vs, diff_tbl$direction),
  FUN = function(x) rank(-x, ties.method = "first")
)

plot_df <- diff_tbl[!is.na(diff_tbl$signif_cat), ]
if (nrow(plot_df) > 0) {
  plot_df <- plot_df[order(plot_df$status_vs, plot_df$direction, plot_df$pair), ]
  plot_df$interaction_name <- factor(plot_df$interaction_name, levels = rev(sort(unique(plot_df$interaction_name))))
  plot_df$pair_status <- factor(plot_df$pair_status, levels = unique(plot_df$pair_status))
  pdf("outputs/CellChat_ILC_adenosine_diff_bubble.pdf", width = 12, height = 6.5)
  print(
    ggplot(plot_df, aes(x = pair_status, y = interaction_name)) +
      geom_point(
        aes(size = signif_cat, fill = prob),
        shape = 21, color = "grey20", alpha = 0.9, stroke = 0.4
      ) +
      scale_fill_gradientn(colors = c("#2E86AB", "#F3DF6C", "#D7263D"), name = "Commun. prob") +
      scale_size_manual(
        values = c("0.01 ≤ p < 0.05" = 3.2, "p < 0.01" = 4.8),
        name = "P value"
      ) +
      facet_grid(status_vs ~ direction, scales = "free", space = "free") +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        panel.grid.major = element_line(color = "grey90"),
        strip.text = element_text(face = "bold", size = 11)
      ) +
      labs(
        x = "发送→接收细胞（含分组信息）",
        y = "Ligand-Receptor",
        title = "ILC腺苷信号差异：加强 vs 减弱",
        subtitle = "只展示p<0.05的交互；颜色=通讯概率，点大小=显著性"
      )
  )
  dev.off()
} else {
  message("p>=0.05，无法生成差异气泡图")
}

write.csv(diff_tbl, "outputs/CellChat_ILC_adenosine_diff_table.csv", row.names = FALSE)

ranked_tbl <- diff_tbl[diff_tbl$rank_within <= 10,
                       c("status_vs", "direction", "interaction_name", "pair",
                         "prob", "prob_Healthy", "prob_diff", "abs_prob_diff",
                         "pval", "signif_cat", "rank_within")]
write.csv(ranked_tbl,
          "outputs/CellChat_ILC_adenosine_diff_top10.csv",
          row.names = FALSE)

message("完成 CellChat 差异气泡图生成")
