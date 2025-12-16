#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(CellChat)
  library(dplyr)
  library(readr)
})

option_list <- list(
  make_option(c("--full_obj"), type = "character", help = "包含全部细胞的Seurat对象（带Tcell_subtype）"),
  make_option(c("--cd73_cells"), type = "character", help = "CD73阴性ILC细胞列表CSV"),
  make_option(c("--output_prefix"), type = "character", default = "outputs/objects/cd73neg_cellchat",
              help = "输出前缀（对象、表格）"),
  make_option(c("--group_col"), type = "character", default = "Site_Status",
              help = "分组列（lesion/non_lesion/healthy）"),
  make_option(c("--conditions"), type = "character", default = "lesion,non_lesion",
              help = "逗号分隔的条件列表"),
  make_option(c("--min_cells"), type = "integer", default = 10,
              help = "CellChat filterCommunication 的最少细胞数"),
  make_option(c("--tag"), type = "character", default = format(Sys.Date(), "%Y%m%d"))
)

opt <- parse_args(OptionParser(option_list = option_list))
required <- c("full_obj", "cd73_cells")
missing <- required[!nzchar(unlist(opt[required]))]
if (length(missing) > 0) {
  stop("缺少参数: ", paste(missing, collapse = ", "))
}
stopifnot(file.exists(opt$full_obj), file.exists(opt$cd73_cells))

conditions <- unique(trimws(strsplit(opt$conditions, ",")[[1]]))
conditions <- conditions[nchar(conditions) > 0]
stopifnot(length(conditions) > 0)

dir.create("outputs/objects", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)

message("读取对象: ", opt$full_obj)
full_obj <- readRDS(opt$full_obj)
meta <- full_obj@meta.data
stopifnot(opt$group_col %in% colnames(meta))
stopifnot("cell_anno" %in% colnames(meta))
stopifnot("Tcell_subtype" %in% colnames(meta))

cd73_cells <- read_csv(opt$cd73_cells, show_col_types = FALSE)$cell
stopifnot(length(cd73_cells) > 0)

CellChatDB <- CellChatDB.human

run_cellchat_condition <- function(condition_label) {
  message("处理条件: ", condition_label)
  cond_idx <- meta[[opt$group_col]] == condition_label
  cond_cells <- rownames(meta)[cond_idx]
  ilc_cells_cond <- intersect(cond_cells, cd73_cells)
  tcell_cells_cond <- cond_cells[meta[cond_cells, "major_annotation"] == "T_cells"]
  keep_cells <- c(ilc_cells_cond, tcell_cells_cond)
  keep_cells <- unique(keep_cells)
  if (length(ilc_cells_cond) < opt$min_cells) {
    warning("条件 ", condition_label, " 的 CD73- ILC 数量 < min_cells，跳过")
    return(NULL)
  }
  if (length(tcell_cells_cond) < opt$min_cells) {
    warning("条件 ", condition_label, " 的 T cell 数量 < min_cells，跳过")
    return(NULL)
  }
  cond_obj <- subset(full_obj, cells = keep_cells)
  cond_meta <- cond_obj@meta.data
  cond_meta$interaction_group <- ifelse(
    rownames(cond_meta) %in% ilc_cells_cond,
    "CD73neg_ILC",
    paste0("T_", ifelse(is.na(cond_meta$Tcell_subtype), "unknown", cond_meta$Tcell_subtype))
  )
  cond_obj@meta.data <- cond_meta
  message("Cell counts - CD73neg_ILC: ", sum(cond_meta$interaction_group == "CD73neg_ILC"),
          ", T cells: ", sum(startsWith(cond_meta$interaction_group, "T_")))
  cellchat <- createCellChat(object = cond_obj, group.by = "interaction_group")
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = opt$min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- rankNet(cellchat, mode = "comparison")
  out_rds <- sprintf("%s_%s_%s.rds", opt$output_prefix, condition_label, opt$tag)
  saveRDS(cellchat, out_rds)
  comm <- subsetCommunication(cellchat)
  out_comm <- sprintf("outputs/tables/cd73neg_cellchat_comm_%s_%s.csv", condition_label, opt$tag)
  write_csv(comm, out_comm)
  message("保存: ", out_rds)
  list(cellchat = cellchat, comm = out_comm, rds = out_rds)
}

results <- lapply(conditions, run_cellchat_condition)
names(results) <- conditions

session_path <- sprintf("outputs/objects/cd73neg_cellchat_sessionInfo_%s.txt", opt$tag)
writeLines(capture.output(sessionInfo()), con = session_path)
message("完成全部条件，SessionInfo: ", session_path)
