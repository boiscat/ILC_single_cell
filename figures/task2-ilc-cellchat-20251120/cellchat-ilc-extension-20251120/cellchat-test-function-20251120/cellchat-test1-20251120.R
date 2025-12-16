#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
})

args <- commandArgs(trailingOnly = TRUE)
tag <- format(Sys.Date(), "%Y%m%d")
if (length(args) > 0) {
  tag_arg <- grep("^--tag=", args, value = TRUE)
  if (length(tag_arg) == 1) {
    tag <- sub("^--tag=", "", tag_arg)
  }
}

base_dir <- file.path("todolist拆分-20251120", "任务2-ILC_CellChat-20251120", "cellchat_ilc扩展-20251120")
input_path <- file.path(base_dir, "inputs", sprintf("full_cellchat_ready_%s.rds", tag))
stopifnot(file.exists(input_path))

obj <- readRDS(input_path)
meta <- obj@meta.data

sel_cells <- function(mask, n = 25) {
  ids <- rownames(meta)[mask]
  if (length(ids) <= n) return(ids)
  sample(ids, n)
}

ilc_cells <- sel_cells(meta$cellchat_celltype == "ILC", n = 25)
tcell_cells <- sel_cells(grepl("^Tcell_", meta$cellchat_celltype), n = 25)
keep <- unique(c(ilc_cells, tcell_cells))
stopifnot(length(keep) >= 50)
obj_test <- subset(obj, cells = keep)

data("CellChatDB.human")

build_cellchat_test <- function(obj) {
  data.use <- GetAssayData(obj, assay = "RNA", slot = "data")
  meta.use <- obj@meta.data
  cc <- createCellChat(object = data.use, meta = meta.use, group.by = "cellchat_celltype")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc@options$population.size <- FALSE
  cc <- computeCommunProb(cc, population.size = FALSE)
  cc <- filterCommunication(cc, min.cells = 1, nonFilter.keep = TRUE)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

message(">>> 运行测试 CellChat，细胞数：", ncol(obj_test))
cc_test <- build_cellchat_test(obj_test)

out_dir <- file.path(base_dir, "cellchat-测试-功能-20251120", "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, sprintf("cellchat_test_subset_%s.rds", tag))
saveRDS(cc_test, out_path)
message("测试对象保存：", out_path)
