#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
})

get_project_root <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  file_arg <- grep(file_flag, args_full, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub(file_flag, "", file_arg))
  } else {
    script_path <- normalizePath(getwd())
  }
  script_dir <- dirname(script_path)
  normalizePath(file.path(script_dir, "..", "..", ".."))
}

project_root <- get_project_root()
pipeline_root <- file.path(project_root, "task6-cd73pos-ilc-interaction-20251120", "cd73pos-ilc-interaction-pipeline-20251120")
obj_path <- file.path(pipeline_root, "outputs", "objects", "CD73_ILC_Tcell_Healthy_20251120.rds")
stopifnot(file.exists(obj_path))

obj <- readRDS(obj_path)
meta_cols <- c("interaction_group", "Site_Status")
stopifnot(all(meta_cols %in% colnames(obj@meta.data)))
stopifnot("RNA" %in% names(obj@assays))

set.seed(20251123)
group_col <- "interaction_group"
meta <- obj@meta.data
split_ids <- split(rownames(meta), meta[[group_col]])
keep <- unlist(lapply(split_ids, function(ids) {
  if (length(ids) <= 20) {
    return(ids)
  }
  sample(ids, 20)
}))
mini_obj <- subset(obj, cells = keep)
stopifnot(ncol(mini_obj) >= 40)
stopifnot(length(unique(mini_obj@meta.data[[group_col]])) >= 2)

data("CellChatDB.human")
cc <- createCellChat(
  object = GetAssayData(mini_obj, assay = "RNA", slot = "data"),
  meta = mini_obj@meta.data,
  group.by = group_col
)
cc@DB <- CellChatDB.human
cc@options$population.size <- FALSE
cc <- subsetData(cc)
cc <- identifyOverExpressedGenes(cc)
cc <- identifyOverExpressedInteractions(cc)
cc <- computeCommunProb(cc, raw.use = TRUE, population.size = FALSE)
cc <- filterCommunication(cc, min.cells = 1)
cc <- computeCommunProbPathway(cc)
cc <- aggregateNet(cc)
df <- subsetCommunication(cc)
stopifnot(nrow(df) > 0)
cat("CD73ILC-cellchat-test1-20251120 completed, interactions:", nrow(df), "\n")
