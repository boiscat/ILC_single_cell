#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(readr)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "全局对象RDS路径", metavar = "file"),
  make_option(c("-o", "--output_prefix"), type = "character", help = "输出前缀（不含扩展名）"),
  make_option(c("-m", "--meta_output"), type = "character", help = "T cell meta输出CSV"),
  make_option(c("-t", "--tag"), type = "character", default = format(Sys.Date(), "%Y%m%d"), help = "日期tag"),
  make_option(c("-c", "--cd73_threshold"), type = "double", default = 0.1, help = "CD73阴性阈值（NT5E表达）"),
  make_option(c("--site_col"), type = "character", default = "Site_Status", help = "病灶分组列名"),
  make_option(c("--celltype_col"), type = "character", default = "celltype", help = "ILC标注列"),
  make_option(c("--tcell_col"), type = "character", default = "Tcell_subtype", help = "T细胞亚型列"),
  make_option(c("--test_fraction"), type = "double", default = 0.02, help = "抽样比例用于测试输出"),
  make_option(c("--test_output"), type = "character", default = NA, help = "测试抽样RDS输出")
)

opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$input), file.exists(opt$input))
stopifnot(!is.null(opt$output_prefix))
stopifnot(!is.null(opt$meta_output))

message("读取对象: ", opt$input)
obj <- readRDS(opt$input)
meta <- obj@meta.data
meta$cell_id <- rownames(meta)
required_cols <- c(opt$site_col, opt$celltype_col, opt$tcell_col)
stopifnot(all(required_cols %in% colnames(meta)))
if (!"NT5E" %in% rownames(obj)) {
  stop("对象中缺少 NT5E 基因，无法根据 CD73 筛选")
}

ilc_labels <- c("ILC", "ILC1", "ILC2", "ILC3", "ILC_CD73pos", "ILC_CD73neg")
ilc_idx <- meta[[opt$celltype_col]] %in% ilc_labels
stopifnot(sum(ilc_idx) > 0)

# 计算 CD73（NT5E）表达量
cd73_expr <- FetchData(obj, vars = "NT5E")[, 1]
meta$CD73_expr <- cd73_expr

cd73_neg_idx <- ilc_idx & cd73_expr <= opt$cd73_threshold
cd73_neg_cells <- colnames(obj)[cd73_neg_idx]
stopifnot(length(cd73_neg_cells) > 0)

cd73neg_obj <- subset(obj, cells = cd73_neg_cells)
output_rds <- paste0(opt$output_prefix, "_cd73neg_ilc_", opt$tag, ".rds")
output_cells <- paste0(opt$output_prefix, "_cd73neg_ilc_cells_", opt$tag, ".csv")
write_rds(cd73neg_obj, output_rds)
write_csv(data.frame(cell = cd73_neg_cells), output_cells)
message("CD73阴性ILC数量: ", length(cd73_neg_cells))

# 输出 T cell meta
tcell_levels <- c("CD4", "CD8", "Treg", "MAIT", "Tex", "T_Naive", "T_CM", "T_EM")
tcell_meta <- meta %>%
  filter(.data[[opt$tcell_col]] %in% tcell_levels) %>%
  transmute(cell = .data$cell_id,
            Tcell_subtype = .data[[opt$tcell_col]],
            Site_Status = .data[[opt$site_col]])
write_csv(tcell_meta, opt$meta_output)

# 可选测试输出
if (!is.na(opt$test_output)) {
  set.seed(2025)
  sample_size <- max(50, ceiling(length(cd73_neg_cells) * opt$test_fraction))
  sample_cells <- sample(cd73_neg_cells, min(sample_size, length(cd73_neg_cells)))
  write_rds(subset(cd73neg_obj, cells = sample_cells), opt$test_output)
}

session_info_path <- paste0(opt$output_prefix, "_sessionInfo_", opt$tag, ".txt")
writeLines(capture.output(sessionInfo()), con = session_info_path)
message("完成：", output_rds)
