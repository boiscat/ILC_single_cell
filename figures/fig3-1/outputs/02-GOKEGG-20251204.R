#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(glue)
  library(tibble)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

option_list <- list(
  make_option("--input", type = "character", default = "CD73_ILC_DE_results.csv",
              help = "差异表达结果csv，至少包含列 gene、avg_log2FC、p_val、p_val_adj"),
  make_option("--output_dir", type = "character", default = "CD73_ILC_GOKEGG_20251204",
              help = "输出文件夹（将自动创建figures和tables子目录）"),
  make_option("--up_log2fc", type = "double", default = 0.5,
              help = "定义上调基因的log2FC下限"),
  make_option("--up_threshold", type = "double", default = 0.05,
              help = "上调基因显著性阈值（可表示padj或p）"),
  make_option("--up_use_pval", type = "logical", default = TRUE,
              help = "TRUE=按p值筛选上调基因，FALSE=按padj"),
  make_option("--down_log2fc", type = "double", default = 0.25,
              help = "定义下调基因的log2FC下限（绝对值）"),
  make_option("--down_p", type = "double", default = 0.05,
              help = "定义下调基因的原始p值阈值（padj为1时退而求其次）"),
  make_option("--top_terms", type = "integer", default = 15,
              help = "富集结果图展示的最大通路数量"),
  make_option("--kegg_use_internal", type = "logical", default = FALSE,
              help = "是否使用clusterProfiler自带的KEGG注释以避免网络超时（需安装KEGG.db）")
)
opt <- parse_args(OptionParser(option_list = option_list))
use_p_for_up <- isTRUE(opt$up_use_pval)
use_internal_kegg <- isTRUE(opt$kegg_use_internal)

message(">>> 读取差异表达结果：", opt$input)
de <- read_csv(opt$input, show_col_types = FALSE)
required_cols <- c("gene", "avg_log2FC", "p_val", "p_val_adj")
stopifnot(all(required_cols %in% names(de)))

format_number <- function(x, digits = 2) {
  formatC(x, digits = digits, format = "f")
}

select_genes <- function(data, direction = c("up", "down"), log2fc, p_cutoff, use_padj = TRUE) {
  direction <- match.arg(direction)
  if (direction == "up") {
    filtered <- data %>%
      filter(avg_log2FC >= log2fc,
             if (use_padj) p_val_adj <= p_cutoff else p_val <= p_cutoff)
  } else {
    filtered <- data %>%
      filter(avg_log2FC <= -log2fc,
             if (use_padj) p_val_adj <= p_cutoff else p_val <= p_cutoff)
  }
  unique(filtered$gene)
}

up_symbols <- select_genes(
  de,
  direction = "up",
  log2fc = opt$up_log2fc,
  p_cutoff = opt$up_threshold,
  use_padj = !use_p_for_up
)

down_symbols <- select_genes(
  de,
  direction = "down",
  log2fc = opt$down_log2fc,
  p_cutoff = opt$down_p,
  use_padj = FALSE
)

message(">>> 上调基因候选：", length(up_symbols))
message(">>> 下调基因候选：", length(down_symbols))

if (length(up_symbols) < 5) {
  stop("上调基因数量过少，无法进行富集分析。")
}
if (length(down_symbols) < 5) {
  warning("下调基因数量低于5，部分富集图可能无法生成。")
}

map2entrez <- function(symbols) {
  if (length(symbols) == 0) {
    return(tibble())
  }
  bitr(
    symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db,
    drop = TRUE
  ) %>%
    distinct(SYMBOL, .keep_all = TRUE)
}

up_genes <- map2entrez(up_symbols)
down_genes <- map2entrez(down_symbols)

stopifnot(nrow(up_genes) >= 5)
if (nrow(down_genes) == 0) {
  warning("没有任何下调基因能映射到ENTREZ ID。")
}

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
plot_dir <- file.path(opt$output_dir, "figures")
table_dir <- file.path(opt$output_dir, "tables")
dir.create(plot_dir, showWarnings = FALSE)
dir.create(table_dir, showWarnings = FALSE)

fc_vec <- de$avg_log2FC
names(fc_vec) <- de$gene

run_go <- function(entrez_df) {
  if (nrow(entrez_df) < 5) {
    return(NULL)
  }
  enrichGO(
    gene = entrez_df$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    keyType = "ENTREZID",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
}

run_kegg <- function(entrez_df) {
  if (nrow(entrez_df) < 5) {
    return(NULL)
  }
  attempt <- function(use_internal) {
    enrichKEGG(
      gene = entrez_df$ENTREZID,
      organism = "hsa",
      pvalueCutoff = 0.1,
      use_internal_data = use_internal
    )
  }
  kk <- tryCatch(
    attempt(use_internal_kegg),
    error = function(e) {
      message(">>> enrichKEGG失败：", conditionMessage(e))
      if (use_internal_kegg) {
        return(NULL)
      }
      message(">>> 尝试改用use_internal_data=TRUE避免网络依赖")
      tryCatch(
        attempt(TRUE),
        error = function(e2) {
          message(">>> enrichKEGG二次失败：", conditionMessage(e2))
          NULL
        }
      )
    }
  )
  if (is.null(kk) || nrow(kk) == 0) {
    return(NULL)
  }
  setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

save_enrich_table <- function(enrich_obj, path) {
  if (is.null(enrich_obj) || nrow(enrich_obj) == 0) {
    message(">>> 无富集结果可写入：", path)
    return(invisible(NULL))
  }
  enrich_obj %>%
    as_tibble() %>%
    arrange(p.adjust) %>%
    write_csv(path)
  message(">>> 已保存：", path)
}

nature_theme <- theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    plot.subtitle = element_text(size = 12, hjust = 0),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.key = element_blank()
  )

plot_go_dot <- function(enrich_obj, title_text, subtitle_text) {
  if (is.null(enrich_obj) || nrow(enrich_obj) == 0) {
    return(NULL)
  }
  show_n <- min(opt$top_terms, nrow(enrich_obj))
  dotplot(enrich_obj, showCategory = show_n, font.size = 12) +
    scale_color_viridis_c(option = "magma", direction = -1, name = "Adjusted p") +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Gene Ratio",
      y = NULL
    ) +
    nature_theme
}

plot_kegg_dot <- function(enrich_obj, title_text, subtitle_text) {
  if (is.null(enrich_obj) || nrow(enrich_obj) == 0) {
    return(NULL)
  }
  show_n <- min(opt$top_terms, nrow(enrich_obj))
  dotplot(enrich_obj, showCategory = show_n, font.size = 12) +
    scale_color_viridis_c(option = "plasma", direction = -1, name = "Adjusted p") +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Gene Ratio",
      y = NULL
    ) +
    nature_theme
}

ggsave_plot <- function(plot_obj, base_path, width = 7, height = 6.5) {
  if (is.null(plot_obj)) {
    message(">>> 跳过未生成的图：", base_path)
    return(invisible(NULL))
  }
  png_path <- paste0(base_path, ".png")
  pdf_path <- paste0(base_path, ".pdf")
  ggsave(
    filename = png_path,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = 600
  )
  ggsave(
    filename = pdf_path,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in"
  )
  message(">>> 已保存：", png_path, " 和 ", pdf_path)
}

up_go <- run_go(up_genes)
down_go <- run_go(down_genes)
up_kegg <- run_kegg(up_genes)
down_kegg <- run_kegg(down_genes)

save_enrich_table(up_go, file.path(table_dir, "upregulated_GO_BP.csv"))
save_enrich_table(down_go, file.path(table_dir, "downregulated_GO_BP.csv"))
save_enrich_table(up_kegg, file.path(table_dir, "upregulated_KEGG.csv"))
save_enrich_table(down_kegg, file.path(table_dir, "downregulated_KEGG.csv"))

sign_label_up <- if (use_p_for_up) {
  glue("p≤{format_number(opt$up_threshold)}")
} else {
  glue("padj≤{format_number(opt$up_threshold)}")
}
sign_label_down <- glue("p≤{format_number(opt$down_p)}")

subtitle_up <- glue(
  "CD73+ vs CD73- ILC; {length(up_symbols)} genes (log2FC≥{format_number(opt$up_log2fc)}, {sign_label_up})"
)
subtitle_down <- glue(
  "CD73+ vs CD73- ILC; {length(down_symbols)} genes (log2FC≤-{format_number(opt$down_log2fc)}, {sign_label_down})"
)

ggsave_plot(
  plot_go_dot(up_go, "Upregulated CD73+ ILC programs (GO BP)", subtitle_up),
  file.path(plot_dir, "fig_up_GO_BP_dotplot")
)
ggsave_plot(
  plot_kegg_dot(up_kegg, "Upregulated CD73+ ILC pathways (KEGG)", subtitle_up),
  file.path(plot_dir, "fig_up_KEGG_dotplot")
)
ggsave_plot(
  plot_go_dot(down_go, "Downregulated CD73+ ILC programs (GO BP)", subtitle_down),
  file.path(plot_dir, "fig_down_GO_BP_dotplot")
)
ggsave_plot(
  plot_kegg_dot(down_kegg, "Downregulated CD73+ ILC pathways (KEGG)", subtitle_down),
  file.path(plot_dir, "fig_down_KEGG_dotplot")
)

metadata <- tibble(
  direction = c("up", "down"),
  genes_input = c(length(up_symbols), length(down_symbols)),
  mapped_genes = c(nrow(up_genes), nrow(down_genes)),
  log2fc_threshold = c(opt$up_log2fc, opt$down_log2fc),
  significance = c(sign_label_up, sign_label_down)
)
write_csv(metadata, file.path(table_dir, "threshold_summary.csv"))
message(">>> 分析完成：所有可用的表格与图保存在 ", opt$output_dir)
