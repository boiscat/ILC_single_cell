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
  make_option("--de_input", type = "character", default = "outputs/CD73_ILC_DE_results.csv",
              help = "差异表达结果csv（来自01脚本），至少包含列 gene、avg_log2FC、p_val、p_val_adj"),
  make_option("--output_dir", type = "character", default = "outputs/CD73_ILC_GO_selected_20251214",
              help = "输出文件夹（将自动创建figures和tables子目录）"),
  make_option("--down_mode", type = "character", default = "all",
              help = "下调基因集选择方式：all=所有下调基因（按阈值过滤）；top=下调TopN（用于和柱状图一致）"),
  make_option("--down_n", type = "integer", default = 30,
              help = "当--down_mode=top时使用：下调TopN"),
  make_option("--down_log2fc", type = "double", default = 0.25,
              help = "当--down_mode=all时使用：|log2FC|阈值（下调为avg_log2FC<=-down_log2fc）"),
  make_option("--down_p", type = "double", default = 0.05,
              help = "当--down_mode=all且--down_use_padj=FALSE时使用：p值阈值"),
  make_option("--down_padj", type = "double", default = 0.05,
              help = "当--down_mode=all且--down_use_padj=TRUE时使用：padj阈值"),
  make_option("--down_use_padj", type = "logical", default = FALSE,
              help = "当--down_mode=all时：TRUE=用padj筛选下调；FALSE=用原始p值筛选"),
  make_option("--down_require_signif", type = "logical", default = TRUE,
              help = "当--down_mode=all时：TRUE=必须满足p/padj阈值；FALSE=仅按avg_log2FC方向筛选（会非常宽松）"),
  make_option("--ont", type = "character", default = "BP",
              help = "GO本体：BP/CC/MF"),
  make_option("--pvalue_cutoff", type = "double", default = 0.05,
              help = "GO富集pvalueCutoff"),
  make_option("--qvalue_cutoff", type = "double", default = 0.2,
              help = "GO富集qvalueCutoff"),
  make_option("--min_genes", type = "integer", default = 5,
              help = "做GO的最小基因数（不足则仅导出基因/映射，不出GO图表）"),
  make_option("--top_terms", type = "integer", default = 15,
              help = "GO dotplot展示的最大条目数"),
  make_option("--use_universe", type = "logical", default = FALSE,
              help = "是否使用DE表中的所有基因作为背景(universe)做GO；FALSE则使用全基因背景（更容易得到结果，但背景更宽）"),
  make_option("--up_genes_file", type = "character", default = NA,
              help = "可选：指定上调基因列表文件（每行一个symbol）；不提供则使用脚本内置30基因")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!opt$ont %in% c("BP", "CC", "MF")) {
  stop("--ont只能是 BP/CC/MF")
}
if (!opt$down_mode %in% c("all", "top")) {
  stop("--down_mode只能是 all 或 top")
}

normalize_gene <- function(x) gsub("[^A-Za-z0-9_.]", "", x)

message(">>> 读取差异结果：", opt$de_input)
de <- read_csv(opt$de_input, show_col_types = FALSE)
required_cols <- c("gene", "avg_log2FC", "p_val", "p_val_adj")
stopifnot(all(required_cols %in% names(de)))

de_tbl <- de %>%
  mutate(gene_norm = normalize_gene(gene)) %>%
  arrange(desc(avg_log2FC))

select_top <- function(tbl, direction, n) {
  is_up <- direction == "Up"
  dir_tbl <- tbl %>%
    filter(if (is_up) avg_log2FC > 0 else avg_log2FC < 0)

  arrange_dir <- function(df) {
    if (is_up) {
      df %>% arrange(desc(avg_log2FC))
    } else {
      df %>% arrange(avg_log2FC)
    }
  }

  sig <- dir_tbl %>%
    filter(p_val_adj < 0.05) %>%
    arrange_dir() %>%
    slice_head(n = n)
  if (nrow(sig) < n) {
    supl <- dir_tbl %>%
      filter(!(gene %in% sig$gene)) %>%
      arrange_dir() %>%
      slice_head(n = n - nrow(sig))
    sig <- bind_rows(sig, supl)
  }
  sig
}

default_up_genes_raw <- c(
  "KLRC3",
  "CD5",
  "CD72",
  "NKRF",
  "TNFRSF9",
  "IL-10RB",
  "LANY",
  "PPARA",
  "UCP2",
  "TRIM5",
  "TET3",
  "GZF1",
  "ADORA2A",
  "TGFBR3L",
  "CD59",
  "BCL2L1",
  "CSK",
  "KLF9",
  "IKZF2",
  "CBL",
  "FKBP1A",
  "STAT5A",
  "SOCS4",
  "SNX9",
  "IL2RG",
  "CD74",
  "TNFSF10",
  "TRIM27",
  "SLAMF1",
  "LY9"
)

read_gene_file <- function(path) {
  lines <- read_lines(path, progress = FALSE)
  lines <- str_trim(lines)
  lines <- lines[lines != ""]
  unique(lines)
}

up_genes_raw <- if (!is.na(opt$up_genes_file) && nzchar(opt$up_genes_file)) {
  message(">>> 使用用户提供的上调基因文件：", opt$up_genes_file)
  read_gene_file(opt$up_genes_file)
} else {
  default_up_genes_raw
}

resolve_requested_genes <- function(requested_raw, de_gene, de_gene_norm) {
  req_norm <- normalize_gene(requested_raw)
  tibble(raw_gene = requested_raw, normalized = req_norm) %>%
    mutate(
      present_in_de = normalized %in% de_gene_norm,
      matched_symbol = ifelse(
        present_in_de,
        vapply(
          normalized,
          function(x) {
            hits <- unique(de_gene[de_gene_norm == x])
            if (length(hits) == 0) "" else hits[[1]]
          },
          character(1)
        ),
        ""
      )
    )
}

up_map_tbl <- resolve_requested_genes(up_genes_raw, de_tbl$gene, de_tbl$gene_norm) %>%
  mutate(
    present_up = matched_symbol != "" &
      matched_symbol %in% de_tbl$gene[de_tbl$avg_log2FC > 0]
  )

up_genes_used <- up_map_tbl %>%
  filter(present_up) %>%
  pull(matched_symbol) %>%
  unique()

select_all_down <- function(tbl) {
  down_tbl <- tbl %>%
    filter(avg_log2FC < 0)
  if (isTRUE(opt$down_require_signif)) {
    if (isTRUE(opt$down_use_padj)) {
      down_tbl <- down_tbl %>% filter(p_val_adj <= opt$down_padj)
    } else {
      down_tbl <- down_tbl %>% filter(p_val <= opt$down_p)
    }
  }
  down_tbl %>%
    filter(avg_log2FC <= -opt$down_log2fc) %>%
    arrange(avg_log2FC)
}

down_tbl <- if (opt$down_mode == "top") {
  select_top(de_tbl, "Down", opt$down_n)
} else {
  select_all_down(de_tbl)
}
down_genes_used <- unique(down_tbl$gene)

all_genes_used <- unique(c(up_genes_used, down_genes_used))

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
plot_dir <- file.path(opt$output_dir, "figures")
table_dir <- file.path(opt$output_dir, "tables")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

write_lines(up_genes_used, file.path(table_dir, "up_selected_genes_used.txt"))
down_tag <- if (opt$down_mode == "top") glue("down_top{opt$down_n}") else "down_all"
write_lines(down_genes_used, file.path(table_dir, glue("{down_tag}_genes_used.txt")))
write_lines(all_genes_used, file.path(table_dir, "all_selected_genes_used.txt"))
write_csv(up_map_tbl, file.path(table_dir, "up_selected_genes_mapping_check.csv"))

message(">>> 上调（用于GO）的基因数：", length(up_genes_used))
if (opt$down_mode == "top") {
  message(">>> 下调Top", opt$down_n, "（用于GO）的基因数：", length(down_genes_used))
} else {
  signif_label <- if (!isTRUE(opt$down_require_signif)) {
    "不筛显著性"
  } else if (isTRUE(opt$down_use_padj)) {
    glue("padj≤{opt$down_padj}")
  } else {
    glue("p≤{opt$down_p}")
  }
  message(">>> 下调全部（用于GO）的基因数：", length(down_genes_used),
          "；阈值：avg_log2FC≤-", opt$down_log2fc, "；", signif_label)
}
message(">>> 合并（用于GO）的基因数：", length(all_genes_used))

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

universe_symbols <- unique(de_tbl$gene)
universe_entrez <- map2entrez(universe_symbols) %>%
  pull(ENTREZID) %>%
  unique()

run_go <- function(entrez_df) {
  if (nrow(entrez_df) < opt$min_genes) {
    return(NULL)
  }
  universe_ids <- if (isTRUE(opt$use_universe)) universe_entrez else NULL
  tryCatch(
    enrichGO(
      gene = entrez_df$ENTREZID,
      universe = universe_ids,
      OrgDb = org.Hs.eg.db,
      ont = opt$ont,
      keyType = "ENTREZID",
      pAdjustMethod = "BH",
      pvalueCutoff = opt$pvalue_cutoff,
      qvalueCutoff = opt$qvalue_cutoff,
      readable = TRUE
    ),
    error = function(e) {
      message(">>> enrichGO失败：", conditionMessage(e))
      NULL
    }
  )
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

plot_go_dot <- function(enrich_obj, title_text, subtitle_text, low_color) {
  if (is.null(enrich_obj) || nrow(enrich_obj) == 0) {
    return(NULL)
  }
  show_n <- min(opt$top_terms, nrow(enrich_obj))
  dotplot(enrich_obj, showCategory = show_n, font.size = 12) +
    scale_color_gradient(low = low_color, high = "#F2F2F2", name = "Adjusted p") +
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
  ggsave(filename = png_path, plot = plot_obj, width = width, height = height, units = "in", dpi = 600)
  ggsave(filename = pdf_path, plot = plot_obj, width = width, height = height, units = "in")
  message(">>> 已保存：", png_path, " 和 ", pdf_path)
}

up_entrez <- map2entrez(up_genes_used)
down_entrez <- map2entrez(down_genes_used)
all_entrez <- map2entrez(all_genes_used)

write_csv(up_entrez, file.path(table_dir, "up_selected_genes_entrez.csv"))
write_csv(down_entrez, file.path(table_dir, glue("{down_tag}_genes_entrez.csv")))
write_csv(all_entrez, file.path(table_dir, "all_selected_genes_entrez.csv"))

up_go <- run_go(up_entrez)
down_go <- run_go(down_entrez)
all_go <- run_go(all_entrez)

save_enrich_table(up_go, file.path(table_dir, glue("up_selected_GO_{opt$ont}.csv")))
save_enrich_table(down_go, file.path(table_dir, glue("{down_tag}_GO_{opt$ont}.csv")))
save_enrich_table(all_go, file.path(table_dir, glue("all_selected_GO_{opt$ont}.csv")))

subtitle_common <- if (isTRUE(opt$use_universe)) {
  glue("CD73+ vs CD73- ILC; background=DE genes (universe={length(universe_entrez)} Entrez)")
} else {
  "CD73+ vs CD73- ILC; background=whole genome"
}
subtitle_up <- glue("{subtitle_common}; up(selected)={length(up_genes_used)} genes")
subtitle_down <- if (opt$down_mode == "top") {
  glue("{subtitle_common}; down(top{opt$down_n})={length(down_genes_used)} genes")
} else {
  glue("{subtitle_common}; down(all)={length(down_genes_used)} genes")
}
subtitle_all <- glue("{subtitle_common}; all(selected)={length(all_genes_used)} genes")

ggsave_plot(
  plot_go_dot(up_go, glue("Up gene set (selected) GO {opt$ont}"), subtitle_up, low_color = "#E64B35"),
  file.path(plot_dir, glue("fig_up_selected_GO_{opt$ont}_dotplot"))
)
ggsave_plot(
  plot_go_dot(
    down_go,
    if (opt$down_mode == "top") glue("Down gene set (top{opt$down_n}) GO {opt$ont}") else glue("Down gene set (all) GO {opt$ont}"),
    subtitle_down,
    low_color = "#4DBBD5"
  ),
  file.path(plot_dir, glue("fig_{down_tag}_GO_{opt$ont}_dotplot"))
)
ggsave_plot(
  plot_go_dot(all_go, glue("All gene set (selected) GO {opt$ont}"), subtitle_all, low_color = "#6A3D9A"),
  file.path(plot_dir, glue("fig_all_selected_GO_{opt$ont}_dotplot"))
)

summary_tbl <- tibble(
  set = c("up_selected", glue("down_top{opt$down_n}"), "all_selected"),
  symbol_n = c(length(up_genes_used), length(down_genes_used), length(all_genes_used)),
  entrez_n = c(nrow(up_entrez), nrow(down_entrez), nrow(all_entrez)),
  go_terms = c(
    if (is.null(up_go)) 0 else nrow(up_go),
    if (is.null(down_go)) 0 else nrow(down_go),
    if (is.null(all_go)) 0 else nrow(all_go)
  )
)
write_csv(summary_tbl, file.path(table_dir, "summary.csv"))
write_lines(c(capture.output(sessionInfo())), file.path(opt$output_dir, "sessionInfo.txt"))

message(">>> 完成：输出目录 ", opt$output_dir)
