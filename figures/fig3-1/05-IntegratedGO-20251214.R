#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(glue)
  library(tibble)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

option_list <- list(
  make_option("--up_genes", type = "character", default = "outputs/CD73_ILC_DE_up_genes.txt",
              help = "上调基因列表文件（每行一个SYMBOL）"),
  make_option("--down_genes", type = "character", default = "outputs/CD73_ILC_DE_down_genes.txt",
              help = "下调基因列表文件（每行一个SYMBOL）"),
  make_option("--output_dir", type = "character", default = "outputs/CD73_ILC_GO_integrated_20251214",
              help = "输出文件夹（将自动创建figures和tables子目录）"),
  make_option("--ont", type = "character", default = "BP",
              help = "GO本体：BP/CC/MF"),
  make_option("--pvalue_cutoff", type = "double", default = 0.05,
              help = "GO富集pvalueCutoff"),
  make_option("--qvalue_cutoff", type = "double", default = 0.2,
              help = "GO富集qvalueCutoff"),
  make_option("--top_terms", type = "integer", default = 15,
              help = "每个方向展示的最大条目数"),
  make_option("--min_genes", type = "integer", default = 5,
              help = "做GO的最小基因数"),
  make_option("--use_universe", type = "logical", default = FALSE,
              help = "是否使用--universe_de_input中的DE基因作为背景(universe)"),
  make_option("--universe_de_input", type = "character", default = "outputs/CD73_ILC_DE_results.csv",
              help = "当--use_universe=TRUE时：从该DE结果csv取gene列作为背景(universe)")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!opt$ont %in% c("BP", "CC", "MF")) stop("--ont只能是 BP/CC/MF")

read_gene_file <- function(path) {
  x <- read_lines(path, progress = FALSE)
  x <- trimws(x)
  x <- x[x != ""]
  unique(x)
}

message(">>> 读取基因文件")
up_symbols <- read_gene_file(opt$up_genes)
down_symbols <- read_gene_file(opt$down_genes)
all_symbols <- unique(c(up_symbols, down_symbols))

message(">>> Up genes: ", length(up_symbols))
message(">>> Down genes: ", length(down_symbols))
message(">>> All (union) genes: ", length(all_symbols))

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
plot_dir <- file.path(opt$output_dir, "figures")
table_dir <- file.path(opt$output_dir, "tables")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

write_lines(up_symbols, file.path(table_dir, "up_genes_input.txt"))
write_lines(down_symbols, file.path(table_dir, "down_genes_input.txt"))
write_lines(all_symbols, file.path(table_dir, "all_genes_union.txt"))

map2entrez <- function(symbols) {
  if (length(symbols) == 0) return(tibble())
  bitr(
    symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db,
    drop = TRUE
  ) %>%
    distinct(SYMBOL, .keep_all = TRUE)
}

universe_entrez <- NULL
subtitle_bg <- "background=whole genome"
if (isTRUE(opt$use_universe)) {
  message(">>> 读取universe(背景)自：", opt$universe_de_input)
  de <- read_csv(opt$universe_de_input, show_col_types = FALSE)
  stopifnot("gene" %in% names(de))
  universe_entrez <- map2entrez(unique(de$gene)) %>%
    pull(ENTREZID) %>%
    unique()
  subtitle_bg <- glue("background=DE genes (universe={length(universe_entrez)} Entrez)")
}

run_go <- function(entrez_df) {
  if (nrow(entrez_df) < opt$min_genes) return(NULL)
  tryCatch(
    enrichGO(
      gene = entrez_df$ENTREZID,
      universe = universe_entrez,
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
  if (is.null(enrich_obj) || nrow(enrich_obj) == 0) return(NULL)
  show_n <- min(opt$top_terms, nrow(enrich_obj))
  dotplot(enrich_obj, showCategory = show_n, font.size = 12) +
    scale_color_gradient(low = low_color, high = "#F2F2F2", name = "Adjusted p") +
    labs(title = title_text, subtitle = subtitle_text, x = "Gene Ratio", y = NULL) +
    nature_theme
}

ggsave_plot <- function(plot_obj, base_path, width = 7, height = 6.5, dpi = 600) {
  if (is.null(plot_obj)) {
    message(">>> 跳过未生成的图：", base_path)
    return(invisible(NULL))
  }
  png_path <- paste0(base_path, ".png")
  pdf_path <- paste0(base_path, ".pdf")
  ggsave(filename = png_path, plot = plot_obj, width = width, height = height, units = "in", dpi = dpi)
  ggsave(filename = pdf_path, plot = plot_obj, width = width, height = height, units = "in")
  message(">>> 已保存：", png_path, " 和 ", pdf_path)
}

up_entrez <- map2entrez(up_symbols)
down_entrez <- map2entrez(down_symbols)
all_entrez <- map2entrez(all_symbols)

write_csv(up_entrez, file.path(table_dir, "up_entrez.csv"))
write_csv(down_entrez, file.path(table_dir, "down_entrez.csv"))
write_csv(all_entrez, file.path(table_dir, "all_entrez.csv"))

up_go <- run_go(up_entrez)
down_go <- run_go(down_entrez)
all_go <- run_go(all_entrez)

save_enrich_table(up_go, file.path(table_dir, glue("up_GO_{opt$ont}.csv")))
save_enrich_table(down_go, file.path(table_dir, glue("down_GO_{opt$ont}.csv")))
save_enrich_table(all_go, file.path(table_dir, glue("all_GO_{opt$ont}.csv")))

subtitle_common <- glue("CD73+ vs CD73- ILC; {subtitle_bg}")
ggsave_plot(
  plot_go_dot(up_go, glue("Up gene set GO {opt$ont}"), glue("{subtitle_common}; up={length(up_symbols)} genes"), low_color = "#E64B35"),
  file.path(plot_dir, glue("fig_up_GO_{opt$ont}_dotplot"))
)
ggsave_plot(
  plot_go_dot(down_go, glue("Down gene set GO {opt$ont}"), glue("{subtitle_common}; down={length(down_symbols)} genes"), low_color = "#4DBBD5"),
  file.path(plot_dir, glue("fig_down_GO_{opt$ont}_dotplot"))
)
ggsave_plot(
  plot_go_dot(all_go, glue("All (union) gene set GO {opt$ont}"), glue("{subtitle_common}; all={length(all_symbols)} genes"), low_color = "#6A3D9A"),
  file.path(plot_dir, glue("fig_all_GO_{opt$ont}_dotplot"))
)

combine_tables <- function(path, direction) {
  if (!file.exists(path)) return(tibble())
  read_csv(path, show_col_types = FALSE) %>%
    mutate(Direction = direction) %>%
    relocate(Direction, .before = 1)
}

integrated_tbl <- bind_rows(
  combine_tables(file.path(table_dir, glue("up_GO_{opt$ont}.csv")), "Up"),
  combine_tables(file.path(table_dir, glue("down_GO_{opt$ont}.csv")), "Down"),
  combine_tables(file.path(table_dir, glue("all_GO_{opt$ont}.csv")), "All")
)
if (nrow(integrated_tbl) > 0) {
  write_csv(integrated_tbl, file.path(table_dir, glue("integrated_GO_{opt$ont}.csv")))
}

summary_tbl <- tibble(
  set = c("up", "down", "all"),
  symbol_n = c(length(up_symbols), length(down_symbols), length(all_symbols)),
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

