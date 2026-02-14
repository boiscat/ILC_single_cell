#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(glue)
  library(ggplot2)
})


validate_hex_color <- function(x, opt_name) {
  if (is.null(x) || !nzchar(x)) stop(opt_name, " 不能为空")
  if (!grepl("^#(?:[0-9A-Fa-f]{6}|[0-9A-Fa-f]{8})$", x)) {
    stop(opt_name, " 不是合法的hex颜色: ", x, "（示例：#2A7F8E）")
  }
  x
}
opt$point_color <- validate_hex_color(opt$point_color, "--point_color")
opt$row_fill <- validate_hex_color(opt$row_fill, "--row_fill")
if (opt$row_alpha < 0 || opt$row_alpha > 1) stop("--row_alpha 必须在 0-1 之间")

if (!file.exists(opt$go_bp_csv)) stop("找不到输入表: ", opt$go_bp_csv)

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
plot_dir <- file.path(opt$output_dir, "figures")
table_dir <- file.path(opt$output_dir, "tables")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

wrap_text <- function(x, width) {
  vapply(
    x,
    FUN = function(s) paste(strwrap(s, width = width), collapse = "\n"),
    FUN.VALUE = character(1)
  )
}

parse_ratio <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  ok <- lengths(parts) == 2
  out <- rep(NA_real_, length(parts))
  out[ok] <- vapply(parts[ok], function(p) as.numeric(p[[1]]) / as.numeric(p[[2]]), numeric(1))
  out
}

message(">>> 读取GO BP结果表：", opt$go_bp_csv)
go <- read_csv(opt$go_bp_csv, show_col_types = FALSE)
required_cols <- c("Description", "GeneRatio", "p.adjust", "Count")
miss <- setdiff(required_cols, names(go))
if (length(miss) > 0) stop("输入表缺少列：", paste(miss, collapse = ", "))

go2 <- go %>%
  mutate(
    GeneRatio_num = parse_ratio(GeneRatio),
    GeneRatio_pct = 100 * GeneRatio_num,
    FDR = as.numeric(p.adjust),
    neglog10_FDR = -log10(FDR)
  ) %>%
  filter(!is.na(GeneRatio_pct), !is.na(neglog10_FDR), is.finite(neglog10_FDR))

if (nrow(go2) == 0) stop("输入表解析后无有效条目（GeneRatio或p.adjust可能有问题）")

normalize_for_match <- function(x) {
  tolower(x) %>%
    gsub("[^a-z0-9]", "", ., perl = TRUE)
}

parse_keyword_vec <- function(x) {
  v <- unlist(strsplit(as.character(x), ";", fixed = TRUE), use.names = FALSE)
  v <- trimws(v)
  v[nzchar(v)]
}

match_keywords <- function(descriptions, keywords_vec) {
  if (length(keywords_vec) == 0) return(rep(FALSE, length(descriptions)))
  desc_norm <- normalize_for_match(descriptions)
  key_norm <- normalize_for_match(keywords_vec)
  hit <- rep(FALSE, length(desc_norm))
  for (k in key_norm) {
    if (!nzchar(k)) next
    hit <- hit | grepl(k, desc_norm, fixed = TRUE)
  }
  hit
}

select_show_terms_dotplot <- function(tbl, top_n, keywords_vec) {
  if (is.null(tbl) || nrow(tbl) == 0) return(character())
  tbl <- tbl %>% arrange(FDR)
  forced <- tbl$Description[match_keywords(tbl$Description, keywords_vec)]
  selected <- character()
  for (d in forced) {
    if (!d %in% selected) selected <- c(selected, d)
  }
  for (d in tbl$Description) {
    if (length(selected) >= top_n) break
    if (!d %in% selected) selected <- c(selected, d)
  }
  selected
}

keep_mode <- nzchar(opt$keep_terms)

if (keep_mode) {
  keep_vec <- parse_keyword_vec(opt$keep_terms)
  if (length(keep_vec) == 0) stop("--keep_terms 解析后为空：请用英文分号;分隔，并确保 term 非空")

  keep_norm <- normalize_for_match(keep_vec)
  go_norm <- normalize_for_match(go2$Description)

  # 按用户给定顺序挑选；每个 term 只取一次
  idx <- integer()
  for (k in keep_norm) {
    hit <- which(go_norm == k)
    if (length(hit) == 0) {
      stop("在筛选结果中找不到 keep term：", keep_vec[which(keep_norm == k)[1]])
    }
    idx <- c(idx, hit[1])
  }
  sel <- go2[idx, , drop = FALSE]
  message(">>> keep_terms 生效：仅保留 ", nrow(sel), " 条 term（按给定顺序）")
} else if (opt$select_mode == "dotplot") {
  message(">>> 选择模式：dotplot（复刻dotplot的showCategory选择规则）")
  keyword_vec <- parse_keyword_vec(opt$dotplot_keywords)
  show_terms <- select_show_terms_dotplot(go2, opt$top_terms, keyword_vec)
  sel <- go2 %>%
    filter(Description %in% show_terms)
} else {
  message(">>> 选择模式：regex（按正则筛选Description）")
  message(">>> regex：", opt$keywords)
  sel <- go2 %>%
    filter(grepl(opt$keywords, Description, ignore.case = TRUE))
}

if (nrow(sel) == 0) {
  stop(
    "筛选后没有条目。请检查 --keywords 是否能匹配到Description。\n",
    "示例：--keywords 'immune|leukocyte|lymphocyte|T cell'"
  )
}

show_n <- if (keep_mode) nrow(sel) else min(opt$top_terms, nrow(sel))
order_by_actual <- opt$order_by
if (order_by_actual == "auto") {
  order_by_actual <- if (keep_mode) "keep" else "GeneRatio_pct"
}

sel <- sel %>% slice_head(n = show_n)
if (!keep_mode || order_by_actual != "keep") {
  if (order_by_actual == "GeneRatio_pct") {
    sel <- sel %>% arrange(desc(GeneRatio_pct), FDR)
  } else if (order_by_actual == "neglog10_FDR") {
    sel <- sel %>% arrange(desc(neglog10_FDR), desc(GeneRatio_pct))
  } else if (order_by_actual != "keep") {
    stop("未知的排序方式：", order_by_actual)
  }
}

sel <- sel %>% mutate(Description_wrapped = wrap_text(Description, opt$wrap_width))

terms_out_df <- sel %>% select(ID, Description, GeneRatio, GeneRatio_pct, p.adjust, FDR, neglog10_FDR, Count)
terms_csv_latest <- file.path(table_dir, "up_selected_GO_BP_terms_used_for_lollipop_dualx.csv")
terms_csv_run <- file.path(table_dir, paste0(opt$out_basename, ".terms_used.csv"))
write_csv(terms_out_df, terms_csv_run)
write_csv(terms_out_df, terms_csv_latest)

levels_desc <- rev(sel$Description_wrapped)

plot_df <- bind_rows(
  sel %>%
    transmute(
      term = Description_wrapped,
      metric = "GeneRatio (%)",
      value = GeneRatio_pct
    ),
  sel %>%
    transmute(
      term = Description_wrapped,
      metric = "-log10\n(FDR)",
      value = neglog10_FDR
    )
) %>%
  mutate(
    term = factor(term, levels = levels_desc),
    metric = factor(metric, levels = c("GeneRatio (%)", "-log10\n(FDR)"))
  )

bg_df <- tibble::tibble(term = factor(levels_desc, levels = levels_desc)) %>%
  mutate(y = as.numeric(term)) %>%
  mutate(row_from_top = length(levels_desc) - y + 1) %>%
  filter(row_from_top %% 2 == 1) %>%
  transmute(
    xmin = -Inf, xmax = Inf,
    ymin = y - 0.5, ymax = y + 0.5
  )

theme_ref <- theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 11),
    legend.position = "none"
  )

p <- ggplot(plot_df, aes(x = value, y = term)) +
  geom_rect(
    data = bg_df,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = opt$row_fill,
    alpha = opt$row_alpha
  ) +
  geom_segment(aes(x = 0, xend = value, yend = term), color = opt$point_color, linewidth = 1.0) +
  geom_point(color = opt$point_color, fill = "white", shape = 21, stroke = 1.0, size = 3.2) +
  facet_grid(. ~ metric, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = opt$set_label,
    subtitle = opt$subtitle,
    x = NULL,
    y = NULL
  ) +
  theme_ref

base_path <- file.path(plot_dir, opt$out_basename)
ggsave(filename = paste0(base_path, ".png"), plot = p, width = opt$width_in, height = opt$height_in, units = "in", dpi = opt$dpi)
ggsave(filename = paste0(base_path, ".pdf"), plot = p, width = opt$width_in, height = opt$height_in, units = "in")

session_latest <- file.path(opt$output_dir, "sessionInfo_immune_Tcell_lollipop_dualx.txt")
session_run <- file.path(opt$output_dir, paste0("sessionInfo_immune_Tcell_lollipop_dualx.", opt$out_basename, ".txt"))
writeLines(c(capture.output(sessionInfo())), session_run)
writeLines(c(capture.output(sessionInfo())), session_latest)
message(">>> 已保存：", paste0(base_path, ".png"))
message(">>> 已保存：", paste0(base_path, ".pdf"))
message(">>> 已保存：", terms_csv_run)
message(">>> 已保存：", session_run)
message(">>> 完成：输出目录 ", opt$output_dir)
