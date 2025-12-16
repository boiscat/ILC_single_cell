#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
stage <- "all"
if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[[i]] == "--stage" && i < length(args)) {
      stage <- args[[i + 1]]
    }
  }
}
if (!stage %in% c("all", "de_only")) {
  stop("stage参数只能是'all'或'de_only'")
}

message("运行阶段: ", stage)

input_path <- "full_obj_saved_1117.rds"
output_dir <- "outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
de_file <- file.path(output_dir, "CD73_ILC_DE_results.csv")
summary_file <- file.path(output_dir, "CD73_ILC_DE_summary.csv")
plot_file <- file.path(output_dir, "CD73_ILC_DE_barplot.png")
plot_file_top30 <- file.path(output_dir, "CD73_ILC_DE_barplot_top30.png")
plot_file_up_top30_highlight <- file.path(output_dir, "CD73_ILC_DE_barplot_up_top30_highlight.png")
highlight_check_file <- file.path(output_dir, "CD73_ILC_DE_up_highlight_check.csv")
plot_file_up_selected30 <- file.path(output_dir, "CD73_ILC_DE_barplot_up_selected30.png")
selected_check_file <- file.path(output_dir, "CD73_ILC_DE_up_selected30_check.csv")
plot_file_down_top30 <- file.path(output_dir, "CD73_ILC_DE_barplot_down_top30.png")

message("读取对象: ", input_path)
obj <- readRDS(input_path)

if (!"cell_anno" %in% colnames(obj@meta.data)) {
  stop("meta.data缺少cell_anno列，无法筛选ILC细胞")
}

ilc <- subset(obj, subset = cell_anno == "ILC")
message("ILC细胞数量: ", ncol(ilc))
if (ncol(ilc) == 0) stop("没有找到ILC细胞")

gene_of_interest <- "NT5E" # CD73
if (!(gene_of_interest %in% rownames(ilc))) {
  stop("在表达矩阵中找不到基因: ", gene_of_interest)
}

nt5e_expr <- FetchData(ilc, vars = gene_of_interest)[, 1]
cd73_threshold <- 0 # 只要在log-normalized数据中表达大于0即视为阳性
ilc$cd73_status <- ifelse(nt5e_expr > cd73_threshold, "CD73_pos", "CD73_neg")

group_counts <- table(ilc$cd73_status)
message("CD73分类: ", paste(names(group_counts), group_counts, sep = "=", collapse = " | "))
if (length(group_counts) < 2) stop("CD73+/-分组不完整，无法差异分析")

Idents(ilc) <- "cd73_status"
set.seed(1234)
de_res <- FindMarkers(
  ilc,
  ident.1 = "CD73_pos",
  ident.2 = "CD73_neg",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)

de_tbl <- de_res %>%
  rownames_to_column("gene") %>%
  arrange(desc(avg_log2FC))

write.csv(de_tbl, de_file, row.names = FALSE)
message("差异结果写入: ", de_file)

summary_tbl <- tibble(
  total_ilc = ncol(ilc),
  cd73_pos = as.integer(group_counts["CD73_pos"]),
  cd73_neg = as.integer(group_counts["CD73_neg"]),
  nt5e_threshold = cd73_threshold,
  logfc_threshold = 0.25,
  min_pct = 0.1
)
write.csv(summary_tbl, summary_file, row.names = FALSE)
message("统计摘要写入: ", summary_file)

if (stage != "de_only") {
  normalize_gene <- function(x) gsub("[^A-Za-z0-9_.]", "", x)

  highlight_genes_raw <- c(
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
  highlight_genes_norm <- unique(normalize_gene(highlight_genes_raw))
  selected_genes_raw <- highlight_genes_raw

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
    mutate(sig, Group = direction)
  }

  create_plot <- function(up_n, down_n, file_path, title_note) {
    plot_data <- bind_rows(
      select_top(de_tbl, "Up", up_n),
      select_top(de_tbl, "Down", down_n)
    ) %>%
      mutate(gene = factor(gene, levels = unique(gene[order(avg_log2FC, decreasing = TRUE)])))

    if (nrow(plot_data) == 0) {
      warning("没有可用基因，跳过绘图: ", title_note)
    } else {
      p <- ggplot(plot_data, aes(x = gene, y = avg_log2FC, fill = Group)) +
        geom_col(width = 0.7, color = "black", linewidth = 0.2) +
        coord_flip() +
        labs(
          title = paste0("Human (", title_note, ")"),
          x = "gene",
          y = "Avg_log2FC"
        ) +
        scale_fill_manual(values = c(Up = "#E64B35", Down = "#4DBBD5")) +
        theme_bw(base_size = 12) +
        theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
        )
      ggsave(file_path, p, width = 4, height = 6, dpi = 300)
      message("绘图保存: ", file_path)
    }
  }

  create_plot_up_only_highlight <- function(up_n, file_path, title_note) {
    plot_data <- select_top(de_tbl, "Up", up_n) %>%
      mutate(
        gene_norm = normalize_gene(gene),
        Highlight = ifelse(gene_norm %in% highlight_genes_norm, "Highlight", "Up"),
        gene = factor(gene, levels = unique(gene[order(avg_log2FC, decreasing = FALSE)]))
      )

    highlight_check <- tibble(raw_gene = highlight_genes_raw) %>%
      mutate(
        normalized = normalize_gene(raw_gene),
        present_in_de = normalized %in% normalize_gene(de_tbl$gene),
        present_in_top = normalized %in% plot_data$gene_norm,
        matched_symbol = ifelse(
          present_in_de,
          vapply(
            normalized,
            function(x) {
              hits <- unique(de_tbl$gene[normalize_gene(de_tbl$gene) == x])
              if (length(hits) == 0) "" else paste(hits, collapse = ";")
            },
            character(1)
          ),
          ""
        )
      )
    write.csv(highlight_check, highlight_check_file, row.names = FALSE)
    message("高亮基因匹配检查写入: ", highlight_check_file)

    if (nrow(plot_data) == 0) {
      warning("没有可用上调基因，跳过绘图: ", title_note)
    } else {
      p <- ggplot(plot_data, aes(x = gene, y = avg_log2FC, fill = Highlight)) +
        geom_col(width = 0.7, color = "black", linewidth = 0.2) +
        coord_flip() +
        labs(
          title = paste0("Human (", title_note, ")"),
          x = "gene",
          y = "Avg_log2FC",
          fill = "Group"
        ) +
        scale_fill_manual(values = c(Up = "#E64B35", Highlight = "#7A1F1F")) +
        theme_bw(base_size = 12) +
        theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
        )
      ggsave(file_path, p, width = 4, height = 6, dpi = 300)
      message("绘图保存: ", file_path)
    }
  }

  create_plot_up_only_selected <- function(file_path, title_note) {
    selected_norm <- normalize_gene(selected_genes_raw)

    plot_data <- de_tbl %>%
      mutate(gene_norm = normalize_gene(gene)) %>%
      filter(gene_norm %in% selected_norm, avg_log2FC > 0) %>%
      group_by(gene_norm) %>%
      slice_max(order_by = avg_log2FC, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(
        gene = factor(gene, levels = unique(gene[order(avg_log2FC, decreasing = FALSE)])),
        Group = "Up"
      )

    selected_check <- tibble(raw_gene = selected_genes_raw) %>%
      mutate(
        normalized = normalize_gene(raw_gene),
        present_in_de = normalized %in% normalize_gene(de_tbl$gene),
        present_up = normalized %in% normalize_gene(de_tbl$gene[de_tbl$avg_log2FC > 0]),
        present_in_plot = normalized %in% plot_data$gene_norm,
        matched_symbol = ifelse(
          present_in_de,
          vapply(
            normalized,
            function(x) {
              hits <- unique(de_tbl$gene[normalize_gene(de_tbl$gene) == x])
              if (length(hits) == 0) "" else paste(hits, collapse = ";")
            },
            character(1)
          ),
          ""
        )
      )
    write.csv(selected_check, selected_check_file, row.names = FALSE)
    message("指定基因匹配检查写入: ", selected_check_file)

    if (nrow(plot_data) == 0) {
      warning("指定基因在上调中没有匹配到任何条目，跳过绘图: ", title_note)
    } else {
      p <- ggplot(plot_data, aes(x = gene, y = avg_log2FC, fill = Group)) +
        geom_col(width = 0.7, color = "black", linewidth = 0.2) +
        coord_flip() +
        labs(
          title = paste0("Human (", title_note, ")"),
          x = "gene",
          y = "Avg_log2FC",
          fill = "Group"
        ) +
        scale_fill_manual(values = c(Up = "#E64B35")) +
        theme_bw(base_size = 12) +
        theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
        )
      ggsave(file_path, p, width = 4, height = 6, dpi = 300)
      message("绘图保存: ", file_path)
    }
  }

  create_plot_down_only_top <- function(down_n, file_path, title_note) {
    plot_data <- select_top(de_tbl, "Down", down_n) %>%
      mutate(
        gene = factor(gene, levels = unique(gene[order(avg_log2FC, decreasing = FALSE)])),
        Group = "Down"
      )

    if (nrow(plot_data) == 0) {
      warning("没有可用下调基因，跳过绘图: ", title_note)
    } else {
      p <- ggplot(plot_data, aes(x = gene, y = avg_log2FC, fill = Group)) +
        geom_col(width = 0.7, color = "black", linewidth = 0.2) +
        coord_flip() +
        labs(
          title = paste0("Human (", title_note, ")"),
          x = "gene",
          y = "Avg_log2FC",
          fill = "Group"
        ) +
        scale_fill_manual(values = c(Down = "#4DBBD5")) +
        theme_bw(base_size = 12) +
        theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
        )
      ggsave(file_path, p, width = 4, height = 6, dpi = 300)
      message("绘图保存: ", file_path)
    }
  }

  create_plot(10, 10, plot_file, "top10 each")
  create_plot(15, 15, plot_file_top30, "top30 total")
  create_plot_up_only_highlight(30, plot_file_up_top30_highlight, "top30 up only")
  create_plot_up_only_selected(plot_file_up_selected30, "selected up genes")
  create_plot_down_only_top(30, plot_file_down_top30, "top30 down only")
}

message("脚本执行完成")
