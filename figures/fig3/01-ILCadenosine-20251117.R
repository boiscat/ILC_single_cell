#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

message(">>> 载入Seurat对象")
seurat_path <- "full_obj_saved_1117.rds"
obj <- readRDS(seurat_path)

meta <- obj@meta.data %>%
  mutate(cell_id = rownames(.))

if (!"Site_Status" %in% colnames(meta)) {
  stop("Site_Status列不存在，无法区分Healthy/lesion/non_lesion")
}

status_levels <- c("Healthy", "non_lesion", "lesion")

# 筛选ILC细胞（cell_anno为ILC或full_clustering以ILC开头）
ilc_cells <- meta %>%
  filter(!is.na(Site_Status)) %>%
  filter(cell_anno == "ILC" | grepl("^ILC", full_clustering))

if (nrow(ilc_cells) == 0) {
  stop("没有识别出ILC细胞，请检查注释列")
}

message(">>> ILC细胞数量: ", nrow(ilc_cells))

ilc_obj <- subset(obj, cells = ilc_cells$cell_id)

adenosine_genes <- c(
  CD39 = "ENTPD1",
  CD73 = "NT5E",
  A1 = "ADORA1",
  A2A = "ADORA2A",
  A2B = "ADORA2B",
  A3 = "ADORA3"
)

present_genes <- intersect(adenosine_genes, rownames(ilc_obj))
missing_genes <- setdiff(adenosine_genes, present_genes)

if (length(present_genes) == 0) {
  stop("腺苷通路基因全部缺失")
}

if (length(missing_genes) > 0) {
  warning("以下基因在对象中未找到，将被跳过: ",
          paste(missing_genes, collapse = ", "))
}

gene_alias_tbl <- tibble(
  alias = names(adenosine_genes),
  gene = adenosine_genes
) %>%
  filter(gene %in% present_genes) %>%
  mutate(alias_label = paste0(alias, " (", gene, ")"))

fetch_vars <- c("Site_Status", present_genes)
expr_df <- FetchData(ilc_obj, vars = fetch_vars) %>%
  mutate(cell_id = rownames(.))

long_expr <- expr_df %>%
  pivot_longer(
    cols = all_of(present_genes),
    names_to = "gene",
    values_to = "expr"
  ) %>%
  left_join(gene_alias_tbl, by = "gene") %>%
  mutate(
    Site_Status = factor(Site_Status, levels = status_levels),
    alias_label = factor(alias_label, levels = unique(gene_alias_tbl$alias_label))
  )

if (any(is.na(long_expr$Site_Status))) {
  stop("存在未定义Site_Status的ILC细胞")
}

summary_df <- long_expr %>%
  group_by(Site_Status, gene, alias_label) %>%
  summarise(
    cells_in_group = n(),
    mean_expr = mean(expr),
    median_expr = median(expr),
    pct_expr = mean(expr > 0) * 100,
    .groups = "drop"
  )

# 简单pairwise比较
pairings <- tibble(
  group1 = c("lesion", "non_lesion", "lesion"),
  group2 = c("Healthy", "Healthy", "non_lesion")
)

pairwise_stats <- expand_grid(gene_alias_tbl, pairings) %>%
  mutate(
    p_value = purrr::pmap_dbl(
      list(gene, group1, group2),
      function(g, g1, g2) {
        v1 <- long_expr %>% filter(gene == g, Site_Status == g1) %>% pull(expr)
        v2 <- long_expr %>% filter(gene == g, Site_Status == g2) %>% pull(expr)
        if (length(v1) < 3 || length(v2) < 3) {
          return(NA_real_)
        }
        suppressWarnings({
          test <- wilcox.test(v1, v2)
          test$p.value
        })
      }
    ),
    log2_fc = purrr::pmap_dbl(
      list(gene, group1, group2),
      function(g, g1, g2) {
        v1 <- long_expr %>% filter(gene == g, Site_Status == g1) %>% pull(expr)
        v2 <- long_expr %>% filter(gene == g, Site_Status == g2) %>% pull(expr)
        mean(v1) - mean(v2)
      }
    ),
    comparison = paste0(group1, "_vs_", group2)
  ) %>%
  select(alias, gene, comparison, log2_fc, p_value)

output_dir <- "outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

summary_path <- file.path(output_dir, "ILC_adenosine_expression_summary.tsv")
pairwise_path <- file.path(output_dir, "ILC_adenosine_pairwise.tsv")

write.table(
  summary_df,
  file = summary_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  pairwise_stats,
  file = pairwise_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(">>> 已写出: ", summary_path)
message(">>> 已写出: ", pairwise_path)

p <- ggplot(summary_df, aes(x = Site_Status, y = alias_label)) +
  geom_point(aes(size = pct_expr, color = mean_expr)) +
  scale_size(name = "%表达细胞", range = c(1, 8)) +
  scale_color_viridis_c(name = "平均表达") +
  theme_bw(base_size = 12) +
  labs(
    x = "部位状态",
    y = "腺苷信号节点",
    title = "ILC腺苷信号通路基因表达（按部位状态）"
  )

plot_path <- file.path(output_dir, "ILC_adenosine_dotplot.pdf")
ggsave(
  filename = plot_path,
  plot = p,
  width = 6,
  height = 4,
  device = cairo_pdf
)
message(">>> 已写出: ", plot_path)

message("分析完成")
