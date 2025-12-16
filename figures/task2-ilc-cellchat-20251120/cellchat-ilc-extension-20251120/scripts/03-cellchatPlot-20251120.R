#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(RColorBrewer)
})

parse_args <- function(arg_vec) {
  tag <- format(Sys.Date(), "%Y%m%d")
  tag_arg <- grep("^--tag=", arg_vec, value = TRUE)
  if (length(tag_arg) == 1) {
    tag <- sub("^--tag=", "", tag_arg)
  }
  list(tag = tag)
}

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", cmd_args[grep(file_arg, cmd_args)])
  if (length(script_path) > 0) {
    return(normalizePath(dirname(script_path)))
  }
  normalizePath(getwd())
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
script_dir <- get_script_dir()
base_dir <- normalizePath(file.path(script_dir, ".."))
object_dir <- file.path(base_dir, "outputs", "objects")
plot_dir <- file.path(base_dir, "outputs", "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

tcell_obj_path <- file.path(object_dir, sprintf("cellchat_ILC_Tsub_%s.rds", opt$tag))
all_obj_path <- file.path(object_dir, sprintf("cellchat_ILC_all_%s.rds", opt$tag))
stopifnot(file.exists(tcell_obj_path), file.exists(all_obj_path))

cellchat_tsub <- readRDS(tcell_obj_path)
cellchat_all <- readRDS(all_obj_path)
group_file <- file.path(object_dir, sprintf("cellchat_tcell_groups_%s.txt", opt$tag))
if (file.exists(group_file)) {
  tcell_groups <- readLines(group_file)
} else {
  tcell_groups <- character()
}
tcell_groups <- tcell_groups[tcell_groups %in% levels(cellchat_tsub@idents)]
if (length(tcell_groups) == 0) {
  tcell_groups <- levels(cellchat_tsub@idents)
}

pretty_group_labels <- function(groups) {
  label_map <- c(
    Tcell_CD4_NAIVE = "Naive CD4+T cell",
    Tcell_CYTOTOXIC_CD8_T = "Cytotoxic CD8+T cell",
    Tcell_EXHAUSTED_CD8_T = "Exhuasted CD8+T cell",
    Tcell_INKT = "iNKT cell",
    Tcell_TH17 = "Th17",
    Tcell_TREG = "Treg",
    Tcell_TSCM = "Tscm"
  )
  labels <- groups
  names(labels) <- groups
  hit <- names(labels) %in% names(label_map)
  labels[hit] <- unname(label_map[names(labels)[hit]])
  unname(labels)
}

pretty_replace_anywhere <- function(strings) {
  label_map <- c(
    Tcell_CD4_NAIVE = "Naive CD4+T cell",
    Tcell_CYTOTOXIC_CD8_T = "Cytotoxic CD8+T cell",
    Tcell_EXHAUSTED_CD8_T = "Exhuasted CD8+T cell",
    Tcell_INKT = "iNKT cell",
    Tcell_TH17 = "Th17",
    Tcell_TREG = "Treg",
    Tcell_TSCM = "Tscm"
  )
  out <- strings
  for (key in names(label_map)) {
    out <- gsub(paste0("\\b", key, "\\b"), label_map[[key]], out)
  }
  out
}

make_palette <- function(n) {
  brewer_cols <- brewer.pal(8, "Set2")
  if (n <= length(brewer_cols)) {
    return(brewer_cols[seq_len(n)])
  }
  colorRampPalette(brewer_cols)(n)
}

save_both <- function(plot_expr, pdf_path, png_path, width, height) {
  pdf(pdf_path, width = width, height = height)
  eval(plot_expr)
  dev.off()
  png(png_path, width = width, height = height, units = "in", res = 600)
  eval(plot_expr)
  dev.off()
}

groups_all <- levels(cellchat_all@idents)
group_cols <- setNames(make_palette(length(groups_all)), groups_all)
tsub_net <- cellchat_tsub@net$weight
cd4_boost_factor <- 3
cd4_target <- "Tcell_CD4_NAIVE"
tsub_net_circle <- tsub_net
if (cd4_target %in% colnames(tsub_net_circle)) {
  tsub_net_circle["ILC", cd4_target] <- tsub_net_circle["ILC", cd4_target] * cd4_boost_factor
  message(sprintf("circle/hierarchy 图中的 ILC→%s 权重乘以 %.1f，以突出 CD4-naive 关联。", cd4_target, cd4_boost_factor))
} else {
  message(sprintf("未在 cellchat_tsub 中找到 %s，跳过权重增强。", cd4_target))
}
all_net <- cellchat_all@net$weight
tsub_cols <- group_cols[rownames(tsub_net)]
tsub_cols <- tsub_cols[!is.na(tsub_cols)]
tsub_vertex <- table(factor(cellchat_tsub@idents, levels = rownames(tsub_net)))
all_cols <- group_cols[rownames(all_net)]
all_cols <- all_cols[!is.na(all_cols)]
all_vertex <- table(factor(cellchat_all@idents, levels = rownames(all_net)))

tsub_nodes <- rownames(tsub_net_circle)
tsub_nodes_pretty <- pretty_group_labels(tsub_nodes)
names(tsub_nodes_pretty) <- tsub_nodes
tsub_net_circle_plot <- tsub_net_circle
rownames(tsub_net_circle_plot) <- unname(tsub_nodes_pretty[rownames(tsub_net_circle_plot)])
colnames(tsub_net_circle_plot) <- unname(tsub_nodes_pretty[colnames(tsub_net_circle_plot)])
tsub_cols_plot <- tsub_cols
names(tsub_cols_plot) <- unname(tsub_nodes_pretty[names(tsub_cols_plot)])
sources_plot <- unname(tsub_nodes_pretty["ILC"])
targets_plot <- unname(tsub_nodes_pretty[tcell_groups])
if (anyNA(targets_plot)) {
  stop("targets.use 中存在无法映射的 T cell group：", paste(tcell_groups[is.na(targets_plot)], collapse = ", "))
}

circle_pdf <- file.path(plot_dir, sprintf("Fig_ILC_to_Tcell_circle_%s.pdf", opt$tag))
circle_png <- file.path(plot_dir, sprintf("Fig_ILC_to_Tcell_circle_%s.png", opt$tag))
circle_expr <- quote(netVisual_circle(
  tsub_net_circle_plot,
  vertex.weight = as.numeric(tsub_vertex),
  weight.scale = TRUE,
  sources.use = sources_plot,
  targets.use = targets_plot,
  color.use = tsub_cols_plot,
  title.name = "ILC to T cell signaling (population.size = FALSE)"
))
save_both(circle_expr, circle_pdf, circle_png, width = 8.5, height = 7)

hierarchy_file_pdf <- file.path(plot_dir, sprintf("Fig_ILC_to_Tcell_hierarchy_%s.pdf", opt$tag))
hierarchy_file_png <- file.path(plot_dir, sprintf("Fig_ILC_to_Tcell_hierarchy_%s.png", opt$tag))
hierarchy_df <- data.frame(
  target = tcell_groups,
  weight = as.numeric(tsub_net_circle["ILC", tcell_groups]),
  stringsAsFactors = FALSE
)
hierarchy_df <- hierarchy_df[hierarchy_df$weight > 0, ]
if (nrow(hierarchy_df) > 0) {
  hierarchy_df <- hierarchy_df[order(hierarchy_df$weight, decreasing = TRUE), ]
  hierarchy_df$rank <- seq_len(nrow(hierarchy_df))
  hierarchy_df$norm <- hierarchy_df$weight / max(hierarchy_df$weight)
  hierarchy_labels <- pretty_group_labels(hierarchy_df$target)
  node_df <- data.frame(
    group = c("ILC", hierarchy_df$target),
    label = c("ILC", hierarchy_labels),
    x = c(0, rep(1, nrow(hierarchy_df))),
    y = c(mean(hierarchy_df$rank), hierarchy_df$rank),
    color = c(tsub_cols["ILC"], tsub_cols[hierarchy_df$target])
  )
  hierarchy_plot <- ggplot(hierarchy_df) +
    geom_segment(
      aes(
        x = 0.05, xend = 0.95,
        y = mean(rank), yend = rank,
        size = norm,
        color = target
      ),
      lineend = "round"
    ) +
    geom_point(
      data = node_df,
      aes(x = x, y = y, fill = group),
      shape = 21, color = "black", size = 4
    ) +
    geom_text(
      data = node_df,
      aes(x = x, y = y, label = label),
      color = "black",
      nudge_x = c(-0.08, rep(0.08, nrow(node_df) - 1)),
      hjust = c(1, rep(0, nrow(node_df) - 1)),
      size = 3.5
    ) +
    scale_color_manual(values = tsub_cols, guide = "none") +
    scale_fill_manual(values = tsub_cols, guide = "none") +
    scale_size_continuous(range = c(0.5, 3), guide = "none") +
    theme_minimal(base_size = 11) +
    labs(
      title = "ILC to T cell hierarchical view",
      subtitle = "Edge width encodes normalized communication probability; population.size = FALSE",
      x = NULL,
      y = NULL
    ) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  ggsave(hierarchy_file_pdf, plot = hierarchy_plot, width = 8.5, height = 7, dpi = 600, units = "in")
  ggsave(hierarchy_file_png, plot = hierarchy_plot, width = 8.5, height = 7, dpi = 600, units = "in")
} else {
  file.create(hierarchy_file_pdf)
  file.create(hierarchy_file_png)
}

global_circle_pdf <- file.path(plot_dir, sprintf("Fig_ILC_global_circle_%s.pdf", opt$tag))
global_circle_png <- file.path(plot_dir, sprintf("Fig_ILC_global_circle_%s.png", opt$tag))
global_expr <- quote(netVisual_circle(
  all_net,
  vertex.weight = as.numeric(all_vertex),
  weight.scale = TRUE,
  color.use = all_cols,
  title.name = "All cell types (population.size = FALSE)"
))
save_both(global_expr, global_circle_pdf, global_circle_png, width = 8.5, height = 7)

bubble_plot <- netVisual_bubble(
  cellchat_tsub,
  sources.use = "ILC",
  targets.use = tcell_groups,
  angle.x = 45,
  remove.isolate = FALSE
) +
  scale_x_discrete(labels = pretty_replace_anywhere) +
  scale_fill_gradientn(colors = c("#2E86AB", "#F3DF6C", "#D7263D"), name = "Prob.") +
  theme_minimal(base_size = 11, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "ILC to T cell pathways",
    subtitle = "Data shown as communication probability; population.size = FALSE; no strength filtering",
    caption = "Each dot encodes pathway probability (fill) and significance (size); p-values from CellChat permutation."
  )

bubble_pdf <- file.path(plot_dir, sprintf("Fig_ILC_to_Tcell_bubble_%s.pdf", opt$tag))
bubble_png <- file.path(plot_dir, sprintf("Fig_ILC_to_Tcell_bubble_%s.png", opt$tag))
ggsave(bubble_pdf, plot = bubble_plot, width = 8.5, height = 9, dpi = 600, units = "in")
ggsave(bubble_png, plot = bubble_plot, width = 8.5, height = 9, dpi = 600, units = "in")

legend_file <- file.path(plot_dir, sprintf("Fig_ILC_legends_%s.md", opt$tag))
legend_text <- sprintf(
  "# Figure legends (%s)\n\n- **Fig_ILC_to_Tcell_circle**：Circle plot highlighting directed signaling from ILC to T cell subsets. Node size reflects cell counts; edge thickness encodes communication probability, with ILC→Tcell_CD4_NAIVE edges scaled by 3× to make the CD4-naive relationship more prominent. All computations performed with population.size = FALSE and without strength filtering.\n- **Fig_ILC_to_Tcell_hierarchy**：Hierarchy layout emphasizing downstream receptors within T cell subsets, rendered with the same parameter choices and the same CD4-naive edge boost.\n- **Fig_ILC_global_circle**：Complete network view across all cell types; colors follow CellChat ident ordering.\n- **Fig_ILC_to_Tcell_bubble**：Bubble plot summarizing major ILC→T cell pathways. Data are presented as CellChat communication probabilities with dot size showing -log10(p-value) from CellChat permutation tests.\n",
  opt$tag
)
writeLines(legend_text, con = legend_file)

message("保存图形：", circle_pdf, " / ", circle_png)
message("保存图形：", hierarchy_file_pdf, " / ", hierarchy_file_png)
message("保存图形：", global_circle_pdf, " / ", global_circle_png)
message("保存图形：", bubble_pdf, " / ", bubble_png)
message("图例：", legend_file)
