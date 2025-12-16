#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(dplyr)
})

parse_args <- function(arg_vec) {
  res <- list(tag = format(Sys.Date(), "%Y%m%d"))
  for (arg in arg_vec) {
    if (startsWith(arg, "--tag=")) {
      res$tag <- sub("^--tag=", "", arg)
    }
    if (startsWith(arg, "--input=")) {
      res$input_rds <- sub("^--input=", "", arg)
    }
    if (startsWith(arg, "--outdir=")) {
      res$outdir <- sub("^--outdir=", "", arg)
    }
    if (startsWith(arg, "--prefix=")) {
      res$prefix <- sub("^--prefix=", "", arg)
    }
  }
  res
}
`%||%` <- function(a, b) if (!is.null(a)) a else b

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
log_dir <- file.path(base_dir, "logs")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

opt$input_rds <- opt$input_rds %||% file.path(object_dir, sprintf("cellchat_ILC_Tsub_%s.rds", opt$tag))
opt$outdir <- opt$outdir %||% plot_dir
opt$prefix <- opt$prefix %||% "Fig_ILC"
opt$log_file <- file.path(log_dir, sprintf("04-cellchatPlotExtensions-%s.log", opt$tag))

log_message <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 paste0(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = opt$log_file, append = TRUE)
}

log_message("读取 CellChat 对象：", opt$input_rds)
stopifnot(file.exists(opt$input_rds))
cellchat <- readRDS(opt$input_rds)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

plot_save <- function(expr, pdf_path, png_path, width = 8, height = 6) {
  pdf(pdf_path, width = width, height = height)
  eval(expr)
  dev.off()
  png(png_path, width = width, height = height, units = "in", res = 600)
  eval(expr)
  dev.off()
  log_message("已输出：", pdf_path, " & ", png_path)
}

circle_expr <- quote(netVisual_circle(cellchat@net$weight,
                                      vertex.weight = apply(cellchat@net$weight, 1, sum),
                                      weight.scale = TRUE,
                                      label.edge = FALSE,
                                      title.name = "ILC/T cell network"))
plot_save(circle_expr,
          file.path(opt$outdir, sprintf("%s_global_circle_%s.pdf", opt$prefix, opt$tag)),
          file.path(opt$outdir, sprintf("%s_global_circle_%s.png", opt$prefix, opt$tag)),
          width = 8, height = 8)

heat_expr <- quote(netVisual_heatmap(cellchat, measure = "weight"))
plot_save(heat_expr,
          file.path(opt$outdir, sprintf("%s_to_Tcell_heatmap_%s.pdf", opt$prefix, opt$tag)),
          file.path(opt$outdir, sprintf("%s_to_Tcell_heatmap_%s.png", opt$prefix, opt$tag)),
          width = 9, height = 6)

scatter_expr <- quote(netAnalysis_signalingRole_scatter(cellchat) +
                        ggtitle("ILC/T cell signaling roles"))
plot_save(scatter_expr,
          file.path(opt$outdir, sprintf("%s_signalingRole_scatter_%s.pdf", opt$prefix, opt$tag)),
          file.path(opt$outdir, sprintf("%s_signalingRole_scatter_%s.png", opt$prefix, opt$tag)),
          width = 7, height = 5)

comm_df <- subsetCommunication(cellchat)
ilc_targets <- comm_df %>%
  filter(source %in% c("ILC")) %>%
  group_by(target) %>%
  summarise(prob_sum = sum(prob), .groups = "drop") %>%
  arrange(desc(prob_sum))
bar_plot <- ggplot(ilc_targets, aes(x = reorder(target, prob_sum), y = prob_sum)) +
  geom_col(fill = "#4DBBD5", color = "black", linewidth = 0.2) +
  coord_flip() +
  labs(x = "Target T cell subtype", y = "Total probability",
       title = "ILC → T cell communication strength") +
  theme_minimal(base_size = 12)
ggsave(file.path(opt$outdir, sprintf("%s_ILC_to_Tcell_strength_%s.pdf", opt$prefix, opt$tag)),
       bar_plot, width = 6, height = 4)
ggsave(file.path(opt$outdir, sprintf("%s_ILC_to_Tcell_strength_%s.png", opt$prefix, opt$tag)),
       bar_plot, width = 6, height = 4, dpi = 600)
log_message("ILC→T cell 强度条形图输出完成。")

log_message("扩展绘图执行结束。sessionInfo：")
for (line in capture.output(sessionInfo())) {
  log_message(line)
}
