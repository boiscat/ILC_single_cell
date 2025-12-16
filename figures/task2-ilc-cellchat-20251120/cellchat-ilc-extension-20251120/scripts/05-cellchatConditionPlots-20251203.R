#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(grid)
  library(png)
  library(circlize)
  library(RColorBrewer)
})

parse_args <- function(arg_vec) {
  res <- list(tag = format(Sys.Date(), "%Y%m%d"), conditions = "lesion,non_lesion")
  i <- 1
  while (i <= length(arg_vec)) {
    token <- arg_vec[[i]]
    if (startsWith(token, "--")) {
      token_clean <- substring(token, 3)
      key <- NA_character_
      value <- NULL
      if (grepl("=", token_clean, fixed = TRUE)) {
        parts <- strsplit(token_clean, "=", fixed = TRUE)[[1]]
        key <- parts[1]
        if (length(parts) > 1 && nzchar(parts[2])) {
          value <- parts[2]
        }
      } else {
        key <- token_clean
        if (i < length(arg_vec)) {
          next_token <- arg_vec[[i + 1]]
          if (!startsWith(next_token, "--")) {
            value <- next_token
            i <- i + 1
          }
        }
      }
      if (!is.na(key)) {
        if (is.null(value)) value <- TRUE
        res[[key]] <- value
      }
    }
    i <- i + 1
  }
  res
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

combine_ggplots <- function(plot_list, pdf_path, png_path, width, height) {
  if (length(plot_list) == 0) {
    warning("没有可用的散点图，跳过：", pdf_path)
    return(invisible(NULL))
  }
  draw_panels <- function(device_fun, ...) {
    device_fun(..., width = width, height = height)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1, length(plot_list))))
    for (i in seq_along(plot_list)) {
      print(plot_list[[i]], vp = viewport(layout.pos.row = 1, layout.pos.col = i))
    }
    dev.off()
  }
  draw_panels(pdf, file = pdf_path)
  draw_panels(png, filename = png_path, units = "in", res = 600)
}

combine_images <- function(image_paths, pdf_path, png_path, width, height) {
  if (length(image_paths) == 0) {
    warning("没有可用的弦图，跳过：", pdf_path)
    return(invisible(NULL))
  }
  imgs <- lapply(image_paths, readPNG)
  draw_panels <- function(device_fun, ...) {
    device_fun(..., width = width, height = height)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1, length(imgs))))
    for (i in seq_along(imgs)) {
      grid.raster(imgs[[i]], vp = viewport(layout.pos.row = 1, layout.pos.col = i))
    }
    dev.off()
  }
  draw_panels(pdf, file = pdf_path)
  draw_panels(png, filename = png_path, units = "in", res = 600)
}

draw_chord_plot <- function(cellchat_obj, title_text, pdf_path, png_path) {
  net_mat <- cellchat_obj@net$weight
  if (is.null(net_mat)) {
    warning("CellChat 对象缺少 net$weight，跳过：", title_text)
    return()
  }
  keep <- rowSums(net_mat) + colSums(net_mat) > 0
  net_mat <- net_mat[keep, keep, drop = FALSE]
  if (nrow(net_mat) == 0) {
    warning("没有可视化的互作：", title_text)
    return()
  }
  groups <- rownames(net_mat)
  pal_fun <- colorRampPalette(brewer.pal(min(8, length(groups)), "Set3"))
  group_cols <- setNames(pal_fun(length(groups)), groups)
  devices <- list(
    list(fun = pdf, args = list(file = pdf_path, width = 6, height = 6)),
    list(fun = png, args = list(filename = png_path, width = 6, height = 6, units = "in", res = 600))
  )
  for (dev in devices) {
    do.call(dev$fun, dev$args)
    circos.clear()
    circos.par(track.margin = c(0.01, 0.01), cell.padding = c(0, 0, 0, 0), start.degree = 90)
    chordDiagram(
      net_mat,
      directional = 1,
      direction.type = c("diffHeight", "arrows"),
      link.arr.type = "big.arrow",
      grid.col = group_cols,
      transparency = 0.35,
      annotationTrack = "grid",
      preAllocateTracks = list(track.height = 0.08)
    )
    circos.track(track.index = 1, panel.fun = function(x, y) {
      sector <- get.cell.meta.data("sector.index")
      circos.text(mean(x), y[1] + 0.15, sector, facing = "clockwise", niceFacing = TRUE, cex = 0.5)
    }, bg.border = NA)
    grid.text(title_text,
              x = 0.5, y = unit(0.98, "npc"),
              gp = gpar(fontface = "bold"))
    grid.text(sprintf("Total cells: %s", nrow(cellchat_obj@meta)),
              x = 0.5, y = unit(0.03, "npc"),
              gp = gpar(fontsize = 8))
    dev.off()
    circos.clear()
  }
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
script_dir <- get_script_dir()
base_dir <- normalizePath(file.path(script_dir, ".."))
object_dir <- file.path(base_dir, "outputs", "objects")
plot_dir <- file.path(base_dir, "outputs", "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
opt$prefix <- opt$prefix %||% "Fig_condition"

condition_vec <- unique(trimws(strsplit(opt$conditions, ",")[[1]]))
condition_vec <- condition_vec[condition_vec != ""]
if (length(condition_vec) == 0) {
  stop("请通过 --conditions 指定至少一个条件")
}

condition_objs <- list()
for (cond in condition_vec) {
  path <- file.path(object_dir, sprintf("cellchat_ILC_all_%s_%s.rds", cond, opt$tag))
  if (!file.exists(path)) {
    warning("缺少 CellChat 对象：", path)
    next
  }
  condition_objs[[cond]] <- readRDS(path)
}

if (length(condition_objs) == 0) {
  stop("未能加载任何条件的 CellChat 对象，请确认 02 脚本已输出对应文件。")
}

scatter_plots <- list()
for (cond in names(condition_objs)) {
  cc <- condition_objs[[cond]]
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  condition_objs[[cond]] <- cc
  total_cells <- nrow(cc@meta)
  scatter_plots[[cond]] <- netAnalysis_signalingRole_scatter(cc) +
    ggtitle(sprintf("%s (n = %s)", cond, total_cells)) +
    theme(
      plot.title = element_text(face = "bold"),
      text = element_text(family = "Helvetica")
    )
}

scatter_width <- max(6, 5 * length(scatter_plots))
scatter_pdf <- file.path(plot_dir, sprintf("%s_signalingRole_%s.pdf", opt$prefix, opt$tag))
scatter_png <- file.path(plot_dir, sprintf("%s_signalingRole_%s.png", opt$prefix, opt$tag))
combine_ggplots(scatter_plots, scatter_pdf, scatter_png, width = scatter_width, height = 5)

chord_panel_paths <- list()
for (cond in names(condition_objs)) {
  cc <- condition_objs[[cond]]
  single_pdf <- file.path(plot_dir, sprintf("%s_chord_%s_%s.pdf", opt$prefix, cond, opt$tag))
  single_png <- file.path(plot_dir, sprintf("%s_chord_%s_%s.png", opt$prefix, cond, opt$tag))
  draw_chord_plot(cc, sprintf("%s network", cond), single_pdf, single_png)
  chord_panel_paths[[cond]] <- single_png
}

combined_width <- max(6, 6 * length(chord_panel_paths))
combined_pdf <- file.path(plot_dir, sprintf("%s_chord_grid_%s.pdf", opt$prefix, opt$tag))
combined_png <- file.path(plot_dir, sprintf("%s_chord_grid_%s.png", opt$prefix, opt$tag))
combine_images(chord_panel_paths, combined_pdf, combined_png, width = combined_width, height = 6)

legend_path <- file.path(plot_dir, sprintf("%s_panels_%s.md", opt$prefix, opt$tag))
legend_lines <- c(
  sprintf("# Condition-specific CellChat panels (%s | prefix = %s)", opt$tag, opt$prefix),
  "",
  "- **Fig_condition_signalingRole**：scatter panel comparing incoming vs outgoing centrality per cell state across指定条件。",
  "- **Fig_condition_chord_*_***：单个条件的弦图（PDF/PNG）。",
  "- **Fig_condition_chord_grid**：自动拼接的弦图多 panel，便于与示例 A/B 对照。"
)
writeLines(legend_lines, con = legend_path)

message("散点图输出：", scatter_pdf, " / ", scatter_png)
message("弦图多 panel 输出：", combined_pdf, " / ", combined_png)
message("图例：", legend_path)
