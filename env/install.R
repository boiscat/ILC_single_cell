#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  "optparse",
  "tidyverse",
  "Seurat",
  "future",
  "patchwork",
  "RColorBrewer",
  "scales",
  "circlize",
  "png",
  "glue",
  "harmony"
)

bioc_pkgs <- c(
  "clusterProfiler",
  "enrichplot",
  "org.Hs.eg.db"
)

install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
}

install_if_missing(cran_pkgs)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

if (!requireNamespace("CellChat", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("sqjin/CellChat")
}


