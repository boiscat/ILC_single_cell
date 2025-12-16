#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- "."
root_arg <- grep("^--root=", args, value = TRUE)
if (length(root_arg) == 1) {
  root <- sub("^--root=", "", root_arg)
}

files <- list.files(root, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("/\\.git/", files)]
files <- sort(unique(files))

if (length(files) == 0) {
  cat("No .R files found under:", root, "\n")
  quit(status = 1)
}

bad <- 0L
for (path in files) {
  parsed <- tryCatch(
    parse(file = path, encoding = "UTF-8"),
    error = function(e) e
  )
  if (inherits(parsed, "error")) {
    bad <- bad + 1L
    cat("[PARSE FAIL]", path, "\n")
    cat("  ", conditionMessage(parsed), "\n", sep = "")
  }
}

if (bad > 0) {
  cat("FAILED:", bad, "file(s) could not be parsed.\n")
  quit(status = 2)
}

cat("OK:", length(files), "file(s) parsed successfully.\n")
