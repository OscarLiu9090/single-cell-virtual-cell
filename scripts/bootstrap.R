# Shared startup; source before loading any analysis package.
if (.Platform$OS.type == "windows") {
  Sys.setenv(LANG = "English_United States.utf8", LC_ALL = "English_United States.utf8")
  suppressWarnings(try(Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent = TRUE))
}
scvp_root <- function() {
  args <- commandArgs(FALSE)
  script <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (!length(script)) stop("Run this entry point with Rscript.")
  normalizePath(file.path(dirname(script[1]), ".."), winslash = "/", mustWork = TRUE)
}
scvp_bootstrap <- function(root) {
  local <- file.path(root, "references", "environment.local.json")
  paths <- strsplit(Sys.getenv("SCVP_R_LIBS"), ";", fixed = TRUE)[[1]]
  # jsonlite is part of the minimum runtime; bootstrap uses base R path hints first.
  hints <- file.path(root, "references", "library_paths.local.txt")
  if (file.exists(hints)) paths <- c(paths, readLines(hints, encoding = "UTF-8", warn = FALSE))
  paths <- paths[nzchar(paths) & dir.exists(paths)]
  .libPaths(unique(c(paths, .libPaths())))
  for (p in c("yaml", "jsonlite", "Matrix", "digest"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p, ". See references/environment.md")
  if (file.exists(local)) return(jsonlite::fromJSON(local))
  list()
}
