args0 <- commandArgs(FALSE)
script0 <- sub("^--file=", "", args0[grepl("^--file=", args0)])[1]
source(file.path(dirname(script0), "bootstrap.R"), encoding = "UTF-8")
root <- scvp_root(); runtime <- scvp_bootstrap(root)
args <- commandArgs(TRUE)
i <- match("--lib", args)
if (is.na(i) || i == length(args)) stop("Usage: Rscript scripts/install_engine.R --lib /path/to/project/library")
lib <- normalizePath(args[i + 1], winslash = "/", mustWork = FALSE)
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
manifest <- jsonlite::fromJSON(file.path(root, "references", "upstream", "manifest.json"))
for (j in seq_len(nrow(manifest))) {
  path <- file.path(root, "references", "upstream", manifest$repo[j], manifest$path[j])
  if (digest::digest(file = path, algo = "sha256") != manifest$sha256[j]) stop("Source hash mismatch: ", path)
}
for (pkg in c("scTenifoldNet", "scTenifoldKnk")) {
  src <- file.path(root, "references", "upstream", pkg)
  desc <- read.dcf(file.path(src, "DESCRIPTION"))
  imports <- trimws(strsplit(gsub("\\([^)]*\\)", "", desc[1, "Imports"]), ",")[[1]])
  missing <- imports[!vapply(imports, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing dependencies in the selected library first: ", paste(missing, collapse = ", "))
  utils::install.packages(src, repos = NULL, type = "source", lib = lib, INSTALL_opts = c("--no-multiarch", "--no-docs", "--no-help"))
  if (!file.exists(file.path(lib, pkg, "DESCRIPTION"))) stop("Installation failed: ", pkg)
}
cat("Installed pinned source packages in", lib, "\nRestart R before use.\n")
