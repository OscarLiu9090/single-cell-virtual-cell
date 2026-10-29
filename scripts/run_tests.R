args0 <- commandArgs(FALSE); script0 <- sub("^--file=", "", args0[grepl("^--file=", args0)])[1]
source(file.path(dirname(script0), "bootstrap.R"), encoding = "UTF-8")
root <- scvp_root(); runtime <- scvp_bootstrap(root)
for (f in c("core.R", "io.R", "network.R", "analysis.R")) source(file.path(root, "scripts/lib", f), encoding = "UTF-8")
check_engine(root)
testthat::test_file(file.path(root, "tests/test_pipeline.R"), reporter = "summary", stop_on_failure = TRUE)
