args0 <- commandArgs(FALSE)
script0 <- sub("^--file=", "", args0[grepl("^--file=", args0)])[1]
source(file.path(dirname(script0), "bootstrap.R"), encoding = "UTF-8")
root <- scvp_root(); runtime <- scvp_bootstrap(root)
for (file in c("core.R", "io.R", "network.R", "analysis.R")) source(file.path(root, "scripts", "lib", file), encoding = "UTF-8")
args <- commandArgs(TRUE)
if ("--help" %in% args || !length(args)) {
  cat("Usage: Rscript scripts/run_perturbation.R --config path.yaml [--preflight]\n")
  quit(status = 0)
}
pos <- match("--config", args)
if (is.na(pos) || pos == length(args)) fail("CONFIG", "--config path.yaml is required.")
config_path <- normalizePath(args[pos + 1], winslash = "/", mustWork = TRUE)
base <- dirname(config_path)
cfg <- yaml::read_yaml(config_path, eval.expr = FALSE)
check_config(cfg)
extra_libs <- cfg$environment$r_libraries %||% character()
.libPaths(unique(c(vapply(extra_libs, resolve_path, character(1), base = base), .libPaths())))
for (p in c("Seurat", "SeuratObject")) if (!requireNamespace(p, quietly = TRUE)) fail("MISSING_PACKAGE", p)
check_engine(root)
out_base <- resolve_path(cfg$output$directory, base)
run_id <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
out <- file.path(out_base, run_id)
if (dir.exists(out)) fail("OUTPUT_EXISTS", "Refusing to overwrite an existing run.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
file.copy(config_path, file.path(out, "config.yaml"))
write_json(cfg, file.path(out, "resolved_config.json"))
logcon <- file(file.path(out, "run.log"), open = "wt", encoding = "UTF-8")
sink(logcon, split = TRUE); sink(logcon, type = "message")
warnings_seen <- character(); statuses <- list(); allresults <- list(); eligibility_all <- list(); status_code <- 0L
append_status <- function(cell_type, condition, holdout, seed, status, detail = "") {
  statuses[[length(statuses) + 1L]] <<- data.frame(cell_type = cell_type, condition = condition, holdout = holdout, seed = seed, status = status, detail = detail)
}
tryCatch(withCallingHandlers({
  cat("SCVP run", run_id, "\n")
  cat("Technical demonstration:", isTRUE(cfg$analysis$technical_demo), "\n")
  data <- load_input(cfg, base, out, root, runtime)
  write_json(list(input_type = cfg$input$type, genes = nrow(data$counts), cells = ncol(data$counts), cell_ids_sha256 = digest::digest(colnames(data$counts), algo = "sha256"), gene_ids_sha256 = digest::digest(rownames(data$counts), algo = "sha256"), counts_sha256 = digest::digest(data$counts, algo = "sha256"), metadata_sha256 = digest::digest(data$meta, algo = "sha256")), file.path(out, "input_summary.json"))
  qc <- qc_counts(data$counts, data$meta, cfg$qc, cfg$analysis$species)
  write_tsv(qc$qc, file.path(out, "cell_qc.tsv"))
  x <- qc$counts; meta <- qc$meta; rm(data, qc); invisible(gc())
  write_tsv(meta, file.path(out, "metadata_after_qc.tsv"))
  if (anyNA(meta$condition) || anyNA(meta$cell_type)) fail("ANNOTATION_REQUIRED", "condition/cell_type missing; annotate and review cell identity before modeling.")
  chosen <- meta$cell_type %in% cfg$analysis$cell_types & meta$condition %in% cfg$analysis$conditions
  if (!any(chosen)) fail("NO_SELECTED_CELLS", "No cells match requested cell types/conditions.")
  x <- x[, chosen, drop = FALSE]; meta <- meta[chosen, , drop = FALSE]
  pb <- donor_expression(x, meta, cfg, out)
  score_signatures(pb, cfg, base, out); rm(pb)
  for (ct in cfg$analysis$cell_types) for (cond in cfg$analysis$conditions) {
    ids <- which(meta$cell_type == ct & meta$condition == cond)
    if (!length(ids)) { append_status(ct, cond, "none", NA, "no_selected_cells"); next }
    xm <- x[, ids, drop = FALSE]; mm <- meta[ids, , drop = FALSE]
    groupdir <- file.path(out, paste0(safe_label(ct), "__", safe_label(cond)))
    dir.create(groupdir, recursive = TRUE, showWarnings = FALSE)
    gm <- gene_metrics(xm, unique(c(cfg$analysis$target_gene, cfg$analysis$readout_genes)))
    gm$status[gm$status == "ok" & gm$detection < cfg$qc$min_target_detection] <- "low_detection"
    write_tsv(gm, file.path(groupdir, "target_readout_eligibility.tsv"))
    eg <- gm; eg$cell_type <- ct; eg$condition <- cond
    eligibility_all[[length(eligibility_all) + 1L]] <- eg
    if (gm$status[gm$gene == cfg$analysis$target_gene] != "ok") {
      append_status(ct, cond, "none", NA, "target_not_evaluable", gm$status[gm$gene == cfg$analysis$target_gene]); next
    }
    for (g in cfg$analysis$readout_genes) if (gm$status[gm$gene == g] != "ok") append_status(ct, cond, "none", NA, paste0("readout_", gm$status[gm$gene == g]), g)
    if (anyNA(mm$donor_id)) { append_status(ct, cond, "none", NA, "missing_donor_id"); next }
    if (ncol(xm) < cfg$qc$min_cells_network) { append_status(ct, cond, "none", NA, "insufficient_cells"); next }
    nd <- length(unique(mm$donor_id))
    if (nd < cfg$statistics$min_donors_network) {
      append_status(ct, cond, "none", NA, "insufficient_donors", paste(nd, "donors; exploratory only"))
      if (!isTRUE(cfg$analysis$allow_exploratory_low_donors)) next
    }
    if ("--preflight" %in% args) { append_status(ct, cond, "none", NA, "preflight_ready"); next }
    holdouts <- c("none", if (isTRUE(cfg$network$leave_one_donor_out) && nd >= 3) sort(unique(mm$donor_id)))
    if (isTRUE(cfg$network$leave_one_donor_out) && nd < 3) append_status(ct, cond, "none", NA, "LODO_unavailable", "Requires >=3 donors.")
    for (holdout in holdouts) {
      use <- if (holdout == "none") rep(TRUE, nrow(mm)) else mm$donor_id != holdout
      xx <- xm[, use, drop = FALSE]; md <- mm[use, , drop = FALSE]
      seeds <- if (holdout == "none") cfg$network$seeds else cfg$network$seeds[1]
      for (seed in seeds) {
        sd <- file.path(groupdir, paste0("holdout_", safe_label(holdout)), paste0("seed_", seed))
        dir.create(sd, recursive = TRUE, showWarnings = FALSE)
        tryCatch({
          # Refit variable genes inside the retained donor set to avoid holdout leakage.
          eligible_readouts <- gm$gene[gm$status == "ok" & gm$gene %in% cfg$analysis$readout_genes]
          genes <- select_genes(xx, unique(c(cfg$analysis$target_gene, eligible_readouts)), cfg$network$n_hvg, cfg$qc$min_cells_gene, cfg$network$exclude_genes_regex)
          controls <- match_controls(xx[genes, , drop = FALSE], cfg$analysis$target_gene, c(cfg$analysis$target_gene, cfg$analysis$readout_genes), cfg$network$n_control_genes)
          write_tsv(data.frame(gene = genes, role = ifelse(genes == cfg$analysis$target_gene, "target", ifelse(genes %in% cfg$analysis$readout_genes, "readout", ifelse(genes %in% controls, "matched_control", "network")))) , file.path(sd, "network_genes.tsv"))
          res <- run_network_seed(xx, md, genes, cfg, seed, sd, controls)
          res$cell_type <- ct; res$condition <- cond; res$holdout <- holdout
          allresults[[length(allresults) + 1L]] <- res
          append_status(ct, cond, holdout, seed, "completed")
        }, error = function(e) append_status(ct, cond, holdout, seed, "failed", conditionMessage(e)))
      }
    }
  }
  results <- if (length(allresults)) do.call(rbind, allresults) else NULL
  if (!is.null(results)) {
    saveRDS(results, file.path(out, "all_results.rds"))
    summarize_stability(results, cfg, out, if (length(eligibility_all)) do.call(rbind, eligibility_all) else NULL)
    run_enrichment(results, cfg, base, out)
  }
  if (!length(allresults) && !"--preflight" %in% args) status_code <- 2L
  if (any(vapply(statuses, function(s) s$status == "failed", logical(1)))) status_code <- 2L
}, warning = function(w) { warnings_seen <<- c(warnings_seen, conditionMessage(w)); invokeRestart("muffleWarning") }), error = function(e) {
  status_code <<- 1L
  append_status("", "", "", NA, "fatal", conditionMessage(e))
})
if (length(statuses)) write_tsv(do.call(rbind, statuses), file.path(out, "run_status.tsv"))
if (length(eligibility_all)) write_tsv(do.call(rbind, eligibility_all), file.path(out, "target_readout_eligibility.tsv"))
writeLines(unique(warnings_seen), file.path(out, "warnings.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"), useBytes = TRUE)
ip <- installed.packages()
write_tsv(data.frame(package = ip[, "Package"], version = ip[, "Version"], library = ip[, "LibPath"]), file.path(out, "installed_packages.tsv"))
write_json(list(exit_code = status_code, technical_demo = isTRUE(cfg$analysis$technical_demo), completed_network_runs = length(allresults), expression_direction = "not_estimated", interpretation = "Network FC is squared manifold distance relative to mean squared non-target distance; not expression fold change."), file.path(out, "run_summary.json"))
report <- c("# 单细胞虚拟扰动运行报告", "", paste0("运行：", run_id), paste0("技术演示：", isTRUE(cfg$analysis$technical_demo)), paste0("完成网络运行数：", length(allresults)), paste0("退出状态：", status_code), "", "网络结果：readout_summary.tsv（如有）、readout_per_run.tsv、stability_all_genes.tsv。", "供体表达与关联：donor_expression.tsv、donor_association.tsv。", "可评估性及失败原因：run_status.tsv、各组 target_readout_eligibility.tsv。", "", "PPARA 等关注基因的表达方向尚不能判断：本流程没有经过适配扰动数据验证的表达响应模型。", "FC、Z 和 distance 衡量网络变化，不代表 mRNA 上下调；不同网络的原始距离不直接比较。", "随机种子和供体留一结果是敏感性检查，不能作为新增生物学重复。", "没有显著基因不等于证明没有效应；未进入网络和不显著必须分别记录。", "", "详见配置、sessionInfo.txt、warnings.txt 及 run.log。")
writeLines(enc2utf8(report), file.path(out, "report_zh.md"), useBytes = TRUE)
sink(type = "message"); sink(); close(logcon)
cat("Output:", out, "\nExit code:", status_code, "\n")
quit(status = status_code)
