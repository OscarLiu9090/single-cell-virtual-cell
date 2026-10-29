perturb_edges <- function(wt, gene, factor) {
  if (!gene %in% rownames(wt) || !identical(rownames(wt), colnames(wt))) fail("NETWORK_ID", "Target absent or adjacency names misaligned.")
  if (!is.finite(factor) || factor < 0) fail("PERTURBATION_FACTOR", "Factor must be nonnegative.")
  perturbed <- wt
  # Verified scTenifoldKnk convention AFTER its transpose: row = source.
  perturbed[gene, ] <- wt[gene, ] * factor
  perturbed
}
null_result <- function(genes, reason) data.frame(gene = genes, distance = 0, Z = NA_real_, FC = NA_real_, p.value = 1, p.adj = 1, result_status = reason)
compare_network <- function(wt, perturbed, gene, cfg, seed) {
  if (identical(wt, perturbed) || max(abs(wt - perturbed)) == 0) return(null_result(rownames(wt), "unchanged_network"))
  set.seed(seed)
  ma <- scTenifoldNet::manifoldAlignment(wt, perturbed, d = cfg$network$ma_dim, nCores = cfg$network$n_cores)
  dr <- getFromNamespace("dRegulation", "scTenifoldKnk")(ma, gKO = gene, empiricalNull = FALSE, direction = NULL)
  dr$result_status <- "ok"
  dr
}
run_network_seed <- function(x, meta, genes, cfg, seed, out, controls) {
  n <- cfg$network; target <- cfg$analysis$target_gene
  cells <- balanced_cells(meta, seed, n$max_cells_per_donor)
  write_tsv(meta[match(cells, meta$cell_id), c("cell_id", "sample_id", "donor_id", "condition", "cell_type")], file.path(out, "sampled_cells.tsv"))
  y <- x[genes, cells, drop = FALSE]
  gm <- gene_metrics(y, genes)
  y <- y[gm$gene[gm$status == "ok"], , drop = FALSE]
  if (!target %in% rownames(y)) fail("TARGET_LOST", "Target has no variation in this balanced sample.")
  if (ncol(y) < cfg$qc$min_cells_network || nrow(y) <= n$n_comp + 2) fail("INSUFFICIENT_CELLS", "Too few cells/genes in balanced sample.")
  estimated_gb <- 8 * nrow(y)^2 * (4 * n$n_net + 20) / 1024^3
  write_json(list(genes = nrow(y), cells = ncol(y), estimated_peak_gb = estimated_gb, budget_gb = n$memory_budget_gb), file.path(out, "memory_estimate.json"))
  if (estimated_gb > n$memory_budget_gb) fail("MEMORY_BUDGET", "Estimated network workspace exceeds configured budget; reduce n_hvg explicitly.")
  # Core package normalizes selected counts to CPM. Use this same WT for every OE/control.
  result <- scTenifoldKnk::scTenifoldKnk(countMatrix = y, gKO = target, qc = FALSE,
    nc_nNet = n$n_net, nc_nCells = min(n$n_cells, ncol(y)), nc_nComp = n$n_comp,
    nc_scaleScores = TRUE, nc_symmetric = FALSE, nc_q = n$q, nc_lambda = 0,
    td_K = n$td_k, td_maxIter = n$td_max_iter, td_maxError = n$td_max_error,
    td_nDecimal = n$td_decimals, ma_nDim = n$ma_dim, ma_method = "manifold",
    dr_empiricalNull = FALSE, dr_direction = FALSE, nCores = n$n_cores, seed = seed)
  wt <- as.matrix(result$tensorNetworks$WT)
  if (any(!is.finite(wt))) fail("NONFINITE_NETWORK", "WT contains non-finite edges.")
  if (isTRUE(cfg$output$save_networks)) saveRDS(wt, file.path(out, "WT.rds"))
  targets <- unique(c(target, controls[controls %in% rownames(wt)]))
  combined <- list(); idx <- 0L
  for (g in targets) {
    factors <- if (g == target) unique(c(0, 1, cfg$analysis$oe_factors)) else 0
    for (f in factors) {
      mode <- if (f == 0) "KO" else if (f == 1) "sham" else "OE_edge_scaling"
      mutated <- perturb_edges(wt, g, f)
      if (g == target && f == 0 && any(wt[g, ] != 0)) {
        dr <- result$diffRegulation; dr$result_status <- "ok"
      } else dr <- compare_network(wt, mutated, g, cfg, seed)
      if (any(!is.finite(dr$distance)) || any(!is.finite(dr$p.adj))) fail("NONFINITE_STATISTICS", "Invalid distance or P value.")
      dr$is_target <- dr$gene == g
      dr$rank <- NA_real_
      informative <- !dr$is_target & dr$result_status == "ok"
      if (any(informative)) dr$rank[informative] <- rank(-dr$distance[informative], ties.method = "min")
      dr$rank_fraction <- dr$rank / sum(!dr$is_target)
      dr$perturbed_gene <- g; dr$mode <- mode; dr$factor <- f; dr$seed <- seed
      dr$is_control <- g != target
      dr$expression_direction <- NA_character_
      dr$direction_status <- "not_estimated_no_validated_perturbation_model"
      dr$n_genes <- nrow(wt); dr$n_cells <- ncol(y); dr$n_donors <- length(unique(meta$donor_id))
      fname <- paste0(safe_label(g), "_", mode, "_", f)
      write_tsv(dr, file.path(out, paste0(fname, "_all.tsv")))
      sig <- dr[!dr$is_target & dr$p.adj < cfg$statistics$fdr, , drop = FALSE]
      write_tsv(sig, file.path(out, paste0(fname, "_significant.tsv")))
      if (g == target && f != 1) plot_network_result(dr, file.path(out, paste0(fname, ".pdf")), cfg$analysis$readout_genes)
      idx <- idx + 1L; combined[[idx]] <- dr
    }
  }
  do.call(rbind, combined)
}
plot_network_result <- function(dr, path, readouts) {
  grDevices::pdf(path, width = 8, height = 4.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  d <- dr[!dr$is_target & is.finite(dr$distance) & is.finite(dr$rank), , drop = FALSE]
  if (!nrow(d)) { graphics::plot.new(); graphics::title("No informative network change"); graphics::text(.5, .5, "Unchanged network / unavailable rank; see result_status.", cex = .8); return(invisible()) }
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  graphics::plot(d$rank, d$distance, pch = 16, cex = 0.5, col = ifelse(d$gene %in% readouts, "#D55E00", "#777777"), xlab = "Network perturbation rank", ylab = "Manifold distance")
  k <- d$gene %in% readouts
  graphics::text(d$rank[k], d$distance[k], labels = d$gene[k], pos = 3, cex = 0.7)
  graphics::plot(d$distance, -log10(pmax(d$p.adj, .Machine$double.xmin)), pch = 16, cex = 0.5, xlab = "Network distance (unsigned)", ylab = "-log10(adjusted P)")
}
