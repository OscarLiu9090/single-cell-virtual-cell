donor_expression <- function(x, meta, cfg, out) {
  keys <- interaction(meta$cell_type, meta$condition, meta$donor_id, drop = TRUE, lex.order = TRUE)
  valid <- !is.na(keys)
  if (!any(valid)) return(invisible(NULL))
  keys <- droplevels(keys[valid])
  design <- Matrix::sparseMatrix(i = seq_along(keys), j = as.integer(keys), x = 1, dims = c(length(keys), nlevels(keys)))
  sums <- x[, valid, drop = FALSE] %*% design
  totals <- Matrix::colSums(sums)
  md <- meta[valid, , drop = FALSE][match(levels(keys), as.character(keys)), c("donor_id", "condition", "cell_type"), drop = FALSE]
  md$n_cells <- tabulate(as.integer(keys), nbins = nlevels(keys)); md$total_counts <- totals
  genes <- unique(c(cfg$analysis$target_gene, cfg$analysis$readout_genes))
  records <- lapply(genes, function(g) {
    found <- g %in% rownames(sums)
    d <- md; d$gene <- g
    d$count_sum <- if (found) as.numeric(sums[g, ]) else NA_real_
    d$log1p_cpm <- if (found) log1p(d$count_sum / pmax(totals, 1) * 1e6) else NA_real_
    d$evidence <- "observational_donor_pseudobulk"
    d
  })
  tab <- do.call(rbind, records)
  write_tsv(tab, file.path(out, "donor_expression.tsv"))
  associations <- list()
  for (ct in unique(md$cell_type)) for (cond in unique(md$condition)) for (readout in cfg$analysis$readout_genes) {
    a <- tab[tab$cell_type == ct & tab$condition == cond & tab$gene == cfg$analysis$target_gene, ]
    b <- tab[tab$cell_type == ct & tab$condition == cond & tab$gene == readout, ]
    a <- a[order(a$donor_id), ]; b <- b[match(a$donor_id, b$donor_id), ]
    ok <- is.finite(a$log1p_cpm) & is.finite(b$log1p_cpm)
    status <- "ok"; rho <- p <- NA_real_
    if (sum(ok) < cfg$statistics$min_donors_association) status <- "insufficient_donors" else if (stats::sd(a$log1p_cpm[ok]) == 0 || stats::sd(b$log1p_cpm[ok]) == 0) status <- "no_variation" else {
      fit <- suppressWarnings(stats::cor.test(a$log1p_cpm[ok], b$log1p_cpm[ok], method = "spearman", exact = FALSE))
      rho <- unname(fit$estimate); p <- fit$p.value
    }
    associations[[length(associations) + 1L]] <- data.frame(cell_type = ct, condition = cond, target = cfg$analysis$target_gene, readout = readout, n_donors = sum(ok), spearman_rho = rho, p_value = p, status = status, evidence = "association_not_intervention")
  }
  ass <- do.call(rbind, associations)
  ass$p_adj <- stats::p.adjust(ass$p_value, "BH")
  write_tsv(ass, file.path(out, "donor_association.tsv"))
  invisible(list(sums = sums, metadata = md))
}
summarize_stability <- function(results, cfg, out, eligibility = NULL) {
  if (is.null(results) || !nrow(results)) return(invisible(NULL))
  readouts <- cfg$analysis$readout_genes
  readout_rows <- list()
  run_groups <- split(results, interaction(results$cell_type, results$condition, results$holdout, results$perturbed_gene, results$mode, results$factor, results$seed, drop = TRUE))
  for (d in run_groups) for (g in readouts) {
    row <- d[d$gene == g, , drop = FALSE]
    if (!nrow(row)) {
      row <- d[1, , drop = FALSE]; row$gene <- g; row$is_target <- g == row$perturbed_gene
      row$result_status <- "not_in_network"
      for (nm in c("distance", "Z", "FC", "p.value", "p.adj", "rank", "rank_fraction")) row[[nm]] <- NA_real_
    }
    if (!is.null(eligibility)) {
      e <- eligibility[eligibility$gene == g & eligibility$cell_type == row$cell_type & eligibility$condition == row$condition, ]
      if (nrow(e) && e$status[1] != "ok") {
        row$result_status <- e$status[1]
        for (nm in c("distance", "Z", "FC", "p.value", "p.adj", "rank", "rank_fraction")) row[[nm]] <- NA_real_
      }
    }
    readout_rows[[length(readout_rows) + 1L]] <- row
  }
  readout_rows <- do.call(rbind, readout_rows)
  write_tsv(readout_rows, file.path(out, "readout_per_run.tsv"))
  results_for_summary <- rbind(results[!results$gene %in% readouts, ], readout_rows)
  # Seeds and leave-one-donor runs are sensitivity analyses, not extra replicates.
  group <- interaction(results_for_summary$cell_type, results_for_summary$condition, results_for_summary$holdout, results_for_summary$perturbed_gene, results_for_summary$mode, results_for_summary$factor, results_for_summary$gene, drop = TRUE, lex.order = TRUE)
  summary <- lapply(split(results_for_summary, group), function(d) {
    z <- d[1, c("cell_type", "condition", "holdout", "perturbed_gene", "mode", "factor", "gene", "is_target", "is_control")]
    z$n_algorithm_runs <- nrow(d)
    z$result_status <- paste(sort(unique(d$result_status)), collapse = ";")
    z$unchanged_network_fraction <- mean(d$result_status == "unchanged_network")
    z$median_rank_fraction <- stats::median(d$rank_fraction, na.rm = TRUE)
    z$rank_fraction_min <- if (all(is.na(d$rank_fraction))) NA_real_ else min(d$rank_fraction, na.rm = TRUE)
    z$rank_fraction_max <- if (all(is.na(d$rank_fraction))) NA_real_ else max(d$rank_fraction, na.rm = TRUE)
    z$median_distance <- stats::median(d$distance)
    z$significant_fraction <- mean(d$p.adj < cfg$statistics$fdr)
    z$expected_seeds <- if (d$holdout[1] == "none") length(cfg$network$seeds) else 1L
    z$complete <- z$n_algorithm_runs == z$expected_seeds
    z$direction <- "not_estimated"
    z
  })
  summary <- do.call(rbind, summary)
  write_tsv(summary, file.path(out, "stability_all_genes.tsv"))
  write_tsv(summary[summary$gene %in% readouts, ], file.path(out, "readout_summary.tsv"))
  # Full-run pairwise ranking agreement; never pool P values.
  primary <- results[results$holdout == "none" & !results$is_target & !results$is_control & results$mode != "sham", ]
  corr <- list()
  if (nrow(primary)) for (d in split(primary, interaction(primary$cell_type, primary$condition, primary$mode, primary$factor, drop = TRUE))) {
    seeds <- unique(d$seed)
    if (length(seeds) < 2) next
    pairs <- utils::combn(seeds, 2)
    for (j in seq_len(ncol(pairs))) {
      a <- d[d$seed == pairs[1, j], ]; b <- d[d$seed == pairs[2, j], ]; k <- match(a$gene, b$gene)
      valid <- is.finite(a$rank_fraction) & is.finite(b$rank_fraction[k])
      rho <- if (sum(valid) >= 3) suppressWarnings(stats::cor(a$rank_fraction[valid], b$rank_fraction[k][valid], method = "spearman")) else NA_real_
      corr[[length(corr) + 1L]] <- data.frame(cell_type = d$cell_type[1], condition = d$condition[1], mode = d$mode[1], factor = d$factor[1], seed_a = pairs[1, j], seed_b = pairs[2, j], spearman_rank = rho)
    }
  }
  if (length(corr)) write_tsv(do.call(rbind, corr), file.path(out, "seed_rank_agreement.tsv"))
  invisible(summary)
}
run_enrichment <- function(results, cfg, base, out) {
  if (!isTRUE(cfg$enrichment$enabled) || is.null(results) || !nrow(results)) return(invisible())
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    write_tsv(data.frame(status = "missing_clusterProfiler"), file.path(out, "enrichment_status.tsv")); return(invisible())
  }
  # Prespecified first seed, all-donor run: prevents choosing the best-looking seed.
  d <- results[results$seed == cfg$network$seeds[1] & results$holdout == "none" & !results$is_control & !results$is_target & results$mode != "sham", ]
  groups <- split(d, interaction(d$cell_type, d$condition, d$mode, d$factor, drop = TRUE))
  statuses <- list()
  for (name in names(groups)) for (db in c("GO", "KEGG")) {
    z <- groups[[name]]; genes <- z$gene[z$p.adj < cfg$statistics$fdr]; universe <- unique(z$gene)
    folder <- file.path(out, "enrichment", safe_label(name), db); dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    tryCatch({
      gmt <- cfg$enrichment[[paste0(tolower(db), "_gmt")]]
      if (!is.null(gmt)) {
        path <- resolve_path(gmt, base)
        sets <- clusterProfiler::read.gmt(path)
        background <- intersect(universe, unique(sets$gene)); selected <- intersect(genes, background)
        write_json(list(source = path, sha256 = digest::digest(file = path, algo = "sha256"), ontology = db), file.path(folder, "source.json"))
        if (!length(selected)) {
          write_tsv(data.frame(gene = background), file.path(folder, "background.tsv"))
          write_tsv(data.frame(gene = selected), file.path(folder, "selected.tsv"))
          statuses[[length(statuses) + 1L]] <- data.frame(group = name, database = db, status = "no_mapped_significant_genes"); next
        }
        fit <- clusterProfiler::enricher(selected, universe = background, TERM2GENE = sets, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = cfg$enrichment$min_set_size, maxGSSize = cfg$enrichment$max_set_size)
      } else if (db == "GO") {
        pkg <- if (cfg$analysis$species == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
        if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing annotation package: ", pkg)
        org <- get(pkg, envir = asNamespace(pkg))
        mapping <- AnnotationDbi::select(org, keys = intersect(universe, AnnotationDbi::keys(org, keytype = "SYMBOL")), columns = c("ENTREZID", "GOALL"), keytype = "SYMBOL")
        mapping <- unique(mapping[!is.na(mapping$ENTREZID) & !is.na(mapping$GOALL), ])
        write_tsv(mapping, file.path(folder, "mapping.tsv"))
        background <- unique(mapping$ENTREZID); selected <- unique(mapping$ENTREZID[mapping$SYMBOL %in% genes])
        write_json(list(package = pkg, version = as.character(utils::packageVersion(pkg)), ontology = "GO_ALL", pool = TRUE), file.path(folder, "source.json"))
        if (!length(selected)) {
          write_tsv(data.frame(gene = background), file.path(folder, "background.tsv"))
          write_tsv(data.frame(gene = selected), file.path(folder, "selected.tsv"))
          statuses[[length(statuses) + 1L]] <- data.frame(group = name, database = db, status = "no_mapped_significant_genes"); next
        }
        fit <- clusterProfiler::enrichGO(selected, universe = background, OrgDb = org, keyType = "ENTREZID", ont = "ALL", pool = TRUE, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = FALSE, minGSSize = cfg$enrichment$min_set_size, maxGSSize = cfg$enrichment$max_set_size)
      } else {
        statuses[[length(statuses) + 1L]] <- data.frame(group = name, database = db, status = "offline_KEGG_GMT_required"); next
      }
      write_tsv(data.frame(gene = background), file.path(folder, "background.tsv"))
      write_tsv(data.frame(gene = selected), file.path(folder, "selected.tsv"))
      tab <- as.data.frame(fit)
      if (!ncol(tab)) tab <- data.frame(ID = character(), Description = character(), pvalue = numeric(), p.adjust = numeric())
      write_tsv(tab, file.path(folder, "all_terms.tsv"))
      write_tsv(tab[!is.na(tab$p.adjust) & tab$p.adjust < cfg$statistics$fdr, , drop = FALSE], file.path(folder, "significant_terms.tsv"))
      statuses[[length(statuses) + 1L]] <- data.frame(group = name, database = db, status = if (nrow(tab)) "ok" else "no_testable_terms")
    }, error = function(e) {
      statuses[[length(statuses) + 1L]] <<- data.frame(group = name, database = db, status = paste0("unavailable: ", conditionMessage(e)))
    })
  }
  if (length(statuses)) write_tsv(do.call(rbind, statuses), file.path(out, "enrichment_status.tsv"))
}
score_signatures <- function(pb, cfg, base, out) {
  if (!isTRUE(cfg$signatures$enabled) || is.null(pb)) return(invisible())
  path <- resolve_path(cfg$signatures$gene_sets, base)
  provenance <- resolve_path(cfg$signatures$provenance, base)
  if (is.null(path) || is.null(provenance) || !file.exists(provenance)) fail("SIGNATURE_PROVENANCE", "Provide gene-set TSV and provenance JSON with source and version.")
  prov <- jsonlite::fromJSON(provenance)
  if (is.null(prov$source) || is.null(prov$version)) fail("SIGNATURE_PROVENANCE", "source and version are required.")
  sets <- read_table(path)
  if (!all(c("set", "gene", "weight", "type") %in% names(sets))) fail("SIGNATURE_SCHEMA", "Require set, gene, weight, type.")
  if (anyDuplicated(paste(sets$set, sets$gene)) || any(!is.finite(sets$weight))) fail("SIGNATURE_SCHEMA", "Duplicate genes per set or invalid weights.")
  l <- log1p(t(t(as.matrix(pb$sums)) / pmax(Matrix::colSums(pb$sums), 1)) * 1e6)
  rows <- list()
  for (ct in unique(pb$metadata$cell_type)) {
    ix <- which(pb$metadata$cell_type == ct)
    z <- t(scale(t(l[, ix, drop = FALSE])))
    for (nm in unique(sets$set)) {
      ss <- sets[sets$set == nm, ]; ss <- ss[!ss$gene %in% c(cfg$analysis$target_gene, cfg$analysis$readout_genes), ]
      usable <- ss$gene[ss$gene %in% rownames(z)]
      usable <- usable[vapply(usable, function(g) all(is.finite(z[g, ])), logical(1))]
      md <- pb$metadata[ix, ]; md$set <- nm; md$type <- paste(unique(ss$type), collapse = ";")
      md$covered_genes <- length(usable); md$total_genes <- nrow(ss)
      md$status <- if (length(usable) < cfg$signatures$min_genes) "insufficient_coverage" else "observational_signature_proxy"
      w <- ss$weight[match(usable, ss$gene)]
      md$score <- if (length(usable) >= cfg$signatures$min_genes && sum(abs(w)) > 0) as.numeric(crossprod(w, z[usable, , drop = FALSE]) / sum(abs(w))) else NA_real_
      rows[[length(rows) + 1L]] <- md
    }
  }
  write_tsv(do.call(rbind, rows), file.path(out, "donor_signature_scores.tsv"))
  write_json(list(provenance = prov, sha256 = digest::digest(file = path, algo = "sha256"), method = "weighted mean of gene-wise z(log1p pseudobulk CPM), across donors/conditions within cell type; excludes perturbed/readout genes"), file.path(out, "signature_provenance.json"))
}
