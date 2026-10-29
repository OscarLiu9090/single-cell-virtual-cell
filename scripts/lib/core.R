`%||%` <- function(x, y) if (is.null(x)) y else x
fail <- function(code, message) stop(structure(list(message = paste0(code, ": ", message), call = NULL, code = code), class = c("scvp_error", "error", "condition")))
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = TRUE, na = "NA", fileEncoding = "UTF-8")
}
write_json <- function(x, path) jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null")
resolve_path <- function(x, base) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  if (!grepl("^([A-Za-z]:|/|\\\\)", x)) x <- file.path(base, x)
  normalizePath(x, winslash = "/", mustWork = FALSE)
}
safe_label <- function(x) gsub("[^[:alnum:]_.-]", "_", x)
check_config <- function(cfg) {
  for (key in c("input", "analysis", "output")) if (is.null(cfg[[key]])) fail("CONFIG", paste("Missing", key))
  a <- cfg$analysis
  if (length(a$target_gene) != 1 || !nzchar(a$target_gene)) fail("CONFIG", "Exactly one target_gene is required.")
  if (!length(a$readout_genes) || !length(a$cell_types) || !length(a$conditions)) fail("CONFIG", "readout_genes, cell_types and conditions are required.")
  for (v in list(a$cell_types, a$conditions)) if (anyDuplicated(safe_label(v))) fail("CONFIG", "Selected labels collide after filesystem sanitization.")
  if (a$target_gene %in% a$readout_genes) fail("CONFIG", "Readout gene must differ from the perturbed gene.")
  if (!a$species %in% c("human", "mouse")) fail("CONFIG", "species must be human or mouse.")
  n <- cfg$network
  for (key in c("n_hvg", "n_net", "n_cells", "n_comp", "td_k", "td_max_iter", "ma_dim", "n_cores")) {
    v <- n[[key]]
    if (length(v) != 1 || !is.numeric(v) || !is.finite(v) || v < 1 || v != as.integer(v)) fail("CONFIG", paste("Invalid network", key))
  }
  if (n$n_comp < 2 || n$n_hvg <= n$n_comp || n$td_k > n$n_net || n$n_cells <= n$n_comp + 1) fail("CONFIG", "Incompatible network dimensions.")
  if (any(!is.finite(a$oe_factors)) || any(a$oe_factors < 1)) fail("CONFIG", "OE factors must be finite and >=1.")
  if (!length(n$seeds) || anyDuplicated(n$seeds) || any(!is.finite(n$seeds))) fail("CONFIG", "seeds must be unique finite integers.")
  if (any(n$seeds != as.integer(n$seeds)) || any(n$seeds < 0)) fail("CONFIG", "Invalid seed.")
  if (!is.numeric(cfg$qc$max_mt_fraction) || cfg$qc$max_mt_fraction < 0 || cfg$qc$max_mt_fraction > 1) fail("CONFIG", "Invalid mitochondrial fraction.")
  if (n$q < 0 || n$q > 1 || cfg$statistics$fdr <= 0 || cfg$statistics$fdr >= 1) fail("CONFIG", "Invalid q/FDR threshold.")
  if (cfg$expression_prediction$enabled %||% FALSE) fail("UNVALIDATED_PREDICTOR", "No validated expression-response adapter is bundled. Supply matched perturbation evidence before implementing one.")
  invisible(cfg)
}
validate_counts <- function(x) {
  if (!is.matrix(x) && !inherits(x, "Matrix")) fail("COUNTS_TYPE", "Expected a numeric genes-by-cells matrix.")
  if (any(dim(x) < 2)) fail("COUNTS_DIM", "At least two genes and cells are required.")
  if (is.null(rownames(x)) || is.null(colnames(x))) fail("COUNTS_NAMES", "Gene and cell identifiers are required.")
  if (anyNA(rownames(x)) || anyNA(colnames(x)) || any(!nzchar(rownames(x))) || any(!nzchar(colnames(x)))) fail("COUNTS_NAMES", "Empty identifiers.")
  if (anyDuplicated(colnames(x))) fail("DUPLICATE_CELL", "Duplicate cell identifiers; add sample prefixes at import.")
  v <- if (inherits(x, "sparseMatrix")) x@x else as.vector(x)
  if (!is.numeric(v) || any(!is.finite(v)) || any(v < 0) || any(abs(v - round(v)) > 1e-6)) fail("RAW_COUNTS_REQUIRED", "Counts must be finite, non-negative integers; normalized/integrated values are not accepted.")
  as(x, "dgCMatrix")
}
collapse_symbols <- function(x, symbols) {
  if (length(symbols) != nrow(x) || anyNA(symbols) || any(!nzchar(symbols))) fail("GENE_ID", "Every feature requires a nonempty symbol; supply a mapping table.")
  uniq <- unique(symbols)
  mapping <- data.frame(original_id = rownames(x), symbol = symbols, collapsed = duplicated(symbols) | duplicated(symbols, fromLast = TRUE))
  if (anyDuplicated(symbols)) {
    aggregation <- Matrix::sparseMatrix(i = match(symbols, uniq), j = seq_along(symbols), x = 1, dims = c(length(uniq), length(symbols)))
    x <- aggregation %*% x
    rownames(x) <- uniq
  } else rownames(x) <- symbols
  list(counts = as(x, "dgCMatrix"), mapping = mapping)
}
align_metadata <- function(meta, cells, mapping = list()) {
  for (dest in names(mapping)) {
    src <- mapping[[dest]]
    if (!src %in% names(meta)) fail("METADATA_COLUMN", paste("Missing", src))
    meta[[dest]] <- meta[[src]]
  }
  if (!"cell_id" %in% names(meta)) fail("METADATA_COLUMN", "cell_id is required.")
  meta$cell_id <- as.character(meta$cell_id)
  if (anyNA(meta$cell_id) || anyDuplicated(meta$cell_id)) fail("DUPLICATE_METADATA", "cell_id must be present and unique.")
  if (!all(cells %in% meta$cell_id)) fail("METADATA_MISMATCH", "Some count-matrix cells have no metadata.")
  extras <- sum(!meta$cell_id %in% cells)
  if (extras) warning(extras, " metadata rows have no counts and are excluded.")
  meta <- meta[match(cells, meta$cell_id), , drop = FALSE]
  for (key in c("sample_id", "donor_id", "condition", "cell_type")) {
    if (!key %in% names(meta)) meta[[key]] <- NA_character_
    meta[[key]] <- as.character(meta[[key]])
    meta[[key]][!nzchar(trimws(meta[[key]]))] <- NA_character_
  }
  rownames(meta) <- cells
  if (any(meta$donor_id == "none", na.rm = TRUE)) fail("DONOR_ID", "donor_id 'none' is reserved for all-donor runs.")
  donors <- unique(stats::na.omit(meta$donor_id))
  if (anyDuplicated(safe_label(donors))) fail("DONOR_ID", "Donor labels collide after filesystem sanitization.")
  meta
}
qc_counts <- function(x, meta, q, species) {
  lib <- Matrix::colSums(x)
  features <- Matrix::colSums(x > 0)
  mt <- grepl(if (species == "human") "^MT-" else "^mt-", rownames(x))
  if (!any(mt)) warning("No mitochondrial symbols detected: verify gene IDs and assay (scRNA/snRNA); fraction is reported as 0.")
  frac <- Matrix::colSums(x[mt, , drop = FALSE]) / pmax(lib, 1)
  keep <- lib >= q$min_counts & features >= q$min_genes & frac <= q$max_mt_fraction
  if (!is.null(q$max_counts)) keep <- keep & lib <= q$max_counts
  if (!is.null(q$max_genes)) keep <- keep & features <= q$max_genes
  if (!is.null(q$doublet_column)) {
    if (!q$doublet_column %in% names(meta)) fail("DOUBLET_COLUMN", "Configured doublet column is absent.")
    calls <- meta[[q$doublet_column]]
    keep <- keep & !is.na(calls) & calls %in% q$singlet_values
  }
  qc <- data.frame(cell_id = colnames(x), counts = lib, genes = features, mt_fraction = frac, retained = keep)
  if (!any(keep)) fail("EMPTY_AFTER_QC", "No cells pass QC.")
  list(counts = x[, keep, drop = FALSE], meta = meta[keep, , drop = FALSE], qc = qc)
}
gene_metrics <- function(x, genes) {
  do.call(rbind, lapply(genes, function(g) {
    if (!g %in% rownames(x)) return(data.frame(gene = g, mean_count = NA_real_, detection = 0, variance = NA_real_, status = "missing_gene"))
    v <- as.numeric(x[g, ])
    data.frame(gene = g, mean_count = mean(v), detection = mean(v > 0), variance = stats::var(v), status = if (!any(v > 0)) "not_detected" else if (stats::var(v) == 0) "no_variation" else "ok")
  }))
}
balanced_cells <- function(meta, seed, max_per_donor) {
  if (anyNA(meta$donor_id)) fail("DONOR_ID_REQUIRED", "Network sampling requires explicit donor IDs. Do not replace donor IDs by condition labels.")
  groups <- split(meta$cell_id, meta$donor_id)
  k <- min(min(lengths(groups)), max_per_donor)
  if (k < 1) fail("NO_CELLS", "No cells to sample.")
  set.seed(seed)
  # sample.int avoids R's length-one numeric sample shortcut.
  unlist(lapply(groups, function(ids) ids[sample.int(length(ids), k)]), use.names = FALSE)
}
select_genes <- function(x, required, n_hvg, min_cells, exclude_pattern) {
  ok <- Matrix::rowSums(x > 0) >= min_cells
  means <- Matrix::rowMeans(x)
  vars <- Matrix::rowMeans(x^2) - means^2
  ok <- ok & is.finite(vars) & vars > 1e-12
  if (!is.null(exclude_pattern) && nzchar(exclude_pattern)) ok <- ok & !grepl(exclude_pattern, rownames(x))
  eligible <- rownames(x)[ok]
  required_ok <- gene_metrics(x, required)
  required_ok <- required_ok$gene[required_ok$status == "ok"]
  candidates <- unique(c(eligible, required_ok))
  if (length(candidates) < 10) fail("INSUFFICIENT_GENES", "Fewer than 10 informative genes.")
  s <- Seurat::CreateSeuratObject(x[candidates, , drop = FALSE], min.cells = 0, min.features = 0)
  s <- Seurat::FindVariableFeatures(s, selection.method = "vst", nfeatures = min(n_hvg, length(candidates)), verbose = FALSE)
  unique(c(SeuratObject::VariableFeatures(s), required_ok))
}
match_controls <- function(x, target, excluded, n = 3) {
  mu <- Matrix::rowMeans(x); pct <- Matrix::rowMeans(x > 0)
  score <- abs(log1p(mu) - log1p(mu[target])) / max(stats::sd(log1p(mu)), 1e-8) + abs(pct - pct[target]) / max(stats::sd(pct), 1e-8)
  candidates <- setdiff(names(sort(score)), excluded)
  head(candidates, n)
}
function_fingerprint <- function(pkg, name) {
  f <- getFromNamespace(name, pkg)
  digest::digest(paste(c(deparse(formals(f), width.cutoff = 500), deparse(body(f), width.cutoff = 500)), collapse = "\n"), algo = "sha256", serialize = FALSE)
}
check_engine <- function(root) {
  lock <- jsonlite::fromJSON(file.path(root, "references", "engine-lock.json"), simplifyVector = FALSE)
  for (pkg in names(lock$packages)) {
    if (!requireNamespace(pkg, quietly = TRUE)) fail("MISSING_ENGINE", paste("Install pinned", pkg, "with scripts/install_engine.R"))
    entry <- lock$packages[[pkg]]
    if (as.character(utils::packageVersion(pkg)) != entry$version) fail("ENGINE_VERSION", paste(pkg, "must be", entry$version))
    for (nm in names(entry$functions))
      if (!identical(function_fingerprint(pkg, nm), entry$functions[[nm]])) fail("ENGINE_FINGERPRINT", paste(pkg, nm, "differs from verified source; revalidate before using."))
  }
  invisible(TRUE)
}
