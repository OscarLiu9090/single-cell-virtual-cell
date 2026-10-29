read_table <- function(path) utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8", quote = "\"", comment.char = "")
read_10x_counts <- function(path, gene_column = 2) {
  find <- function(names) {
    paths <- file.path(path, names)
    found <- paths[file.exists(paths)]
    if (!length(found)) fail("10X_FILE", paste("Missing", paste(names, collapse = " or ")))
    found[1]
  }
  mf <- find(c("matrix.mtx.gz", "matrix.mtx"))
  ff <- find(c("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"))
  bf <- find(c("barcodes.tsv.gz", "barcodes.tsv"))
  con <- if (grepl("[.]gz$", mf)) gzfile(mf, "rt") else file(mf, "rt")
  on.exit(close(con), add = TRUE)
  m <- Matrix::readMM(con)
  features <- utils::read.delim(ff, header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  bars <- utils::read.delim(bf, header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (nrow(m) != nrow(features) || ncol(m) != nrow(bars) || gene_column > ncol(features)) fail("10X_DIM", "Feature/barcode dimensions do not match matrix.")
  rownames(m) <- features[[gene_column]]; colnames(m) <- bars[[1]]
  if (ncol(features) >= 3) m <- m[features[[3]] == "Gene Expression", , drop = FALSE]
  as(m, "dgCMatrix")
}
seurat_counts <- function(object, assay = "RNA") {
  if (!inherits(object, "Seurat")) fail("SEURAT_OBJECT", "Selected object is not a Seurat object.")
  if (!assay %in% names(object@assays)) fail("RNA_ASSAY", paste("Missing assay", assay))
  if (assay != "RNA") fail("RNA_ASSAY", "Use the original RNA count assay, not SCT/integrated residuals.")
  layers <- SeuratObject::Layers(object[[assay]], search = "^counts($|[.])")
  if (!length(layers)) fail("RAW_COUNTS_REQUIRED", "No counts layer found.")
  if (length(layers) > 1) {
    layer_cells <- unlist(lapply(layers, function(layer) colnames(SeuratObject::LayerData(object[[assay]], layer = layer))), use.names = FALSE)
    if (anyDuplicated(layer_cells)) fail("OVERLAPPING_COUNT_LAYERS", "Multiple count layers share cells; select the intended raw layer explicitly before import.")
    object <- SeuratObject::JoinLayers(object, assay = assay, layers = "counts", new = "counts")
  }
  counts <- SeuratObject::LayerData(object, assay = assay, layer = if (length(layers) == 1) layers else "counts")
  meta <- object[[]]
  meta$cell_id <- rownames(meta)
  list(counts = counts, meta = meta)
}
load_input <- function(cfg, base, out, root, runtime) {
  inp <- cfg$input
  type <- inp$type
  meta <- NULL
  if (type %in% c("seurat_rds", "seurat_rdata")) {
    path <- resolve_path(inp$path, base)
    if (type == "seurat_rds") obj <- readRDS(path) else {
      env <- new.env(parent = emptyenv())
      loaded <- load(path, envir = env)
      if (is.null(inp$object_name)) {
        found <- loaded[vapply(loaded, function(n) inherits(env[[n]], "Seurat"), logical(1))]
        if (length(found) != 1) fail("AMBIGUOUS_RDATA", "Specify input.object_name for the intended Seurat object.")
        obj <- env[[found]]
      } else {
        if (!inp$object_name %in% loaded) fail("RDATA_OBJECT", "Configured object was not found.")
        obj <- env[[inp$object_name]]
      }
    }
    data <- seurat_counts(obj, inp$assay %||% "RNA")
    x <- data$counts; meta <- data$meta
    rm(obj, data); invisible(gc())
  } else if (type == "10x") {
    samples <- read_table(resolve_path(inp$samples_file, base))
    if (!all(c("sample_id", "path") %in% names(samples)) || anyDuplicated(samples$sample_id)) fail("SAMPLE_MANIFEST", "Unique sample_id and path columns are required.")
    mats <- metas <- list()
    sample_base <- dirname(resolve_path(inp$samples_file, base))
    for (i in seq_len(nrow(samples))) {
      s <- samples[i, ]; path <- resolve_path(s$path, sample_base)
      m <- read_10x_counts(path, gene_column = inp$gene_column %||% 2)
      # Collapse genuine duplicate symbols, never silently make.unique().
      collapsed <- collapse_symbols(m, rownames(m)); m <- collapsed$counts
      write_tsv(collapsed$mapping, file.path(out, paste0("gene_mapping_", i, ".tsv")))
      colnames(m) <- paste0(s$sample_id, "_", colnames(m))
      md <- data.frame(cell_id = colnames(m), sample_id = s$sample_id)
      for (key in intersect(c("donor_id", "condition", "cell_type"), names(s))) md[[key]] <- as.character(s[[key]])
      mats[[i]] <- m; metas[[i]] <- md
    }
    # Feature union preserves counts and permits sample-specific feature lists.
    genes <- unique(unlist(lapply(mats, rownames), use.names = FALSE))
    mats <- lapply(mats, function(m) {
      map <- Matrix::sparseMatrix(i = match(rownames(m), genes), j = seq_len(nrow(m)), x = 1, dims = c(length(genes), nrow(m)))
      o <- map %*% m; rownames(o) <- genes; o
    })
    x <- do.call(cbind, mats)
    allcols <- unique(unlist(lapply(metas, names)))
    metas <- lapply(metas, function(m) { for (nm in setdiff(allcols, names(m))) m[[nm]] <- NA_character_; m[, allcols, drop = FALSE] })
    meta <- do.call(rbind, metas)
  } else if (type == "h5ad") {
    python <- cfg$environment$python %||% runtime$python
    if (is.null(python) || !file.exists(python)) fail("PYTHON", "Configure an interpreter with h5py, scipy and pandas.")
    bridge <- file.path(out, "h5ad_import"); dir.create(bridge, recursive = TRUE)
    args <- c(shQuote(file.path(root, "scripts", "h5ad_to_mtx.py")), "--input", shQuote(resolve_path(inp$path, base)), "--output", shQuote(bridge), "--layer", shQuote(inp$counts_layer %||% "counts"))
    if (!is.null(inp$gene_symbol_column)) args <- c(args, "--gene-symbol-column", shQuote(inp$gene_symbol_column))
    status <- system2(python, args, stdout = file.path(out, "h5ad_import.log"), stderr = file.path(out, "h5ad_import.stderr.log"))
    if (status != 0) fail("H5AD_IMPORT", "Conversion failed; inspect h5ad_import logs.")
    x <- Matrix::readMM(file.path(bridge, "matrix.mtx"))
    rownames(x) <- readLines(file.path(bridge, "genes.txt"), encoding = "UTF-8")
    colnames(x) <- readLines(file.path(bridge, "cells.txt"), encoding = "UTF-8")
    meta <- read_table(file.path(bridge, "metadata.tsv"))
  } else if (type == "matrix_rds") {
    x <- readRDS(resolve_path(inp$path, base))
    meta <- data.frame(cell_id = colnames(x))
  } else fail("INPUT_TYPE", "Supported: seurat_rds, seurat_rdata, 10x, h5ad, matrix_rds.")
  x <- validate_counts(x)
  symbols <- rownames(x)
  if (!is.null(inp$gene_mapping)) {
    gm <- read_table(resolve_path(inp$gene_mapping, base))
    if (!all(c("gene_id", "symbol") %in% names(gm)) || anyDuplicated(gm$gene_id)) fail("GENE_MAPPING", "Require unique gene_id and symbol columns.")
    symbols <- gm$symbol[match(rownames(x), gm$gene_id)]
  } else if (any(grepl("^ENS(MUS)?G[0-9]", symbols))) fail("GENE_ID", "Ensembl IDs found: supply an explicit gene_id-to-symbol mapping.")
  collapsed <- collapse_symbols(x, symbols); x <- collapsed$counts
  write_tsv(collapsed$mapping, file.path(out, "gene_mapping.tsv"))
  if (!is.null(inp$metadata)) {
    ext <- read_table(resolve_path(inp$metadata, base))
    maps <- inp$metadata_columns %||% list()
    for (dest in names(maps)) {
      if (!maps[[dest]] %in% names(ext)) fail("METADATA_COLUMN", paste("Missing", maps[[dest]]))
      ext[[dest]] <- ext[[maps[[dest]]]]
    }
    if (!"cell_id" %in% names(ext) || anyDuplicated(ext$cell_id)) fail("DUPLICATE_METADATA", "External metadata needs unique cell_id.")
    if (!all(colnames(x) %in% ext$cell_id)) fail("METADATA_MISMATCH", "External metadata must cover all input cells.")
    ext <- ext[match(meta$cell_id, ext$cell_id), , drop = FALSE]
    for (nm in setdiff(names(ext), "cell_id")) meta[[nm]] <- ext[[nm]]
  } else for (dest in names(inp$metadata_columns %||% list())) {
    src <- inp$metadata_columns[[dest]]
    if (!src %in% names(meta)) fail("METADATA_COLUMN", paste("Missing", src))
    meta[[dest]] <- meta[[src]]
  }
  meta <- align_metadata(meta, colnames(x))
  list(counts = x, meta = meta)
}
