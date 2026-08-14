#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
})

parse_args <- function(args) {
  out <- list(
    project_dir = normalizePath(Sys.getenv("PROJECT", "."), mustWork = FALSE),
    bridge_manifest = "results/task_c_cellchat_input_bridge/cellchat_input_bridge_manifest.tsv",
    sample_id = "GSM9060732_AdjIII-0019",
    output_prefix = "cellchat_full_gsm9060732",
    input_dir = "",
    meta_path = "",
    out_dir = "results/task_c_cellchat_full_sample_baseline",
    min_cells = 10,
    trim = 0.1,
    seed = 20260712
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop(paste("Invalid argument sequence near", key))
    }
    name <- sub("^--", "", key)
    if (!name %in% names(out)) {
      stop(paste("Unknown argument", key))
    }
    out[[name]] <- args[[i + 1]]
    i <- i + 2
  }
  out$min_cells <- as.integer(out$min_cells)
  out$trim <- as.numeric(out$trim)
  out$seed <- as.integer(out$seed)
  out
}

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

write_json_simple <- function(path, values) {
  scalar <- function(value) {
    if (is.logical(value)) return(ifelse(value, "true", "false"))
    if (is.numeric(value)) return(as.character(value))
    paste0('"', gsub('"', '\\"', as.character(value), fixed = TRUE), '"')
  }
  keys <- names(values)
  lines <- c("{")
  for (idx in seq_along(keys)) {
    suffix <- if (idx < length(keys)) "," else ""
    lines <- c(lines, paste0('  "', keys[[idx]], '": ', scalar(values[[idx]]), suffix))
  }
  writeLines(c(lines, "}"), path)
}

standardize_cellchat <- function(df, dataset, sample_id) {
  if (nrow(df) == 0) {
    return(data.frame(
      method = character(), dataset = character(), sample_id = character(),
      sender = character(), receiver = character(), ligand = character(),
      receptor = character(), pathway = character(), probability = numeric(),
      p_value = numeric(), interaction_name = character(), notes = character(),
      stringsAsFactors = FALSE
    ))
  }
  pick <- function(names_vec) {
    hit <- names_vec[names_vec %in% colnames(df)]
    if (length(hit) == 0) rep(NA_character_, nrow(df)) else as.character(df[[hit[[1]]]])
  }
  data.frame(
    method = "cellchat_full_sample",
    dataset = dataset,
    sample_id = sample_id,
    sender = pick(c("source", "sources")),
    receiver = pick(c("target", "targets")),
    ligand = pick(c("ligand", "ligands")),
    receptor = pick(c("receptor", "receptors")),
    pathway = pick(c("pathway_name", "pathway")),
    probability = suppressWarnings(as.numeric(pick(c("prob")))),
    p_value = suppressWarnings(as.numeric(pick(c("pval", "p.value")))),
    interaction_name = pick(c("interaction_name", "interaction_name_2")),
    notes = "CellChat full-expression single-sample resource calibration; computational compartments",
    stringsAsFactors = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
set.seed(args$seed)
project_dir <- normalizePath(args$project_dir, mustWork = TRUE)
manifest <- read_tsv(file.path(project_dir, args$bridge_manifest))
manifest <- manifest[manifest$sample_id == args$sample_id, , drop = FALSE]
if (nrow(manifest) != 1) {
  stop(paste("Expected one bridge row for", args$sample_id, "found", nrow(manifest)))
}
record <- manifest[1, ]
input_dir <- if (nzchar(args$input_dir)) {
  normalizePath(args$input_dir, mustWork = TRUE)
} else {
  file.path(project_dir, "data/processed/minimal_inputs", record$dataset, record$sample_id)
}
matrix_path <- file.path(input_dir, "matrix.mtx.gz")
features_path <- file.path(input_dir, "features.tsv.gz")
barcodes_path <- file.path(input_dir, "barcodes.tsv.gz")
meta_path <- if (nzchar(args$meta_path)) args$meta_path else file.path(project_dir, record$meta_path)
required <- c(matrix_path, features_path, barcodes_path, meta_path)
if (!all(file.exists(required))) stop(paste("Missing input:", paste(required[!file.exists(required)], collapse = ", ")))

features <- read.delim(gzfile(features_path), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
barcodes <- read.delim(gzfile(barcodes_path), header = FALSE, sep = "\t", stringsAsFactors = FALSE)[[1]]
counts <- as(Matrix::readMM(gzfile(matrix_path)), "CsparseMatrix")
if (nrow(counts) != nrow(features) || ncol(counts) != length(barcodes)) {
  stop("Matrix dimensions do not match features and barcodes")
}
gene_type <- if (ncol(features) >= 3) features[[3]] else rep("Gene Expression", nrow(features))
keep_features <- gene_type == "Gene Expression" & nzchar(features[[2]])
counts <- counts[keep_features, , drop = FALSE]
gene_symbols <- features[[2]][keep_features]
duplicate_symbols <- sum(duplicated(gene_symbols))
if (duplicate_symbols > 0) {
  unique_symbols <- unique(gene_symbols)
  aggregation <- Matrix::sparseMatrix(
    i = match(gene_symbols, unique_symbols),
    j = seq_along(gene_symbols),
    x = 1,
    dims = c(length(unique_symbols), length(gene_symbols))
  )
  counts <- aggregation %*% counts
  gene_symbols <- unique_symbols
}
rownames(counts) <- gene_symbols
colnames(counts) <- barcodes
nonzero_genes <- Matrix::rowSums(counts) > 0
counts <- counts[nonzero_genes, , drop = FALSE]

meta <- read_tsv(meta_path)
if (!"labels" %in% colnames(meta)) {
  if (!"spot_compartment" %in% colnames(meta)) stop("Metadata lacks labels and spot_compartment columns")
  meta$labels <- meta$spot_compartment
}
meta$labels <- as.character(meta$labels)
label_counts <- table(meta$labels)
rare_labels <- names(label_counts[label_counts < args$min_cells])
if (length(rare_labels) > 0) {
  collapse_target <- if ("ambiguous_or_low_signal" %in% names(label_counts)) {
    "ambiguous_or_low_signal"
  } else {
    names(label_counts)[which.max(label_counts)]
  }
  meta$labels[meta$labels %in% rare_labels] <- collapse_target
}
column_order <- match(meta$barcode, colnames(counts))
if (any(is.na(column_order))) stop("Metadata contains barcodes absent from the count matrix")
counts <- counts[, column_order, drop = FALSE]
rownames(meta) <- meta$barcode
meta$labels <- factor(meta$labels)
meta$samples <- factor(rep(as.character(record$sample_id), nrow(meta)))
if (nlevels(meta$labels) < 2) stop("CellChat requires at least two metadata groups")

library_size <- Matrix::colSums(counts)
scale_factors <- 10000 / pmax(library_size, 1)
data_input <- counts %*% Matrix::Diagonal(x = scale_factors)
data_input@x <- log1p(data_input@x)
rownames(data_input) <- rownames(counts)
colnames(data_input) <- colnames(counts)

cellchat <- createCellChat(object = data_input, meta = meta, group.by = "labels")
data("CellChatDB.human", package = "CellChat")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = args$trim)
cellchat <- filterCommunication(cellchat, min.cells = args$min_cells)
raw_df <- subsetCommunication(cellchat)
std_df <- standardize_cellchat(raw_df, record$dataset, record$sample_id)

out_dir <- if (grepl("^/", args$out_dir)) args$out_dir else file.path(project_dir, args$out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
raw_path <- file.path(out_dir, paste0(args$output_prefix, "_raw_subsetCommunication.tsv"))
std_path <- file.path(out_dir, paste0(args$output_prefix, "_candidates.tsv"))
summary_path <- file.path(out_dir, paste0(args$output_prefix, "_summary.json"))
write.table(raw_df, raw_path, sep = "\t", quote = FALSE, row.names = FALSE, eol = "\n")
write.table(std_df, std_path, sep = "\t", quote = FALSE, row.names = FALSE, eol = "\n")
write_json_simple(summary_path, list(
  run_id = Sys.getenv("RUN_ID", "manual"), dataset = record$dataset, sample_id = record$sample_id,
  spots = ncol(data_input), input_gene_symbols = length(gene_symbols), nonzero_genes = nrow(data_input),
  duplicate_symbols_aggregated = duplicate_symbols, groups = nlevels(meta$labels),
  rare_labels_collapsed = paste(rare_labels, collapse = ";"),
  raw_rows = nrow(raw_df), standardized_rows = nrow(std_df), database_mode = "CellChatDB.human_all_categories",
  min_cells = args$min_cells, trim = args$trim, seed = args$seed,
  output_prefix = args$output_prefix, cellchat_version = as.character(utils::packageVersion("CellChat")),
  nmf_version = as.character(utils::packageVersion("NMF")), status = "completed"
))
cat("CellChat full-sample calibration completed\n")
cat("raw_rows=", nrow(raw_df), "\n", sep = "")
cat("standardized_rows=", nrow(std_df), "\n", sep = "")
