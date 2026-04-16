#!/usr/bin/env Rscript

project_root <- "/Users/hai/CPSC-445-Project"
local_rlib <- file.path(project_root, ".Rlibs")
if (dir.exists(local_rlib)) {
  .libPaths(c(local_rlib, .libPaths()))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(topicmodels)
  library(slam)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

log_step <- function(...) {
  msg <- paste(...)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
}

parse_args <- function(x) {
  if (length(x) == 0) {
    return(list())
  }
  parsed <- strsplit(sub("^--", "", x), "=", fixed = TRUE)
  keys <- vapply(parsed, `[`, character(1), 1)
  vals <- vapply(parsed, function(y) if (length(y) > 1) y[2] else "", character(1))
  as.list(stats::setNames(vals, keys))
}

clean_colnames <- function(df) {
  nm <- names(df)
  nm <- gsub("^\\ufeff", "", nm, useBytes = TRUE)
  nm <- gsub("^X\\.\\.\\.", "", nm)
  if (length(nm) > 0 && grepl("cluster$", nm[1])) {
    nm[1] <- "cluster"
  }
  if (length(nm) > 0 && grepl("Gene$", nm[1])) {
    nm[1] <- "Gene"
  }
  names(df) <- nm
  df
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

if (!length(args) || "--help" %in% commandArgs(trailingOnly = TRUE)) {
  cat(
    "Usage:\n",
    "  Rscript lda_escc_programs.R \\\n",
    "    --input=/path/to/scESCC_map.rds \\\n",
    "    --signatures=/path/to/BS_BK_DK_sig.csv \\\n",
    "    --output=/path/to/output_dir \\\n",
    "    --k_grid=2,3,4,5,6,7,8\n\n",
    "Default scope:\n",
    "  Fits LDA to squamous epithelium cells only and compares LDA mixture\n",
    "  entropy to the original BS/BK/DK-based CCI score.\n",
    sep = ""
  )
  quit(save = "no")
}

input_path <- if (!is.null(args$input)) args$input else stop("--input is required")
signature_path <- if (!is.null(args$signatures)) args$signatures else stop("--signatures is required")
output_dir <- if (!is.null(args$output)) args$output else stop("--output is required")
k_grid <- if (!is.null(args$k_grid) && nzchar(args$k_grid)) as.integer(strsplit(args$k_grid, ",", fixed = TRUE)[[1]]) else 2:8
max_cells <- if (!is.null(args$max_cells) && nzchar(args$max_cells)) as.integer(args$max_cells) else NA_integer_
n_hvg <- if (!is.null(args$n_hvg) && nzchar(args$n_hvg)) as.integer(args$n_hvg) else 2000L

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_path)
log_step("Loaded Seurat object from", input_path)

if (!all(c("new_anno2", "new_group2") %in% colnames(obj@meta.data))) {
  stop("Expected metadata columns new_anno2 and new_group2 are missing from the Seurat object")
}

subset_cells <- rownames(subset(obj@meta.data, new_anno2 == "Squamous Epithelium"))
if (length(subset_cells) == 0) {
  stop("No squamous epithelium cells found in new_anno2")
}
if (!is.na(max_cells) && length(subset_cells) > max_cells) {
  set.seed(1)
  subset_cells <- sample(subset_cells, max_cells)
}
log_step("Selected", length(subset_cells), "squamous epithelium cells for LDA input")

counts_all <- GetAssayData(obj, layer = "counts")[, subset_cells, drop = FALSE]
meta_sub <- obj@meta.data[subset_cells, , drop = FALSE]
obj_sub <- CreateSeuratObject(counts = counts_all, meta.data = meta_sub)
obj_sub <- NormalizeData(obj_sub, verbose = FALSE)
log_step("Created lightweight Seurat object and normalized data")

# Use highly variable genes to keep the LDA model interpretable and computationally manageable.
obj_sub <- FindVariableFeatures(obj_sub, selection.method = "vst", nfeatures = n_hvg, verbose = FALSE)
hvg <- VariableFeatures(obj_sub)
log_step("Selected", length(hvg), "highly variable genes")

counts <- GetAssayData(obj_sub, layer = "counts")[hvg, , drop = FALSE]
keep_cells <- Matrix::colSums(counts) >= 100
keep_genes <- Matrix::rowSums(counts > 0) >= 20
counts <- counts[keep_genes, keep_cells, drop = FALSE]
obj_sub <- CreateSeuratObject(counts = counts, meta.data = meta_sub[colnames(counts), , drop = FALSE])
obj_sub <- NormalizeData(obj_sub, verbose = FALSE)
log_step("Filtered matrix to", nrow(counts), "genes and", ncol(counts), "cells")

if (ncol(counts) < 100) {
  stop("Too few cells remain after filtering for stable LDA fitting")
}

set.seed(1)
cell_ids <- colnames(counts)
train_ids <- sample(cell_ids, size = floor(0.8 * length(cell_ids)))
test_ids <- setdiff(cell_ids, train_ids)

train_dtm <- slam::as.simple_triplet_matrix(Matrix::t(counts[, train_ids, drop = FALSE]))
test_dtm <- slam::as.simple_triplet_matrix(Matrix::t(counts[, test_ids, drop = FALSE]))
full_dtm <- slam::as.simple_triplet_matrix(Matrix::t(counts))
log_step("Prepared train/test/full document-term matrices")

fit_one_k <- function(k) {
  log_step("Starting LDA model fit for K =", k)
  model <- topicmodels::LDA(
    x = train_dtm,
    k = k,
    method = "VEM",
    control = list(seed = 1, estimate.alpha = TRUE, verbose = 0)
  )
  data.frame(
    k = k,
    log_likelihood = as.numeric(logLik(model)),
    perplexity = perplexity(model, newdata = test_dtm)
  )
}

model_selection <- do.call(rbind, lapply(k_grid, fit_one_k))
write.csv(model_selection, file.path(output_dir, "tables", "lda_model_selection.csv"), row.names = FALSE)
log_step("Finished model selection grid:", paste(k_grid, collapse = ", "))

best_k <- model_selection$k[which.min(model_selection$perplexity)]
log_step("Selected best K =", best_k, "based on minimum held-out perplexity")

final_model <- topicmodels::LDA(
  x = full_dtm,
  k = best_k,
  method = "VEM",
  control = list(seed = 1, estimate.alpha = TRUE, verbose = 0)
)
log_step("Finished final LDA fit on full matrix")

post <- posterior(final_model)
theta <- as.data.frame(post$topics)
colnames(theta) <- paste0("topic_", seq_len(ncol(theta)))
theta$cell_id <- rownames(post$topics)

beta <- as.data.frame(post$terms)
beta$topic <- paste0("topic_", seq_len(nrow(beta)))

topic_top_genes <- do.call(rbind, lapply(seq_len(nrow(post$terms)), function(i) {
  ord <- order(post$terms[i, ], decreasing = TRUE)[1:20]
  data.frame(
    topic = paste0("topic_", i),
    rank = seq_along(ord),
    gene = colnames(post$terms)[ord],
    beta = post$terms[i, ord]
  )
}))

write.csv(theta, file.path(output_dir, "tables", "lda_topic_weights_per_cell.csv"), row.names = FALSE)
write.csv(topic_top_genes, file.path(output_dir, "tables", "lda_top_genes_per_topic.csv"), row.names = FALSE)
saveRDS(final_model, file.path(output_dir, "lda_final_model.rds"))
log_step("Saved final model and topic tables")

# LDA-derived "identity mixture" summary.
theta_mat <- as.matrix(theta[, grepl("^topic_", colnames(theta)), drop = FALSE])
entropy <- -rowSums(theta_mat * log(theta_mat + 1e-12)) / log(ncol(theta_mat))
dominance <- apply(theta_mat, 1, max)
dominant_topic <- colnames(theta_mat)[max.col(theta_mat, ties.method = "first")]

lda_metrics <- data.frame(
  cell_id = theta$cell_id,
  lda_entropy = entropy,
  lda_mixture_score = 1 - dominance,
  dominant_topic = dominant_topic
)

# Recompute the original CCI score on the same cells for direct comparison.
sig <- clean_colnames(read.csv(signature_path, check.names = FALSE))

sig_dk <- subset(sig, cluster == "Differentiated_keratinocyte")$gene
sig_bk <- subset(sig, cluster == "Basel_keratinocyte")$gene
sig_bs <- subset(sig, cluster == "Basel_Stem_Cells")$gene

expr_data <- GetAssayData(obj, layer = "data")[, colnames(counts), drop = FALSE]
cci <- data.frame(
  cell_id = colnames(expr_data),
  DK = Matrix::colMeans(expr_data[intersect(sig_dk, rownames(expr_data)), , drop = FALSE]),
  BK = Matrix::colMeans(expr_data[intersect(sig_bk, rownames(expr_data)), , drop = FALSE]),
  BS = Matrix::colMeans(expr_data[intersect(sig_bs, rownames(expr_data)), , drop = FALSE])
)

a <- cci$BK
b <- cci$DK
c_ <- cci$BS
den1 <- 2 * a * sqrt(a^2 + b^2 + c_^2)
den2 <- 2 * b * sqrt(a^2 + b^2 + c_^2)
den3 <- 2 * c_ * sqrt(a^2 + b^2 + c_^2)

cci$cos1 <- ifelse(den1 == 0, 2, (a^2 + (a^2 + b^2 + c_^2) - (b^2 + c_^2)) / den1)
cci$cos2 <- ifelse(den2 == 0, 2, (b^2 + (a^2 + b^2 + c_^2) - (a^2 + c_^2)) / den2)
cci$cos3 <- ifelse(den3 == 0, 2, (c_^2 + (a^2 + b^2 + c_^2) - (a^2 + b^2)) / den3)
cci$CCI_score <- 10 / apply(cci[, c("cos1", "cos2", "cos3")], 1, sd)

meta <- meta_sub[colnames(counts), , drop = FALSE] %>%
  tibble::rownames_to_column("cell_id") %>%
  select(cell_id, new_group2, new_group, new_anno2)

comparison_df <- meta %>%
  inner_join(lda_metrics, by = "cell_id") %>%
  inner_join(cci, by = "cell_id")

write.csv(comparison_df, file.path(output_dir, "tables", "lda_vs_original_cci.csv"), row.names = FALSE)

cor_test <- suppressWarnings(cor.test(comparison_df$lda_entropy, comparison_df$CCI_score, method = "spearman"))
cor_summary <- data.frame(
  best_k = best_k,
  spearman_rho = unname(cor_test$estimate),
  p_value = cor_test$p.value,
  n_cells = nrow(comparison_df)
)
write.csv(cor_summary, file.path(output_dir, "tables", "lda_vs_cci_correlation.csv"), row.names = FALSE)
log_step("Computed LDA-vs-CCI comparison metrics")

# Plots
p_k <- ggplot(model_selection, aes(k, perplexity)) +
  geom_line() +
  geom_point() +
  theme_classic(base_size = 12) +
  labs(title = "LDA model selection", x = "Number of latent programs (K)", y = "Held-out perplexity")
ggsave(file.path(output_dir, "01_lda_model_selection.pdf"), p_k, width = 5, height = 4)

p_entropy <- ggplot(comparison_df, aes(new_group2, lda_entropy, fill = new_group2)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = quantile(comparison_df$lda_entropy, c(0.01, 0.99), na.rm = TRUE)) +
  theme_classic(base_size = 12) +
  labs(title = "LDA mixture entropy by group", x = NULL, y = "Normalized topic entropy")
ggsave(file.path(output_dir, "02_lda_entropy_by_group.pdf"), p_entropy, width = 5, height = 4)

p_compare <- ggplot(comparison_df, aes(CCI_score, lda_entropy, color = new_group2)) +
  geom_point(alpha = 0.35, size = 0.5) +
  theme_classic(base_size = 12) +
  labs(
    title = paste0("Original CCI vs LDA entropy (Spearman rho = ", round(cor_summary$spearman_rho, 3), ")"),
    x = "Original CCI score",
    y = "LDA normalized entropy"
  )
ggsave(file.path(output_dir, "03_lda_entropy_vs_original_cci.pdf"), p_compare, width = 5.5, height = 4.5)

if ("tsne" %in% names(obj@reductions)) {
  emb <- as.data.frame(Embeddings(obj, "tsne")[colnames(counts), 1:2, drop = FALSE])
  emb$cell_id <- rownames(emb)
  plot_df <- inner_join(emb, theta, by = "cell_id")
  topic_cols <- colnames(theta)[grepl("^topic_", colnames(theta))]
  n_plot_topics <- min(length(topic_cols), 6L)
  plot_df_long <- tidyr::pivot_longer(plot_df, cols = all_of(topic_cols[1:n_plot_topics]), names_to = "topic", values_to = "weight")
  p_topics <- ggplot(plot_df_long, aes(tSNE_1, tSNE_2, color = weight)) +
    geom_point(size = 0.2) +
    facet_wrap(~topic, ncol = 3) +
    scale_color_gradient(low = "#d9d9d9", high = "#b30000") +
    theme_classic(base_size = 10) +
    labs(title = "Top LDA topic weights on tSNE", color = "Weight")
  ggsave(file.path(output_dir, "04_lda_topics_on_tsne.pdf"), p_topics, width = 9, height = 6)
}

log_step("LDA ESCC program analysis complete:", output_dir)
