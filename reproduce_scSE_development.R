#!/usr/bin/env Rscript

project_root <- "/Users/hai/CPSC-445-Project"
local_rlib <- file.path(project_root, ".Rlibs")
local_cache <- file.path(project_root, ".cache")
dir.create(local_cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = local_cache)
if (dir.exists(local_rlib)) {
  .libPaths(c(local_rlib, .libPaths()))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(Matrix)
  library(pheatmap)
  library(RColorBrewer)
  library(GSVA)
  library(msigdbr)
})

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(x) {
  if (length(x) == 0) {
    return(list())
  }
  parsed <- strsplit(sub("^--", "", x), "=", fixed = TRUE)
  keys <- vapply(parsed, `[`, character(1), 1)
  vals <- vapply(parsed, function(y) if (length(y) > 1) y[2] else "", character(1))
  as.list(stats::setNames(vals, keys))
}
args <- parse_args(args)

if (!length(args) || "--help" %in% commandArgs(trailingOnly = TRUE)) {
  cat(
    "Usage:\n",
    "  Rscript reproduce_scSE_development.R \\\n",
    "    --normal=/path/to/normalSE_development_map.rds \\\n",
    "    --integrated=/path/to/scHCA_SE_development.rds \\\n",
    "    --monocle=/path/to/scHCA_SE_development_monocle.rds \\\n",
    "    --output=/path/to/output_dir\n\n",
    "The script recreates the scSE_development.md workflow as closely as possible.\n",
    "Because monocle3 is unavailable in this environment, the trajectory panel is\n",
    "reproduced as a PHATE/pseudotime proxy using the integrated Seurat object.\n",
    sep = ""
  )
  quit(save = "no")
}

required_args <- c("normal", "integrated", "monocle", "output")
missing_args <- required_args[!required_args %in% names(args) | vapply(required_args, function(k) args[[k]] == "", logical(1))]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

output_dir <- args$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "intermediates"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "summary_tables"), recursive = TRUE, showWarnings = FALSE)

normal_obj <- readRDS(args$normal)
integrated_obj <- readRDS(args$integrated)
monocle_obj <- readRDS(args$monocle)

cell_types <- c("Basel_keratinocyte", "Basel_Stem_Cells", "Differentiated_keratinocyte")
cell_type_colors <- c(
  "Basel_keratinocyte" = "#e76f51",
  "Basel_Stem_Cells" = "#2a9d8f",
  "Differentiated_keratinocyte" = "#264653"
)

feature_palette <- c("#007BBF", "#FFF485", "#FF0000")

pseudo_bulk_mean <- function(seurat_obj, group_labels) {
  expr <- GetAssayData(seurat_obj, layer = "data")
  split_cells <- split(colnames(seurat_obj), group_labels)
  split_cells <- split_cells[vapply(split_cells, length, integer(1)) > 0]
  out <- lapply(split_cells, function(cells) {
    Matrix::rowMeans(expr[, cells, drop = FALSE])
  })
  out <- do.call(cbind, out)
  colnames(out) <- names(split_cells)
  out
}

pseudo_bulk_split_mean <- function(seurat_obj, cells, nsplit, prefix) {
  cells <- sample(cells)
  bucket <- cut(seq_along(cells), breaks = min(nsplit, length(cells)), labels = FALSE)
  expr <- GetAssayData(seurat_obj, layer = "data")[, cells, drop = FALSE]
  out <- lapply(sort(unique(bucket)), function(i) {
    Matrix::rowMeans(expr[, bucket == i, drop = FALSE])
  })
  out <- do.call(cbind, out)
  colnames(out) <- paste0(prefix, "_", seq_len(ncol(out)))
  out
}

save_feature_grid <- function(obj, features, reduction, filename, ncol = 2) {
  features <- features[features %in% rownames(obj) | features %in% colnames(obj@meta.data)]
  p <- FeaturePlot(
    object = obj,
    features = features,
    reduction = reduction,
    cols = feature_palette,
    ncol = ncol,
    pt.size = 0.2,
    combine = TRUE
  ) & NoLegend()
  ggsave(filename, plot = p, width = 7.5, height = 5)
}

# 1. Normal sample PHATE by annotated subtype
normal_obj$new_anno3 <- factor(as.character(normal_obj$new_anno3), levels = cell_types)
p1 <- DimPlot(
  object = normal_obj,
  reduction = "phate",
  group.by = "new_anno3",
  cols = cell_type_colors,
  pt.size = 0.3
) +
  labs(title = "Normal esophagus epithelium") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(output_dir, "01_normal_epi_phate_by_cell_type.pdf"), p1, width = 8.5, height = 5)

# 2. Marker feature plots on PHATE
save_feature_grid(
  obj = normal_obj,
  features = c("TP63", "SOX2", "KRT13", "IVL"),
  reduction = "phate",
  filename = file.path(output_dir, "02_normal_epi_featureplot_TP63_SOX2_KRT13_IVL.pdf"),
  ncol = 2
)

# 3. Marker identification
Idents(normal_obj) <- normal_obj$new_anno3
marker_rds_path <- file.path(output_dir, "intermediates", "Heso6_normal_only_epi_hetero.marker.rds")
marker_csv_path <- file.path(output_dir, "intermediates", "Heso6_normal_only_epi_hetero.marker.csv")
if (file.exists(marker_rds_path)) {
  markers <- readRDS(marker_rds_path)
} else {
  markers <- FindAllMarkers(normal_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers <- subset(markers, p_val_adj < 0.05)
  write.csv(markers, marker_csv_path, row.names = FALSE)
  saveRDS(markers, marker_rds_path)
}

markers_dk <- subset(markers, cluster == "Differentiated_keratinocyte" & pct.2 < 0.4)
markers_bk <- subset(markers, cluster == "Basel_keratinocyte" & pct.2 < 0.4)
markers_bs <- subset(markers, cluster == "Basel_Stem_Cells" & pct.2 < 0.4)

# 4. HCA/ours integrated cell-cycle scatter and sample similarity
integrated_obj$predicted.id <- factor(as.character(integrated_obj$predicted.id), levels = c("Basel_Stem_Cells", "Basel_keratinocyte", "Differentiated_keratinocyte"))

cell_cycle_df <- FetchData(integrated_obj, vars = c("G2M.Score", "S.Score", "predicted.id", "Phase"))
p2 <- ggplot(cell_cycle_df, aes(x = G2M.Score, y = S.Score, color = Phase)) +
  geom_point(alpha = 0.12, size = 0.4) +
  facet_wrap(~predicted.id) +
  theme_classic(base_size = 12) +
  labs(title = "Cell-cycle state by normal SE subtype")
ggsave(file.path(output_dir, "03_integrated_cell_cycle_scatter.pdf"), p2, width = 15, height = 5)

integrated_obj$new_group <- as.character(integrated_obj$donor_time)
integrated_obj$new_group[is.na(integrated_obj$new_group)] <- "Heso_Ours"
group_means <- pseudo_bulk_mean(integrated_obj, integrated_obj$new_group)
group_cor <- cor(group_means)
write.csv(group_cor, file.path(output_dir, "summary_tables", "sample_similarity_correlation_matrix.csv"))
pheatmap(
  group_cor,
  border_color = NA,
  color = colorRampPalette(c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026"))(100),
  filename = file.path(output_dir, "04_sample_similarity_heatmap.pdf"),
  width = 10,
  height = 8
)

# 5. Integrated UMAP by predicted subtype
p3 <- DimPlot(
  object = integrated_obj,
  reduction = "umap",
  group.by = "predicted.id",
  cols = cell_type_colors,
  pt.size = 0.15
) +
  labs(title = "Integrated normal SE cell types") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(output_dir, "05_integrated_umap_by_predicted_id.pdf"), p3, width = 8.5, height = 5)

# 6. Marker dot plot on integrated object
integrated_obj$predicted.id <- factor(as.character(integrated_obj$predicted.id), levels = c("Differentiated_keratinocyte", "Basel_keratinocyte", "Basel_Stem_Cells"))
Idents(integrated_obj) <- integrated_obj$predicted.id
marker_features <- c("SOX2", "TP63", "KRT15", "KRT5", "KRT4", "KRT13", "IVL")
marker_features <- marker_features[marker_features %in% rownames(integrated_obj)]
p4 <- DotPlot(
  integrated_obj,
  features = marker_features,
  cols = c("#ffffff", "#B30000"),
  scale = TRUE,
  col.min = 0,
  col.max = 5
) + RotatedAxis() + theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(output_dir, "06_integrated_marker_dotplot.pdf"), p4, width = 8, height = 5)

# 7. DK/BK/BS signature scores on UMAP
score_signature <- function(obj, genes, field_name) {
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) {
    obj[[field_name]] <- 0
    return(obj)
  }
  expr <- GetAssayData(obj, layer = "data")[genes, , drop = FALSE]
  obj[[field_name]] <- Matrix::colMeans(expr)
  obj
}
integrated_obj <- score_signature(integrated_obj, markers_dk$gene, "Differentiated_keratinocyte")
integrated_obj <- score_signature(integrated_obj, markers_bk$gene, "Basel_keratinocyte")
integrated_obj <- score_signature(integrated_obj, markers_bs$gene, "Basel_Stem_Cells")
save_feature_grid(
  obj = integrated_obj,
  features = c("Differentiated_keratinocyte", "Basel_keratinocyte", "Basel_Stem_Cells"),
  reduction = "umap",
  filename = file.path(output_dir, "07_integrated_DK_BK_BS_signature_featureplots.pdf"),
  ncol = 3
)

# 8. Trajectory proxy panel
phate_embed <- Embeddings(integrated_obj, "phate")[, 1:2]
trajectory_df <- data.frame(
  PHATE_1 = phate_embed[, 1],
  PHATE_2 = phate_embed[, 2],
  predicted.id = integrated_obj$predicted.id,
  monocle3_pseudotime = integrated_obj$monocle3_pseudotime
)
trajectory_df <- trajectory_df[complete.cases(trajectory_df), ]

quantile_path <- trajectory_df %>%
  mutate(bin = dplyr::ntile(monocle3_pseudotime, 30)) %>%
  group_by(bin) %>%
  summarize(PHATE_1 = median(PHATE_1), PHATE_2 = median(PHATE_2), .groups = "drop")

p5 <- ggplot(trajectory_df, aes(PHATE_1, PHATE_2, color = predicted.id)) +
  geom_point(size = 0.2, alpha = 0.25) +
  geom_path(data = quantile_path, aes(PHATE_1, PHATE_2), inherit.aes = FALSE, color = "black", linewidth = 0.7) +
  scale_color_manual(values = cell_type_colors) +
  theme_classic(base_size = 12) +
  labs(title = "Trajectory proxy on PHATE", subtitle = "Approximated from stored pseudotime in the integrated Seurat object")
ggsave(file.path(output_dir, "08_trajectory_proxy_on_phate.pdf"), p5, width = 8.5, height = 5)

# 9. GSVA pathway heatmap
Idents(integrated_obj) <- integrated_obj$predicted.id
bulk_dk <- pseudo_bulk_split_mean(integrated_obj, WhichCells(integrated_obj, idents = "Differentiated_keratinocyte"), 30, "DK")
bulk_bk <- pseudo_bulk_split_mean(integrated_obj, WhichCells(integrated_obj, idents = "Basel_keratinocyte"), 30, "BK")
bulk_bs <- pseudo_bulk_split_mean(integrated_obj, WhichCells(integrated_obj, idents = "Basel_Stem_Cells"), 30, "BS")
subset_matrix <- cbind(bulk_dk, bulk_bk, bulk_bs)

target_sets <- c(
  "GO_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT",
  "GO_EPITHELIAL_CELL_FATE_COMMITMENT",
  "KEGG_RIBOSOME",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "KEGG_SPLICEOSOME",
  "KEGG_DNA_REPLICATION",
  "GO_CONDENSIN_COMPLEX",
  "GO_MAINTENANCE_OF_DNA_METHYLATION",
  "GO_EPITHELIAL_CELL_CELL_ADHESION",
  "GO_CORNIFICATION",
  "GO_KERATINIZATION"
)

msig_df <- msigdbr(species = "Homo sapiens")
gs_list <- split(msig_df$gene_symbol, msig_df$gs_name)
gs_list <- gs_list[target_sets[target_sets %in% names(gs_list)]]
gs_list <- lapply(gs_list, intersect, rownames(subset_matrix))
gs_list <- gs_list[vapply(gs_list, length, integer(1)) >= 5]

gsva_param <- gsvaParam(
  exprData = as.matrix(subset_matrix),
  geneSets = gs_list,
  minSize = 5,
  maxSize = 500,
  maxDiff = TRUE,
  verbose = FALSE
)
gsva_scores <- gsva(gsva_param, verbose = FALSE)
gsva_scores <- gsva_scores[target_sets[target_sets %in% rownames(gsva_scores)], , drop = FALSE]

gsva_meta <- data.frame(
  group = c(rep("DK", ncol(bulk_dk)), rep("BK", ncol(bulk_bk)), rep("BS", ncol(bulk_bs))),
  row.names = colnames(gsva_scores)
)
all_gsva_seurat <- CreateSeuratObject(counts = gsva_scores, assay = "RNA", project = "GSVA", min.cells = 0, meta.data = gsva_meta)
saveRDS(all_gsva_seurat, file.path(output_dir, "intermediates", "all_gsva_seurat_GSVA_seurat.rds"))

ordered_cols <- c(paste0("BS_", seq_len(ncol(bulk_bs))), paste0("BK_", seq_len(ncol(bulk_bk))), paste0("DK_", seq_len(ncol(bulk_dk))))
ordered_cols <- ordered_cols[ordered_cols %in% colnames(gsva_scores)]
gsva_plot_mat <- gsva_scores[, ordered_cols, drop = FALSE]
gsva_plot_mat <- t(scale(t(gsva_plot_mat)))
gsva_plot_mat[is.na(gsva_plot_mat)] <- 0

annotation_col <- data.frame(group = factor(c(rep("BS", ncol(bulk_bs)), rep("BK", ncol(bulk_bk)), rep("DK", ncol(bulk_dk))), levels = c("BS", "BK", "DK")))
rownames(annotation_col) <- ordered_cols

pheatmap(
  gsva_plot_mat,
  annotation_col = annotation_col,
  color = colorRampPalette(brewer.pal(10, "RdBu"))(101),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  border_color = NA,
  filename = file.path(output_dir, "09_gsva_heatmap.pdf"),
  width = 11,
  height = 5
)

write.csv(gsva_scores, file.path(output_dir, "summary_tables", "gsva_scores.csv"))

# 10. Cell-cycle pie charts
phase_table <- as.data.frame(table(integrated_obj$Phase, integrated_obj$predicted.id), stringsAsFactors = FALSE)
colnames(phase_table) <- c("Phase", "predicted.id", "Freq")
phase_table <- phase_table %>%
  group_by(predicted.id) %>%
  mutate(normal_ratio = 100 * Freq / sum(Freq)) %>%
  ungroup()
write.csv(phase_table, file.path(output_dir, "summary_tables", "cell_cycle_phase_distribution.csv"), row.names = FALSE)

make_pie <- function(df, subtype) {
  ggplot(subset(df, predicted.id == subtype), aes(x = "", y = normal_ratio, fill = Phase)) +
    geom_col(width = 1) +
    coord_polar("y", start = 0) +
    theme_void(base_size = 12) +
    labs(title = subtype)
}
p6 <- cowplot::plot_grid(
  make_pie(phase_table, "Basel_Stem_Cells"),
  make_pie(phase_table, "Basel_keratinocyte"),
  make_pie(phase_table, "Differentiated_keratinocyte"),
  nrow = 1
)
ggsave(file.path(output_dir, "10_cell_cycle_pies.pdf"), p6, width = 8, height = 5)

# 11. Summary of inputs
summary_df <- data.frame(
  object = c("normalSE_development_map", "scHCA_SE_development", "scHCA_SE_development_monocle"),
  class = c(class(normal_obj)[1], class(integrated_obj)[1], class(monocle_obj)[1]),
  cells = c(ncol(normal_obj), ncol(integrated_obj), NA),
  genes = c(nrow(normal_obj), nrow(integrated_obj), NA)
)
write.csv(summary_df, file.path(output_dir, "summary_tables", "object_summary.csv"), row.names = FALSE)

message("scSE_development reproduction complete: ", output_dir)
