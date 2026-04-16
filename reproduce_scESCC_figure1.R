#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
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
    "  Rscript reproduce_scESCC_figure1.R \\\n",
    "    --input=/path/to/scESCC_map.rds \\\n",
    "    --output=/path/to/output_dir\n\n",
    "Outputs:\n",
    "  01_tsne_by_cell_type.pdf\n",
    "  02_tsne_by_patient.pdf\n",
    "  03_celltype_composition_by_patient.pdf\n",
    "  04_marker_dotplot.pdf\n",
    "  05_nFeature_RNA_SE_boxplot.pdf\n",
    "  summary_tables/*.csv\n",
    sep = ""
  )
  quit(save = "no")
}

input_path <- args$input
output_dir <- args$output

if (is.null(input_path) || input_path == "") {
  stop("--input is required")
}
if (is.null(output_dir) || output_dir == "") {
  stop("--output is required")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "summary_tables"), recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_path)

required_meta <- c("new_anno2", "new_group", "new_group2", "nFeature_RNA")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
}
if (!"tsne" %in% names(obj@reductions)) {
  stop("This object does not contain a tsne reduction")
}

cell_type_levels <- c(
  "Squamous Epithelium",
  "Glandular epithelium",
  "Fibroblast",
  "T cells",
  "Vascular endothelial",
  "Macrophages",
  "Neutrophil",
  "Plasma cells",
  "B cells",
  "Smooth muscle",
  "Lymphatic Vascular endothelial",
  "Mast cells",
  "Plasmacytoid dendritic cells"
)

cell_type_colors <- c(
  "#f000e6", "#21A0A0", "#E0607E", "#00db00",
  "#2f3eff", "#ff0098", "#3A6EA5", "#942093", "#00b0f0", "#45503B",
  "#00eef3", "#f7f85e", "#36C9C6"
)
names(cell_type_colors) <- cell_type_levels

patient_levels <- c("Patient 2", "Patient 3", "Patient 1", "Patient 4", "Patient 5", "Patient 6")
patient_colors <- cell_type_colors[seq_along(patient_levels)]
names(patient_colors) <- patient_levels

obj$new_anno2 <- factor(as.character(obj$new_anno2), levels = cell_type_levels)
obj$new_group <- factor(as.character(obj$new_group), levels = patient_levels)
obj$new_group2 <- factor(as.character(obj$new_group2), levels = c("Normal", "ESCC"))

theme_paper <- function() {
  theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank()
    )
}

# 1. tSNE by cell type
p_celltype <- DimPlot(
  object = obj,
  reduction = "tsne",
  group.by = "new_anno2",
  cols = cell_type_colors,
  pt.size = 0.35
) +
  labs(title = "Cell type") +
  theme_paper()

ggsave(
  filename = file.path(output_dir, "01_tsne_by_cell_type.pdf"),
  plot = p_celltype,
  width = 8,
  height = 6
)

# 2. tSNE by patient
p_patient <- DimPlot(
  object = obj,
  reduction = "tsne",
  group.by = "new_group",
  cols = patient_colors,
  pt.size = 0.35
) +
  labs(title = "Patient") +
  theme_paper()

ggsave(
  filename = file.path(output_dir, "02_tsne_by_patient.pdf"),
  plot = p_patient,
  width = 8,
  height = 6
)

# 3. Cell-type composition by patient
composition <- as.data.frame(table(obj$new_group, obj$new_anno2), stringsAsFactors = FALSE)
colnames(composition) <- c("patient", "cell_type", "count")
composition <- composition %>%
  filter(!is.na(patient), !is.na(cell_type)) %>%
  group_by(patient) %>%
  mutate(percent = 100 * count / sum(count)) %>%
  ungroup()

write.csv(
  composition,
  file = file.path(output_dir, "summary_tables", "celltype_composition_by_patient.csv"),
  row.names = FALSE
)

p_comp <- ggplot(composition, aes(x = patient, y = percent, fill = cell_type)) +
  geom_col() +
  coord_flip() +
  geom_hline(yintercept = 50, linetype = 2, color = "grey45") +
  scale_fill_manual(values = cell_type_colors, drop = FALSE) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Percent", title = "Cell-type composition by patient") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(
  filename = file.path(output_dir, "03_celltype_composition_by_patient.pdf"),
  plot = p_comp,
  width = 8,
  height = 6
)

# 4. Marker dot plot
marker_features <- c("CD79A", "THY1", "KRT7", "PDPN", "CD14", "KIT", "S100A8", "IGJ", "IL3RA", "ACTA2", "KRT5", "CD3E", "ICAM1")
marker_features <- marker_features[marker_features %in% rownames(obj)]

Idents(obj) <- obj$new_anno2
p_dot <- DotPlot(
  object = obj,
  features = marker_features,
  cols = c("#ffffff", "#B30000"),
  scale = TRUE,
  col.min = 0,
  col.max = 5
) +
  RotatedAxis() +
  labs(title = "Marker genes across cell types") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(
  filename = file.path(output_dir, "04_marker_dotplot.pdf"),
  plot = p_dot,
  width = 10,
  height = 5
)

# 5. nFeature_RNA in squamous epithelium
se_cells <- obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  filter(new_anno2 == "Squamous Epithelium", !is.na(new_group2))

write.csv(
  se_cells,
  file = file.path(output_dir, "summary_tables", "squamous_epithelium_metadata.csv"),
  row.names = FALSE
)

p_box <- ggboxplot(
  data = se_cells,
  x = "new_group2",
  y = "nFeature_RNA",
  fill = "new_group2",
  legend = "none",
  outlier.shape = NA,
  notch = TRUE
) +
  stat_compare_means(
    comparisons = list(c("Normal", "ESCC")),
    method = "t.test",
    label = "p.signif"
  ) +
  labs(title = "Number of genes in squamous epithelium", x = NULL, y = "nFeature_RNA") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(
  filename = file.path(output_dir, "05_nFeature_RNA_SE_boxplot.pdf"),
  plot = p_box,
  width = 5,
  height = 5
)

write.csv(
  data.frame(
    metric = c("cells", "genes"),
    value = c(ncol(obj), nrow(obj))
  ),
  file = file.path(output_dir, "summary_tables", "object_dimensions.csv"),
  row.names = FALSE
)

message("Figure 1-style reproduction complete: ", output_dir)
