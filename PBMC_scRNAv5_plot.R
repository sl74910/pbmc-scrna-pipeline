#!/usr/bin/env Rscript

# PBMC v5 summary figure
# Reads the saved Harmony checkpoint and saved ScType cluster labels.  No
# integration or cell-type annotation is rerun here.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(scales)
})
source(file.path(getwd(), "R", "pbmc_helpers.R"))

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
checkpoint_file <- file.path(
  project_dir, "outputs", "PBMC_scRNAv3_harmony", "checkpoints",
  "PBMC_v3_harmony_integrated.RData"
)
label_file <- file.path(
  project_dir, "outputs", "PBMC_scRNAv4_ScType",
  "PBMC_v4_ScType_cluster_labels.csv"
)
out_dir <- file.path(project_dir, "outputs", "PBMC_scRNAv5_plot")
figure_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")
pdf_file <- file.path(figure_dir, "PBMC_v5_summary_figure.pdf")
day_levels <- c("Day0", "Day1", "Day3", "Day5")
random_seed <- 1234L

pbmc_validate_parameters(list(
  checkpoint_file = checkpoint_file, label_file = label_file,
  out_dir = out_dir, day_levels = day_levels, random_seed = random_seed
))
pbmc_make_dirs(c(out_dir, figure_dir, table_dir))
if (!file.exists(checkpoint_file)) stop("Missing Harmony checkpoint: ", checkpoint_file)
if (!file.exists(label_file)) stop("Missing ScType label table: ", label_file)
set.seed(random_seed)

message("Loading saved Harmony object (this may take a minute)...")
pbmc <- pbmc_load_rdata_object(checkpoint_file, "pbmc")
if (!inherits(pbmc, "Seurat")) stop("Checkpoint object `pbmc` is not a Seurat object.")
pbmc_require_metadata(pbmc, c("sample_id", "day_group", "measurement_id"))
pbmc_require_reductions(pbmc, "umap.harmony")
if (!"RNA" %in% Assays(pbmc)) stop("The checkpoint must contain an RNA assay.")
DefaultAssay(pbmc) <- "RNA"

labels <- read.csv(label_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("cluster", "type") %in% colnames(labels))) {
  stop("ScType label table must contain `cluster` and `type` columns.")
}
labels$cluster <- as.character(labels$cluster)
labels$type <- trimws(as.character(labels$type))
labels <- labels[!duplicated(labels$cluster), , drop = FALSE]
clusters <- as.character(pbmc[["clusters.harmony"]][, 1L])
missing_clusters <- setdiff(unique(clusters), labels$cluster)
if (length(missing_clusters)) {
  stop("ScType labels are missing clusters: ", paste(missing_clusters, collapse = ", "))
}
celltype <- unname(setNames(labels$type, labels$cluster)[clusters])
type_levels <- c(
  "Naive B cells", "Memory CD4+ T cells", "Memory CD8+ T cells",
  "CD4+ NKT-like cells", "CD8+ NKT-like cells", "Natural killer  cells",
  "Non-classical monocytes", "Intermediate monocytes", "Neutrophils", "Granulocytes"
)
type_levels <- c(intersect(type_levels, unique(celltype)), setdiff(unique(celltype), type_levels))
short_names <- c(
  "Naive B cells" = "Naive B cells",
  "Memory CD4+ T cells" = "CD4 T cells",
  "Memory CD8+ T cells" = "CD8 T cells",
  "CD4+ NKT-like cells" = "CD4+ NKT",
  "CD8+ NKT-like cells" = "CD8+ NKT",
  "Natural killer  cells" = "NK cells",
  "Non-classical monocytes" = "Non-classical monocytes",
  "Intermediate monocytes" = "Intermediate monocytes",
  "Neutrophils" = "Neutrophils",
  "Granulocytes" = "Granulocytes"
)
short_names <- short_names[type_levels]
short_names[is.na(short_names)] <- type_levels[is.na(short_names)]
pbmc$sc_type_label <- factor(celltype, levels = type_levels)
pbmc$sc_type_short <- factor(unname(short_names[celltype]), levels = unname(short_names))
pbmc$sc_type_dot <- factor(unname(short_names[celltype]), levels = rev(unname(short_names)))
pbmc$day_group <- factor(as.character(pbmc$day_group), levels = day_levels)
measurement_numbers <- sort(unique(suppressWarnings(as.integer(sub("^R", "", as.character(pbmc$measurement_id))))))
measurement_levels <- paste0("R", measurement_numbers[is.finite(measurement_numbers)])
pbmc$measurement_id <- factor(as.character(pbmc$measurement_id), levels = measurement_levels)

# Keep a compact cell-level table for provenance and downstream reuse.
metadata <- pbmc[[]] %>%
  transmute(cell = rownames(pbmc[[]]), sample_id = as.character(sample_id),
            day_group = as.character(day_group), measurement_id = as.character(measurement_id),
            cluster = clusters, celltype = as.character(pbmc$sc_type_label),
            celltype_short = as.character(pbmc$sc_type_short))
write.csv(metadata, file.path(table_dir, "PBMC_v5_cell_metadata.csv"), row.names = FALSE)

# A restrained, publication-style palette shared by all panels.
type_palette <- c(
  "Naive B cells" = "#D89000",
  "CD4 T cells" = "#E07A5F",
  "CD8 T cells" = "#C06C9B",
  "CD4+ NKT" = "#F2B134",
  "CD8+ NKT" = "#56A7D8",
  "NK cells" = "#009E73",
  "Non-classical monocytes" = "#0072B2",
  "Intermediate monocytes" = "#4C78A8",
  "Neutrophils" = "#7A5195",
  "Granulocytes" = "#E5C93D"
)
type_palette <- type_palette[unname(short_names)]
names(type_palette) <- unname(short_names)

# Panel a: Harmony UMAP with saved ScType labels.
umap_plot <- DimPlot(
  pbmc, reduction = "umap.harmony", group.by = "sc_type_short",
  label = TRUE, repel = TRUE, label.size = 2.2, raster = TRUE,
  raster.dpi = c(240, 240), cols = type_palette
) +
  labs(title = "a", x = "UMAP 1", y = "UMAP 2", colour = NULL) +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 13),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 7), legend.key.height = unit(0.32, "cm"),
        plot.margin = margin(3, 3, 3, 3))

# Panel b: composition for every saved library, ordered by day and measurement.
composition <- metadata %>%
  mutate(
    day_group = factor(day_group, levels = day_levels),
    measurement_id = factor(measurement_id, levels = measurement_levels),
    celltype_short = factor(celltype_short, levels = unname(short_names))
  ) %>%
  count(day_group, measurement_id, celltype_short, name = "n", .drop = FALSE) %>%
  group_by(day_group, measurement_id) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()
write.csv(composition, file.path(table_dir, "PBMC_v5_celltype_composition.csv"), row.names = FALSE)

# Facets create four clean day blocks while retaining all five measurements in
# each block.  The same cell-type palette is decoded by the UMAP legend.
composition_plot <- ggplot(composition, aes(measurement_id, proportion, fill = celltype_short)) +
  geom_col(width = 0.88, colour = "white", linewidth = 0.08) +
  facet_grid(. ~ day_group, scales = "free_x", space = "free", switch = "x") +
  scale_y_continuous(labels = percent_format(accuracy = 25), breaks = c(0, .25, .5, .75, 1),
                     expand = c(0, 0)) +
  scale_fill_manual(values = type_palette, drop = FALSE) +
  labs(title = "b", x = NULL, y = "Proportion of cells", fill = NULL) +
  theme_classic(base_size = 8) +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 13),
        axis.text.x = element_text(size = 6), axis.title.y = element_text(size = 8),
        legend.position = "none", strip.background = element_blank(),
        strip.placement = "outside", strip.text = element_text(face = "bold", size = 8),
        panel.spacing.x = unit(0.55, "lines"), panel.border = element_rect(
          colour = "grey80", fill = NA, linewidth = 0.35),
        plot.margin = margin(3, 3, 3, 3))

# Panel d: average marker expression by annotated cell type.  Markers are
# intentionally broad and mouse-compatible; absent genes are silently omitted.
marker_sets <- list(
  "B" = c("Cd79a", "Ms4a1", "Cd74", "H2-Ab1", "Cd37"),
  "Granulocytes" = c("S100a8", "S100a9", "Il1b", "Cxcl2", "Retnlg"),
  "Monocytes" = c("Lyz2", "Lst1", "Ctss", "Fcgr3", "Ccr2"),
  "Neutrophils" = c("Ly6g", "Mpo", "Camp", "Ltf", "Ngp"),
  "NK/NKT" = c("Nkg7", "Klrd1", "Ncr1", "Prf1", "Gzmb"),
  "T cells" = c("Cd3d", "Cd3e", "Trbc1", "Cd4", "Cd8a")
)
marker_genes <- unique(intersect(unlist(marker_sets, use.names = FALSE), rownames(pbmc)))
if (length(marker_genes) < 5L) stop("Too few marker genes match the object.")
if (length(grep("^counts\\.", Layers(pbmc[["RNA"]]))) > 1L) {
  pbmc <- JoinLayers(pbmc, assay = "RNA")
}
if (!"data" %in% Layers(pbmc[["RNA"]])) {
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000,
                        verbose = FALSE)
}
marker_to_group <- setNames(rep(names(marker_sets), lengths(marker_sets)), unlist(marker_sets))
marker_group <- unname(marker_to_group[marker_genes])
dot_features <- split(marker_genes, marker_group)
dot_features <- dot_features[vapply(dot_features, length, integer(1)) > 0L]
dot_plot <- DotPlot(
  pbmc, features = dot_features, group.by = "sc_type_dot",
  cols = c("#F3F3F3", "#D55E00"), dot.scale = 4.8
) +
  RotatedAxis() +
  labs(title = "c", x = NULL, y = NULL, colour = "Average expression", size = "% expressed") +
  theme_classic(base_size = 8) +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 13),
        axis.text.x = element_text(angle = 50, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        strip.text.x = element_text(size = 6, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
        legend.position = "right", plot.margin = margin(3, 3, 3, 3))

# Assemble the large figure: the evidence-rich dot plot occupies the full
# right column, matching the visual hierarchy of a main-text figure.
figure <- patchwork::wrap_plots(
  A = umap_plot, B = dot_plot, C = composition_plot,
  design = "AB\nCB"
) +
  plot_layout(widths = c(1.06, 1.24), heights = c(1.03, 0.97)) +
  plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))

pbmc_save_pdf(figure, pdf_file, width = 16, height = 13)
ggsave(file.path(figure_dir, "PBMC_v5_summary_figure.png"), figure,
       width = 16, height = 13, dpi = 220, bg = "white")
capture.output(sessionInfo(), file = file.path(table_dir, "sessionInfo.txt"))
message("PBMC v5 summary figure written to: ", pdf_file)
