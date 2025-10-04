# ================================
# 📦 Load Required Libraries
# ================================
library(phyloseq)
library(tidyverse)
library(vegan)
library(microbiome)
library(pheatmap)
library(ggpubr)
library(readr)
library(ggplot2)
library(viridis)
library(FSA)
library(qiime2R)
library(dplyr)
library(tidyr)
library(ggsignif)
library(patchwork)

# ================================
# 📥 Load and Preprocess Data
# ================================
# Load ASV Table
asv_table <- read_qza("table.qza")$data
asv_table <- otu_table(asv_table, taxa_are_rows = TRUE)

# Load Taxonomy Table
taxonomy <- read_qza("silva_taxonomy.qza")$data %>%
  as.data.frame() %>%
  column_to_rownames("Feature.ID")
taxonomy <- tax_table(as.matrix(taxonomy))

# Load Metadata
metadata <- read.delim("Metadata.tsv", sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(
    SampleName = trimws(SampleName),
    Condition = trimws(Condition),
    ExperimentLocation = trimws(ExperimentLocation),
    Plant = trimws(Plant)
  ) %>%
  column_to_rownames("SampleID")

metadata$SampleID <- tolower(rownames(metadata))
rownames(metadata) <- tolower(rownames(metadata))
colnames(asv_table) <- tolower(colnames(asv_table))

# Load Tree
tree <- read_qza("rooted-tree.qza")$data
tree <- phy_tree(tree)

# Create phyloseq object
asv_table <- asv_table[, rownames(metadata), drop = FALSE]
metadata <- metadata[colnames(asv_table), , drop = FALSE]
taxa_names(asv_table) <- rownames(asv_table)
taxa_names(taxonomy) <- rownames(asv_table)
physeq <- phyloseq(asv_table, taxonomy, sample_data(metadata), tree)

# ================================
# 🧬 Ordination + Merge Metadata
# ================================
unifrac_dist <- phyloseq::distance(physeq, method = "unifrac")
ordination <- ordinate(physeq, method = "PCoA", distance = unifrac_dist)

ordination_df <- data.frame(ordination$vectors[, 1:3])
ordination_df$SampleID <- tolower(rownames(ordination_df))
colnames(ordination_df)[1:3] <- c("Axis1", "Axis2", "Axis3")

metadata_df <- as(sample_data(physeq), "data.frame")
metadata_df$SampleID <- tolower(rownames(metadata_df))

df <- inner_join(ordination_df, metadata_df, by = "SampleID")

# ================================
# 📈 Compute % variance for axes (for labels)
# ================================
rel_eig <- ordination$values$Relative_eig
if (is.null(rel_eig) && "Eigenvalues" %in% names(ordination$values)) {
  rel_eig <- ordination$values$Eigenvalues / sum(ordination$values$Eigenvalues)
}
pct1 <- round(100 * rel_eig[1], 1)
pct2 <- round(100 * rel_eig[2], 1)
pct3 <- round(100 * rel_eig[3], 1)

lab_x1 <- paste0("Axis 1 (", pct1, "%)")
lab_y2 <- paste0("Axis 2 (", pct2, "%)")
lab_y3 <- paste0("Axis 3 (", pct3, "%)")

# ================================
# 🌿 Define Highlight Colors & Shapes
# ================================
highlight_samples <- c(
  "Arabidopsis thaliana", "Solanum lycopersicum", "Panicum virgatum",
  "Phaseolus vulgaris", "Oryza sativa", "Zea mays",
  "Sorghum bicolor", "Medicago sativa L."
)
highlight_colors <- c(
  "Arabidopsis thaliana" = "#D62728",
  "Solanum lycopersicum" = "#9467BD",
  "Panicum virgatum" = "#2CA02C",
  "Phaseolus vulgaris" = "#FF7F0E",
  "Oryza sativa" = "#1F77B4",
  "Zea mays" = "#8C564B",
  "Sorghum bicolor" = "#E377C2",
  "Medicago sativa L." = "#17BECF",
  "Other" = "grey80"
)
df$Highlight <- ifelse(df$SampleName %in% highlight_samples, df$SampleName, "Other")
df$Highlight <- factor(df$Highlight, levels = names(highlight_colors))
df$Plant <- factor(df$Plant, levels = c("Monocot", "Dicot"))

# ================================
# ✏️ Theme (bold axis titles)
# ================================
base_theme <- theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10)
  )

# ================================
# 📂 Output folder
# ================================
out_dir <- "beta-diversity"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ================================
# 📊 Create 6 PCoA Panels (axis labels without 'PCoA', with %)
# ================================
# a) Axis1 vs Axis2 - Color only
p1 <- ggplot(df, aes(Axis1, Axis2)) +
  geom_point(aes(color = Highlight), size = 3.5, alpha = 0.85) +
  scale_color_manual(values = highlight_colors) +
  labs(x = lab_x1, y = lab_y2) +
  ggtitle("a") +
  base_theme

# b) Axis1 vs Axis3 - Color only
p2 <- ggplot(df, aes(Axis1, Axis3)) +
  geom_point(aes(color = Highlight), size = 3.5, alpha = 0.85) +
  scale_color_manual(values = highlight_colors) +
  labs(x = lab_x1, y = lab_y3) +
  ggtitle("b") +
  base_theme +
  theme(legend.position = "none")

# c) Axis1 vs Axis2 - with Ellipses
p3 <- ggplot(df, aes(Axis1, Axis2)) +
  stat_ellipse(aes(group = Condition, color = Condition), size = 1, alpha = 0.5) +
  geom_point(aes(color = Highlight), size = 3.5, alpha = 0.85) +
  scale_color_manual(values = c("Control" = "#1f77b4", "Drought" = "#e6550d", highlight_colors)) +
  labs(x = lab_x1, y = lab_y2) +
  ggtitle("c") +
  base_theme

# d) Axis1 vs Axis3 - with Ellipses
p4 <- ggplot(df, aes(Axis1, Axis3)) +
  stat_ellipse(aes(group = Condition, color = Condition), size = 1, alpha = 0.5) +
  geom_point(aes(color = Highlight), size = 3.5, alpha = 0.85) +
  scale_color_manual(values = c("Control" = "#1f77b4", "Drought" = "#e6550d", highlight_colors)) +
  labs(x = lab_x1, y = lab_y3) +
  ggtitle("d") +
  base_theme +
  theme(legend.position = "none")

# e) Axis1 vs Axis2 - Shape + Ellipses
p5 <- ggplot(df, aes(Axis1, Axis2)) +
  stat_ellipse(aes(group = Condition, color = Condition), size = 1, alpha = 0.5) +
  geom_point(aes(color = Highlight, shape = Plant), size = 3.5, alpha = 0.85) +
  scale_shape_manual(values = c("Monocot" = 17, "Dicot" = 19)) +
  scale_color_manual(values = c("Control" = "#1f77b4", "Drought" = "#e6550d", highlight_colors)) +
  labs(x = lab_x1, y = lab_y2) +
  ggtitle("e") +
  base_theme

# f) Axis1 vs Axis3 - Shape + Ellipses
p6 <- ggplot(df, aes(Axis1, Axis3)) +
  stat_ellipse(aes(group = Condition, color = Condition), size = 1, alpha = 0.5) +
  geom_point(aes(color = Highlight, shape = Plant), size = 3.5, alpha = 0.85) +
  scale_shape_manual(values = c("Monocot" = 17, "Dicot" = 19)) +
  scale_color_manual(values = c("Control" = "#1f77b4", "Drought" = "#e6550d", highlight_colors)) +
  labs(x = lab_x1, y = lab_y3) +
  ggtitle("f") +
  base_theme +
  theme(legend.position = "none")

# ================================
# 🧩 Compose final and grouped figures
# ================================
final_plot <- (p1 + p2) / (p3 + p4) / (p5 + p6)
pair_ab <- p1 + p2
pair_cd <- p3 + p4
pair_ef <- p5 + p6

# ================================
# 💾 Save TIFFs (600 dpi, LZW) to beta-diversity/
# ================================
ggsave(file.path(out_dir, "figure4.tiff"), plot = final_plot,
       width = 174 / 25.4, height = 234 / 25.4, dpi = 600,
       units = "in", device = "tiff", compression = "lzw", bg = "white")

ggsave(file.path(out_dir, "panel_ab.tiff"), plot = pair_ab,
       width = 12, height = 5.5, dpi = 600,
       units = "in", device = "tiff", compression = "lzw", bg = "white")

ggsave(file.path(out_dir, "panel_cd.tiff"), plot = pair_cd,
       width = 12, height = 5.5, dpi = 600,
       units = "in", device = "tiff", compression = "lzw", bg = "white")

ggsave(file.path(out_dir, "panel_ef.tiff"), plot = pair_ef,
       width = 12, height = 5.5, dpi = 600,
       units = "in", device = "tiff", compression = "lzw", bg = "white")

# ================================
# 💾 Save PDFs (vector) to beta-diversity/
# ================================
ggsave(file.path(out_dir, "figure4.pdf"), plot = final_plot,
       width = 174 / 25.4, height = 234 / 25.4,
       units = "in", device = cairo_pdf, bg = "white")

ggsave(file.path(out_dir, "panel_ab.pdf"), plot = pair_ab,
       width = 12, height = 5.5,
       units = "in", device = cairo_pdf, bg = "white")

ggsave(file.path(out_dir, "panel_cd.pdf"), plot = pair_cd,
       width = 12, height = 5.5,
       units = "in", device = cairo_pdf, bg = "white")

ggsave(file.path(out_dir, "panel_ef.pdf"), plot = pair_ef,
       width = 12, height = 5.5,
       units = "in", device = cairo_pdf, bg = "white")

message("✅ Saved TIFFs (600 dpi) and PDFs to 'beta-diversity/' with axis labels like 'Axis 1 (xx.x%)'.")
