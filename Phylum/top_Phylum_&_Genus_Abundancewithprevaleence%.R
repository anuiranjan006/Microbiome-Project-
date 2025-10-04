# ================================
# 📦 Load Required Libraries
# ================================
library(phyloseq)
library(tidyverse)
library(vegan)
library(microbiome)
library(ggplot2)
library(ggrepel)
library(qiime2R)
library(rstatix)
library(magick)   # for rasterized 600-dpi PDF export
library(stringr)  # for str_remove()

# ================================
# 📥 Load Input Data
# ================================
asv_table <- read_qza("table.qza")$data %>% otu_table(taxa_are_rows = TRUE)

taxonomy <- read_qza("silva_taxonomy.qza")$data %>%
  as.data.frame() %>%
  column_to_rownames("Feature.ID") %>%
  mutate(Taxon = as.character(Taxon)) %>%
  separate(
    Taxon,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = ";", fill = "right", extra = "merge"
  ) %>%
  mutate(across(everything(), ~str_remove(., "^[a-zA-Z]__")))
taxonomy <- tax_table(as.matrix(taxonomy))

metadata <- read.delim("Metadata.tsv", sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(SampleID = tolower(trimws(SampleID)))
rownames(metadata) <- metadata$SampleID

colnames(asv_table) <- tolower(trimws(colnames(asv_table)))
common_samples <- intersect(colnames(asv_table), rownames(metadata))
asv_table <- asv_table[, common_samples]
metadata  <- metadata[common_samples, ]

tree  <- read_qza("rooted-tree.qza")$data %>% phy_tree()
physeq <- phyloseq(asv_table, taxonomy, sample_data(metadata), tree)

# ================================
# 🚫 Remove Archaea globally (before normalization)
# ================================
physeq_no_archaea <- tryCatch(
  subset_taxa(physeq, is.na(Kingdom) | Kingdom != "Archaea"),
  error = function(e) physeq
)

# ================================
# 📁 Output Directory Setup
# ================================
output_dir <- "Figure8"
plot_dir   <- file.path(output_dir, "Plots")
table_dir  <- file.path(output_dir, "Tables")
stat_dir   <- file.path(output_dir, "SampleStats")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, showWarnings = FALSE)
dir.create(stat_dir,  showWarnings = FALSE)

# ================================
# 🔄 Normalize Data (after Archaea removal)
# ================================
physeq_rel <- transform_sample_counts(physeq_no_archaea, function(x) x / sum(x))

remove_unclassified <- function(df, rank_col) {
  df %>%
    filter(!grepl("uncultured|unclassified|metagenome|_sp|_bacterium|[0-9]{2,}",
                  !!sym(rank_col), ignore.case = TRUE)) %>%
    filter(!is.na(!!sym(rank_col)), !!sym(rank_col) != "")
}

# ================================
# 🔹 Top 10 Phyla
# ================================
phy_phylum <- tax_glom(physeq_rel, taxrank = "Phylum")
df_phy <- psmelt(phy_phylum) %>%
  left_join(metadata, by = c("Sample" = "SampleID")) %>%
  mutate(
    SampleName = coalesce(SampleName.x, SampleName.y),
    Condition  = coalesce(Condition.x, Condition.y)
  ) %>%
  select(SampleName, Condition, Phylum, Abundance)
df_phy <- remove_unclassified(df_phy, "Phylum")

top_phyla <- df_phy %>%
  group_by(Phylum) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  slice_max(mean_abund, n = 10) %>%
  pull(Phylum)

df_phy_stats <- df_phy %>%
  filter(Phylum %in% top_phyla) %>%
  group_by(Condition, Phylum) %>%
  summarise(
    PrevalencePercent = 100 * mean(Abundance > 0),
    MeanAbundance     = mean(Abundance) * 100,
    SDAbundance       = sd(Abundance) * 100,
    .groups = "drop"
  )

df_sample_phy <- df_phy %>%
  filter(Phylum %in% top_phyla) %>%
  group_by(SampleName, Condition, Phylum) %>%
  summarise(Abundance = sum(Abundance) * 100, .groups = "drop") %>%
  left_join(df_phy_stats, by = c("Condition", "Phylum"))

# ================================
# 🔹 Top 10 Genera
# ================================
phy_genus <- tax_glom(physeq_rel, taxrank = "Genus")
df_genus <- psmelt(phy_genus) %>%
  left_join(metadata, by = c("Sample" = "SampleID")) %>%
  mutate(
    SampleName = coalesce(SampleName.x, SampleName.y),
    Condition  = coalesce(Condition.x, Condition.y)
  ) %>%
  select(SampleName, Condition, Genus, Abundance)
df_genus <- remove_unclassified(df_genus, "Genus")

top_genera <- df_genus %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  slice_max(mean_abund, n = 10) %>%
  pull(Genus)

df_genus_stats <- df_genus %>%
  filter(Genus %in% top_genera) %>%
  group_by(Condition, Genus) %>%
  summarise(
    PrevalencePercent = 100 * mean(Abundance > 0),
    MeanAbundance     = mean(Abundance) * 100,
    SDAbundance       = sd(Abundance) * 100,
    .groups = "drop"
  )

df_sample_genus <- df_genus %>%
  filter(Genus %in% top_genera) %>%
  group_by(SampleName, Condition, Genus) %>%
  summarise(Abundance = sum(Abundance) * 100, .groups = "drop") %>%
  left_join(df_genus_stats, by = c("Condition", "Genus"))

# ================================
# 💾 Save All Final Tables
# ================================
write.csv(df_phy_stats,    file.path(table_dir, "Top10_Phylum_Stats.csv"),                           row.names = FALSE)
write.csv(df_genus_stats,  file.path(table_dir, "Top10_Genus_Stats.csv"),                            row.names = FALSE)
write.csv(df_sample_phy,   file.path(stat_dir,  "SampleWise_Top10_Phylum_Abundance_Prevalence.csv"), row.names = FALSE)
write.csv(df_sample_genus, file.path(stat_dir,  "SampleWise_Top10_Genus_Abundance_Prevalence.csv"),  row.names = FALSE)

# ================================
# 📊 Plot Function (bold axis titles + 600 dpi + PDFs)
# ================================
make_plot <- function(df, tax, cond, filename) {
  p <- ggplot(df %>% filter(Condition == cond),
              aes(x = MeanAbundance, y = PrevalencePercent, label = !!sym(tax))) +
    geom_point(size = 2.5, color = ifelse(cond == "Control", "#D73027", "#4575B4")) +
    geom_text_repel(size = 3, box.padding = 0.4, family = "Arial") +
    labs(
      title = paste(tax, "-", cond),
      x = "Relative abundance (%)",
      y = "Prevalence (%)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      text = element_text(family = "Arial"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 14, face = "bold"),  # bold x-axis title
      axis.title.y = element_text(size = 14, face = "bold"),  # bold y-axis title
      axis.text    = element_text(size = 10),
      plot.title   = element_text(size = 11, face = "bold", hjust = 0.5)
    )
  
  # --- Paths ---
  tiff_path    <- file.path(plot_dir, filename)
  pdf_vec_path <- file.path(plot_dir, sub("\\.tiff$", ".pdf", filename))
  pdf_600_path <- file.path(plot_dir, sub("\\.tiff$", "_600dpi.pdf", filename))
  
  # --- Save TIFF @ 600 dpi ---
  ggsave(
    tiff_path, plot = p,
    width = 174, height = 130, units = "mm",
    dpi = 600, device = "tiff", compression = "lzw", bg = "white"
  )
  
  # --- Save vector PDF (resolution-independent) ---
  ggsave(
    pdf_vec_path, plot = p,
    width = 174, height = 130, units = "mm",
    device = cairo_pdf, bg = "white"
  )
  
  # --- Save rasterized 600-dpi PDF by embedding the 600-dpi TIFF ---
  img <- magick::image_read(tiff_path)         # reads the 600-dpi TIFF
  magick::image_write(img, path = pdf_600_path, format = "pdf")
}

# ================================
# 🖼 Generate and Save All Plots
# ================================
make_plot(df_phy_stats,   "Phylum", "Control", "Phylum_AbundanceVsPrevalence_Control.tiff")
make_plot(df_phy_stats,   "Phylum", "Drought", "Phylum_AbundanceVsPrevalence_Drought.tiff")
make_plot(df_genus_stats, "Genus",  "Control", "Genus_AbundanceVsPrevalence_Control.tiff")
make_plot(df_genus_stats, "Genus",  "Drought", "Genus_AbundanceVsPrevalence_Drought.tiff")

cat("✅ Saved TIFF (600 dpi), vector PDF, and rasterized 600-dpi PDF under 'Figure8/Plots/'.\n")
