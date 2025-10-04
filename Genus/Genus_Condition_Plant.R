# ========================
# 📦 Load Required Libraries
# ========================
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
library(rstatix)
library(RColorBrewer)
library(patchwork)
library(stringr)

# ========================
# 📥 Load Input Data
# ========================
asv_table <- read_qza("table.qza")$data
asv_table <- otu_table(asv_table, taxa_are_rows = TRUE)

taxonomy <- read_qza("silva_taxonomy.qza")$data %>%
  as.data.frame() %>%
  column_to_rownames("Feature.ID") %>%
  mutate(Taxon = as.character(Taxon)) %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", extra = "merge") %>%
  mutate(across(everything(), ~str_remove(., "^[a-zA-Z]__")))
taxonomy <- tax_table(as.matrix(taxonomy))

metadata <- read.delim("Metadata.tsv", sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(SampleID = tolower(trimws(SampleID)))
rownames(metadata) <- metadata$SampleID

colnames(asv_table) <- tolower(trimws(colnames(asv_table)))
common_samples <- intersect(colnames(asv_table), rownames(metadata))
asv_table <- asv_table[, common_samples]
metadata <- metadata[common_samples, ]
tree <- read_qza("rooted-tree.qza")$data
tree <- phy_tree(tree)

physeq <- phyloseq(asv_table, taxonomy, sample_data(metadata), tree)

# ========================
# 🔄 Relative Abundance + Agglomerate to Genus
# ========================
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))
physeq_genus <- tax_glom(physeq_rel, taxrank = "Genus")

df_genus <- psmelt(physeq_genus) %>%
  left_join(metadata, by = c("Sample" = "SampleID")) %>%
  mutate(
    Condition = coalesce(Condition.x, Condition.y),
    Plant = coalesce(Plant.x, Plant.y)
  ) %>%
  select(OTU, Sample, Abundance, Kingdom, Genus, Condition, Plant)

# ========================
# ❌ Remove WD2101_soil_group
# ========================
df_genus <- df_genus %>% filter(Genus != "WD2101_soil_group")

# ========================
# 🔝 Top 10 Genera (After Removal)
# ========================
top10_genera <- df_genus %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

df_genus <- df_genus %>% filter(Genus %in% top10_genera)

# 🎨 Assign Colors
set.seed(123)
genus_colors <- brewer.pal(n = 10, name = "Paired")
names(genus_colors) <- top10_genera

# ========================
# 🌱 Map Plant Type
# ========================
df_genus <- df_genus %>%
  mutate(
    Plant_clean = str_to_lower(gsub("[^a-zA-Z ]", "", str_trim(Plant))),
    Type = case_when(
      str_detect(Plant_clean, "oryza sativa|panicum virgatum|sorghum bicolor|zea mays") ~ "Monocot",
      str_detect(Plant_clean, "arabidopsis thaliana|medicago sativa|phaseolus vulgaris|solanum lycopersicum") ~ "Dicot",
      str_detect(Plant_clean, "monocot") ~ "Monocot",
      str_detect(Plant_clean, "dicot") ~ "Dicot",
      TRUE ~ NA_character_
    )
  )

df_genus <- df_genus %>% filter(!is.na(Type))

# ========================
# 📊 Pie Chart Data (Genus Level)
# ========================
df_cond_type <- df_genus %>%
  group_by(Condition, Type, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Condition, Type) %>%
  mutate(
    Percentage = Abundance / sum(Abundance) * 100,
    Label = paste0(round(Percentage, 1), "%"),
    LabelText = ifelse(Percentage >= 5, Label, "")  # Show only if >= 5%
  )

# ========================
# 🧩 Function for Clean Pie Chart
# ========================
get_plot <- function(cond, type) {
  df_sub <- df_cond_type %>% filter(Condition == cond, Type == type)
  
  ggplot(df_sub, aes(x = "", y = Percentage, fill = Genus)) +
    geom_col(width = 1, color = "white") +
    geom_text(aes(label = LabelText),
              position = position_stack(vjust = 0.5),
              size = 4.2, color = "black", fontface = "bold") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = genus_colors, na.translate = FALSE) +
    theme_void(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13)
    ) +
    labs(title = paste(cond, "-", type))
}

# ========================
# 🖼️ Build 2x2 Panel
# ========================
p1 <- get_plot("Control", "Monocot")
p2 <- get_plot("Control", "Dicot")
p3 <- get_plot("Drought", "Monocot")
p4 <- get_plot("Drought", "Dicot")

combined <- (p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# ========================
# 💾 Save Output
# ========================
ggsave("Refined_PieCharts_Genus.png", combined, width = 14, height = 12, dpi = 600, bg = "white")
ggsave("Refined_PieCharts_Genus.pdf", combined, width = 14, height = 12, bg = "white")

cat("✅ Refined pie charts saved with label clarity and WD2101 removed.\n")
