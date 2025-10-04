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
library(rstatix)
library(RColorBrewer)
library(patchwork)
library(stringr)
library(ggrepel)
library(webshot2)
library(magick)
library(egg)

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
    into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
    fill = "right", sep = ";"
  ) %>%
  as.matrix() %>%
  tax_table()

metadata <- read.delim("Metadata.tsv", sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(across(everything(), trimws)) %>%
  column_to_rownames("SampleID")
metadata$SampleID <- rownames(metadata)

tree <- read_qza("rooted-tree.qza")$data %>% phy_tree()

# ================================
# 🌿 Create Phyloseq Object
# ================================
physeq <- phyloseq(asv_table, taxonomy, sample_data(metadata), tree)

# ================================
# 🧮 Calculate Alpha Diversity: Chao1 Only
# ================================
alpha_df <- microbiome::alpha(physeq, index = "chao1")
alpha_df$SampleID <- rownames(alpha_df)

# Merge with metadata
alpha_df <- inner_join(alpha_df, metadata, by = "SampleID") %>%
  rename(Chao1 = chao1)

# Ensure Condition factor
alpha_df$Condition <- factor(alpha_df$Condition, levels = c("Control","Drought"))

# ================================
# 📊 Wilcoxon Test (numeric p-values + stars)
# ================================
stats_df <- alpha_df %>%
  group_by(SampleName) %>%
  wilcox_test(Chao1 ~ Condition) %>%
  adjust_pvalue(method = "bonferroni") %>%
  mutate(
    stars = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE ~ ""
    ),
    p.label = paste0("p = ", signif(p.adj, 3), " ", stars)
  )

# y-positions for brackets per facet
max_y_vals <- alpha_df %>%
  group_by(SampleName) %>%
  summarise(y.position = max(Chao1, na.rm = TRUE) * 1.05)

stats_df <- left_join(stats_df, max_y_vals, by = "SampleName")

# ================================
# 📈 Plot Chao1 (Control = #FF5E0E, Drought = #2E8B57)
# ================================
p <- ggplot(alpha_df, aes(x = Condition, y = Chao1, fill = Condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.9) +
  facet_wrap(~ SampleName, scales = "free_y", ncol = 4) +
  stat_pvalue_manual(
    stats_df,
    label = "p.label",
    y.position = "y.position",
    xmin = "group1",
    xmax = "group2",
    tip.length = 0.01,
    size = 3,
    bracket.size = 0.6,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c("Control" = "#FF5E0E", "Drought" = "#2E8B57")) +
  labs(
    x = expression(bold("Condition")),
    y = expression(bold("Chao1 Index"))
  ) +
  theme_bw(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 8),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# ================================
# 💾 Save as TIFF and PDF in "Figure3" folder
# ================================
if (!dir.exists("Figure3")) dir.create("Figure3")

# TIFF
ggsave(
  "Figure3/Chao1_AlphaDiversity_4x2_Final_Journal.tiff",
  p, width = 6.85, height = 9.2, units = "in", dpi = 600,
  compression = "lzw", bg = "white"
)

# PDF
ggsave(
  "Figure3/Chao1_AlphaDiversity_4x2_Final_Journal.pdf",
  p, width = 6.85, height = 9.2, units = "in",
  device = cairo_pdf, bg = "white"
)

message("✅ Saved TIFF and PDF with Control=#FF5E0E (orange), Drought=#2E8B57 (green) in 'Figure3/'")
