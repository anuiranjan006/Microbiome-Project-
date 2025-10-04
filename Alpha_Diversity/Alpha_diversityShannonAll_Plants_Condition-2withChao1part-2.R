# ================================
# ✅ Shannon Plot (Exactly like Chao1)
# ================================
shannon_df <- microbiome::alpha(physeq, index = "shannon")
shannon_df$SampleID <- rownames(shannon_df)

# Merge with metadata
shannon_df <- inner_join(shannon_df, metadata, by = "SampleID") %>%
  rename(Shannon = diversity_shannon)

# Make sure Condition order matches Chao1
shannon_df$Condition <- factor(shannon_df$Condition, levels = c("Control","Drought"))

# 📊 Wilcoxon Test
stats_df_shannon <- shannon_df %>%
  group_by(SampleName) %>%
  wilcox_test(Shannon ~ Condition) %>%
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

# Annotation positions
max_y_vals_shannon <- shannon_df %>%
  group_by(SampleName) %>%
  summarise(y.position = max(Shannon, na.rm = TRUE) * 1.05)

stats_df_shannon <- left_join(stats_df_shannon, max_y_vals_shannon, by = "SampleName")

# 🎨 Plot Shannon
p_shannon <- ggplot(shannon_df, aes(x = Condition, y = Shannon, fill = Condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.9) +
  facet_wrap(~ SampleName, scales = "free_y", ncol = 4) +
  stat_pvalue_manual(
    stats_df_shannon,
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
    y = expression(bold("Shannon Index"))
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

# Ensure output folder exists
if (!dir.exists("Figure3")) dir.create("Figure3")

# 💾 Save as TIFF
ggsave("Figure3/Shannon_AlphaDiversity_4x2_Final_Journal.tiff", p_shannon,
       width = 6.85, height = 9.2, units = "in", dpi = 600,
       compression = "lzw", bg = "white")

# 💾 Save as PDF (to match Chao1)
ggsave("Figure3/Shannon_AlphaDiversity_4x2_Final_Journal.pdf", p_shannon,
       width = 6.85, height = 9.2, units = "in",
       device = cairo_pdf, bg = "white")

message("✅ Figures saved as TIFF and PDF in 'Figure3/' for Shannon (Control = #FF5E0E, Drought = #2E8B57)")
