# Load necessary libraries
# Load necessary libraries
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

# Ensure SampleID is present as a column in metadata and row names are properly cleaned
metadata$SampleID <- rownames(metadata)
metadata$SampleID <- trimws(tolower(metadata$SampleID))  # Ensure consistent case and no extra spaces

# Standardize Sample Names
rownames(metadata) <- tolower(trimws(rownames(metadata)))
colnames(asv_table) <- tolower(trimws(colnames(asv_table)))

# Load Phylogenetic Tree
tree <- read_qza("rooted-tree.qza")$data
tree <- phy_tree(tree)

# Create Phyloseq Object
asv_table <- asv_table[, rownames(metadata), drop = FALSE]
metadata <- metadata[colnames(asv_table), , drop = FALSE]
taxa_names(asv_table) <- rownames(asv_table)
taxa_names(taxonomy) <- rownames(asv_table)

physeq <- phyloseq(asv_table, taxonomy, sample_data(metadata), tree)

# Compute alpha diversity using microbiome package
alpha_div <- microbiome::alpha(physeq, index = c("shannon", "simpson", "chao1"))

# Ensure SampleID is added to alpha_div
alpha_div$SampleID <- rownames(alpha_div)

# Merge with metadata for statistical analysis
alpha_div <- inner_join(alpha_div, metadata, by = "SampleID")

# Rename columns for clarity
alpha_div <- alpha_div %>% rename(
  Shannon = diversity_shannon,
  Simpson = evenness_simpson,
  Chao1 = chao1
)

# Define categorical variables for analysis
categorical_vars <- c("SampleName", "Plant", "Condition", "ExperimentLocation")

# Initialize results data frame for statistical tests
stat_results <- data.frame()

# Loop through diversity metrics and categorical variables for statistical tests
for (metric in c("Shannon", "Simpson", "Chao1")) {
  for (var in categorical_vars) {
    n_groups <- length(unique(alpha_div[[var]]))
    
    if (n_groups < 2) {
      next  # Skip if not enough groups
    }
    
    # Non-normal distribution assumed for all
    if (n_groups == 2) {
      test <- wilcox.test(alpha_div[[metric]] ~ alpha_div[[var]])
      test_type <- "Wilcoxon"
    } else {
      test <- kruskal.test(alpha_div[[metric]] ~ alpha_div[[var]])
      test_type <- "Kruskal-Wallis"
    }
    
    stat_results <- rbind(stat_results, data.frame(
      Test_Type = test_type,
      Scenario = var,
      Diversity_Metric = metric,
      p_value = test$p.value
    ))
  }
}

# Save statistical results to CSV
write.csv(stat_results, "alpha_div_stats_results.csv", row.names = FALSE)
print("✅ Statistical analysis completed! Results saved as alpha_diversity_statistical_tests.csv.")

# Reshape for visualization
alpha_div_long <- pivot_longer(alpha_div, cols = c("Shannon", "Simpson", "Chao1"), names_to = "Metric", values_to = "Value")

# Create plots
p_alpha <- ggplot(alpha_div_long, aes(x = Condition, y = Value, fill = SampleName, shape = Plant)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_bw() +
  labs(title = "Alpha Diversity by Condition", x = "Condition", y = "Diversity Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Load necessary packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Custom dark color palette
fill_colors <- c(
  "Arabidopsis thaliana" = "#b22222",
  "Medicago sativa L." = "#1e90ff",
  "Oryza sativa" = "#228b22",
  "Panicum virgatum" = "#8a2be2",
  "Phaseolus vulgaris" = "#ff8c00",
  "Solanum lycopersicum" = "#ffd700",
  "Sorghum bicolor" = "#a0522d",
  "Zea mays" = "#da70d6"
)

# Function to create plot for a given metric
plot_alpha_diversity <- function(metric_name) {
  df_metric <- alpha_div_long %>% filter(Metric == metric_name)
  
  # Calculate shape points for median values
  alpha_shape_points <- df_metric %>%
    group_by(SampleName, Condition, Metric, Plant) %>%
    summarise(Value = median(Value), .groups = "drop")
  
  # Plot
  p <- ggplot(df_metric, aes(x = Condition, y = Value, fill = SampleName)) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.8,
      width = 0.6,
      color = "black",
      aes(group = interaction(Condition, SampleName))
    ) +
    geom_point(
      data = alpha_shape_points,
      aes(x = Condition, y = Value, shape = Plant, group = SampleName, fill = SampleName),
      size = 3,
      color = "black",
      position = position_dodge2(width = 0.6, preserve = "single"),
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = fill_colors) +
    scale_shape_manual(values = c("Dicot" = 16, "Monocot" = 17)) +
    labs(
      title = paste(metric_name, "Index"),
      x = "Condition",
      y = paste(metric_name, "Value"),
      fill = "Sample Name",
      shape = "Plant Type"
    ) +
    guides(
      fill = guide_legend(override.aes = list(shape = NA)),
      shape = guide_legend(override.aes = list(fill = "black"))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 12),
      legend.position = "right",
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", fill = NA),
      plot.title = element_text(face = "bold")
    )
  
  return(p)
}

# Generate plots for each metric
p1 <- plot_alpha_diversity("Chao1")
p2 <- plot_alpha_diversity("Shannon")
p3 <- plot_alpha_diversity("Simpson")

# Print plots
print(p1)
print(p2)
print(p3)

# Save plots
ggsave("chao1_plot.png", p1, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("shannon_plot.png", p2, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("simpson_plot.png", p3, width = 10, height = 6, dpi = 300, bg = "white")


# Load required libraries
library(dplyr)

# ========================
# 🎨 Custom color palette for SampleNames
# ========================
sample_colors <- c(
  "Arabidopsis thaliana"   = "#b22222",
  "Medicago sativa L."     = "#1f78b4",
  "Oryza sativa"           = "#33a02c",
  "Panicum virgatum"       = "#ff7f00",
  "Phaseolus vulgaris"     = "#6a3d9a",
  "Solanum lycopersicum"   = "#a6cee3",
  "Sorghum bicolor"        = "#b2df8a",
  "Zea mays"               = "#fb9a99"
)

# ========================
# ⚙️ Function: Run pairwise Wilcoxon for a variable and grouping factor
# ========================
run_pairwise_wilcox <- function(data, variable, group) {
  tryCatch({
    result <- pairwise.wilcox.test(data[[variable]], data[[group]], p.adjust.method = "bonferroni")
    df <- as.data.frame(as.table(result$p.value)) %>%
      rename(Group1 = Var1, Group2 = Var2, p.value = Freq) %>%
      mutate(
        Variable = variable,
        Grouping = group,
        Test = "Pairwise Wilcoxon"
      )
    return(df)
  }, error = function(e) {
    data.frame(
      Variable = variable,
      Grouping = group,
      Group1 = NA,
      Group2 = NA,
      p.value = NA,
      Test = "Pairwise Wilcoxon",
      Error = e$message
    )
  })
}

# ========================
# 🔁 Run tests for all alpha metrics & groupings
# ========================
alpha_metrics <- c("Shannon", "Simpson", "Chao1")
grouping_factors <- c("SampleName", "Condition", "ExperimentLocation", "Plant")

# Run all tests
pairwise_all_results <- do.call(rbind, lapply(alpha_metrics, function(var) {
  do.call(rbind, lapply(grouping_factors, function(grp) {
    run_pairwise_wilcox(alpha_div, var, grp)
  }))
}))

# ========================
# ✳️ Add significance stars
# ========================
pairwise_all_results <- pairwise_all_results %>%
  mutate(Significance = case_when(
    is.na(p.value) ~ NA_character_,
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

# ========================
# 📊 Summary of significant results
# ========================
sig_summary <- pairwise_all_results %>%
  filter(Significance %in% c("*", "**", "***")) %>%
  group_by(Variable, Grouping) %>%
  summarise(Significant_Pairs = n(), .groups = "drop")

print("📌 Summary of significant pairwise comparisons (p < 0.05):")
print(sig_summary)

# ========================
# 💾 Save all results
# ========================
write.csv(pairwise_all_results, "pairwise_wilcoxon_all_groups_with_significance.csv", row.names = FALSE)
write.csv(sig_summary, "pairwise_wilcoxon_significance_summary.csv", row.names = FALSE)

# ========================
# ✅ Done
# ========================
print("✅ All pairwise Wilcoxon tests completed across SampleName, Condition, ExperimentLocation, and Plant.")
print("📁 Full results: pairwise_wilcoxon_all_groups_with_significance.csv")
print("📁 Summary: pairwise_wilcoxon_significance_summary.csv")
print("🎨 SampleName color mapping ready for ggplot or visualizations.")


# Compute beta diversity distances
bray_dist <- distance(physeq, method = "bray")

# Perform ordination
ordination_pcoa <- ordinate(physeq, method = "PCoA", distance = "bray")
ordination_nmds <- ordinate(physeq, method = "NMDS", distance = "bray")

# Convert ordination results to dataframe
ordination_df <- data.frame(ordination_pcoa$vectors[, 1:2])
ordination_df$SampleID <- rownames(ordination_df)

# Merge with metadata
ordination_df <- inner_join(ordination_df, metadata, by = "SampleID")

# Rename columns
colnames(ordination_df)[1:2] <- c("PCoA1", "PCoA2")

# Plot PCoA
p_pcoa <- ggplot(ordination_df, aes(x = PCoA1, y = PCoA2, color = Condition)) +
  geom_point(size = 4) + theme_bw() +
  labs(title = "PCoA - Bray-Curtis")

ggsave("beta_diversity_pcoa.png", p_pcoa, width = 8, height = 5, dpi = 300)

# PERMANOVA tests
permanova_condition <- adonis2(bray_dist ~ Condition, data = metadata, permutations = 999)
permanova_plant <- adonis2(bray_dist ~ Plant, data = metadata, permutations = 999)
permanova_experiment <- adonis2(bray_dist ~ ExperimentLocation, data = metadata, permutations = 999)

# Save PERMANOVA results
permanova_results <- data.frame(
  Factor = c("Condition", "Plant", "ExperimentLocation"),
  R2 = c(permanova_condition$R2[1], permanova_plant$R2[1], permanova_experiment$R2[1]),
  p_value = c(permanova_condition$`Pr(>F)`[1], permanova_plant$`Pr(>F)`[1], permanova_experiment$`Pr(>F)`[1])
)
write.csv(permanova_results, "permanova_results.csv", row.names = FALSE)

print("✅ Microbiome diversity analysis completed. Results saved.")
