# Load required libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(compositions)
library(zCompositions)
library(betapart)
library(umap)
library(Rtsne)
library(microbiome)
library(ggpubr)
library(dplyr)
library(pheatmap)

# Load metadata
metadata <- read_tsv("Metadata.tsv")
metadata$SampleID <- trimws(as.character(metadata$SampleID))
metadata$Condition <- as.factor(metadata$Condition)
metadata$ExperimentLocation <- as.factor(metadata$ExperimentLocation)
metadata$Plant <- as.factor(metadata$Plant)

# Ensure metadata and phyloseq object compatibility
sample_data(physeq) <- sample_data(metadata)

# Filter low-abundance taxa (5% prevalence threshold)
physeq <- filter_taxa(physeq, function(x) sum(x > 0) > 0.05 * nsamples(physeq), TRUE)

# Compute distance matrices
methods <- list(
  bray = "bray",
  unifrac = "unifrac",
  wunifrac = "wunifrac",
  aitchison = "aitchison",
  w.jaccard = "weighted_jaccard"
)

distance_matrices <- list()
for (method in names(methods)) {
  tryCatch({
    if (method == "aitchison") {
      otu_comp <- cmultRepl(as.matrix(otu_table(physeq)), method = "CZM", output = "p-counts")
      distance_matrices[[method]] <- dist(clr(otu_comp))
    } else if (method == "w.jaccard") {
      beta_abund <- beta.pair.abund(as.matrix(otu_table(physeq)), index.family = "jaccard")
      distance_matrices[[method]] <- beta_abund$beta.jac
    } else {
      distance_matrices[[method]] <- distance(physeq, method = methods[[method]])
    }
  }, error = function(e) {
    message("Distance calculation failed for ", method, ": ", e$message)
  })
}

# Perform ordinations
ordination_methods <- list(PCoA = c("bray", "unifrac", "wunifrac", "aitchison"),
                           NMDS = c("bray", "w.jaccard"),
                           tSNE = c("aitchison", "bray"),
                           UMAP = c("aitchison", "bray"))

ordination_results <- list()
for (ord in names(ordination_methods)) {
  for (dist_name in ordination_methods[[ord]]) {
    try({
      dist_mat <- as.matrix(distance_matrices[[dist_name]])
      rownames(dist_mat) <- colnames(dist_mat) <- sample_names(physeq)
      
      ord_result <- switch(ord,
                           tSNE = Rtsne(dist_mat, is_distance = TRUE, perplexity = 5, verbose = FALSE)$Y,
                           UMAP = umap(dist_mat, umap.defaults)$layout,
                           ordinate(physeq, method = ord, distance = dist_name)
      )
      
      ordination_results[[paste(ord, dist_name, sep = "_")]] <- ord_result
    }, silent = FALSE)
  }
}

# PERMANOVA Analysis
permanova_formula <- ~ Condition + ExperimentLocation + Plant
all_results <- list()
for (dist_name in names(distance_matrices)) {
  try({
    result <- adonis2(distance_matrices[[dist_name]] ~ Condition + ExperimentLocation + Plant,
                      data = metadata, permutations = 9999, by = "margin")
    result$OmegaSq <- (result$SumOfSqs - (result$R2 * result$SumOfSqs[1])) / sum(result$SumOfSqs)
    all_results[[dist_name]] <- result
  })
}

# Compile results with FDR correction
results_df <- do.call(rbind, lapply(names(all_results), function(dist_name) {
  res <- all_results[[dist_name]]
  if (!is.null(res)) {
    data.frame(Distance = dist_name, Term = rownames(res), R2 = res$R2,
               OmegaSq = res$OmegaSq, F = res$F, p.value = res$`Pr(>F)`,
               stringsAsFactors = FALSE)
  }
}))

results_df$p.adj <- p.adjust(results_df$p.value, method = "fdr")

# Filter significant results
significant_results <- filter(results_df, p.adj < 0.05 & OmegaSq > 0.01 & !Term %in% c("Residual", "Total"))

# Generate plots if significant results exist
if (nrow(significant_results) > 0) {
  for (pair in names(ordination_results)) {
    print(plot_ordination(physeq, ordination_results[[pair]], color = "Condition") +
            geom_point(size = 3) + stat_ellipse(level = 0.9) +
            ggtitle(paste("Ordination:", pair)) + theme_bw())
  }
} else {
  message("No significant associations found. Generating diagnostic plots...")
  physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))
  print(plot_bar(physeq_rel, fill = "Phylum") + facet_wrap(~Condition, scales = "free_x"))
  alpha_df <- estimate_richness(physeq)
  alpha_df$Condition <- metadata$Condition
  print(ggboxplot(alpha_df, x = "Condition", y = "Shannon", color = "Condition") +
          stat_compare_means())
}

# Save PERMANOVA and dispersion results
write.csv(results_df, "permanova_results.csv", row.names = FALSE)

# Print final summary
if (exists("significant_results")) {
  print(significant_results[order(-significant_results$OmegaSq), ])
} else {
  message("No significant associations found. Consider increasing sample size or checking for confounders.")
}
