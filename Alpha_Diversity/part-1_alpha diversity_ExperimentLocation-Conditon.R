# ================================
# ✅ Required Libraries (as before)
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
library(plotly)
library(ggsignif)
library(htmlwidgets)
library(webshot2)
library(gridExtra)
library(png)
library(grid)

# ================================
# ✅ Preprocess Data (same as before)
# ================================
asv_table <- read_qza("table.qza")$data
asv_table <- otu_table(asv_table, taxa_are_rows = TRUE)

taxonomy <- read_qza("silva_taxonomy.qza")$data %>%
  as.data.frame() %>%
  column_to_rownames("Feature.ID")
taxonomy <- tax_table(as.matrix(taxonomy))

metadata <- read.delim("Metadata.tsv", sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(across(everything(), trimws)) %>%
  column_to_rownames("SampleID")

metadata$SampleID <- rownames(metadata)
metadata$SampleID <- tolower(trimws(metadata$SampleID))
rownames(metadata) <- metadata$SampleID
colnames(asv_table) <- tolower(trimws(colnames(asv_table)))

tree <- read_qza("rooted-tree.qza")$data
tree <- phy_tree(tree)

asv_table <- asv_table[, rownames(metadata), drop = FALSE]
metadata <- metadata[colnames(asv_table), , drop = FALSE]

taxa_names(asv_table) <- rownames(asv_table)
taxa_names(taxonomy) <- rownames(asv_table)

physeq <- phyloseq(asv_table, taxonomy, sample_data(metadata), tree)

# ================================
# ✅ Alpha Diversity Calculation
# ================================
alpha_div <- microbiome::alpha(physeq, index = c("shannon", "simpson", "chao1"))
alpha_div$SampleID <- rownames(alpha_div)
alpha_div <- inner_join(alpha_div, metadata, by = "SampleID") %>%
  rename(
    Shannon = diversity_shannon,
    Simpson = evenness_simpson,
    Chao1 = chao1
  )

alpha_div$Blank <- ""

# ================================
# ✅ Split Violin Plot Function
# ================================
split_violin_plot_save <- function(df, group_col, metric, title_text, filename, folder = "ExpLoc_Output") {
  group_levels <- unique(sort(df[[group_col]]))
  if (length(group_levels) != 2) {
    message("Skipping plot: ", filename, " (requires exactly 2 groups)")
    return(NULL)
  }
  
  pval <- kruskal.test(as.formula(paste(metric, "~", group_col)), data = df)$p.value
  signif_star <- ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))
  p_text <- paste0("p = ", formatC(pval, digits = 4), " ", signif_star)
  
  y_max <- max(df[[metric]], na.rm = TRUE)
  y_sd <- sd(df[[metric]], na.rm = TRUE)
  y_min <- 0
  bracket_y <- y_max + 1.5 * y_sd
  bracket_y2 <- y_max + 1.2 * y_sd
  pval_y <- bracket_y + 0.3 * y_sd
  yaxis_range <- c(y_min, pval_y + 0.2 * y_sd)
  
  control_color <- "#1f77b4"
  drought_color <- "#e41a1c"
  
  p <- plot_ly(df, type = 'violin') %>%
    add_trace(
      x = df$Blank[df[[group_col]] == group_levels[1]],
      y = df[[metric]][df[[group_col]] == group_levels[1]],
      name = group_levels[1],
      side = "negative",
      box = list(visible = TRUE, fillcolor = 'white', line = list(color = 'black')),
      color = I(control_color),
      points = FALSE
    ) %>%
    add_trace(
      x = df$Blank[df[[group_col]] == group_levels[2]],
      y = df[[metric]][df[[group_col]] == group_levels[2]],
      name = group_levels[2],
      side = "positive",
      box = list(visible = TRUE, fillcolor = 'white', line = list(color = 'black')),
      color = I(drought_color),
      points = FALSE
    ) %>%
    layout(
      title = list(text = paste0("<b>", title_text, "</b>"), x = 0.5, font = list(size = 28)),
      yaxis = list(title = metric, range = yaxis_range, titlefont = list(size = 24), tickfont = list(size = 22)),
      xaxis = list(title = "", tickfont = list(size = 20), zeroline = FALSE),
      legend = list(font = list(size = 18), orientation = "v", x = 1.05, y = 0.9),
      margin = list(t = 100, b = 80, r = 120),
      shapes = list(
        list(type = "line", x0 = -0.1, x1 = 0.1, y0 = bracket_y, y1 = bracket_y, line = list(color = "black", width = 1)),
        list(type = "line", x0 = -0.1, x1 = -0.1, y0 = bracket_y2, y1 = bracket_y, line = list(color = "black", width = 1)),
        list(type = "line", x0 = 0.1, x1 = 0.1, y0 = bracket_y2, y1 = bracket_y, line = list(color = "black", width = 1))
      ),
      annotations = list(
        list(x = 0, y = pval_y, text = p_text, xref = "x", yref = "y", showarrow = FALSE,
             font = list(size = 24, color = "black"))
      )
    )
  
  htmlwidgets::saveWidget(as_widget(p), file = "temp_plot.html", selfcontained = TRUE)
  webshot2::webshot("temp_plot.html", file = file.path(folder, paste0(filename, ".png")), vwidth = 1400, vheight = 1100)
}

# ================================
# ✅ Create Output Folder
# ================================
dir.create("ExpLoc_Output", showWarnings = FALSE)

# ================================
# ✅ Generate Plots by ExperimentLocation + Condition
# ================================
plot_list <- list()
metrics <- c("Chao1", "Shannon", "Simpson")
exp_locations <- unique(alpha_div$ExperimentLocation)

for (loc in exp_locations) {
  df_loc <- alpha_div %>%
    filter(ExperimentLocation == loc, Condition %in% c("Control", "Drought"))
  
  if (nrow(df_loc) == 0) next
  
  for (metric in metrics) {
    filename <- paste0(gsub(" ", "_", loc), "_", metric)
    title <- paste(loc, "-", metric)
    
    split_violin_plot_save(
      df = df_loc,
      group_col = "Condition",
      metric = metric,
      title_text = title,
      filename = filename,
      folder = "ExpLoc_Output"
    )
    plot_list <- append(plot_list, filename)
  }
}

# ================================
# ✅ Compile to PDF
# ================================
grobs <- lapply(plot_list, function(name) {
  img <- readPNG(file.path("ExpLoc_Output", paste0(name, ".png")))
  rasterGrob(img, interpolate = TRUE)
})

pdf_grid <- marrangeGrob(grobs = grobs, nrow = 3, ncol = 2, top = "Alpha Diversity by Experiment Location and Condition")
ggsave("ExperimentLocation_AlphaDiversity.pdf", pdf_grid, width = 14, height = 18, dpi = 600)
