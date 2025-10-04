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
library(plotly)
library(ggsignif)
library(magick)
library(htmlwidgets)
library(webshot2)

# ================================
# 🎨 Match example colors (semi-transparent)
# ================================
col_left_fill  <- "rgba(27,158,119,0.85)"  # teal fill (α=0.85)
col_right_fill <- "rgba(217,95,2,0.85)"    # orange fill (α=0.85)
col_left_line  <- "rgba(16,107,83,1)"      # darker teal edge
col_right_line <- "rgba(173,73,1,1)"       # darker orange edge
col_left_pts   <- "rgba(27,158,119,0.85)"  # point color
col_right_pts  <- "rgba(217,95,2,0.85)"

# ================================
# 📊 Load and Preprocess Data
# ================================
asv_table <- read_qza("table.qza")$data
asv_table <- otu_table(asv_table, taxa_are_rows = TRUE)

taxonomy <- read_qza("silva_taxonomy.qza")$data %>%
  as.data.frame() %>% column_to_rownames("Feature.ID")
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
metadata  <- metadata[colnames(asv_table), , drop = FALSE]

taxa_names(asv_table) <- rownames(asv_table)
taxa_names(taxonomy)  <- rownames(asv_table)

physeq <- phyloseq(asv_table, taxonomy, sample_data(metadata), tree)

# ================================
# 🧬 Compute Alpha Diversity
# ================================
alpha_div <- microbiome::alpha(physeq, index = c("shannon", "simpson", "chao1"))
alpha_div$SampleID <- rownames(alpha_div)
alpha_div <- inner_join(alpha_div, metadata, by = "SampleID") %>%
  rename(
    Shannon = diversity_shannon,
    Simpson = evenness_simpson,
    Chao1   = chao1
  )
alpha_div$Blank <- ""

# ================================
# 🎨 Plotly Split Violin Plot Function (legend inside, full labels)
# ================================
split_violin_plot <- function(df, group_col, metric, flipped = FALSE, force_y_min = NULL){
  neg_num <- 1; pos_num <- 2
  if (flipped) { neg_num <- 2; pos_num <- 1 }
  group_levels <- sort(unique(as.character(df[[group_col]])))  # ensure plain text labels
  
  # stats text
  pval <- kruskal.test(as.formula(paste(metric, "~", group_col)), data = df)$p.value
  signif_star <- ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "")))
  p_text <- paste0("p = ", formatC(pval, digits = 4), " ", signif_star)
  
  y_max <- max(df[[metric]], na.rm = TRUE)
  bracket_y  <- y_max * 1.15
  bracket_y2 <- y_max * 1.10
  y_text     <- bracket_y + 0.03 * y_max
  y_top_range <- y_text + 0.05 * y_max
  mid_x <- 0
  main_title <- paste(metric, ": ", group_levels[neg_num], " vs ", group_levels[pos_num])
  
  # y-axis (optionally force min)
  yaxis_obj <- list(
    title    = list(text = paste0("<b>", metric, "</b>"), font = list(size = 20, family = "Arial")),
    tickfont = list(size = 16, family = "Arial")
  )
  if (!is.null(force_y_min)) yaxis_obj$range <- c(force_y_min, y_top_range)
  
  plot_ly(df, type = "violin") %>%
    add_trace(
      x = df$Blank[df[[group_col]] == group_levels[neg_num]],
      y = df[[metric]][df[[group_col]] == group_levels[neg_num]],
      name = group_levels[neg_num],               # legend label
      legendgroup = group_levels[neg_num],
      showlegend = TRUE,
      side = "negative",
      fillcolor = col_left_fill,
      line = list(color = col_left_line, width = 3),
      box = list(visible = TRUE, fillcolor = "darkblue", line = list(color = "black")),
      meanline = list(visible = TRUE),
      color = I(col_left_line),
      points = "all", pointpos = -0.5, jitter = 0.5,
      marker = list(size = 4, opacity = 0.85, color = col_left_pts)
    ) %>%
    add_trace(
      x = df$Blank[df[[group_col]] == group_levels[pos_num]],
      y = df[[metric]][df[[group_col]] == group_levels[pos_num]],
      name = group_levels[pos_num],               # legend label
      legendgroup = group_levels[pos_num],
      showlegend = TRUE,
      side = "positive",
      fillcolor = col_right_fill,
      line = list(color = col_right_line, width = 3),
      box = list(visible = TRUE, fillcolor = "darkblue", line = list(color = "black")),
      meanline = list(visible = TRUE),
      color = I(col_right_line),
      points = "all", pointpos = 0.5, jitter = 0.5,
      marker = list(size = 4, opacity = 0.85, color = col_right_pts)
    ) %>%
    layout(
      title = list(
        text = paste0("<b>", main_title, "</b>"),
        x = 0.5,
        font = list(size = 22, family = "Arial"),
        pad  = list(t = 8, b = 8)
      ),
      # legend INSIDE the plotting canvas to avoid cropping in PNG/PDF
      margin = list(t = 120, r = 100, b = 70, l = 80),
      yaxis = yaxis_obj,
      xaxis = list(
        title = list(text = paste0("<b>", group_col, "</b>"), font = list(size = 20, family = "Arial")),
        zeroline = FALSE,
        tickfont = list(size = 12, family = "Arial")
      ),
      font = list(family = "Arial"),
      legend = list(
        title = list(text = sprintf("<b>%s</b>", group_col)),  # legend title (“Plant” or “Condition”)
        x = 0.90, y = 1, xanchor = "right", yanchor = "top",   # inside, top-right; room for full text
        bgcolor = "rgba(255,255,255,0.95)", bordercolor = "rgba(0,0,0,0.2)", borderwidth = 1,
        orientation = "v",
        font = list(size = 16, family = "Arial")
      ),
      shapes = list(
        list(type = "line", x0 = -0.1, x1 = 0.1,  y0 = bracket_y,  y1 = bracket_y,  line = list(color = "black", width = 1)),
        list(type = "line", x0 = -0.1, x1 = -0.1, y0 = bracket_y2, y1 = bracket_y,   line = list(color = "black", width = 1)),
        list(type = "line", x0 = 0.1,  x1 = 0.1,  y0 = bracket_y2, y1 = bracket_y,   line = list(color = "black", width = 1))
      ),
      annotations = list(
        list(
          x = mid_x, y = y_text, text = p_text, xref = "x", yref = "y",
          showarrow = FALSE, font = list(size = 20, family = "Arial")
        )
      )
    )
}

# ================================
# 🪴 Create Split Violin Plots
# ================================
df_mono_dico <- alpha_div %>% filter(Plant %in% c("Monocot", "Dicot"))
df_cd        <- alpha_div %>% filter(Condition %in% c("Control", "Drought"))

p_md_shannon <- split_violin_plot(df_mono_dico, "Plant", "Shannon", flipped = TRUE)
p_md_chao1   <- split_violin_plot(df_mono_dico, "Plant", "Chao1",   flipped = TRUE, force_y_min = 0)
p_md_simpson <- split_violin_plot(df_mono_dico, "Plant", "Simpson", flipped = TRUE)

p_cd_shannon <- split_violin_plot(df_cd, "Condition", "Shannon", flipped = TRUE, force_y_min = 0)
p_cd_chao1   <- split_violin_plot(df_cd, "Condition", "Chao1",   flipped = TRUE)
p_cd_simpson <- split_violin_plot(df_cd, "Condition", "Simpson", flipped = TRUE)

# ================================
# 💾 Save HTML and base PNG using webshot2
# ================================
dir.create("Plots_HTML", showWarnings = FALSE)
dir.create("Plots_PNG",  showWarnings = FALSE)

save_plot_as_png <- function(plot, html_file, png_file) {
  saveWidget(plot, file = html_file, selfcontained = TRUE)
  webshot2::webshot(html_file, file = png_file, vwidth = 1200, vheight = 900, zoom = 3)
}

save_plot_as_png(p_cd_chao1,   "Plots_HTML/cd_chao1.html",   "Plots_PNG/Control_vs_Drought_Chao1.png")
save_plot_as_png(p_cd_shannon, "Plots_HTML/cd_shannon.html", "Plots_PNG/Control_vs_Drought_Shannon.png")
save_plot_as_png(p_cd_simpson, "Plots_HTML/cd_simpson.html", "Plots_PNG/Control_vs_Drought_Simpson.png")
save_plot_as_png(p_md_chao1,   "Plots_HTML/md_chao1.html",   "Plots_PNG/Monocot_vs_Dicot_Chao1.png")
save_plot_as_png(p_md_shannon, "Plots_HTML/md_shannon.html", "Plots_PNG/Monocot_vs_Dicot_Shannon.png")
save_plot_as_png(p_md_simpson, "Plots_HTML/md_simpson.html", "Plots_PNG/Monocot_vs_Dicot_Simpson.png")

# ================================
# 🖨 600-dpi Exports (PNG/TIFF/PDF)
# ================================
png_files <- c(
  "Plots_PNG/Control_vs_Drought_Chao1.png",
  "Plots_PNG/Control_vs_Drought_Shannon.png",
  "Plots_PNG/Control_vs_Drought_Simpson.png",
  "Plots_PNG/Monocot_vs_Dicot_Chao1.png",
  "Plots_PNG/Monocot_vs_Dicot_Shannon.png",
  "Plots_PNG/Monocot_vs_Dicot_Simpson.png"
)
out_stems <- c("Fig1","Fig2","Fig3","Fig4","Fig5","Fig6")

dir.create("Final_PNG_600dpi",   showWarnings = FALSE)
dir.create("Final_TIFFs_600dpi", showWarnings = FALSE)
dir.create("Final_PDF_600dpi",   showWarnings = FALSE)

write_600dpi <- function(in_png, stem){
  img <- image_read(in_png) |> image_convert(colorspace = "sRGB")
  image_write(img, path = file.path("Final_PNG_600dpi",  paste0(stem, ".png")),
              format = "png",  density = "600x600")
  tiff_path <- file.path("Final_TIFFs_600dpi", paste0(stem, ".tiff"))
  tryCatch(
    image_write(img, path = tiff_path, format = "tiff", density = "600x600", compression = "lzw"),
    error = function(e) image_write(img, path = tiff_path, format = "tiff", density = "600x600", compression = "zip")
  )
  image_write(img, path = file.path("Final_PDF_600dpi",   paste0(stem, ".pdf")),
              format = "pdf", density = "600x600")
}
for (i in seq_along(png_files)) write_600dpi(png_files[i], out_stems[i])

# ================================
# 💾 Extra: per-figure PDFs directly from plotly (same size as PNGs)
# ================================
dir.create("Final_HTML", showWarnings = FALSE)
dir.create("Final_PDF",  showWarnings = FALSE)

save_plot_as_pdf <- function(p, stem, vwidth = 1200, vheight = 900, zoom = 3) {
  html_file <- file.path("Final_HTML", paste0(stem, ".html"))
  pdf_file  <- file.path("Final_PDF",  paste0(stem, ".pdf"))
  saveWidget(p, file = html_file, selfcontained = TRUE)
  webshot2::webshot(html_file, file = pdf_file, vwidth = vwidth, vheight = vheight, zoom = zoom)
}

save_plot_as_pdf(p_cd_chao1,   "Control_vs_Drought_Chao1")
save_plot_as_pdf(p_cd_shannon, "Control_vs_Drought_Shannon")
save_plot_as_pdf(p_cd_simpson, "Control_vs_Drought_Simpson")
save_plot_as_pdf(p_md_chao1,   "Monocot_vs_Dicot_Chao1")
save_plot_as_pdf(p_md_shannon, "Monocot_vs_Dicot_Shannon")
save_plot_as_pdf(p_md_simpson, "Monocot_vs_Dicot_Simpson")

message("✅ Legend title + full labels shown inside (no clipping). Outputs match your sizes.")
