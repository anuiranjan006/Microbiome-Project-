# ================================
# 📦 Load Required Libraries
# ================================
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
library(extrafont)  # Optional, for Arial font
# loadfonts(device = "win")  # Uncomment for Windows if needed

# ================================
# 📁 Prepare Cleaned Data
# ================================
alpha_div_clean <- alpha_div %>%
  select(
    Chao1, Shannon, Simpson,
    Condition,
    ExperimentLocation
  ) %>%
  filter(!is.na(Condition), !is.na(ExperimentLocation)) %>%
  mutate(
    Condition = factor(Condition, levels = c("Control", "Drought")),
    ExperimentLocation = factor(ExperimentLocation, levels = c("Field", "Green House", "Growth Chamber")),
    Group = paste(ExperimentLocation, Condition, sep = "_")
  )

# ================================
# 🗂️ Output folder
# ================================
out_dir <- "Experiment"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ================================
# 🎨 Function to Generate TIFF + PDF with Bold Labels
# ================================
plot_all_comparisons <- function(metric, panel_letter) {
  df <- alpha_div_clean %>% mutate(Group = factor(Group))
  
  # Statistical test
  stat_test <- df %>%
    wilcox_test(as.formula(paste0(metric, " ~ Group"))) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    mutate(
      p.adj.signif = case_when(
        p.adj < 0.001 ~ "***",
        p.adj < 0.01  ~ "**",
        p.adj < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    )
  
  y.max <- max(df[[metric]], na.rm = TRUE)
  stat_test$y.position <- seq(
    from = y.max * 1.05,
    by = y.max * 0.07,
    length.out = nrow(stat_test)
  )
  
  condition_colors <- c("Control" = "#0072B2", "Drought" = "#D55E00")
  
  # Bold axis labels (Experiment Location + Condition on x; metric on y)
  x_lab <- expression(bold("Experiment Location") * " × " * bold("Condition"))
  y_lab <- bquote(bold(.(metric)))
  
  # ✨ Final ggplot
  p <- ggplot() +
    geom_violin(
      data = df,
      aes(x = Group, y = .data[[metric]], fill = Condition),
      alpha = 0.4, width = 0.9, trim = TRUE, color = NA
    ) +
    geom_boxplot(
      data = df,
      aes(x = Group, y = .data[[metric]]),
      width = 0.08, fill = NA, color = "black", outlier.shape = NA, size = 0.4
    ) +
    geom_jitter(
      data = df,
      aes(x = Group, y = .data[[metric]], color = Condition),
      width = 0.15, size = 1.2, alpha = 0.5
    ) +
    stat_summary(
      data = df,
      aes(x = Group, y = .data[[metric]]),
      fun = mean, geom = "point", shape = 21, size = 2.5,
      fill = "white", color = "black"
    ) +
    stat_pvalue_manual(
      data = stat_test,
      label = "p.adj.signif",
      y.position = "y.position",
      xmin = "group1",
      xmax = "group2",
      tip.length = 0.01,
      bracket.size = 0.5,
      size = 4,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = condition_colors) +
    scale_color_manual(values = condition_colors) +
    labs(
      title = paste0(metric, " Diversity: All Group Comparisons"),
      x = x_lab,
      y = y_lab
    ) +
    annotate("text", x = -Inf, y = Inf, label = paste0("(", panel_letter, ")"),
             hjust = -0.5, vjust = 2, fontface = "bold", size = 6) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),  # metric/condition bolding handled via expressions above
      axis.title.y = element_text(size = 12),
      legend.position = "top",
      legend.title = element_blank(),
      text = element_text(family = "Arial")
    )
  
  # Paths
  tiff_path <- file.path(out_dir, paste0("Figure2", panel_letter, ".tiff"))
  pdf_path  <- file.path(out_dir, paste0("Figure2", panel_letter, ".pdf"))
  
  # 💾 Save as TIFF (RGB, 600 dpi)
  ggsave(
    filename = tiff_path,
    plot = p,
    width = 174 / 25.4,   # 174 mm
    height = 120 / 25.4,  # ~120 mm
    dpi = 600,
    units = "in",
    device = "tiff",
    compression = "lzw",
    bg = "white"
  )
  
  # 💾 Save as PDF (vector)
  ggsave(
    filename = pdf_path,
    plot = p,
    width = 174 / 25.4,
    height = 120 / 25.4,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
  
  return(p)
}

# ================================
# 📊 Generate and Save (TIFF + PDF) in Experiment/
# ================================
p_chao1   <- plot_all_comparisons("Chao1", "a")
p_shannon <- plot_all_comparisons("Shannon", "b")
p_simpson <- plot_all_comparisons("Simpson", "c")

# ✅ Inline preview (optional)
print(p_chao1); print(p_shannon); print(p_simpson)
