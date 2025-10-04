# Microbiome Project — QIIME2 + R

This repository contains R scripts to analyze 16S rRNA QIIME2 outputs and generate publication-ready figures.

## 📂 Repository Layout
- `Alpha_Diversity/` — Scripts/outputs for Shannon, Chao1, Simpson (split violin, p-values).
- `Beta_Diversity/` — PCoA / UniFrac and related plots.
- `Phylum/` — Phylum-level relative abundance and summary.
- `Genus/` — Genus-level relative abundance, prevalence, top taxa.

> Tip: Keep raw QIIME2 files in a `data/` folder (not tracked or via Git LFS).

## 🔧 Requirements (R ≥ 4.0)
```r
install.packages(c("tidyverse","vegan","pheatmap","ggpubr","ggplot2","viridis",
                   "FSA","plotly","ggsignif","magick","htmlwidgets","webshot2"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq","microbiome"))
if(!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("jbisanz/qiime2R")
