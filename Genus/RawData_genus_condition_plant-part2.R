dir.create("RAW_CSV", showWarnings = FALSE)

# 1) Data that directly drives the pies (percentages per Condition × Type × Genus) — long
df_cond_type %>%
  arrange(Condition, Type, desc(Percentage)) %>%
  write_csv("RAW_CSV/pie_input_percentage_long.csv")

# 2) Same percentages — wide matrix (one row per Condition+Type)
df_cond_type %>%
  select(Condition, Type, Genus, Percentage) %>%
  pivot_wider(names_from = Genus, values_from = Percentage, values_fill = 0) %>%
  arrange(Condition, Type) %>%
  write_csv("RAW_CSV/pie_input_percentage_wide.csv")

# 3) Per-sample relative abundance at Genus level (post WD2101 removal & Top 10 filter) — long
df_genus %>%
  write_csv("RAW_CSV/relative_abundance_top10_long.csv")

# 4) Top 10 genera list (order preserved)
tibble(Top10_Genus = top10_genera) %>%
  write_csv("RAW_CSV/top10_genera.csv")

# 5) Color key used in the plot
tibble(Genus = names(genus_colors), Color = unname(genus_colors)) %>%
  write_csv("RAW_CSV/genus_color_map.csv")

# 6) RAW COUNTS version at Genus level (pre-normalization), aligned to your filters
physeq_genus_counts <- tax_glom(physeq, taxrank = "Genus")

df_genus_counts <- psmelt(physeq_genus_counts) %>%
  left_join(metadata, by = c("Sample" = "SampleID")) %>%
  mutate(
    Condition = coalesce(Condition.x, Condition.y),
    Plant     = coalesce(Plant.x, Plant.y),
    Plant_clean = str_to_lower(gsub("[^a-zA-Z ]", "", str_trim(Plant))),
    Type = case_when(
      str_detect(Plant_clean, "oryza sativa|panicum virgatum|sorghum bicolor|zea mays") ~ "Monocot",
      str_detect(Plant_clean, "arabidopsis thaliana|medicago sativa|phaseolus vulgaris|solanum lycopersicum") ~ "Dicot",
      str_detect(Plant_clean, "monocot") ~ "Monocot",
      str_detect(Plant_clean, "dicot") ~ "Dicot",
      TRUE ~ NA_character_
    )
  ) %>%
  select(OTU, Sample, Abundance, Kingdom, Genus, Condition, Plant, Type) %>%
  filter(!is.na(Type)) %>%
  filter(Genus != "WD2101_soil_group") %>%
  filter(Genus %in% top10_genera)

write_csv(df_genus_counts, "RAW_CSV/raw_counts_top10_long.csv")

# 7) RAW COUNTS summarized by Condition × Type × Genus (+ percentages from counts)
df_cond_type_counts <- df_genus_counts %>%
  group_by(Condition, Type, Genus) %>%
  summarise(Counts = sum(Abundance), .groups = "drop") %>%
  group_by(Condition, Type) %>%
  mutate(Percentage = Counts / sum(Counts) * 100)

write_csv(df_cond_type_counts, "RAW_CSV/raw_counts_summary_by_condition_type.csv")

cat("✅ RAW CSVs written to ./RAW_CSV/\n")