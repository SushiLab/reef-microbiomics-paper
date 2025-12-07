# This script creates an upset plot showing the microbial species intersection between the coral groups

# Load packages -----
require(data.table)
require(tidyverse)
require(openxlsx)
require(ComplexUpset)
require(ggplot2)

# Load data -----

metagenomes_tpac = openxlsx::read.xlsx("https://zenodo.org/records/17829964/files/supplementary-table-1.xlsx?download=1", sheet = "Sheet 1")

metagenomes_ext = openxlsx::read.xlsx("https://zenodo.org/records/17829964/files/supplementary-table-2.xlsx?download=1", sheet = "Sheet 2") %>%
  filter(profiled == "Yes")

summary_genomes = fread("https://zenodo.org/records/17844029/files/genomes-aggregated-summary.tsv?download=1", header = TRUE, data.table = FALSE) %>%
  left_join(select(metagenomes_tpac, biosample, biome, biome_group), by = "biosample") %>%
  left_join(select(metagenomes_ext, biosample, biome, biome_group), by = "biosample") %>%
  mutate(biome = ifelse(!is.na(biome.x), biome.x, biome.y)) %>%
  mutate(biome_group = ifelse(!is.na(biome_group.x), biome_group.x, biome_group.y)) %>%
  select(-c(biome.x, biome.y, biome_group.x, biome_group.y))

# Prep data -----
toplot = summary_genomes %>%
  group_by(secondary_cluster) %>%
  mutate(fire_coral = ifelse(any(grepl("fire", biome_group)), TRUE, FALSE)) %>%
  mutate(blue_coral = ifelse(any(grepl("blue", biome_group)), TRUE, FALSE)) %>%
  mutate(soft_coral = ifelse(any(grepl("soft", biome_group)), TRUE, FALSE)) %>%
  mutate(stony_coral = ifelse(any(grepl("stony", biome_group)), TRUE, FALSE)) %>%
  mutate(sponge = ifelse(any(grepl("sponge", biome_group)), TRUE, FALSE)) %>%
  ungroup

# Dereplicate (at 95 % ANI) dataframe
derep = distinct(toplot, pick(contains("secondary_cluster_representative")), .keep_all = TRUE)

# Generate plots -----

# Pull source names
source = c("soft_coral", "fire_coral", "stony_coral")

# Generate upset plot
upset(derep, source,
      width_ratio = 0.3, min_degree = 1,
      # Add novelty information to top barplot
      base_annotations = list(
        "Number of species (intersection)" = intersection_size(counts = TRUE)),
      # Add novelty information to right barplot
      set_sizes = (upset_set_size(
        geom = geom_bar(width = 0.8),
        position = "right") +
          ylab("Number of species (set)")),
      labeller = ggplot2::as_labeller(c(
        "stony_coral" = "stony coral MAGs",
        "soft_coral" = "soft coral MAGs",
        "fire_coral" = "fire coral MAGs")),
      stripes = "transparent",
      name = NULL,
      sort_sets = FALSE,
      guides = "over")

# Save plot as pdf file
ggsave(filename = "figure-2-b-upset-coral-groups.pdf", width = 180, height = 110, units = "mm")

# Generate plot indicating number of samples to overlay on set size
df = data.frame(
  # Source
  source = c("Tara Pacific", "Tara Pacific", "Tara Pacific",
             "Ext coral", "Ext coral", "Ext coral"),
  # Group
  group = c("Stony coral", "Fire coral", "Soft coral",
            "Stony coral", "Fire coral", "Soft coral"),
  # Number of samples
  n_samples = c(metagenomes_tpac %>% filter(grepl("stony", biome_group)) %>% nrow(),
                metagenomes_tpac %>% filter(grepl("fire", biome_group)) %>% nrow(),
                metagenomes_tpac %>% filter(grepl("soft", biome_group)) %>% nrow(),
                metagenomes_ext %>% filter(grepl("stony", biome_group)) %>% nrow(),
                metagenomes_ext %>% filter(grepl("fire", biome_group)) %>% nrow(),
                metagenomes_ext %>% filter(grepl("soft", biome_group)) %>% nrow())) %>%
  # Determine order
  mutate(group = factor(group, levels = c("Soft coral", "Fire coral", "Stony coral")))

ggplot(df, aes(x = group, y = n_samples, fill = source)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = c("#7F7F7F", "#E9540D")) +
  scale_y_continuous(breaks = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), labels = c("100", "200", "300", "400", "500", "600", "700", "800", "900", "1000"), minor_breaks = 0) +
  ylab("Number of samples") +
  theme_bw() +
  theme(rect = element_blank(),
        plot.margin = margin(0, 4, 0, 0, "mm"),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 12, margin = margin(0, 0, 1, 0, "mm")),
        axis.title.x = element_text(size = 12, margin = margin(1, 0, 0, 0, "mm")),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank())

# Save plot as pdf file
ggsave(filename = "figure-2-b-samples.pdf", width = 180, height = 110, units = "mm")