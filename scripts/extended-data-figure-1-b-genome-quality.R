# This script creates a plot showcasing the genome quality

# Load packages ----

require(data.table)
require(tidyverse)
require(openxlsx)
require(patchwork)

# Load sources ----

source("https://raw.github.com/SushiLab/reef-microbiomics-paper/main/resources/color_palettes.R")

# Load data ----

summary_genomes = fread("https://zenodo.org/records/10182967/files/genomes-aggregated-summary.tsv?download=1", header = TRUE, data.table = FALSE) %>%
  filter(!grepl("ISOG", data_type))

metagenomes_tpac = openxlsx::read.xlsx("https://zenodo.org/records/10182966/files/supplementary-table-1.xlsx?download=1", sheet = "Sheet 1")

metagenomes_ext = openxlsx::read.xlsx("https://zenodo.org/records/10182966/files/supplementary-table-2.xlsx?download=1", sheet = "Sheet 2") %>%
  filter(profiled == "Yes") %>%
  filter(!is.na(biome_group))

# Append biome details ----
summary_biomes = summary_genomes %>%
  filter(biosample %in% c(metagenomes_tpac$biosample, metagenomes_ext$biosample)) %>% #keep only biosamples that occur in either TPAC or external metagenomes
  left_join(select(metagenomes_ext,  biosample, biome_group), by = "biosample") %>% #join biome_group from both sources
  left_join(select(metagenomes_tpac, biosample, biome_group), by = "biosample") %>%
  mutate(biome_group = coalesce(biome_group.x, biome_group.y)) %>% #resolve biome_group
  mutate(biome_details = case_when(grepl("TPAC", data_type) ~ "coral", #assign biome_details
                                   grepl("EXTERNAL-MAGS", data_type) & grepl("coral", biome_group, ignore.case = TRUE) ~ "coral",
                                   grepl("EXTERNAL-MAGS", data_type) & grepl("sponge", biome_group, ignore.case = TRUE) ~ "sponge",
                                   TRUE ~ NA_character_),
         source = case_when(grepl("TPAC", data_type) ~ "Tara Pacific",
                            grepl("EXTERNAL-MAGS", data_type) & grepl("coral", biome_group, ignore.case = TRUE) ~ "External coral",
                            grepl("EXTERNAL-MAGS", data_type) & grepl("sponge", biome_group, ignore.case = TRUE) ~ "External sponge",
                            TRUE ~ NA_character_)) %>%
  select(-c(biome_group.x, biome_group.y)) #clean up temporary columns

# Compute mean completeness and contamination scores ----

summary_scores = summary_biomes %>%
  mutate(mcpl = rowMeans(select(., ANVIO_COMPLETENESS, CHECKM_COMPLETENESS)),
         mctn = rowMeans(select(., ANVIO_CONTAMINATION, CHECKM_CONTAMINATION)))

# Classify genomes as of high, good, medium, or fair quality ----

summary_quality = summary_scores %>%
  mutate(quality = ifelse(mcpl >= 50 & mctn<= 10, "Medium quality", "Fair quality")) %>%
  mutate(quality = ifelse(mcpl >= 70 & mctn<= 10, "Good quality", quality)) %>%
  mutate(quality = ifelse(mcpl >= 90 & mctn <= 5, "High quality", quality))

table(summary_quality$quality, summary_quality$data_type)

# Determine order of genome sources
toplot = summary_quality %>%
  mutate(source = factor(source, levels = c("External sponge", "External coral", "Tara Pacific"))) %>%
  mutate(quality = factor(quality, levels = c("Fair quality", "Medium quality", "Good quality", "High quality")))

# Plot
ggplot(toplot, aes(x = source, fill = source, alpha = quality)) +
  geom_bar(position = "fill") +
  coord_flip() +
  scale_fill_manual(values = tpac_colours) +
  scale_alpha_manual(values = c(0.4, 0.6, 0.8, 1)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0 %", "25 %", "50 %", "75 %", "100 %")) +
  theme_bw() +
  guides(fill = "none") +
  theme(rect = element_blank(),
        plot.margin = margin(0, 4, 0, 0, "mm"),
        legend.position = "right",
        legend.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.text = element_text(size = 12, margin = margin(0, 0, 1, 0, "mm")),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.background = element_rect(fill = "white"))

# Save plot
ggsave(filename = "figure-1_genome-quality.pdf", width = 180, height = 110, units = "mm")
