# Script to plot the mOTUs abundances of coral MAGs

# Libraries -----

library(tidyverse)
library(patchwork)

# Resources -----

source("https://raw.github.com/SushiLab/reef-microbiomics-paper/main/resources/color_palettes.R")

genus_dict = c(
  "TPAC_CORAL_POR" = "Porites", 
  "TPAC_CORAL_POC" = "Pocillopora", 
  "TPAC_CORAL_MIL" = "Millepora", 
  "TPAC_CORAL_HEL" = "Heliopora"
)

# Load data -----

tpac_coral_motus = read_tsv("https://zenodo.org/records/10182967/files/RMD-motus-profiles-TPAC-corals.tsv.gz")
motu_membership = read_tsv("https://zenodo.org/records/10182967/files/RMD-motus-membership.tsv")
aggregated_summary = read_tsv("https://zenodo.org/records/10182967/files/genomes-aggregated-summary.tsv")
tpac_metadata = googlesheets4::read_sheet(ss = "1ukivQCmnq6PePM5vl6r3PdRJXSwZ_iChEN1LRbGD0EQ", sheet = "Sheet 1")
tpac_metadata_water = googlesheets4::read_sheet(ss = "1ukivQCmnq6PePM5vl6r3PdRJXSwZ_iChEN1LRbGD0EQ", sheet = "Sheet 2")

# Compare mOTUs and dRep -----

motu_vs_drep = aggregated_summary %>% left_join(motu_membership) %>% select(secondary_cluster, motu) %>% filter(!is.na(motu))
length(unique(motu_vs_drep$secondary_cluster))
length(unique(motu_vs_drep$motu))
sabre::vmeasure(motu_vs_drep$secondary_cluster, motu_vs_drep$motu)

# Process data ------

# prep metadata
metadata = rbind(
  tpac_metadata %>% filter(profiled == "Yes") %>% select(sample, type = biome_group, PAN_GENOMES_FOLDER) %>% mutate(biosample = gsub("TARA_|_METAG", "", sample)),
  tpac_metadata_water %>% filter(profiled == "Yes") %>% select(sample, type = biome_group, PAN_GENOMES_FOLDER) %>% mutate(biosample = gsub("TARA_|_METAG", "", sample))
)

# we want to select the mOTUs clusters from the genomes reconstructed in Tara Pacific metagenomes
tpac_coral_motus_membership = aggregated_summary %>%
  select(genome) %>%
  filter(grepl("^TARA", genome)) %>%
  left_join(motu_membership) %>%
  filter(!is.na(motu))
tpac_coral_motus_membership %>% pull(motu) %>% unique() %>% gsub("_.*", "", .) %>% table

# Just some safety checks
assertthat::assert_that(all(tpac_coral_motus_membership$motu %in% tpac_coral_motus$motu))
assertthat::assert_that(all(tpac_coral_motus$motu %>% unique %in% tpac_coral_motus_membership$motu))

tpac_coral_motus = tpac_coral_motus %>%
  mutate(relabd = ifelse(is.na(relabd), 0, relabd)) %>%
  left_join(metadata) %>%
  mutate(type = ifelse(grepl("CORAL", PAN_GENOMES_FOLDER), genus_dict[PAN_GENOMES_FOLDER], type)) %>%
  filter(!is.na(type) & type != "blue coral") # we want to exclude samples that are not from Tara Pacific & Heliopora

# For log transformation we need to add a pseudocount, which we define as an order of magnitude below the minimal values
pseudocount = min(tpac_coral_motus %>% filter(relabd > 0) %>% pull(relabd)) / 10

tpac_coral_motus_for_plot = tpac_coral_motus %>%
  mutate(logrelabd = log(relabd + pseudocount)) %>%
  mutate(colour_pastel = paste(type, "pastel")) %>% 
  mutate(type = factor(type, levels = c("Millepora", "Porites", "Pocillopora", "Coral-surrounding water", "Reef water", "Open ocean water")))

# Plot figure -----

tpac_coral_motus_for_plot %>%
  filter(abd > 0) %>%
  group_by(type) %>%
  summarize(n = length(unique(motu)))

tpac_coral_motus_for_plot %>%
  filter(abd > 0) %>%
  ggplot() +
  geom_jitter(aes(x = type, y = logrelabd, colour = colour_pastel), size = .2, alpha = .5) + 
  geom_violin(aes(x = type, y = logrelabd, colour = type), fill = NA) +
  scale_colour_manual(values = c(tpac_colours, "Tara Pacific corals" = tpac_colours[["Tara Pacific"]])) + 
  ylim(NA, 0) +
  ylab("log(Relative abundance)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.title.x = element_blank())

