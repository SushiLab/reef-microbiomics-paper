# The goal of this script is to plot the 16S distribution of Acidos superproducers 

# Libraries -----

library(tidyverse)

# Load data -----

aggregated_summary = read_tsv("https://zenodo.org/records/10182967/files/genomes-aggregated-summary.tsv")
tpac_metadata = googlesheets4::read_sheet(ss = "1ukivQCmnq6PePM5vl6r3PdRJXSwZ_iChEN1LRbGD0EQ", sheet = "Sheet 1")
otu = read_tsv("https://zenodo.org/record/4451892/files/TARA_PACIFIC_16S_otu.tsv.gz")

# Prep data -----

coral_internal_samples = tpac_metadata %>% filter(profiled == "Yes") %>% pull(internal) %>% gsub("_METAG$", "", .)
coral_samples = tpac_metadata %>% filter(profiled == "Yes") %>% pull(sample)
assertthat::assert_that(length(coral_internal_samples) == 820)

otu_colnames = names(otu)
sum(coral_internal_samples %in% otu_colnames)

otu_abd = otu %>% 
  filter(euk_prok == "PROK") %>% 
  select(any_of(c("otu", "taxonomy_full", coral_internal_samples))) %>%
  gather(key = sample, value = abd, -otu, -taxonomy_full) %>%
  group_by(sample) %>%
  mutate(relabd = abd / sum(abd)) %>%
  ungroup()

# Combining the different taxonomic levels -----

# Based on taxonomic matching analysis (first, do we have a 16S in our data, second, does any genome at GTDB from said lineage has a 16S? Then get the 16S taxonomic annotation, check for consistency)
# | phylum;Acidobacteriota
# |- class;Subgroup 21 / d__Bacteria;p__Acidobacteriota;c__UBA6911;o__RPQK01;f__;g__;s__
# |--- family;Vicinamibacteraceae / d__Bacteria;p__Acidobacteriota;c__Vicinamibacteria;o__Vicinamibacterales;f__UBA8438;g__WTFV01;s__
# |---- family;Thermoanaerobaculaceae|genus;Subgroup 10 / d__Bacteria;p__Acidobacteriota;c__Thermoanaerobaculia;o__UBA5704;f__UBA5704;g__;s__
# |--- family;Acanthopleuribacteraceae| / d__Bacteria;p__Acidobacteriota;c__Holophagae;o__Acanthopleuribacterales;f__Acanthopleuribacteraceae;g__;s__
#    |- A.pedis (tpac_16S_v1_otu0004009, based on 100% 16S rRNA match)
#    |- S.corallicola (tpac_16S_v1_otu0009664, based on >99% 16S rRNA match)


acidos_otus = otu_abd %>%
  filter(grepl("phylum;Acidobacteriota", taxonomy_full))

acidos_otus_summary = acidos_otus %>%
  group_by(sample) %>%
  summarize(relabd = sum(relabd))

sub21_otus_summary = acidos_otus %>%
  filter(grepl("phylum;Acidobacteriota.*class;Subgroup 21", taxonomy_full)) %>%
  group_by(sample) %>%
  summarize(relabd = sum(relabd))

vicina_otus_summary = acidos_otus %>%
  filter(grepl("phylum;Acidobacteriota.*family;Vicinamibacteraceae", taxonomy_full)) %>%
  group_by(sample) %>%
  summarize(relabd = sum(relabd))

thermoana_otus_summary = acidos_otus %>%
  filter(grepl("phylum;Acidobacteriota.*family;Thermoanaerobaculaceae.*genus;Subgroup 10", taxonomy_full)) %>%
  group_by(sample) %>%
  summarize(relabd = sum(relabd))

acantho_otus_summary = acidos_otus %>%
  filter(grepl("phylum;Acidobacteriota.*family;Acanthopleuribacteraceae", taxonomy_full)) %>%
  group_by(sample) %>%
  summarize(relabd = sum(relabd))

a.pedis_otus_summary = otu_abd %>%
  filter(otu == "tpac_16S_v1_otu0009664") %>%
  group_by(sample) %>%
  summarize(relabd = sum(relabd))

s.corallicola_otus_summary = otu_abd %>%
  filter(otu == "tpac_16S_v1_otu0004009") %>%
  group_by(sample) %>%
  summarize(relabd = sum(relabd))

acidos_combined_otus_summary = rbind(
  acidos_otus_summary %>% mutate(taxo = "Acidos"),
  sub21_otus_summary %>% mutate(taxo = "Sub21"),
  vicina_otus_summary %>% mutate(taxo = "Vicinas"),
  thermoana_otus_summary %>% mutate(taxo = "Thermos"),
  acantho_otus_summary %>% mutate(taxo = "Acanthos"),
  s.corallicola_otus_summary %>% mutate(taxo = "S.corallicola"),
  a.pedis_otus_summary %>% mutate(taxo = "A.pedis")
) %>%
  mutate(taxo = factor(taxo, levels = rev(c("Acidos", "Sub21", "Vicinas", "Thermos", "Acanthos", "S.corallicola", "A.pedis")))) %>%
  left_join(tpac_metadata %>% select(sample = internal, coral = PAN_GENOMES_FOLDER) %>% mutate(sample = gsub("_METAG$", "", sample)))

# Plot figure -----

acidos_combined_otus_summary %>%
  filter(!is.na(coral)) %>%
  ggplot() +
  geom_tile(aes(x = sample, y = taxo, fill = log10(relabd)), height = .9) +
  scale_fill_viridis_c(na.value = "grey90") +
  facet_grid(. ~ coral, scales = "free", space = "free", switch = "y") +
  theme_minimal() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_blank(),
        strip.text = element_text(size = 6),
        axis.text.x = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 0, vjust = 0.5),
        panel.grid = element_blank(),
        strip.placement = "outside",
        legend.position = "none")
