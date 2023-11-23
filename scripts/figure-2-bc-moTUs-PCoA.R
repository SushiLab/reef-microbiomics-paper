# The goal of this script is to plot the distances based on mOTUs taxonomic profiles of the Reef Microbiome Data

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

rmd_motus = read_tsv("https://zenodo.org/records/10182967/files/RMD-motus-profiles.tsv.gz")
motu_membership = read_tsv("https://zenodo.org/records/10182967/files/RMD-motus-membership.tsv")
aggregated_summary = read_tsv("https://zenodo.org/records/10182967/files/genomes-aggregated-summary.tsv")
tpac_metadata = googlesheets4::read_sheet(ss = "1ukivQCmnq6PePM5vl6r3PdRJXSwZ_iChEN1LRbGD0EQ", sheet = "Sheet 1")
tpac_metadata_water = googlesheets4::read_sheet(ss = "1ukivQCmnq6PePM5vl6r3PdRJXSwZ_iChEN1LRbGD0EQ", sheet = "Sheet 2")
external_metadata = googlesheets4::read_sheet(ss = "1XkqDDo7OaUN7zrRv96ps0KZEMdKYK4xfyutJhIp09kw", sheet = "Sheet 2")

# Prep data & compute distances -----

# prep metadata
metadata = rbind(
  tpac_metadata %>% filter(profiled == "Yes") %>% select(sample, type = biome_group, PAN_GENOMES_FOLDER) %>% mutate(biosample = gsub("TARA_|_METAG", "", sample)),
  tpac_metadata_water %>% filter(profiled == "Yes") %>% select(sample, type = biome_group, PAN_GENOMES_FOLDER) %>% mutate(biosample = gsub("TARA_|_METAG", "", sample)),
  external_metadata %>% filter(profiled == "Yes") %>% select(sample, biosample, type = biome_group) %>% mutate(PAN_GENOMES_FOLDER = NA)
)

# note that we are only interested in the mOTUs that have genomes in the RMD
# then we want to minimize noise and sparsity in the profiles, for that

# first, let's only keep the mOTUs detected in enough samples (out of the 2,074 total)
rmd_motus_occurences = rmd_motus %>%
  filter(motu %in% motu_membership$motu) %>%
  filter(motu != "unassigned") %>% # remove the unassigned fraction (shouldn't be part of the distance computation, only to get the right proportions)
  group_by(motu) %>%
  summarize(occ = sum(abd > 0))

rmd_motus_occurences %>%
  ggplot() +
  geom_density(aes(x = log10(occ))) +
  geom_vline(xintercept = log10(3)) +
  geom_vline(xintercept = log10(10)) +
  theme_minimal()

motus_to_keep = rmd_motus_occurences %>% # filter out motus with not enough observations
  filter(occ >= 10) %>%
  pull(motu)

# second, we want to remove samples with very low depth
motus_count_per_sample = rmd_motus %>%
  filter(motu %in% motus_to_keep) %>%
  group_by(sample) %>%
  summarize(totabd = sum(abd))

motus_count_per_sample %>%
  ggplot() +
  geom_density(aes(x = log10(totabd))) +
  geom_vline(xintercept = log10(3)) +
  theme_minimal()

samples_to_remove = motus_count_per_sample %>% # filter out samples with not enough observations
  filter(totabd <= 3) %>%
  pull(sample)

rmd_motus_filtered = rmd_motus %>%
  filter(sample %in% metadata$sample & !sample %in% samples_to_remove) %>%
  filter(motu %in% motus_to_keep)

rmd_motus_filtered_spread = rmd_motus_filtered %>%
  select(sample, motu, relabd) %>%
  spread(motu, relabd)

# PCoA with Tara Corals and water -----

rmd_motus_filtered_spread_tara = rmd_motus_filtered_spread %>%
  filter(grepl("TARA_", sample))

# Here we are interested in presence/absence to compute distances, i.e. Jaccard
rmd_motus_filtered_jaccard_dist_tara = stats::dist(select(rmd_motus_filtered_spread_tara, -sample), method = "binary")

# PCoA:
rmd_motus_filtered_pcoa_tara = ape::pcoa(rmd_motus_filtered_jaccard_dist_tara)

rmd_motus_filtered_pcoa_tara_tbl = 
  as_tibble(rmd_motus_filtered_pcoa_tara$vectors) %>% select(`Axis 1` = Axis.1, `Axis 2` = Axis.2) %>%
  mutate(sample = rmd_motus_filtered_spread_tara$sample) %>%
  left_join(metadata) %>%
  mutate(type = ifelse(grepl("CORAL", PAN_GENOMES_FOLDER), genus_dict[PAN_GENOMES_FOLDER], type)) %>%
  filter(!is.na(type) & type != "Heliopora")

rmd_motus_filtered_pcoa_tara_tbl %>%
  ggplot(aes(x = `Axis 1`, y = `Axis 2`, color = type)) +
  geom_point(size = 1, alpha = .8) +
  stat_ellipse(linewidth = .7) + 
  scale_color_manual(values = tpac_colours) +
  theme_minimal()

# get the statistical significance with PERMANOVA
table(rmd_motus_filtered_pcoa_tara_tbl$PAN_GENOMES_FOLDER)
vegan::adonis2(rmd_motus_filtered_jaccard_dist_tara ~ rmd_motus_filtered_pcoa_tara_tbl$PAN_GENOMES_FOLDER)

# how much variance is explained by each axis? (raw value, rounded, without negative eigenvalues)
round(rmd_motus_filtered_pcoa_tara$values$Relative_eig[1]*100, digits = 1);rmd_motus_filtered_pcoa_tara$values$Relative_eig[1];rmd_motus_filtered_pcoa_tara$values$Eigenvalues[1] / sum(rmd_motus_filtered_pcoa_tara$values$Eigenvalues[rmd_motus_filtered_pcoa_tara$values$Eigenvalues > 0])
round(rmd_motus_filtered_pcoa_tara$values$Relative_eig[2]*100, digits = 1);rmd_motus_filtered_pcoa_tara$values$Relative_eig[2];rmd_motus_filtered_pcoa_tara$values$Eigenvalues[2] / sum(rmd_motus_filtered_pcoa_tara$values$Eigenvalues[rmd_motus_filtered_pcoa_tara$values$Eigenvalues > 0])

# how many of the mOTUs are detected in the water?
rmd_motus_filtered_detect_tara = rmd_motus %>%
  filter(motu %in% motu_membership$motu) %>%
  filter(motu != "unassigned") %>%
  filter(grepl("TARA_", sample)) %>%
  left_join(metadata) %>%
  mutate(type = ifelse(grepl("CORAL", PAN_GENOMES_FOLDER), "Coral", "Water")) %>%
  group_by(motu, type) %>%
  summarize(detection = any(abd > 0)) %>%
  ungroup()

rmd_motus_filtered_detect_tara %>%
  spread(key = type, value = detection) %>%
  filter(Coral) %>%
  gather(key = type, value = detection, -motu) %>%
  group_by(motu) %>%
  summarize(detection = paste0(type[detection], collapse = ";")) %>%
  pull(detection) %>%
  table

# PCoA with Tara Corals only -----

rmd_motus_filtered_spread_tara_corals = rmd_motus_filtered_spread_tara %>%
  filter(sample %in% (metadata %>% filter(grepl("TPAC_CORAL_[MIL|POR|POC]", PAN_GENOMES_FOLDER)) %>% pull(sample)))

# Here we are interested in presence/absence to compute distances, i.e. Jaccard
rmd_motus_filtered_jaccard_dist_tara_corals = stats::dist(select(rmd_motus_filtered_spread_tara_corals, -sample), method = "binary")

# PCoA:
rmd_motus_filtered_pcoa_tara_corals = ape::pcoa(rmd_motus_filtered_jaccard_dist_tara_corals)

rmd_motus_filtered_pcoa_tara_corals_tbl = 
  as_tibble(rmd_motus_filtered_pcoa_tara_corals$vectors) %>% select(`Axis 1` = Axis.1, `Axis 2` = Axis.2) %>%
  mutate(sample = rmd_motus_filtered_spread_tara_corals$sample) %>%
  left_join(metadata) %>%
  mutate(type = ifelse(grepl("CORAL", PAN_GENOMES_FOLDER), genus_dict[PAN_GENOMES_FOLDER], type))

rmd_motus_filtered_pcoa_tara_corals_tbl %>%
  ggplot(aes(x = `Axis 1`, y = `Axis 2`, color = type)) +
  geom_point(size = 1, alpha = .8) +
  stat_ellipse(linewidth = .7) + 
  scale_color_manual(values = tpac_colours) +
  theme_minimal()

# get the statistical significance with PERMANOVA
table(rmd_motus_filtered_pcoa_tara_corals_tbl$PAN_GENOMES_FOLDER)
vegan::adonis2(rmd_motus_filtered_jaccard_dist_tara_corals ~ rmd_motus_filtered_pcoa_tara_corals_tbl$PAN_GENOMES_FOLDER)

# how much variance is explained by each axis? (raw value, rounded, without negative eigenvalues)
round(rmd_motus_filtered_pcoa_tara_corals$values$Relative_eig[1]*100, digits = 1);rmd_motus_filtered_pcoa_tara_corals$values$Relative_eig[1];rmd_motus_filtered_pcoa_tara_corals$values$Eigenvalues[1] / sum(rmd_motus_filtered_pcoa_tara_corals$values$Eigenvalues[rmd_motus_filtered_pcoa_tara_corals$values$Eigenvalues > 0])
round(rmd_motus_filtered_pcoa_tara_corals$values$Relative_eig[2]*100, digits = 1);rmd_motus_filtered_pcoa_tara_corals$values$Relative_eig[2];rmd_motus_filtered_pcoa_tara_corals$values$Eigenvalues[2] / sum(rmd_motus_filtered_pcoa_tara_corals$values$Eigenvalues[rmd_motus_filtered_pcoa_tara_corals$values$Eigenvalues > 0])

# PCoA with External Corals -----

coral_studies = external_metadata %>% filter(biome == "Coral") %>% pull(dataset) %>% unique
coral_studies = c(coral_studies, "WEBE19-1") # WEBE19-1 is a coral study but only has water samples so needs to be added back in manually

rmd_motus_filtered_spread_ext_corals = rmd_motus_filtered_spread %>%
  mutate(prj = gsub("_.*", "", sample)) %>%
  filter(prj %in% coral_studies) %>%
  select(-prj)

# Here we are interested in presence/absence to compute distances, i.e. Jaccard
rmd_motus_filtered_jaccard_dist_ext_corals = stats::dist(select(rmd_motus_filtered_spread_ext_corals, -sample), method = "binary")

# PCoA:
rmd_motus_filtered_pcoa_ext_corals = ape::pcoa(rmd_motus_filtered_jaccard_dist_ext_corals)

rmd_motus_filtered_pcoa_ext_corals_tbl = 
  as_tibble(rmd_motus_filtered_pcoa_ext_corals$vectors) %>% select(`Axis 1` = Axis.1, `Axis 2` = Axis.2) %>%
  mutate(sample = rmd_motus_filtered_spread_ext_corals$sample) %>%
  left_join(external_metadata %>% select(sample, biome))

rmd_motus_filtered_pcoa_ext_corals_tbl %>%
  ggplot(aes(x = `Axis 1`, y = `Axis 2`, color = biome)) +
  geom_point(size = 4) + 
  scale_color_manual(values = c("Coral" = "#E377C2", "Water" = "#AEC7E8")) +
  stat_ellipse() +
  theme_minimal()

# Test significance of the clustering (PERMANOVA):
vegan::adonis2(rmd_motus_filtered_jaccard_dist_ext_corals ~ rmd_motus_filtered_pcoa_ext_corals_tbl$biome)

# PERMANOVA is sensitive to group imbalance
table(rmd_motus_filtered_pcoa_ext_corals_tbl$biome)

# So to provide a better estimate of the R2, we subsample the larger group and repeat the estimation:
corals_r2 = NULL
for (i in 1:100) {
  corals_subset = c(
    rmd_motus_filtered_pcoa_ext_corals_tbl %>% filter(biome == "Coral") %>% pull(sample) %>% sample(30),
    rmd_motus_filtered_pcoa_ext_corals_tbl %>% filter(biome == "Water") %>% pull(sample)
  )
  corals_subset_filter = rmd_motus_filtered_pcoa_ext_corals_tbl$sample %in% corals_subset
  corals_subset_distances = as.dist(as.matrix(rmd_motus_filtered_jaccard_dist_ext_corals)[corals_subset_filter, corals_subset_filter])
  corals_subset_permanova = vegan::adonis2(corals_subset_distances ~ rmd_motus_filtered_pcoa_ext_corals_tbl$biome[corals_subset_filter])
  corals_r2 = c(corals_r2, corals_subset_permanova$R2[1])
}
ggplot() + geom_density(aes(x = corals_r2)) + theme_minimal()
summary(corals_r2)

# PCoA with External Sponges -----

sponge_studies = external_metadata %>% filter(biome == "Sponge") %>% pull(dataset) %>% unique

rmd_motus_filtered_spread_ext_sponges = rmd_motus_filtered_spread %>%
  mutate(prj = gsub("_.*", "", sample)) %>%
  filter(prj %in% sponge_studies) %>%
  select(-prj)

# Here we are interested in presence/absence to compute distances, i.e. Jaccard
rmd_motus_filtered_jaccard_dist_ext_sponges = stats::dist(select(rmd_motus_filtered_spread_ext_sponges, -sample), method = "binary")

# PCoA:
rmd_motus_filtered_pcoa_ext_sponges = ape::pcoa(rmd_motus_filtered_jaccard_dist_ext_sponges)

rmd_motus_filtered_pcoa_ext_sponges_tbl = 
  as_tibble(rmd_motus_filtered_pcoa_ext_sponges$vectors) %>% select(`Axis 1` = Axis.1, `Axis 2` = Axis.2) %>%
  mutate(sample = rmd_motus_filtered_spread_ext_sponges$sample) %>%
  left_join(external_metadata %>% select(sample, biome))

rmd_motus_filtered_pcoa_ext_sponges_tbl %>%
  ggplot(aes(x = `Axis 1`, y = `Axis 2`, color = biome)) +
  geom_point(size = 4) + 
  scale_color_manual(values = c("Sponge" = "#9467BD", "Water" = "#AEC7E8")) +
  stat_ellipse() +
  theme_minimal()

# Test significance of the clustering (PERMANOVA):
vegan::adonis2(rmd_motus_filtered_jaccard_dist_ext_sponges ~ rmd_motus_filtered_pcoa_ext_sponges_tbl$biome)

# PERMANOVA is sensitive to group imbalance
table(rmd_motus_filtered_pcoa_ext_sponges_tbl$biome)

# So to provide a better estimate of the R2, we subsample the larger group and repeat the estimation:
sponges_r2 = NULL
for (i in 1:100) {
  sponges_subset = c(
    rmd_motus_filtered_pcoa_ext_sponges_tbl %>% filter(biome == "Sponge") %>% pull(sample) %>% sample(50),
    rmd_motus_filtered_pcoa_ext_sponges_tbl %>% filter(biome == "Water") %>% pull(sample)
  )
  sponges_subset_filter = rmd_motus_filtered_pcoa_ext_sponges_tbl$sample %in% sponges_subset
  sponges_subset_distances = as.dist(as.matrix(rmd_motus_filtered_jaccard_dist_ext_sponges)[sponges_subset_filter, sponges_subset_filter])
  sponges_subset_permanova = vegan::adonis2(sponges_subset_distances ~ rmd_motus_filtered_pcoa_ext_sponges_tbl$biome[sponges_subset_filter])
  sponges_r2 = c(sponges_r2, sponges_subset_permanova$R2[1])
}
ggplot() + geom_density(aes(x = sponges_r2)) + theme_minimal()
summary(sponges_r2)
