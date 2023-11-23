# The goal of this script is to plot the data about BGC-rich lineages

# Libraries -----

library(tidyverse)
library(patchwork)

# Resources -----

taxonomy_colors = c(
  "Acidos" = "#C9458E",
  "Cyanos" = "#55BEC5",
  "Others" = "grey50"
)

# Load data -----

aggregated_summary = read_tsv("https://zenodo.org/records/10182967/files/genomes-aggregated-summary.tsv")
tpac_metadata = googlesheets4::read_sheet(ss = "1ukivQCmnq6PePM5vl6r3PdRJXSwZ_iChEN1LRbGD0EQ", sheet = "Sheet 1")
external_metadata = googlesheets4::read_sheet(ss = "1XkqDDo7OaUN7zrRv96ps0KZEMdKYK4xfyutJhIp09kw", sheet = "Sheet 2")

# Prep data -----

# prep metadata
metadata = rbind(
  tpac_metadata %>% filter(profiled == "Yes") %>% select(sample, type = biome_group, PAN_GENOMES_FOLDER) %>% mutate(biosample = gsub("TARA_|_METAG", "", sample)),
  external_metadata %>% filter(profiled == "Yes") %>% select(sample, biosample, type = biome_group) %>% mutate(PAN_GENOMES_FOLDER = NA)
)

# We define BGC-rich lineages as having at least 15 BGCs (based on previous assessments)
candidate_bgc_rich_sp = aggregated_summary %>% filter(n_bgcs >= 15) %>% pull(secondary_cluster) %>% unique

# However, due to high genome fragmentation (low N50, high number of scaffolds), some may BGCs may be spread across multiple genomic fragments and inflate the number of BGCs
aggregated_summary %>% filter(secondary_cluster %in% candidate_bgc_rich_sp) %>% View()
# After manual inspection of the data (aggregated_summary), we identified 6 false positives:
false_positives = c("9_1", "344_1", "7_0", "2086_0", "2105_0", "1175_0")

superprod_pre_summary = aggregated_summary %>%
  filter(secondary_cluster %in% candidate_bgc_rich_sp & !secondary_cluster %in% false_positives) %>%
  select(genome, biosample, species = secondary_cluster, GTDBTK_TAXONOMY, n_bgcs, n_bgcs_complete, N50, QSCORE, COMPLETENESS, CONTAMINATION) %>%
  mutate(Taxonomy = ifelse(grepl("p__Acido", GTDBTK_TAXONOMY), "Acidos", ifelse(grepl("p__Cyano", GTDBTK_TAXONOMY), "Cyanos", "Others"))) %>%
  left_join(metadata)

# We want to select a representative genome per species
# for that we use a composite index to balance the fragmentation and number of BGCs
superprod_pre_summary_prep = superprod_pre_summary %>%
  mutate(n_bgcs_incomplete = n_bgcs - n_bgcs_complete) %>% 
  group_by(species) %>%
  mutate(norm_n50 = log(N50)/max(log(N50))*100,
         norm_bgcs = log(n_bgcs)/max(log(n_bgcs))*100) %>%
  mutate(index = QSCORE + norm_n50 + norm_bgcs) %>%
  filter(index == max(index)) %>%
  ungroup() %>% 
  mutate(taxo = paste0(GTDBTK_TAXONOMY, " | ", species),
         phylum = gsub(".*;p__|;c__.*", "", GTDBTK_TAXONOMY))

# We want to reorder the data table based on taxonomy first, number of BGC second
sp_factor_levels = superprod_pre_summary_prep %>% arrange(desc(n_bgcs)) %>% pull(taxo)
phylum_factor_levels = superprod_pre_summary_prep %>% arrange(desc(n_bgcs)) %>% pull(phylum) %>% unique

superprod_pre_summary_plot = superprod_pre_summary_prep %>%
  select(phylum, taxo, n_bgcs_complete, n_bgcs_incomplete) %>%
  gather(key = complete, value = n, -taxo, -phylum) %>%
  mutate(taxo = factor(taxo, levels = sp_factor_levels),
         phylum = factor(phylum, levels = phylum_factor_levels)) 

# Plot figure of superprod taxonomy -----

superprod_pre_summary_plot %>%
  mutate(complete = factor(complete, levels = c("n_bgcs_incomplete", "n_bgcs_complete"))) %>%
  ggplot() +
  geom_col(aes(x = taxo, y = n, fill = complete)) +
  scale_fill_manual(values = c("grey65", "grey30")) +
  facet_grid(.~phylum, scales = "free_x", space = "free_x") +
  theme_minimal() +
  geom_vline(xintercept = 15) +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size = 6),
        legend.position = "none")

# Prep & plot distribution of superprod lineages across hosts -----

superprod_summary = superprod_pre_summary %>%
  group_by(Taxonomy, type) %>%
  reframe(n = length(unique(species))) %>%
  rbind(tribble(~Taxonomy, ~type, ~n, "Others", "soft coral", 0)) %>%
  mutate(Taxonomy = factor(Taxonomy, levels = c("Others", "Cyanos", "Acidos")),
         type = factor(type, levels = rev(c("stony coral", "fire coral", "soft coral", "sponge"))))

# Raw number of species
superprod_summary %>%
  mutate() %>%
  ggplot() + 
  geom_col(aes(y = type, x = n, fill = Taxonomy)) +
  xlab("# BGC-rich species reconstructed") +
  scale_fill_manual(values = taxonomy_colors) +
  scale_x_continuous(position = "top") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.title.y = element_blank())

# normalised per samples
superprod_summary %>%
  left_join(metadata %>% group_by(type) %>% summarize(n_samples = n())) %>%
  mutate(n_superprod_per_sample = n / n_samples) %>%
  mutate(type = factor(type, levels = rev(c("stony coral", "fire coral", "soft coral", "sponge")))) %>%
  ggplot() + 
  geom_col(aes(y = type, x = n_superprod_per_sample, fill = Taxonomy)) +
  xlab("# BGC-rich species / sample") +
  scale_fill_manual(values = taxonomy_colors) +
  scale_x_continuous(position = "top") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())

# normalized per total species recovered
biome_type_n_sp = aggregated_summary %>%
  select(genome, biosample, species = secondary_cluster) %>%
  left_join(metadata) %>%
  group_by(type) %>%
  summarize(n_sp = length(unique(species)))

superprod_summary %>%
  left_join(biome_type_n_sp) %>%
  mutate(n_superprod_per_sp = n / n_sp)%>%
  mutate(type = factor(type, levels = rev(c("stony coral", "fire coral", "soft coral", "sponge")))) %>%
  ggplot() + 
  geom_col(aes(y = type, x = n_superprod_per_sp*100, fill = Taxonomy)) +
  xlab("% of BGC-rich species") +
  scale_fill_manual(values = taxonomy_colors) +
  scale_x_continuous(position = "top") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
