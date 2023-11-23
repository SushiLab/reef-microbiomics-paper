# This script plots the distance of RMD GCFs to RefSeq

# Libraries -----

library(tidyverse)

# Resources -----

source("https://raw.github.com/SushiLab/reef-microbiomics-paper/main/resources/color_palettes.R")

# Load data -----

rmd_bgcs = read_tsv("https://zenodo.org/records/10182967/files/RMD-biosynthetic-genes.tsv") %>% mutate(region_clustomatic = gsub("_", "-", region))
rmd_gcfs = read_tsv("https://zenodo.org/records/10182967/files/RMD-biosynthetic-novelty-light.clustomatic.out", col_names = c("bgc", "gcf")) %>%
  mutate(dataset = ifelse(grepl("^MIBIG", bgc), "MIBIG", "OMD")) %>%
  mutate(dataset = ifelse(bgc %in% rmd_bgcs$region_clustomatic & grepl("^TARA", bgc), "RMD_internal", dataset)) %>%
  mutate(dataset = ifelse(bgc %in% rmd_bgcs$region_clustomatic & !grepl("^TARA", bgc), "RMD_external", dataset))
aggregated_summary = read_tsv("https://zenodo.org/records/10182967/files/genomes-aggregated-summary.tsv")
antismash_results = read_tsv("https://zenodo.org/records/10182967/files/genomes-antismash-summary.tsv")
antismash_category = read_tsv("https://raw.github.com/SushiLab/reef-microbiomics-paper/main/resources/bgc_category.tsv")
antismash_category_dict = antismash_category$Summary
names(antismash_category_dict) = antismash_category$Antismash
bgcs_dists = read_tsv("https://zenodo.org/records/10182967/files/RMD-bigslice-gcf-distances.tsv")
bgcs_ids = read_tsv("https://zenodo.org/records/10182967/files/RMD-bigslice-bgc-ids.tsv")

# Prep data -----

bgcs_dists_proc = antismash_results %>%
  mutate(bgc = gsub("_", "-", region)) %>%
  left_join(rmd_gcfs) %>%
  left_join(bgcs_ids, by = c("genbank file" = "orig_filename")) %>%
  left_join(bgcs_dists, by = c("id" = "bgc_id")) %>%
  left_join(aggregated_summary %>% select(genome, secondary_cluster))

# We want only 1 BGC per species per GCF to minimize sampling bias
dereplicated_bgcs = bgcs_dists_proc %>%
  group_by(gcf, secondary_cluster) %>%
  filter(length == max(length)) %>%
  filter(id == id[1]) %>%
  select(`genbank file`, length, id, gcf_id, membership_value, rank, secondary_cluster) %>%
  ungroup()

gcfs_dists_summary =  bgcs_dists_proc %>%
  filter(id %in% dereplicated_bgcs$id) %>%
  group_by(gcf) %>%
  summarize(dist_to_refseq = mean(membership_value))

# Based on the initial BiG-FAM/BiG-SLICE papers, the GCF threshold in eucl. distances is of 900
# See https://doi.org/10.1093/nar/gkaa812 
gcfs_dists_summary %>% mutate(new = dist_to_refseq > 900) %>% pull(new) %>% table / nrow(gcfs_dists_summary)

gcf_classes = bgcs_dists_proc %>%
  filter(id %in% dereplicated_bgcs$id) %>%
  separate_rows(products, sep = ";") %>%
  group_by(gcf) %>%
  mutate(n_bgcs = n()) %>%
  ungroup() %>%
  group_by(gcf, products) %>%
  reframe(r_products = n() / unique(n_bgcs),
          class = antismash_category_dict[products])

# Plot figure -----

gcfs_dists_summary %>%
  left_join(gcf_classes) %>%
  ggplot() +
  geom_histogram(aes(x = dist_to_refseq, fill = class, color = class), color = "white", bins = 40, size =.3) +
  geom_vline(xintercept = 900) +
  xlab("Avg. distance to BiG-FAM")+
  ylab("GCF count") +
  scale_fill_manual(values = bgc_colors) +
  scale_color_manual(values = bgc_colors) +
  theme_minimal() +
  theme(text = element_text(size = 6),
        rect = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        legend.position = 'none',
        plot.margin = margin(0,0,0,0, 'mm'),
        panel.grid.minor = element_blank())
