# This script plots the proportion NP classes across the RMD and OMD

# Libraries -----

library(tidyverse)

# Resources -----

source("https://raw.github.com/SushiLab/reef-microbiomics-paper/main/resources/color_palettes.R")

# Load data -----

rmd_bgcs = read_tsv("https://zenodo.org/records/10182967/files/RMD-biosynthetic-genes.tsv") %>% mutate(region_clustomatic = gsub("_", "-", region))
clst = read_tsv("https://zenodo.org/records/10182967/files/RMD-biosynthetic-novelty-light.clustomatic.out", col_names = c("bgc", "gcf")) %>%
  mutate(dataset = ifelse(grepl("^MIBIG", bgc), "MIBIG", "OMD")) %>%
  mutate(dataset = ifelse(bgc %in% rmd_bgcs$region_clustomatic & grepl("^TARA", bgc), "RMD_internal", dataset)) %>%
  mutate(dataset = ifelse(bgc %in% rmd_bgcs$region_clustomatic & !grepl("^TARA", bgc), "RMD_external", dataset))
aggregated_summary = read_tsv("https://zenodo.org/records/10182967/files/genomes-aggregated-summary.tsv")
antismash_results = read_tsv("https://zenodo.org/records/10182967/files/genomes-antismash-summary.tsv")
antismash_results_omd = read_tsv("https://sunagawalab.ethz.ch/share/microbiomics/ocean/suppl_data/antismash-bgcs-genomes-filtered-summary.tsv.gz")
antismash_category = read_tsv("https://raw.github.com/SushiLab/reef-microbiomics-paper/main/resources/bgc_category.tsv")
antismash_category_dict = antismash_category$Summary
names(antismash_category_dict) = antismash_category$Antismash

# Prep data -----

gcf_results = antismash_results %>%
  mutate(bgc = gsub("_", "-", region)) %>%
  left_join(clst) %>%
  group_by(gcf) %>%
  summarize(products = paste0(products, collapse = ";"),
            n_bgcs = n()) %>%
  separate_rows(products, sep = ";") %>%
  group_by(gcf) %>%
  mutate(n_products = n()) %>%
  ungroup()

# Check that we have the NP class for all products
unique(gcf_results$products)[!(unique(gcf_results$products) %in% names(antismash_category_dict))]

gcf_classes = gcf_results %>% 
  mutate(np_class = antismash_category_dict[products]) %>%
  group_by(gcf, n_bgcs, np_class) %>%
  summarize(n_products = unique(n_products),
            n_np_class = n(),
            r_np_class = n_np_class/n_products) %>%
  arrange(desc(r_np_class)) %>%
  group_by(gcf) %>%
  #summarize(np_class = np_class[1])
  mutate(r = n_np_class / sum(n_np_class)) %>%
  ungroup()

gcf_classes %>% group_by(np_class) %>% summarize(freq = sum(r)) %>% mutate(freq = freq / sum(freq))

# And for OMD?
omd_gcf_results = antismash_results_omd %>%
  mutate(bgc = gsub("_", "-", region)) %>%
  left_join(clst) %>%
  group_by(gcf) %>%
  summarize(products = paste0(products, collapse = ";"),
            n_bgcs = n()) %>%
  separate_rows(products, sep = ";") %>%
  group_by(gcf) %>%
  mutate(n_products = n()) %>%
  ungroup()

# Check that we have the NP class for all products
unique(omd_gcf_results$products)[!(unique(omd_gcf_results$products) %in% names(antismash_category_dict))]

omd_gcf_classes = omd_gcf_results %>% 
  mutate(np_class = antismash_category_dict[products]) %>%
  group_by(gcf, n_bgcs, np_class) %>%
  summarize(n_products = unique(n_products),
            n_np_class = n(),
            r_np_class = n_np_class/n_products) %>%
  arrange(desc(r_np_class)) %>%
  group_by(gcf) %>%
  #summarize(np_class = np_class[1])
  mutate(r = n_np_class / sum(n_np_class)) %>%
  ungroup()

omd_gcf_classes %>% group_by(np_class) %>% summarize(freq = sum(r)) %>% mutate(freq = freq / sum(freq))

# Comparison for specific cases:
gcf_results %>% filter(grepl("arylpolyene", products)) %>% pull(gcf) %>% unique() %>% length / omd_gcf_results %>% filter(grepl("arylpolyene", products)) %>% pull(gcf) %>% unique() %>% length

gcf_results %>% filter(grepl("betalactone", products)) %>% pull(gcf) %>% unique() %>% length / omd_gcf_results %>% filter(grepl("betalactone", products)) %>% pull(gcf) %>% unique() %>% length

gcf_results %>% filter(grepl("hserlactone", products)) %>% pull(gcf) %>% unique() %>% length / omd_gcf_results %>% filter(grepl("hserlactone", products)) %>% pull(gcf) %>% unique() %>% length

gcf_results %>% filter(grepl("ectoine", products)) %>% pull(gcf) %>% unique() %>% length / omd_gcf_results %>% filter(grepl("ectoine", products)) %>% pull(gcf) %>% unique() %>% length

gcf_results %>% filter(grepl("phosphonate", products)) %>% pull(gcf) %>% unique() %>% length / omd_gcf_results %>% filter(grepl("phosphonate", products)) %>% pull(gcf) %>% unique() %>% length

gcf_results %>% filter(grepl("siderophore", products)) %>% pull(gcf) %>% unique() %>% length / omd_gcf_results %>% filter(grepl("siderophore", products)) %>% pull(gcf) %>% unique() %>% length

# Zoom in on the different host types
clst_rmd_meta = clst %>%
  mutate(genome = gsub("-", "_", gsub("-scaffold.*", "", bgc))) %>%
  left_join(aggregated_summary %>% select(genome, biosample) %>% left_join(metadata) %>% mutate(genome = gsub("-", "_", genome)))

corals_type_gcfs = clst_rmd_meta %>%
  filter(!is.na(type) & type != "NA") %>%
  select(gcf, type) %>%
  unique() %>%
  left_join(gcf_classes) %>%
  group_by(type, np_class) %>%
  summarize(n = n()) %>% #View()
  group_by(type) %>%
  mutate(r = n / sum(n),
         n_tot = sum(n)) %>% #View()
  ungroup() %>%
  select(-n, -n_tot) %>%
  spread(type, r)
View(corals_type_gcfs)
sum(corals_type_gcfs$`fire coral`,na.rm=T)
sum(corals_type_gcfs$`soft coral`,na.rm=T)
