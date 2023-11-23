# This script plots the distribution of GCFs across datasets

# Libraries -----

library(tidyverse)
library(UpSetR)

# Load data -----

rmd_bgcs = read_tsv("https://zenodo.org/records/10182967/files/RMD-biosynthetic-genes.tsv") %>% mutate(region_clustomatic = gsub("_", "-", region))
clst = read_tsv("https://zenodo.org/records/10182967/files/RMD-biosynthetic-novelty-light.clustomatic.out", col_names = c("bgc", "gcf")) %>%
  mutate(dataset = ifelse(grepl("^MIBIG", bgc), "MIBIG", "OMD")) %>%
  mutate(dataset = ifelse(bgc %in% rmd_bgcs$region_clustomatic & grepl("^TARA", bgc), "RMD_internal", dataset)) %>%
  mutate(dataset = ifelse(bgc %in% rmd_bgcs$region_clustomatic & !grepl("^TARA", bgc), "RMD_external", dataset))
aggregated_summary = read_tsv("https://zenodo.org/records/10182967/files/genomes-aggregated-summary.tsv")
tpac_metadata = googlesheets4::read_sheet(ss = "1ukivQCmnq6PePM5vl6r3PdRJXSwZ_iChEN1LRbGD0EQ", sheet = "Sheet 1")
external_metadata = googlesheets4::read_sheet(ss = "1XkqDDo7OaUN7zrRv96ps0KZEMdKYK4xfyutJhIp09kw", sheet = "Sheet 2")

# Prep data -----

clst_summary = clst %>%
  mutate(value = 1) %>%
  spread(key = dataset, value = value) %>%
  group_by(gcf) %>%
  summarise(n_bgcs = n(),
            n_datasets = sum(any(MIBIG == 1), any(RMD_internal == 1), any(RMD_external == 1), any(OMD == 1), na.rm = T),
            RMD_internal = sum(RMD_internal == 1, na.rm = T),
            RMD_external = sum(RMD_external == 1, na.rm = T),
            OMD = sum(OMD == 1, na.rm = T),
            MIBIG = sum(MIBIG == 1, na.rm = T))

# prep metadata
metadata = rbind(
  external_metadata %>% filter(profiled == "Yes") %>% select(sample, biosample, type = biome_group),
  tpac_metadata %>% filter(profiled == "Yes") %>% select(sample, type = biome_group) %>% mutate(biosample = gsub("TARA_|_METAG", "", sample))
)

# Zoom in on the different host types
clst_rmd_meta = clst %>%
  mutate(genome = gsub("-", "_", gsub("-scaffold.*", "", bgc))) %>%
  left_join(aggregated_summary %>% select(genome, biosample) %>% left_join(metadata) %>% mutate(genome = gsub("-", "_", genome)))

# Get some stats -----

rmd_bgcs$region %>% unique %>% length

table(clst$dataset)
clst$gcf %>% unique %>% length

clst %>% filter(dataset == "MIBIG") %>% pull(gcf) %>% unique %>% length
clst %>% filter(dataset == "RMD_internal") %>% pull(gcf) %>% unique %>% length
clst %>% filter(dataset == "RMD_external") %>% pull(gcf) %>% unique %>% length
clst %>% filter(grepl("RMD", dataset)) %>% pull(gcf) %>% unique %>% length
clst %>% filter(dataset == "OMD") %>% pull(gcf) %>% unique %>% length

clst_summary %>% filter(RMD_internal > 0 & RMD_external == 0 & OMD == 0 & MIBIG == 0) %>% nrow()
clst_summary %>% filter(RMD_internal > 0 & RMD_external == 0 & OMD == 0 & MIBIG == 0) %>% nrow() / clst %>% filter(dataset == "RMD_internal") %>% pull(gcf) %>% unique %>% length
clst_summary %>% filter(RMD_internal > 0 & OMD == 0 & MIBIG == 0) %>% nrow() / clst %>% filter(dataset == "RMD_internal") %>% pull(gcf) %>% unique %>% length
clst_summary %>% filter((RMD_internal > 0 | RMD_external > 0) & OMD == 0 & MIBIG == 0) %>% nrow()
clst_summary %>% filter((RMD_internal > 0 | RMD_external > 0) & OMD == 0 & MIBIG == 0) %>% nrow() / (clst %>% filter(dataset == "RMD_internal") %>% pull(gcf) %>% unique %>% length + clst %>% filter(dataset == "RMD_external") %>% pull(gcf) %>% unique %>% length)
clst_summary %>% filter(RMD_internal == 0 & RMD_external > 0 & OMD == 0 & MIBIG == 0) %>% nrow()
clst_summary %>% filter(RMD_internal == 0 & RMD_external > 0 & OMD == 0 & MIBIG == 0) %>% nrow() / clst %>% filter(dataset == "RMD_external") %>% pull(gcf) %>% unique %>% length
clst_summary %>% filter(RMD_external > 0 & OMD == 0 & MIBIG == 0) %>% nrow() / clst %>% filter(dataset == "RMD_external") %>% pull(gcf) %>% unique %>% length

clst_rmd_meta %>% filter(type == "soft coral") %>% pull(gcf) %>% unique() %>% length()
56/26
clst_rmd_meta %>% filter(type == "stony coral") %>% pull(gcf) %>% unique() %>% length()
1403/516
clst_rmd_meta %>% filter(type == "fire coral") %>% pull(gcf) %>% unique() %>% length()
1630/410
table(clst_rmd_meta$type)

# Plot figures -----

panel_b = upset(
  fromList(
    list(
      "MIBIG" = clst %>% filter(dataset == "MIBIG") %>% pull(gcf) %>% unique,
      "RMD" = clst %>% filter(grepl("^RMD", dataset)) %>% pull(gcf) %>% unique,
      "OMD" = clst %>% filter(dataset == "OMD") %>% pull(gcf) %>% unique
    )),
  intersections = list(
    list("MIBIG"),
    list("OMD"),
    list("RMD"),
    list("RMD", "MIBIG"),
    list("RMD", "OMD"),
    list("MIBIG", "OMD"),
    list("RMD", "MIBIG", "OMD")
  ),
  keep.order = T,
  empty.intersections = "on",
  text.scale = 0.8,
  line.size = NA,
  mb.ratio = c(.6, .4),
  sets.x.label = "# GCFs",
  mainbar.y.label = "# GCFs", 
)
panel_b

panel_e = upset(
  fromList(
    list(
      "Stony_corals" = clst_rmd_meta %>% filter(type == "stony coral") %>% pull(gcf) %>% unique,
      "Fire_corals" = clst_rmd_meta %>% filter(type == "fire coral") %>% pull(gcf) %>% unique,
      "Soft_corals" = clst_rmd_meta %>% filter(type == "soft coral") %>% pull(gcf) %>% unique,
      "Sponges" = clst_rmd_meta %>% filter(type == "sponge") %>% pull(gcf) %>% unique
    )),
  intersections = list(
    list("Sponges"),
    list("Soft_corals"),
    list("Fire_corals"),
    list("Stony_corals"),
    list("Soft_corals", "Sponges"),
    list("Fire_corals", "Sponges"),
    list("Stony_corals", "Sponges"),
    list("Fire_corals", "Soft_corals"),
    list("Stony_corals", "Soft_corals"),
    list("Stony_corals", "Fire_corals"),
    list("Fire_corals", "Soft_corals", "Sponges"),
    list("Stony_corals", "Soft_corals", "Sponges"),
    list("Stony_corals", "Fire_corals", "Sponges"),
    list("Stony_corals", "Fire_corals", "Soft_corals"),
    list("Stony_corals", "Fire_corals", "Soft_corals", "Sponges")
  ),
  keep.order = T,
  empty.intersections = "on",
  text.scale = 0.8,
  line.size = NA,
  mb.ratio = c(.6, .4),
  sets.x.label = "# GCFs",
  mainbar.y.label = "# GCFs", 
)
panel_e
