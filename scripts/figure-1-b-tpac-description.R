# This script plots a summary of the metagenomes collected during the Tara Pacific expedition

# Load packages -----

require(data.table)
require(tidyverse)
require(openxlsx)
require(stringr)
require(ComplexHeatmap)
require(grDevices)

# Load data -----
metagenomes_tpac_coral = openxlsx::read.xlsx("https://zenodo.org/records/17829964/files/supplementary-table-1.xlsx?download=1", sheet = "Sheet 1") %>%
  mutate(biome_group = word(type, 2))

metagenomes_tpac_water = openxlsx::read.xlsx("https://zenodo.org/records/17829964/files/supplementary-table-1.xlsx?download=1", sheet = "Sheet 2")

metadat_tpac = fread("https://zenodo.org/record/6299409/files/TARA-PACIFIC_samples-provenance_20220131d.tsv?download=1", header = TRUE, data.table = FALSE, skip = 1) %>%
  rename(biosample = `sample-id_biosamples`)

# Prep data -----

metagenomes_tpac = metadat_tpac %>%
  select(biosample, `sampling-design_label`) %>%
  filter(biosample %in% append(metagenomes_tpac_coral$biosample, metagenomes_tpac_water$biosample)) %>%
  left_join(select(metagenomes_tpac_water, biosample, biome, biome_group), by = "biosample") %>%
  left_join(select(metagenomes_tpac_coral, biosample, biome, biome_group), by = "biosample") %>%
  mutate(biome = ifelse(!is.na(biome.x), biome.x, biome.y)) %>%
  mutate(biome_group = ifelse(!is.na(biome_group.x), biome_group.x, biome_group.y)) %>%
  select(-c(biome.x, biome.y, biome_group.x, biome_group.y)) %>%
  mutate(biome_group = ifelse(is.na(biome_group), biome,biome_group)) %>%
  separate(`sampling-design_label`, c(NA, "island", "site", "colony"), "-")
  
# Summarise number of samples per host/water type and site
toplot = metagenomes_tpac %>%
  select(biome_group, island, site, colony) %>%
  pivot_wider(names_from = c(biome_group, site), values_from = colony, values_fn = length, values_fill = 0) %>%
  mutate("W1_C" = rowSums(across(contains("Coral-surrounding"))), .keep = "unused") %>%
  mutate("W2_R" = rowSums(across(contains("Reef"))), .keep = "unused") %>%
  rename("W3_O" = contains("Open")) %>%
  arrange(island) %>%
  select(order(colnames(.))) %>%
  column_to_rownames("island")

# Group sites with more than 4 samples together
toplot[toplot >= 4] = "4+"

# Produce plot -----

# Generate the colour palette
colour_palette = colorRampPalette(c("#FFFFFF", "#E9540D"))
colours = colour_palette(5)

# Produce heatmap
p = Heatmap(toplot, col = colours, rect_gp = gpar(col = "grey"),
        row_names_side = "left", column_names_side = "top",
        column_labels = sub(".*_", "", colnames(toplot)),
        # Add splits
        column_split = c(rep("0 - Heliopora (2)", 1), rep("1 - Millepora (201)", 5), rep("2 - Pocillopora (308)", 5), rep("3 - Porites (309)", 5), rep("5 - Seawater (387)", 3)),
        column_title_rot = 60,
        heatmap_legend_param = list(title = "Number of samples"),
        width = ncol(toplot) * unit(5, "mm"), 
        height = nrow(toplot) * unit(5, "mm"))

# Save plot as pdf file
pdf("figure-1-b-tpac-description.pdf", width = 7, height = 10)
p
dev.off()
