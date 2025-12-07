# This script creates barplots to describe the data used in the RMD

# Load packages -----

require(data.table)
require(tidyverse)
require(openxlsx)
require(patchwork)
require(ggbreak)

# Resources -----

source("https://raw.github.com/SushiLab/reef-microbiomics-paper/main/resources/color_palettes.R")

# Load data -----

metagenomes_tpac = openxlsx::read.xlsx("https://zenodo.org/records/17829964/files/supplementary-table-1.xlsx?download=1", sheet = "Sheet 1")

metagenomes_ext = openxlsx::read.xlsx("https://zenodo.org/records/17829964/files/supplementary-table-2.xlsx?download=1", sheet = "Sheet 2") %>%
  filter(profiled == "Yes")

isogenomes_ext = openxlsx::read.xlsx("https://zenodo.org/records/17829964/files/supplementary-table-4.xlsx?download=1", sheet = "Sheet 2")

summary_genomes = fread("https://zenodo.org/records/10182967/files/genomes-aggregated-summary.tsv?download=1", header = TRUE, data.table = FALSE) %>%
  left_join(select(metagenomes_tpac, biosample, biome, biome_group), by = "biosample") %>%
  left_join(select(metagenomes_ext, biosample, biome, biome_group), by = "biosample") %>%
  mutate(biome = ifelse(!is.na(biome.x), biome.x, biome.y)) %>%
  mutate(biome_group = ifelse(!is.na(biome_group.x), biome_group.x, biome_group.y)) %>%
  select(-c(biome.x, biome.y, biome_group.x, biome_group.y))

# Prep data -----

# Tara Pacific: number of included metagenomes
n_samples_tpac = metagenomes_tpac %>%
  nrow()

# Tara Pacific: number of included MAGs
n_genomes_tpac = summary_genomes %>%
  filter(grepl("TPAC", data_type)) %>%
  nrow()

# External METAG studies: number of included studies
n_studies_sponge = metagenomes_ext %>%
  filter(grepl("Sponge", biome)) %>%
  pull(dataset) %>% unique() %>%
  length()
n_studies_coral = metagenomes_ext %>%
  filter(grepl("Coral", biome)) %>%
  pull(dataset) %>% unique() %>%
  length()

# External METAG studies: number of included metagenomes
n_samples_sponge = metagenomes_ext %>%
  filter(grepl("Sponge", biome)) %>%
  pull(dataset) %>%
  length()
n_samples_coral = metagenomes_ext %>%
  filter(grepl("Coral", biome)) %>%
  pull(dataset) %>%
  length()

# External METAG studies: number of included MAGs
n_genomes_sponge = summary_genomes %>%
  filter(grepl("Sponge", biome) & grepl("EXTERNAL-MAGS", data_type)) %>%
  nrow()
n_genomes_coral = summary_genomes %>%
  filter(grepl("Coral", biome) & grepl("EXTERNAL-MAGS", data_type)) %>%
  nrow()

# External ISOG studies: number of included studies
n_studies_isog = isogenomes_ext %>%
  pull(dataset) %>% unique() %>%
  length()

# External ISOG studies: number of included genomes
n_genomes_isog = summary_genomes %>%
  filter(grepl("ISOG", data_type)) %>%
  nrow()

# Create data frame
df<-data.frame(
  # Source
  source = c("Tara Pacific",
             "External coral",
             "External sponge",
             "Isolate genome"),
  # Number of studies
  n_studies = c(1,
                n_studies_coral,
                n_studies_sponge,
                n_studies_isog),
  # Number of samples
  n_samples = c(n_samples_tpac,
                n_samples_coral,
                n_samples_sponge,
                n_genomes_isog),
  # Number of genomes
    n_genomes = c(n_genomes_tpac,
                  n_genomes_coral,
                  n_genomes_sponge,
                  n_genomes_isog)) %>%
  # Determine order of genome sources
  mutate(source = factor(source,levels = c("Isolate genome",
                                           "External sponge",
                                           "External coral",
                                           "Tara Pacific")))

# Produce plot -----

# Plot number of studies
plot_studies = ggplot(df, aes(x = source, y = n_studies, fill = source)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = tpac_colours) +
  scale_y_continuous(breaks = c(5, 10, 15, 20, 25), labels = c("5", "10", "15", "20", "25"), minor_breaks = 0) +
  ylab("Number of studies") +
  theme_bw() +
  theme(rect = element_blank(),
        plot.margin = margin(0, 4, 0, 0, "mm"),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 12, margin = margin(0, 0, 1, 0, "mm")),
        axis.title.x = element_text(size = 12, margin = margin(1, 0, 0, 0, "mm")),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank())

# Plot number of samples
plot_samples = ggplot(df, aes(x = source, y = n_samples, fill = source)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = tpac_colours) +
  scale_y_continuous(breaks = c(200, 400, 600, 800), labels = c("200", "400", "600", "800"), minor_breaks = 0) +
  ylab("Number of samples") +
  theme_bw() +
  theme(rect = element_blank(),
        plot.margin = margin(0, 4, 0, 0, "mm"),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 12, margin = margin(0, 0, 1, 0, "mm")),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 12, margin = margin(0, 0, 0, 0, "mm")),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank())

# Plot number of genomes
ggplot(df, aes(x = source, y = n_genomes, fill = source)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = tpac_colours) +
  scale_y_continuous(breaks =c(500, 1000, 1500, 2000), minor_breaks = c(0)) +
  scale_y_break(breaks = c(2000, 10000), space = 0.5, ticklabels = c(10000, 10500, 11000, 11500)) +
  ylab("Number of genomes") +
  theme_bw() +
  theme(rect = element_blank(),
        plot.margin = margin(0, 4, 0, 0, "mm"),
        axis.ticks = element_blank(),
        axis.text.x.top = element_blank(),
        axis.text.x.bottom = element_text(size = 12, margin = margin(0, 0, 1, 0, "mm")),
        axis.text.y = element_text(size = 12, margin = margin(0, 0, 1, 0, "mm")),
        axis.title.x = element_text(size = 12, margin = margin(1, 0, 0, 0, "mm")),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank())

# Save plots as pdf files
ggsave(filename = "figure-1-e.pdf", width = 180, height = 110, units = "mm")
plots<-list(plot_studies, plot_samples)
patchwork::wrap_plots(plots, guides = "collect", ncol = 2)
ggsave(filename = "figure-1-cd.pdf", width = 180, height = 110, units = "mm")
