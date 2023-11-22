# This script creates an upset plot showing the microbial species intersection between the different datasets along with their novelty

# Load packages -----
require(data.table)
require(tidyverse)
require(googlesheets4)
require(ComplexUpset)
require(ggplot2)

# Load data -----
metagenomes_tpac<-read_sheet("https://docs.google.com/spreadsheets/d/1ukivQCmnq6PePM5vl6r3PdRJXSwZ_iChEN1LRbGD0EQ/edit#gid=0",sheet="Sheet 1") %>%
  separate(sample,into=c(NA,"biosample",NA),sep="_",remove=FALSE)

metadat_tpac<-fread("https://zenodo.org/record/6299409/files/TARA-PACIFIC_samples-provenance_20220131d.tsv?download=1",header=TRUE,data.table=FALSE,skip=1) %>%
  rename(biosample=`sample-id_biosamples`) %>%
  filter(biosample %in% metagenomes_tpac$biosample)

metagenomes_ext<-read_sheet("https://docs.google.com/spreadsheets/d/1XkqDDo7OaUN7zrRv96ps0KZEMdKYK4xfyutJhIp09kw/edit#gid=0",sheet="Sheet 2") %>%
  filter(profiled=="Yes")

isogenomes_ext<-read_sheet("https://docs.google.com/spreadsheets/d/1wWs7pNcjAg5tJTeFonR4tO--01Ap30zlCxORs464Juw/edit#gid=0",sheet="Sheet 2")

summary_genomes<-fread("https://zenodo.org/records/10182967/files/genomes-aggregated-summary.tsv?download=1",header=TRUE,data.table=FALSE) %>%
  left_join(select(metagenomes_tpac,biosample,biome,biome_group),by="biosample") %>%
  left_join(select(metagenomes_ext,biosample,biome,biome_group),by="biosample") %>%
  mutate(biome=ifelse(!is.na(biome.x),biome.x,biome.y)) %>%
  mutate(biome_group=ifelse(!is.na(biome_group.x),biome_group.x,biome_group.y)) %>%
  select(-c(biome.x,biome.y,biome_group.x,biome_group.y))

# Prep data -----

summary_novelty<-summary_genomes %>%
  mutate(novelty_category=ifelse(!grepl(";s__$",GTDBTK_TAXONOMY),"Known Species",NA)) %>%
  mutate(novelty_category=ifelse(grepl(";s__$",GTDBTK_TAXONOMY),"New Species",novelty_category)) %>%
  mutate(novelty_category=ifelse(grepl(";g__;s__",GTDBTK_TAXONOMY),"New Clade",novelty_category)) %>%
  group_by(secondary_cluster) %>%
  mutate(tara_pacific=ifelse(any(grepl("TPAC",data_type)),TRUE,FALSE)) %>%
  mutate(external_coral=ifelse(any(grepl("EXTERNAL-MAGS",data_type)&grepl("Coral",biome)),TRUE,FALSE)) %>%
  mutate(external_sponge=ifelse(any(grepl("EXTERNAL-MAGS",data_type)&grepl("Sponge",biome)),TRUE,FALSE)) %>%
  mutate(external_isolate=ifelse(any(grepl("EXTERNAL-ISOG",data_type)),TRUE,FALSE)) %>%
  ungroup()

# Dereplicate (at 95 % ANI) dataframe
derep<-distinct(summary_novelty,pick(contains("secondary_cluster_representative")),.keep_all=TRUE) %>%
  select(secondary_cluster_representative,novelty_category,tara_pacific,external_coral,external_sponge,external_isolate) %>%
  mutate(novelty_category=factor(novelty_category,levels=c("Known Species","New Species","New Clade"))) #change variable to factor

# Generate plot -----

# Pull source names
source<-rev(colnames(derep)[3:6])

# Generate upset plot
upset(derep,source,
      width_ratio=0.3,
      # Add novelty information to top barplot
      base_annotations=list(
        "Number of species (intersection)"=intersection_size(counts=TRUE,mapping=aes(alpha=novelty_category)) +
          scale_alpha_manual(values=c(0.2,0.6,1),guide="none")),
      # Add novelty information to right barplot
      set_sizes=(upset_set_size(
        geom=geom_bar(aes(alpha=novelty_category,x=group),width=0.8),
        position="right") +
          scale_alpha_manual(values=c(0.2,0.6,1)) +
          ylab("Number of species (set)")),
      labeller=ggplot2::as_labeller(c(
        "tara_pacific"="Tara Pacific MAGs",
        "external_coral"="Coral MAGs",
        "external_sponge"="Sponge MAGs",
        "external_isolate"="Isolate genomes")),
      stripes="transparent",
      name=NULL,
      sort_sets=FALSE,
      guides="over")

# Save plot as pdf file
ggsave(filename="figure-2-a-upset-databases.pdf",width=180,height=110,units="mm")
