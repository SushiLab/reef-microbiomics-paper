## Project: Reef Microbiomics Database
## Script purpose: Code used to generate figure 1 - genome quality
## Author: Fabienne Wiederkehr
## Date: 21/06/2023

# Load packages
require(data.table)
require(tidyverse)
require(googlesheets4)
require(patchwork)

# Load data
summary_genomes<-fread("../../../data/processed/genome-summaries/genomes-aggregated-summary_biome-taxnovelty.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  filter(!grepl("ISOG",data_type))

# Compute mean completeness and contamination scores
summary_genomes<-summary_genomes %>%
  mutate(mcpl=rowMeans(select(.,ANVIO_COMPLETENESS,CHECKM_COMPLETENESS)),
         mctn=rowMeans(select(.,ANVIO_CONTAMINATION,CHECKM_CONTAMINATION)))

# Classify genomes as of high, good, medium, or fair quality
summary_genomes_quality<-summary_genomes %>%
  mutate(quality=ifelse(mcpl>=50&mctn<=10,"Medium quality","Fair quality")) %>%
  mutate(quality=ifelse(mcpl>=70&mctn<=10,"Good quality",quality)) %>%
  mutate(quality=ifelse(mcpl>=90&mctn<=5,"High quality",quality))

table(summary_genomes_quality$quality,summary_genomes_quality$data_type)

# Determine order of genome sources
toplot<-summary_genomes_quality %>%
  mutate(source=ifelse(grepl("TPAC",data_type),"Tara Pacific",NA)) %>%
  mutate(source=ifelse(grepl("EXTERNAL-MAGS",data_type)&grepl("coral",biome_details),"External coral",source)) %>%
  mutate(source=ifelse(grepl("EXTERNAL-MAGS",data_type)&grepl("sponge",biome_details),"External sponge",source)) %>%
  mutate(source=factor(source,levels=c("External sponge","External coral","Tara Pacific"))) %>%
  mutate(quality=factor(quality,levels=c("Fair quality","Medium quality","Good quality","High quality")))

# Set colour palettes
source_colours<-readRDS("../../../data/processed/figures/tpac_colours.RDS")

# Plot
ggplot(toplot,aes(x=source,fill=source,alpha=quality)) +
  geom_bar(position="fill") +
  coord_flip() +
  scale_fill_manual(values=source_colours) +
  scale_alpha_manual(values=c(0.4,0.6,0.8,1)) +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0 %","25 %","50 %","75 %","100 %")) +
  theme_bw() +
  guides(fill="none") +
  theme(rect=element_blank(),
        plot.margin=margin(0,4,0,0,"mm"),
        legend.position="right",
        legend.background=element_rect(fill="white"),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        axis.title=element_blank(),
        axis.text=element_text(size=12,margin=margin(0,0,1,0,"mm")),
        axis.ticks=element_blank(),
        panel.grid.major.y=element_blank(),
        plot.background=element_rect(fill="white"))

#ggsave(filename="../../../data/processed/figures/figure-1/figure-1_genome-quality.png",width=10,height=5,dpi=500)
ggsave(filename="../../../data/processed/figures/figure-1/figure-1_genome-quality.pdf",width=180,height=110,units="mm")
