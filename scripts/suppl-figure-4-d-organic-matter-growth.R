# This script plots the growth curves for Sulfidibacter corallicola grown with coral-derived dissolved organic matter

# Load packages
require(data.table)
require(tidyverse)
require(ggforce)

# Load data
metadat<-fread("organic-matter-growth_metadata.csv",header=TRUE,data.table=FALSE) %>%
  filter(!timepoint=="t8") %>% #abort rate high
  filter(!organic_matter=="DOM") #not used in paper (coral-exuded dissolved organic matter)
flowdat<-fread("organic-matter-growth_cell-counts.csv",header=TRUE,data.table=FALSE) %>%
  filter(id %in% metadat$id)
  
# Append metadata to flow cytometry data
dat<-flowdat %>%
  left_join(metadat,by="id")

# Plot growth curves
ggplot(data=dat %>% filter(bacteria=="bac"),aes(x=hours,y=counts)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~organic_matter) +
  scale_y_log10() +
  theme_bw()
ggsave(filename="organic-matter-growth_growth-curves.png",width=12,height=6,dpi=500)
