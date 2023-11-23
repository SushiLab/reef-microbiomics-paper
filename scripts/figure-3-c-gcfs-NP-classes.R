# This script plots the proportion NP classes across the RMD and OMD

# Libraries -----

library(tidyverse)

# Resources -----

source("https://raw.github.com/SushiLab/reef-microbiomics-paper/main/resources/color_palettes.R")

# Load data -----

gcfs_prop = googlesheets4::read_sheet(ss = "1tDBQ0sZ346AFExoFCaeBDh4HMQMIMq8OJ1Ka7mImcfE", sheet = "Sheet 3") %>%
  filter(!is.na(Class)) %>%
  select(Class, OMD = `Proportion of OMD GCFs (seq. identity)`, RMD = `Proportion of RMD GCFs (seq. identity)`) %>%
  mutate(OMD = OMD / sum(OMD),
         RMD = RMD / sum(RMD))

# Plot figure -----

gcfs_prop %>%
  gather(key = db, value = frequency, - Class) %>%
  mutate(db = factor(db, levels = c("RMD", "OMD"))) %>%
  ggplot() + 
  geom_col(aes(x = db, y = frequency, fill = Class)) +
  ylab("Proportion") +
  scale_fill_manual(values = bgc_colors_fullname) +
  theme_minimal() +
  theme(text = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.title.x = element_blank(),
        legend.position = "none")
