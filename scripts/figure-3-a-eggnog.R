# This script plots the proportion of genes annotated in the RMD and OMD

# Libraries -----

library(tidyverse)

# Resources -----

bar_colors = c(
  "RMD" = "#C73C74",
  "OMD" = "#6AABC6"
)

# Load data -----

# We can get the numbers from:
# RMD annotated: https://sunagawalab.ethz.ch/share/microbiomics/reef/suppl_data/RMDv1-gene-catalog-annotations-eggnogv2.1.7-5.0.2.tsv.gz
# RMD total: https://sunagawalab.ethz.ch/share/microbiomics/reef/suppl_data/RMDv1-gene-catalog-membership.tsv.gz
# OMD annotated: https://sunagawalab.ethz.ch/share/microbiomics/reef/suppl_data/OMDv1-gene-catalog-annotations-eggnogv2.1.7-5.0.2.tsv.gz
# OMD total: https://sunagawalab.ethz.ch/share/microbiomics/ocean/suppl_data/gene-catalog-membership.tsv.gz

eggnog = tribble(~proportion, ~dataset,
                0.66, "RMD", # 10,806,377 / 16,293,006
                0.84, "OMD") %>% # 14,948,783 / 17,732,801
  mutate(dataset = factor(dataset, levels = dataset))

# Plot figure -----

eggnog %>%
  ggplot() +
  geom_col(aes(x = dataset, y = proportion*100, fill = dataset)) +
  scale_fill_manual(values = bar_colors) +
  ylab("% of genes annotated (eggNOG)") +
  ylim(0,100) +
  theme_minimal() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")