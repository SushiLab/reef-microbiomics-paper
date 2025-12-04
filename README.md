# reef-microbiomics-paper

Repository associated with the manuscript 'Coral microbiomes as reservoirs of novel genomic and biosynthetic diversity' and the Reef Microbiomics Database.

It contains the code to generate the analyses and figures reported in the manuscript based on the data available at [Zenodo](https://zenodo.org/doi/10.5281/zenodo.10182966).

Link to manuscript to come.

## Structure

```
├── resources/ --> small tables, icons or color palettes used by the scripts
└── scripts/ --> the scripts used to reproduce analyses, figures and tables
```

## Installation

To reproduce those analysis, you need a local installation of git and R.

You can clone the repository with: `git clone git@github.com:SushiLab/reef-microbiomics-paper.git`

Install the following R packages:
- tidyverse
- googlesheets4
- UpSetR
- patchwork
- sabre
- vegan
- htmlwidgets
- ape

You can then start running the R scripts.
