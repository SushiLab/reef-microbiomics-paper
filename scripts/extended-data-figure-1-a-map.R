# This script plots a world map showing the sampling sites of Tara Pacific and the included external coral and sponge studies

# Load packages -----

require(data.table)
require(tidyverse)
require(leaflet)
require(openxlsx)

# Load functions -----

wrapper_save_html = function(html_widget, path, self = TRUE){
  current = getwd()
  dir = dirname(path)
  file = basename(path)
  setwd(dir)
  htmlwidgets::saveWidget(html_widget, file = file, selfcontained = self)
  setwd(current)
}

# Load data -----

tpac_metagenomes = openxlsx::read.xlsx("https://zenodo.org/records/10182966/files/supplementary-table-1.xlsx?download=1", sheet = "Sheet 1")
ext_metagenomes = openxlsx::read.xlsx("https://zenodo.org/records/10182966/files/supplementary-table-2.xlsx?download=1", sheet = "Sheet 2")
metadat_tpac = fread("https://zenodo.org/record/6299409/files/TARA-PACIFIC_samples-provenance_20220131d.tsv?download=1", header = TRUE, data.table = FALSE, skip = 1)

# Modify map data ----

# Tara Pacific data
mapdat_tpac = metadat_tpac %>%
  rename(biosample = `sample-id_biosamples`) %>%
  filter(biosample %in% tpac_metagenomes$biosample) %>%
  transmute(longitude = as.numeric(`sampling-event_longitude_start_ddd.dddddd`),
            latitude = as.numeric(`sampling-event_latitude_start_dd.dddddd`),
            sample_label = substr(`sampling-design_label`, 1, 9)) %>%
  group_by(sample_label) %>%
  summarise(longitude = mean(longitude, na.rm = TRUE),
            latitude  = mean(latitude,  na.rm = TRUE),
            .groups = "drop") %>%
  mutate(dataset = "Tara Pacific",
         label = substr(sample_label, nchar(sample_label) - 2, nchar(sample_label)))

# External studies data
mapdat_ext <- ext_metagenomes %>%
  filter(biome_group != "NA") %>%
  mutate(biome_group = ifelse(grepl("sponges", biome_group), "sponge", biome_group),
         longitude = as.numeric(lon),
         latitude = as.numeric(lat),
         sample_label = sample,
         dataset = case_when(grepl("sponge", biome_group) ~ "External sponge",
                             grepl("coral", biome_group) ~ "External coral",
                             TRUE ~ NA_character_)) %>%
  select(sample_label, dataset, longitude, latitude)

# Finalise the data ----

sample_marker = mapdat_ext %>%
  bind_rows(mapdat_tpac) %>%
  mutate(text = paste0("<b>Sample label:</b> ", sample_label, "<br>"),
         # shift eastern longitudes to negative Pacific view
         longitude = ifelse(longitude > -10, longitude - 360, longitude),
         # assign colours based on dataset
         colour = case_when(dataset == "External sponge" ~ "#9467BD",
                            dataset == "External coral" ~ "#E377C2",
                            dataset == "Tara Pacific" ~ "#E9540D"))         

# Produce map ----

leaflet_map <- leaflet(data = sample_marker,
                       options = leafletOptions(zoomControl = TRUE, dragging = TRUE)) %>%
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addMapPane("tara", zIndex = 615) %>%
  addCircleMarkers(lng = ~longitude,
                   lat = ~latitude,
                   radius = 6,
                   stroke = FALSE,
                   fillColor = ~colour,
                   fillOpacity = 0.9,
                   popup = ~text) %>%
  setMaxBounds(lng1 = -360, lat1 = -90, lng2 = 0, lat2 = 90) %>%
  setView(lng = -180, lat = 20, zoom = 2)

# Save map as html file
wrapper_save_html(leaflet_map, "extended-data-figure-1-a.html")
