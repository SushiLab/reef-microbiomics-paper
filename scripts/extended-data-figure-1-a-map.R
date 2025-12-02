## Project: Reef Microbiomics Database
## Script purpose: Code used to generate figure 1 - map
## Author: Fabienne Wiederkehr
## Date: 30/08/2023

# Load packages
require(data.table)
require(tidyverse)
require(leaflet)
require(googlesheets4)

# Load functions
wrapper_save_html<-function(html_widget,path,self=TRUE){
  current=getwd()
  dir=dirname(path)
  file=basename(path)
  setwd(dir)
  htmlwidgets::saveWidget(html_widget,file=file,selfcontained=self)
  setwd(current)
}

# Load data
tpac_biosample_list<-fread("../../../data/processed/genome-summaries/tpac-biosample-list.csv",header=FALSE,data.table=FALSE)

# Tara Pacific metadata
metadat_tpac<-fread("https://zenodo.org/record/6299409/files/TARA-PACIFIC_samples-provenance_20220131d.tsv?download=1",header=TRUE,data.table=FALSE,skip=1) %>%
  rename(biosample=`sample-id_biosamples`) %>%
  filter(biosample %in% tpac_biosample_list$V1) %>%
  select(longitude=`sampling-event_longitude_start_ddd.dddddd`,
         latitude=`sampling-event_latitude_start_dd.dddddd`,
         sample_label=`sampling-design_label`) %>%
  mutate(longitude=as.numeric(longitude),
         latitude=as.numeric(latitude),
         sample_label=substr(sample_label,0,9)) %>%
  group_by(sample_label) %>%
  summarise_at(vars(longitude:latitude),mean) %>%
  ungroup() %>%
  mutate(dataset="Tara Pacific") %>%
  select(sample_label,dataset,longitude,latitude) %>%
  mutate(label=substr(sample_label,nchar(sample_label)-2,nchar(sample_label)))

# External studies metadata
metadat_metag<-read_sheet("https://docs.google.com/spreadsheets/d/14sEHDVQ7gqqxKCHvC1kCouP-6WBcirEzNNCsZI9K9XE/edit#gid=0",sheet="Samples for metagenomic projects") %>%
  filter(biome_group!="NA") %>%
  mutate(biome_details=ifelse(grepl("sponges",biome_details),"sponge",biome_details)) %>%
  mutate(longitude=as.numeric(lon),
         latitude=as.numeric(lat),
         sample_label=sample) %>%
  mutate(dataset=ifelse(grepl("sponge",biome_details),"External sponge",NA)) %>%
  mutate(dataset=ifelse(grepl("coral",biome_details),"External coral",dataset)) %>%
  select(sample_label,dataset,longitude,latitude)

# Finalise the data
sample_marker<-metadat_metag %>%
  bind_rows(metadat_tpac) %>%
  mutate(text=paste0("<b>Sample label:</b> ",sample_label,"<br>")) %>%
  mutate(longitude=ifelse(longitude>-10,longitude-360,longitude)) 

# Set icons  
custom_icons<-icons(
  iconUrl=ifelse(sample_marker$dataset=="External sponge","../../../data/processed/figures/figure-1/figure-1_circle_ext_sponge.png",
                 ifelse(sample_marker$dataset=="External coral",
                 "../../../data/processed/figures/figure-1/figure-1_circle_ext_coral.png",
                 "../../../data/processed/figures/figure-1/figure-1_circle_tpac.png")),
  iconWidth=30,iconHeight=30)

# Generate map
leaflet_map<-leaflet(data=sample_marker,options=leafletOptions(zoomControl=TRUE,dragging=TRUE)) %>%
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addMapPane("tara",zIndex=615) %>%
  addMarkers(~longitude,~latitude,icon=custom_icons) %>%
  setMaxBounds(lng1=-360,lat1=-90,lng2=0,lat2=90) %>%
  setView(lng=-180,lat=20,zoom=2) 
wrapper_save_html(leaflet_map,"../../../data/processed/figures/figure-1/figure-1_map_ext.html")
