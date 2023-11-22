# This script plots a world map showing the sampling sites of Tara Pacific

# Load packages -----

require(leaflet)
require(data.table)
require(tidyverse)
require(googlesheets4)

# Load functions -----

wrapper_save_html<-function(html_widget,path,self=TRUE){
  current=getwd()
  dir=dirname(path)
  file=basename(path)
  setwd(dir)
  htmlwidgets::saveWidget(html_widget,file=file,selfcontained=self)
  setwd(current)
}

# Load data -----

metagenomes_tpac<-read_sheet("https://docs.google.com/spreadsheets/d/1ukivQCmnq6PePM5vl6r3PdRJXSwZ_iChEN1LRbGD0EQ/edit#gid=0",sheet="Sheet 1") %>%
  separate(sample,into=c(NA,"biosample",NA),sep="_",remove=FALSE)

metadat_tpac<-fread("https://zenodo.org/record/6299409/files/TARA-PACIFIC_samples-provenance_20220131d.tsv?download=1",header=TRUE,data.table=FALSE,skip=1) %>%
  rename(biosample=`sample-id_biosamples`) %>%
  filter(biosample %in% metagenomes_tpac$biosample) %>%
  select(longitude=`sampling-event_longitude_start_ddd.dddddd`,
         latitude=`sampling-event_latitude_start_dd.dddddd`,
         sample_label=`sampling-design_label`) %>%
  mutate(longitude=as.numeric(longitude),
         latitude=as.numeric(latitude),
         sample_label=substr(sample_label,0,9)) %>%
  group_by(sample_label) %>%
  summarise_at(vars(longitude:latitude),mean) %>%
  ungroup()

# Finalise the data -----

sample_marker<-metadat_tpac %>%
  mutate(groups="Tara Pacific",
         dataset="Tara Pacific") %>%
  mutate(text=paste0("<b>Sample label:</b> ",sample_label,"<br>")) %>%
  mutate(longitude=ifelse(longitude>0,longitude-360,longitude)) %>%
  mutate(label=substr(sample_label,nchar(sample_label)-2,nchar(sample_label)))

# Produce map -----

# Define custom icons
custom_icons<-iconList("Tara Pacific"=makeIcon("https://raw.github.com/SushiLab/reef-microbiomics-paper/main/resources/figure-1_circle_TPac.png",iconWidth=30,iconHeight=30))

# Produce leaflet
leaflet_map<-sample_marker %>%
  leaflet(data=.,options=leafletOptions(zoomControl=TRUE,dragging=TRUE)) %>%
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addMapPane("tara",zIndex=615) %>% 
  addMarkers(data=sample_marker %>% filter(dataset=="Tara Pacific"),~longitude,~latitude,popup=~text,label=~label,
             icon=~custom_icons[dataset],group=~groups,options=pathOptions(pane="tara"),labelOptions=labelOptions(permanent=TRUE,textsize="20px")) %>%
  setMaxBounds(lng1=-360,lat1=-90,lng2=0,lat2=90) %>%
  setView(lng=-180,lat=20,zoom=2)

# Save map as html file
wrapper_save_html(leaflet_map,"figure-1-a-tpac-map.html")
