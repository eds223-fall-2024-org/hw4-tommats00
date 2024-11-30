
library(here)
library(terra)
library(sf)
library(tidyverse)
library(tmap)


## What I need in my function ## 
# species, max & min depth, max & min sst for suitable aquaculture zone
# 


aqua_fun <- function(species, ){
  
  # Load eez data
  eez <- read_sf(here("data", "wc_regions_clean.shp"))
  
  # Load depth data
  depth <- rast(here("data", "depth.tif"))
  
  # Load sst data, use file path
  files <- list.files(
    here("data"), pattern = "*sst", 
    full.names = TRUE)
  
  sst <- c(rast(files))
  
  # Set CRS all the same
  eez <- st_transform(eez, crs("EPSG:4326"))
  depth <- project(depth, crs("EPSG:4326"))
  sst <- project(sst, crs("EPSG:4326"))
  
 
}