
library(here)
library(terra)
library(sf)
library(tidyverse)
library(tmap)


## What I need in my function ## 
# species, max & min depth, max & min sst for suitable aquaculture zone
# 


aqua_fun <- function(species, max_sst, min_sst, max_depth, min_depth ){
  
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
  
  # Find mean SST from 2008-2012
  mean_sst <- app(sst, fun = mean, na.rm = TRUE)
  
  # Convert average SST from Kelvin to Celsius
  sst_c <- mean_sst - 273.15
  
  # Crop depth raster to match extent of SST raster
  depth_crop <- crop(depth, sst_c)
  
  # The resolutions of the SST and depth data do not match
  # Resample the depth data to match the resolution of the SST data using the nearest neighbor approach
  depth_resample <- resample(depth_crop, sst_c, method = "near")
  
  # Set CRS the same
  # depth_resample <- project(depth_resample, crs(sst))
  
  # Check that depth and SST match in resolution, extent and CRS
  # sst_depth <- c(sst_c, depth_resample)
  
  
  
  # Create reclassification matrix for SST
  sst_rcl <- matrix(c(-Inf, min_sst, 0,
                      min_sst, max_sst, 1,
                      max_sst, Inf, 0),
                    ncol = 3, byrow = TRUE)
  
  # Use reclassification matrix to reclassify sst raster
  sst_reclassify <- classify(sst_c, rcl = sst_rcl)
  
  # Check preliminary map
  # plot(sst_reclassify)
  
  # Create classification matrix for Depth
  depth_rcl <- matrix(c(-Inf, -max_depth, 0,
                        -max_depth, -min_depth, 1,
                        -min_depth, Inf, 0),
                      ncol = 3, byrow = TRUE)
  
  # Reclassify Depth raster
  depth_reclassify <- classify(depth_resample, rcl = depth_rcl)
  
  # Check preliminary map
  # plot(depth_reclassify)
  
  # Find areas that satistfy both SST and Depth
  
  # Use lapp() to multiple values of rasters together
  sst_depth_condition <- lapp(c(sst_reclassify, depth_reclassify), fun = "*")
  
  # sst_depth_condition <- sst_reclassify * depth_reclassify
  
  # Plot both variables together to find suitable areas
  # plot(sst_depth_condition)
  
  crs(sst_depth_condition) == crs(eez)
  
  crs(sst_depth_condition) <- crs(eez)
  
  crs(sst_depth_condition) == crs(eez)
  
  sst_depth_NA <- sst_depth_condition
  sst_depth_NA[sst_depth_NA == 0] <- NA
  
  # Find suitable area
  suitable_area <- cellSize(sst_depth_NA,
                            mask = TRUE,
                            unit = "km")
  
  plot(suitable_area)
  
  # Rasterize
  eez_raster <- rasterize(eez, suitable_area, "rgn")
  #plot(eez_raster)
  
  eez_suitable <- zonal(x = suitable_area, 
                        z = eez_raster, 
                        fun = "sum", 
                        na.rm = TRUE)
  
  
  
  print(paste("Suitable area for", species, "in", eez_suitable$rgn[1], "is", eez_suitable$area[1]))
  print(paste("Suitable area for", species, "in", eez_suitable$rgn[2], "is", eez_suitable$area[2]))
  print(paste("Suitable area for", species, "in", eez_suitable$rgn[3], "is", eez_suitable$area[3]))
  print(paste("Suitable area for", species, "in", eez_suitable$rgn[4], "is", eez_suitable$area[4]))
  print(paste("Suitable area for", species, "in", eez_suitable$rgn[5], "is", eez_suitable$area[5]))
  
  # Plot
  tm_shape(depth_crop) +
    tm_raster(legend.show = FALSE) +
    tm_shape(sst_depth_condition) +
    tm_raster(legend.show = FALSE) +
    tm_shape(eez_raster) +
    tm_raster(palette = "viridis",
              title = "Region") + 
    tm_shape(sst_depth_condition) +
    tm_raster(alpha = 0.7,
              legend.show = FALSE) +
    tm_layout(frame = FALSE,
              legend.position = c(0.05,0.05)) +
    tm_compass(position = c(0, 0.85)) +
    tm_scale_bar(position = c(0.615, 0.90)) 
}
