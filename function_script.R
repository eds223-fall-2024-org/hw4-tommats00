
library(here)
library(terra)
library(sf)
library(tidyverse)
library(tmap)



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
  
  # Checking that the CRSs reprojected to EPSG:4326
  if (st_crs(eez) == st_crs(depth) && st_crs(depth) == st_crs(sst)){
    print("All CRSs are the same")
  } else {
    warning("CRSs are NOT the same")
  }
  
  # Find mean SST from 2008-2012
  mean_sst <- app(sst, fun = mean, na.rm = TRUE)
  
  # Convert average SST from Kelvin to Celsius
  sst_c <- mean_sst - 273.15
  
  # Crop depth raster to match extent of SST raster
  depth_crop <- crop(depth, sst_c)
  
  # The resolutions of the SST and depth data do not match
  # Resample the depth data to match the resolution of the SST data using the nearest neighbor approach
  depth_resample <- resample(depth_crop, sst_c, method = "near")
  
  
  # Check that depth and SST match in resolution, extent and CRS
  sst_depth <- c(sst_c, depth_resample)
  
  
  
  # Create reclassification matrix for SST
  sst_rcl <- matrix(c(-Inf, min_sst, 0,
                      min_sst, max_sst, 1,
                      max_sst, Inf, 0),
                    ncol = 3, byrow = TRUE)
  
  # Use reclassification matrix to reclassify sst raster
  sst_reclassify <- classify(sst_c, rcl = sst_rcl)

  # Create classification matrix for Depth
  depth_rcl <- matrix(c(-Inf, -max_depth, 0,
                        -max_depth, -min_depth, 1,
                        -min_depth, Inf, 0),
                      ncol = 3, byrow = TRUE)
  
  # Reclassify Depth raster
  depth_reclassify <- classify(depth_resample, rcl = depth_rcl)
  
  # Find areas that satistfy both SST and Depth
  # Use lapp() to multiple values of rasters together
  sst_depth_condition <- lapp(c(sst_reclassify, depth_reclassify), fun = "*")
  
  # Assign raster values of 0 to NA for cellSize calculation
  sst_depth_NA <- sst_depth_condition
  sst_depth_NA[sst_depth_NA == 0] <- NA
  
  # Find suitable area
  suitable_area <- cellSize(sst_depth_NA,
                            mask = TRUE,
                            unit = "km")
  
  
  # Rasterize the eez data
  eez_raster <- rasterize(eez, suitable_area, "rgn")
  
  # Sum suitable areas for each eez 
  eez_suitable <- zonal(x = suitable_area, 
                        z = eez_raster, 
                        fun = "sum", 
                        na.rm = TRUE)
  
  
  # Create a table of suitable area for the species in each eez 
  suitable_area_table <- kableExtra::kable(eez_suitable, col.names = c("Region", "Suitable Area (km\u00B2)"), align = "c")
  
  # Join eez_suitable data with the original eez data to map
  eez_join <- left_join(eez, eez_suitable, by = "rgn")
  
  
  # Map of suitable areas 
  map <- tm_shape(depth_crop) +
    tm_raster(palette = paletteer_c("ggthemes::Green-Blue Diverging", 7, direction = -1),
              title = "Depth (m)",
              midpoint = NA) +
    tm_shape(eez_join) +
    tm_polygons(col = "area",
                palette = paletteer_c("ggthemes::Classic Area-Brown", 5),
                title = "Suitable area (km\u00B2)",
                alpha = 0.9) +
    tm_text("rgn", 
            size = 0.45,
            col = "black") + 
    tm_layout(frame = FALSE,
              legend.outside = TRUE,
              legend.text.size = 0.65,
              main.title = paste("Suitable aquaculture area (km\u00B2) for", species, "along \nthe West Coast Exclusive Economic Zone (EEZ)"),
              main.title.size = 1.00) +
    tm_compass(position = c(0, 0.125),
               size = 1.5) +
    tm_scale_bar(position = c(0.05, 0.025)) 
  
  # Print both outputs, Kable table and Map 
  print(suitable_area_table)
  
  print(map)
}
