#install.packages("biomod2")
#install.packages("readr")
#install.packages("terra")
#install.packages("spThin")
#install.packages("tidyterra")
#install.packages("ggtext")
install.packages("rio")
setwd("C:/Users/Jannis/OneDrive - Scientific Network South Tyrol/Documents/Master - EMMA/3. Semester/Southtyrol-hunting data")

library(biomod2)
library(readr)
library(terra)
library(dplyr)
library(spThin)
library(ggplot2)
library(tidyterra)
library(ggtext)
library(rio) 

#skip to line 121 from here

#### Worlclim raster layer creation ####
#create precipitation sum per year raster for southtyrol
rastersprecip <- rast(paste0("precip_raster/precip_clip_", 1:12, ".tif"))
annual_precip <- sum(rastersprecip)
#writeRaster(annual_precip, "annual_precip_southtyrol.tif", overwrite = TRUE)

#create precipitation sum WINTER raster for southtyrol
rastersprecip_WINTER <- rast(paste0("precip_raster/precip_clip_", c(12, 1, 2), ".tif"))
winter_precip <- sum(rastersprecip_WINTER)
#writeRaster(winter_precip, "precip_raster/winter_precip_southtyrol.tif", overwrite = TRUE)

#create mean annual temp raster for southtyrol
rasters <- rast(paste0("climate_raster/temp_clip_mean_", 1:12, ".tif"))
mean_annual_temp <- sum(rasters) / 12
#writeRaster(mean_annual_temp, "mean_annual_temp_southtyrol.tif", overwrite = TRUE)

#create mean WINTER (Dez, jan, feb) temp raster for southtyrol
rasters_WINTER <- rast(paste0("climate_raster/temp_clip_mean_", c(12, 1, 2), ".tif"))
mean_WINTER_temp <- sum(rasters_WINTER) / 3
#writeRaster(mean_WINTER_temp, "climate_raster/mean_winter_temp_southtyrol.tif", overwrite = TRUE)

#create Tmax raster for southtyrol
rasters_tmax <- rast(paste0("Tmax_raster/Tmax_mean_clipped_", 1:12, ".tif"))
Tmax_annual_temp <- sum(rasters_tmax) / 12
#writeRaster(Tmax_annual_temp, "Tmax_raster/Tmax_annual_temp_southtyrol.tif", overwrite = TRUE)

#create Tmin raster for southtyrol
rasters_tmin <- rast(paste0("Tmin_raster/Tmin_mean_clipped_", 1:12, ".tif"))
Tmin_annual_temp <- sum(rasters_tmin) / 12
#writeRaster(Tmin_annual_temp, "Tmin_raster/Tmin_annual_temp_southtyrol.tif", overwrite = TRUE)

#### set CRS ####
dem <- rast("dem_alps_100m.tif")
aspect <- rast("Layer/aspect.tif")
slope <- rast("Layer/slope.tif")
bio1 <-rast("Layer/bio1.tif")
bio2 <-rast("Layer/bio2.tif")
bio11 <-rast("Layer/bio11.tif")
bio12 <-rast("Layer/bio12.tif")
bio19 <-rast("Layer/bio19.tif")
forest_cover <- rast("Layer/FTY_2018_010m_03035_V1_0.tif")
grassland_cover <- rast("Layer/GRA_2018_010m_03035_V1_0.tif")

#set EPSG:25832 as target crs into which we want to reproject all other layers
target_crs <- crs(dem)
#dem, aspect and slope are already in EPSG:25832
bio1 <- project(bio1, target_crs)
bio2 <- project(bio2, target_crs)
bio11 <- project(bio11, target_crs)
bio12 <- project(bio12, target_crs)
bio19 <- project(bio19, target_crs)
forest_cover <- project(forest_cover, target_crs)
grassland_cover <- project(grassland_cover, target_crs)

#### Match resolutions of the layers ####
#we start with downscaling the climate variables from 1km2 to 100m2
target_res <- 100
template_raster <- rast(dem)  # Copy extent and CRS from the reference raster (dem)
res(template_raster) <- target_res  # Ensure resolution is set to 100x100 meters

# Step 3: Resample bio1 to match the template raster
bio1_100m <- resample(bio1, template_raster, method = "bilinear")
bio2_100m <- resample(bio2, template_raster, method = "bilinear")
bio11_100m <- resample(bio11, template_raster, method = "bilinear")
bio12_100m <- resample(bio12, template_raster, method = "bilinear")
bio19_100m <- resample(bio19, template_raster, method = "bilinear")

forest_cover_100m <- aggregate(forest_cover, fact = 10, fun = mean) #factor 10 due to grid size before is 10m
forest_cover_100m_aligned <- resample(forest_cover_100m, template_raster, method = "bilinear")

grassland_cover_100m <- aggregate(grassland_cover, fact = 10, fun = mean) #factor 10 due to grid size before is 10m
grassland_cover_100m_aligned <- resample(grassland_cover_100m, template_raster, method = "bilinear")

#### Cropping all layers to the same extent ####
#start with setting the vector layer of south tyrols border as extent we want to crop the other layers to
border_southtyrol <- vect("Layer/border_southTyrol_withoutNP.shp")

dem_crop <- mask(crop (dem, border_southtyrol), border_southtyrol)
aspect_crop <- mask(crop(aspect, border_southtyrol), border_southtyrol)
slope_crop <- mask(crop(slope, border_southtyrol), border_southtyrol)
bio1_crop <- mask(crop(bio1_100m, border_southtyrol), border_southtyrol)
bio2_crop <- mask(crop(bio2_100m, border_southtyrol), border_southtyrol)
bio11_crop <- mask(crop(bio11_100m, border_southtyrol), border_southtyrol)
bio12_crop <- mask(crop(bio12_100m, border_southtyrol), border_southtyrol)
bio19_crop <- mask(crop(bio19_100m, border_southtyrol), border_southtyrol)
forest_cover_crop <- mask(crop(forest_cover_100m_aligned, border_southtyrol), border_southtyrol)
grassland_cover_crop <- mask(crop(grassland_cover_100m_aligned, border_southtyrol), border_southtyrol)

#### Save aligned rasters ####
writeRaster(dem_crop, "aligned_rasters/dem_100m.tif", overwrite = TRUE)
writeRaster(aspect_crop, "aligned_rasters/aspect_100m.tif", overwrite = TRUE)
writeRaster(slope_crop, "aligned_rasters/slope_100m.tif", overwrite = TRUE)
writeRaster(bio1_crop, "aligned_rasters/bio1_100m.tif", overwrite = TRUE)
writeRaster(bio2_crop, "aligned_rasters/bio2_100m.tif", overwrite = TRUE)
writeRaster(bio11_crop, "aligned_rasters/bio11_100m.tif", overwrite = TRUE)
writeRaster(bio12_crop, "aligned_rasters/bio12_100m.tif", overwrite = TRUE)
writeRaster(bio19_crop, "aligned_rasters/bio19_100m.tif", overwrite = TRUE)
writeRaster(forest_cover_crop, "aligned_rasters/forest_cover_100m.tif", overwrite = TRUE)
writeRaster(grassland_cover_crop, "aligned_rasters/grassland_cover_100m.tif", overwrite = TRUE)

#### Read in Aligned rasters ####
dem <- rast("aligned_rasters/dem_100m.tif")
aspect <- rast("aligned_rasters/aspect_100m.tif")
slope <- rast("aligned_rasters/slope_100m.tif")
bio1 <-rast("aligned_rasters/bio1_100m.tif")
bio2 <-rast("aligned_rasters/bio2_100m.tif")
bio11 <-rast("aligned_rasters/bio11_100m.tif")
bio12 <-rast("aligned_rasters/bio12_100m.tif")
bio19 <-rast("aligned_rasters/bio19_100m.tif")
forest_cover <- rast("aligned_rasters/forest_cover_100m.tif")
grassland_cover <- rast("aligned_rasters/grassland_cover_100m.tif")

env_stack <- c(dem, aspect, slope, bio1, bio2, bio11, bio19, forest_cover, grassland_cover)

#PRÜFEN OB DIE VARIAblen abhängig voneinander sind, dann braucht man sie nciht?!

#### read and thin presence data ####
presence_data <- import("Rehe_Unfall_und_Abschuss.csv", header = TRUE)

colnames(presence_data)[colnames(presence_data) == "GPS (lat)"] <- "lat"
colnames(presence_data)[colnames(presence_data) == "GPS (lng)"] <- "lon"
colnames(presence_data)

#change from character to numeric the coordinates
presence_data$lat <- as.numeric(gsub(",", ".", presence_data$lat))
presence_data$lon <- as.numeric(gsub(",", ".", presence_data$lon))

#removing Nas
presence_data <- presence_data %>%
  filter(!is.na(lat), !is.na(lon))

#Change from old WGS 84 to EPSG:35832
#create  Vector from the old dataframe
presence_vect <- vect(presence_data, geom = c("lon", "lat"), crs = "EPSG:4326")  
presence_vect_25832 <- project(presence_vect, "EPSG:25832") # Transform to EPSG:25832
coords <- crds(presence_vect_25832) # extract transformed coordinates
# Add the transformed coordinates back to the data frame
presence_data$X_25832 <- coords[, 1]  # X coordinate
presence_data$y_25832 <- coords[, 2]  # Y coordinate

#Remove all the points that are outside our southtyrol boundary
border_southtyrol <- vect("Layer/border_southTyrol_withoutNP.shp")
presence_vect_25832 <- vect(presence_data, geom = c("X_25832", "y_25832"), crs = "EPSG:25832")  
presence_25832_within_southtyrol <- presence_vect_25832[border_southtyrol, ]
#crs(presence_25832_within_southtyrol)
clipped_data <- as.data.frame(presence_25832_within_southtyrol)
#add new coords back to the dataframe
# Ensure the new coordinates are included in the data frame
clipped_data$X_25832 <- crds(presence_25832_within_southtyrol)[, 1]  # X coordinate (EPSG:25832)
clipped_data$y_25832 <- crds(presence_25832_within_southtyrol)[, 2]  # Y coordinate (EPSG:25832)
export(clipped_data, "Roedeer_within_ST_25832.csv")

#create table for biomod2 with just x,y coordinates, the presence indicator and the species name
presence_data_biomod <- clipped_data[, c("y_25832", "X_25832")]
presence_data_biomod$presence <- 1
presence_data_biomod$species <- "Roe_deer"
colnames(presence_data_biomod)[colnames(presence_data_biomod) == "y_25832"] <- "y"
colnames(presence_data_biomod)[colnames(presence_data_biomod) == "X_25832"] <- "x"

#Thin Dataset down from the 52.900 occurences
r <- rast(ext(border_southtyrol), resolution = 1000, crs = "EPSG:25832") #Set the extent and resolution for the raster grid
r_points <- rasterize(presence_25832_within_southtyrol, r, fun = "first", background = NA) #Rasterize the points (assign each point to a grid cell)
unique_points <- as.points(r_points, na.rm = TRUE) #Extract unique points based on the raster cells
writeVector(unique_points, "thinned_occurence_data_1km.shp", overwrite = TRUE)

#### Format data and generate pseudo-absences ####
Roedeer_data <- BIOMOD_FormatingData(
  resp.var = occurence_data_thinned$presence,           # Presence data
  expl.var = env_stack,                        # Environmental variables
  resp.xy = occurence_data_thinned[, c("x", "y")],      # Coordinates of presences
  resp.name = "species",                      # Name of the species
  PA.nb.rep = 2,                               # Number of pseudo-absence replicates
  PA.nb.absences = 5000,                       # Number of pseudo-absences
  PA.strategy = "random"                       # Strategy for generating pseudo-absences
)

Roedeer_data
