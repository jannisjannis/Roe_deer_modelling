#install.packages("biomod2")
#install.packages("readr")
#install.packages("terra")
#install.packages("spThin")
#install.packages("tidyterra")
#install.packages("ggtext")
setwd("C:/Users/Jannis/OneDrive - Scientific Network South Tyrol/Documents/Master - EMMA/3. Semester/Southtyrol-hunting data")

library(biomod2)
library(readr)
library(terra)
library(dplyr)
library(spThin)
library(ggplot2)
library(tidyterra)
library(ggtext)

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

env_stack <- c(dem, aspect, slope, bio1, bio2, bio11, bio19, forest_cover) #grassland_cover)

#### read and thin presence data ####
presence_data <- read.csv("Rehwilddaten_punktverortet.csv", header = TRUE, sep = ",")

#change from character to numeric the coordinates
presence_data$lat <- as.numeric(gsub(",", ".", presence_data$GPS..lat.))
presence_data$lon <- as.numeric(gsub(",", ".", presence_data$GPS..lng.))

presence_data_biomod <- presence_data[, c("lon", "lat")]
presence_data_biomod$presence <- 1
presence_data_biomod$species <- "Roe_deer"

colnames(presence_data_biomod)[colnames(presence_data_biomod) == "lon"] <- "x"
colnames(presence_data_biomod)[colnames(presence_data_biomod) == "lat"] <- "y"

#Thin Dataset down from the 46.600 occurences
# Define the extent based on the range of your presence data
xmin <- min(presence_data_biomod$x)
xmax <- max(presence_data_biomod$x)
ymin <- min(presence_data_biomod$y)
ymax <- max(presence_data_biomod$y)
extent_obj <- ext(xmin, xmax, ymin, ymax)

# Create a raster grid at the desired resolution
r <- rast(extent_obj, resolution = 0.01)  # 0.01 degrees (~1 km)

coords <- as.matrix(presence_data_biomod[, c("x", "y")])
grid_cells <- cellFromXY(r, coords)

# Select unique cells and corresponding points
unique_cells <- unique(grid_cells)
occurence_data_thinned <- presence_data_biomod[match(unique_cells, grid_cells), ]

#save thinned dataframe
write.csv(occurence_data_thinned, "thinned_occurence_data.csv", row.names = FALSE)

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
