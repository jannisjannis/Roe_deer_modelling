#-------------------------------------------------------------------#
# Course:
# Topic: Road kill risk Roe deer South Tyrol
# Script for glm
# Author: Carla Behringer
# Date: 01.12.2024
#-------------------------------------------------------------------#

# clear workspace from any previous entry
rm(list = ls()) 

# set working directory
setwd("/Users/carlabehringer/iCloud Drive (Archive)/Documents/Documents – Carlas MacBook Air/dokumente/Uni/Master/second_year/year2_sem1/Com_EnvMan/Projectstudy/data")  

# install + load packages
#install.packages("raster") # attention number of bands
#install.packages("terra")
library(raster)
library(terra)


# Filepaths for data
roadKillFilePath = "GIS/roadKills_presence_res100.tif"     # response variable
roadKillCountsFilePath = "GIS/roadKills_count_raster_res100.tif"
## predictor variables
  # broad-leaved forest --> band 1
BLforestFilePath = "variables_glm_clipped_ST_minus_NP/forest.tif"
  # coniferous forest --> band 2
CforestFilePath = "variables_glm_clipped_ST_minus_NP/forest.tif"
  # mixed forest --> band 313
MixforestFilePath = "variables_glm_clipped_ST_minus_NP/landuse.tif"
transitionalLandFilePath = "forGLM/rastersized_data/smallWoodyFeatures_EPSG25832.tif"
  #pastures --> band 231
pasturesFilePath = "variables_glm_clipped_ST_minus_NP/landuse.tif"
grasslandsFilePath = "variables_glm_clipped_ST_minus_NP/grassland_cover_100m.tif"
# !!!!!! sdmFilePath = "/forGLM/
watercoursesFilePath = "variables_glm_clipped_ST_minus_NP/watercourses.tif"
lakesFilePath = "variables_glm_clipped_ST_minus_NP/lakes.tif"
roadTypeFilePath = "variables_glm_clipped_ST_minus_NP/RoadType.tif"
humansFilePath = "variables_glm_clipped_ST_minus_NP/human_influence.tif"


ST_borderFilePath = "GIS/border_southTyrol_withoutNP.shp"

# # # Load all the files # # #
 # subds allows to load specific band
roadKills <- rast(roadKillFilePath)
roadKillsCount <- rast(roadKillCountsFilePath)
broadleave <- rast(BLforestFilePath, subds = 1)
coniferous <- rast(CforestFilePath, subds = 2)
mixedForest <- rast(MixforestFilePath, subds = 25)
transitionalLand <- rast(transitionalLandFilePath)
pastures <- rast(pasturesFilePath, subds = 18)
grasslands <- rast(grasslandsFilePath)
roadType <- rast(roadTypeFilePath)
#roadNetwork <- rast(roadTypeFilePath)
humanInfluence <- rast(humansFilePath)
watercourses <- rast(watercoursesFilePath)
lakes<- rast(lakesFilePath)
#sdm <- rast(sdmFilePath)

# list of all variables             !!!!!!! add sdm !!!!!!!
all_variables <- list(roadKills, roadKillsCount, broadleave, coniferous,
                      mixedForest, transitionalLand, pastures, grasslands,
                      roadType, humanInfluence, watercourses, lakes)
names(all_variables) <- c("roadKills", "roadKillsCount", "broadleave","coniferous",
                          "mixedForest", "transitionalLand", "pastures", "grasslands",
                          "roadType", "humanInfluence", "watercourses", "lakes")

# get an overview over raster properties - function crs(raster) gives 
# coordinate reference system, function res(r) gives resolution
print(grasslands)

# ----------------------------------------------------------------------------
# BEFORE checking assumptions and running the GLM
# ----------------------------------------------------------------------------
# Make sure that...
# a. all layers have the same CRS (can be made sure of in R) --> EPSG 25832
# b. all layers have the same resolution (can be done in R) --> 100m by 100m
# c. all layers have the same spatial extent (can be done in R, except if minimum 
#     extent is not given) --> border South Tyrol (not the one from NUTS!)

# ----------------------------------------------------------------------------
#                               (A)  - CRS

# Define the target CRS 
target_crs <- "EPSG:25832"

# Reproject all raster variables to the target CRS
for (i in seq_along(all_variables)) {
  obj <- all_variables[[i]]
  
  # Reproject to target CRS
  all_variables[[i]] <- project(obj, target_crs)  # Reproject to target CRS
  cat("Reprojected raster", i, "to the target CRS.\n")
}


# ----------------------------------------------------------------------------
#                               (B)  - Resolution

# Define the target resolution (100 meters)
target_res <- 100  # Target resolution in meters
template_raster <- rast(roadKillFilePath)  # Copy extent and CRS from the reference raster (roadKills)
res(template_raster) <- target_res  # Ensure resolution is set to 100x100 meters

# Loop through all rasters and resample to the target resolution
for (i in seq_along(all_variables)) {
  obj <- all_variables[[i]]
  
  # Resample to the target resolution (100 meters)
  all_variables[[i]] <- resample(x=obj, y = template_raster, method = "near")  # 'bilinear' for continuous, 'near' for categorical
  cat("Resampled raster", i, "to 100 meters resolution.\n")
}


# ----------------------------------------------------------------------------
#                               (C)  - extent

#start with setting the vector layer of south tyrols border as extent we want to crop the other layers to
ST_border <- vect(ST_borderFilePath)
crs(ST_border) <- target_crs

roadKills_crop <- mask(crop (roadKills, ST_border), ST_border)
roadKillsCount_crop <- mask(crop (roadKillsCount, ST_border), ST_border)
broadleave_crop <- mask(crop (broadleave, ST_border), ST_border)
coniferous_crop <- mask(crop (coniferous, ST_border), ST_border)
mixedForest_crop <- mask(crop (mixedForest, ST_border), ST_border)
transitionalLand_crop <- mask(crop (transitionalLand, ST_border), ST_border)
pastures_crop <- mask(crop (pastures, ST_border), ST_border)
grasslands_crop <- mask(crop (grasslands, ST_border), ST_border)
roadType_crop <- mask(crop (roadType, ST_border), ST_border)
humanInfluence_crop <- mask(crop (humanInfluence, ST_border), ST_border)
watercourses_crop <- mask(crop (watercourses, ST_border), ST_border)
lakes_crop <- mask(crop (lakes, ST_border), ST_border)
#sdm_crop <- mask(crop (sdm_crop, ST_border), ST_border)

# ----------------------------------------------------------------------------
# Visualization of the data
# ----------------------------------------------------------------------------

library(ggplot2) # needed for some graphs

# ----------------- (1) Extract values from rasters
# Extract values and raster names

roadKills_values <- values(roadKills_crop)
roadKills_binary <- ifelse(roadKills_values == 1, "Yes", "No") #binary data
roadKillsCount_values <- values(roadKillsCount_crop)
broadleave_values <- values(broadleave_crop)
broadleave_binary <- ifelse(broadleave_values == 1, "Yes", "No") #binary data
coniferous_values <- values(coniferous_crop)
coniferous_binary <- ifelse(coniferous_values == 1, "Yes", "No") #binary data
mixedForest_values <- values(mixedForest_crop)
mixedForest_binary <- ifelse(mixedForest_values == 1, "Yes", "No") #binary data
transitionalLand_values <- values(transitionalLand_crop)
transitionalLand_binary <- ifelse(transitionalLand_values == 1, "Yes", "No") #binary data
pastures_values <- values(pastures_crop)
pastures_binary <- ifelse(pastures_values == 1, "Yes", "No") #binary data
grasslands_values <- values(grasslands_crop)
grasslands_binary <- ifelse(grasslands_values == 1, "Yes", "No") #binary data
roadType_values <- values(roadType_crop)
humanInfluence_values <- values(humanInfluence_crop)
watercourses_values <- values(watercourses_crop)
watercourses_binary <- ifelse(watercourses_values == 1, "Yes", "No") #binary data
lakes_values <- values(lakes_crop)
lakes_binary <- ifelse(lakes_values == 1, "Yes", "No") #binary data

df <- data.frame(roadKills_binary , roadKillsCount_values, broadleave_binary,
                 coniferous_binary,mixedForest_binary,transitionalLand_binary,
                 transitionalLand_binary,pastures_binary, grasslands_binary,
                 roadType_values, humanInfluence_values, watercourses_binary,
                 lakes_binary)
head(df)

# Beispiel für einen Boxplot
boxplot(roadKillsCount_values, horizontal = TRUE)

# # Assuming roadnetwork is already loaded as a raster
# # Convert roadnetwork raster to "yes" or "no"
# roadTypes_values <- values(roadType_crop)
# #roadTypes_factor <- ifelse(roadTypes_values == 1, "Yes", "No")
# 
# roadKills_values <- values(roadKillsCount_crop)
# #roadKills_factor <- ifelse(roadKills_values == 1, "Yes", "No")
# 
# # Combine both datasets into a data frame
# df <- data.frame(road = roadTypes_values, roadkill = roadKills_values)
# 
# # Ensure that the number of roadkill and road network values match
# df <- na.omit(df)  # Remove any NA values
# 
# # Check the resulting data frame
# head(df)
# 
# 
# # Filter out values equal to 0
# filtered_df <- df[df$roadKills_count_raster_res100 > 0, ]
# head(filtered_df)

# Beispiel für einen Boxplot
boxplot(filtered_df$roadKills_count_raster_res100, horizontal = TRUE)

ggplot(filtered_df, aes(x = factor(roadType), y = roadKills_count_raster_res100, fill = factor(roadType))) +
  geom_boxplot() +
  labs(title = "Road Presence (0-8 Scale) vs. Roadkill Occurrence",
       x = "Road Presence (0-8)",
       y = "Roadkill Occurrence") +
  theme_minimal()

# ----------------------------------------------------------------------------
#                                   glm
# ----------------------------------------------------------------------------
# Activate the libraries needed for this 
library(car) # needed for the Anova() function, contains function vif(), which offers an easy alternative way to check predictor independence.

