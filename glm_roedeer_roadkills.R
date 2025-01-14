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
BLforestFilePath = "forGLM/landcover/Results/FTY_2018_010m_03035_V1_0/FTY_2018_010m_03035_V1_0/FTY_2018_010m_03035_V1_0_ST.tif"
  # coniferous forest --> band 2
CforestFilePath = "forGLM/landcover/Results/FTY_2018_010m_03035_V1_0/FTY_2018_010m_03035_V1_0/FTY_2018_010m_03035_V1_0_ST.tif"
  # mixed forest --> band 313
MixforestFilePath = "forGLM/landcover/U2018_CLC2018_extent_ST.tif"
transitionalLandFilePath = "forGLM/rastersized_data/smallWoodyFeatures_EPSG25832.tif"
  #pastures --> band 231
pasturesFilePath = "forGLM/landcover/U2018_CLC2018_extent_ST.tif"
grasslandsFilePath = "forGLM/landcover/Results/GRA_2018_010m_03035_V1_0/GRA_2018_010m_03035_V1_0/GRA_2018_010m_03035_V1_0.tif"
# !!!!!! sdmFilePath = "/forGLM/
watercoursesFilePath = "variables_glm_clipped_ST_minus_NP/watercourses.tif"
lakesFilePath = "variables_glm_clipped_ST_minus_NP/lakes.tif"
roadTypeFilePath = "variables_glm_clipped_ST_minus_NP/relevant_roads.tif"
humansFilePath = "variables_glm_clipped_ST_minus_NP/human_influence.tif"
# !!!!!!!! roadnetworkFilePath =

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

# list of all variables             !!!!!!! ad sdm !!!!!!!
all_variables <- list(roadKills, roadKillsCount, broadleave, coniferous,
                      mixedForest, transitionalLand, pastures, grasslands,
                      roadType, humanInfluence, watercourses, lakes)

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
#broadleave_crop <- mask(crop (broadleave, ST_border), ST_border)
#coniferous_crop <- mask(crop (coniferous, ST_border), ST_border)
#mixedForest_crop <- mask(crop (mixedForest, ST_border), ST_border)
transitionalLand_crop <- mask(crop (transitionalLand, ST_border), ST_border)
#pastures_crop <- mask(crop (pastures, ST_border), ST_border)
#grasslands_crop <- mask(crop (grasslands, ST_border), ST_border)
roadType_crop <- mask(crop (roadType, ST_border), ST_border)
#roadNetwork_crop <- mask(crop (roadNetwork, ST_border), ST_border)
humanInfluence_crop <- mask(crop (humanInfluence, ST_border), ST_border)
watercourses_crop <- mask(crop (watercourses, ST_border), ST_border)
lakes_crop <- mask(crop (lakes, ST_border), ST_border)
#sdm_crop <- mask(crop (sdm_crop, ST_border), ST_border)

# ----------------------------------------------------------------------------
# Visualization of the data
# ----------------------------------------------------------------------------

# Activate the libraries needed for this 
library(car) # needed for the Anova() function, contains function vif(), which offers an easy alternative way to check predictor independence.
library(ggplot2) # needed for some graphs

# Assuming roadnetwork is already loaded as a raster
# Convert roadnetwork raster to "yes" or "no"
roadTypes_values <- getValues(roadType_crop)
roadTypesfactor <- ifelse(roadTypes_values == 1, "Yes", "No")

roadKills_values <- getValues(roadKills)
roadKills_factor <- ifelse(roadKills_values == 1, "Yes", "No")

# Combine both datasets into a data frame
df <- data.frame(road = roadnetwork_factor, roadkill = roadkill_factor)

# Ensure that the number of roadkill and road network values match
df <- na.omit(df)  # Remove any NA values

# Check the resulting data frame
head(df)

# Plot a boxplot comparing road presence and roadkill occurrence
ggplot(df, aes(x = roadnetwork_factor, y = roadkill_factor, fill = roadnetwork_factor)) +
  geom_boxplot() +
  labs(title = "Road Network vs. Roadkill Occurrence",
       x = "Road Presence",
       y = "Roadkill Occurrence") +
  theme_minimal() +
  scale_y_discrete(labels = c("No", "Yes"))

