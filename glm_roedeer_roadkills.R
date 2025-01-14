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
humanInfluence <- rast(humansFilePath)
watercourses <- rast(watercoursesFilePath)
lakes<- rast(lakesFilePath)
#sdmRaster <- raster(sdmFilePath)

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

# Loop through all rasters and resample to the target resolution
for (i in seq_along(all_variables)) {
  obj <- all_variables[[i]]
  
  # Resample to the target resolution (100 meters)
  all_variables[[i]] <- resample(obj, res = target_res, method = "near")  # 'bilinear' for continuous, 'near' for categorical
  cat("Resampled raster", i, "to 100 meters resolution.\n")
}


# ----------------------------------------------------------------------------
#                               (C)  - extent

#start with setting the vector layer of south tyrols border as extent we want to crop the other layers to
ST_border <- vect(ST_borderFilePath)
ST_border <- ST_border(target_crs)

# Loop through all rasters and apply crop and mask
for (i in seq_along(all_variables)) {
  # Get the name and raster of the current object
  var_name <- names(all_variables)[i]
  obj <- all_variables[[i]]
  
  # Crop the raster to the boundary of South Tyrol
  cropped_raster <- crop(obj, ST_border)
  
  # Mask the cropped raster to keep only values within South Tyrol
  masked_raster <- mask(cropped_raster, ST_border)
  
  # Dynamically create a new variable name and assign the processed raster to it
  assign(paste0(var_name, "_crop"), masked_raster)
  
  # Print a message to indicate the process is done for the current raster
  cat("Processed raster", var_name, "to create", paste0(var_name, "_crop"), "\n")
}

# ----------------------------------------------------------------------------
# GLM set up
# ----------------------------------------------------------------------------

