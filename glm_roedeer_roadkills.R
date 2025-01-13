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
install.packages("raster") # attention number of bands
install.packages("terra")
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

# # # Load all the files # # #
 # subds allows to load specific band
roadKills <- raster(roadKillFilePath)
roadKillsCount <- raster(roadKillCountsFilePath)
broadleave <- raster(BLforestFilePath, subds = 1)
coniferous <- raster(CforestFilePath, subds = 2)
mixedForest <- raster(MixforestFilePath, subds = 25)
transitionalLand <- raster(transitionalLandFilePath)
pastures <- raster(pasturesFilePath, subds = 18)
grasslands <- raster(grasslandsFilePath)
roadType <- raster(roadTypeFilePath)
humanInfluence <- raster(humansFilePath)
watercourses <- raster(watercoursesFilePath)
lakes<- raster(lakesFilePath)
#sdmRaster <- raster(sdmFilePath)

# get an overview over raster properties - function crs(raster) gives 
# coordinate reference system, function res(r) gives resolution
print(grasslands)

# ----------------------------------------------------------------------------
# BEFORE checking assumptions and running the GLM
# ----------------------------------------------------------------------------
# Make sure that...
# a. 	all layers have the same resolution (can be done in R) --> 100m by 100m
# b. all layers have the same spatial extent (can be done in R, except if minimum 
#     extent is not given) --> border South Tyrol (not the one from NUTS!)
# c. all layers have the same CRS (can be made sure of in R) --> EPSG 25832




