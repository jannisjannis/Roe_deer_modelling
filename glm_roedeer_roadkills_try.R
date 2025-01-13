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
roadKillFilePath = "GIS/roadKills_presence.tif"     # response variable
## predictor variables
  # broad-leaved forest --> band 1
BLforestFilePath = "forGLM/landcover/Results/FTY_2018_010m_03035_V1_0/FTY_2018_010m_03035_V1_0/FTY_2018_010m_03035_V1_0.tif"
  # coniferous forest --> band 2
CforestFilePath = "forGLM/landcover/Results/FTY_2018_010m_03035_V1_0/FTY_2018_010m_03035_V1_0/FTY_2018_010m_03035_V1_0.tif"
  # mixed forest --> band 313
MixforestFilePath = "forGLM/landcover/Results/U2018_CLC2018_V2020_20u1/U2018_CLC2018_V2020_20u1/U2018_CLC2018_V2020_20u1.tif"
transitionalLandFilePath = "forGLM/rastersized_data/smallWoodyFeatures_EPSG25832.tif"
  #pastures --> band 231
pasturesFilePath = "forGLM/landcover/Results/U2018_CLC2018_V2020_20u1/U2018_CLC2018_V2020_20u1/U2018_CLC2018_V2020_20u1.tif"
grasslandsFilePath = "forGLM/landcover/Results/GRA_2018_010m_03035_V1_0/GRA_2018_010m_03035_V1_0/GRA_2018_010m_03035_V1_0.tif"
# sdmFilePath = "/forGLM/
RwaterFilePath = "forGLM/Fliessgewaesser/DownloadService/Watercourses_line.shp"
lakesFilePath = "forGLM/Seen/DownloadService/Lakes_polygon.shp"
roadTypeFilePath = "forGLM/relevant_roads_raster_10by10.tif"
humansFilePath = "forGLM/PopDensity/popDens_clippedTo_villages.shp"

# # # Load all the files # # #
# ----- vector -----
runningWaterVect <- vect(RwaterFilePath)
lakesVect <- vect(lakesFilePath)
humanInfluenceVect <- vect(humansFilePath)

# ----- raster ----- # subds allows to load specific band
roadKillsRaster <- raster(roadKillFilePath)
broadleaveRaster <- raster(BLforestFilePath, subds = 1)
coniferousRaster <- raster(CforestFilePath, subds = 2)
mixedForestRaster <- raster(MixforestFilePath, subds = 313)
transitionalLandRaster <- raster(transitionalLandFilePath)
pasturesRaster <- raster(pasturesFilePath, subds = 231)
grasslandsRaster <- raster(grasslandsFilePath)
roadTypeRaster <- vect(roadTypeFilePath)
#sdmRaster <- raster(sdmFilePath)

# get an overview over raster properties - function crs(raster) gives 
# coordinate reference system, function res(r) gives resolution
print(grasslandsRaster)




