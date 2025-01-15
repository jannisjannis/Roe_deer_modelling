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

###############################################################################
# Filepaths to export data
###############################################################################
rasterExportPath <- "Export_R/raster"

######
#############       SKIPPING OF CODE - S T A R T                    ##############                                         
######

# ###############################################################################
# # Filepaths for data
# ###############################################################################
# roadKillFilePath = "GIS/roadKills_presence_res100.tif"     # response variable
# roadKillCountsFilePath = "GIS/roadKills_count_raster_res100.tif"
# ## predictor variables
#   # broad-leaved forest --> band 1
# BLforestFilePath = "variables_glm_clipped_ST_minus_NP/forest.tif"
#   # coniferous forest --> band 2
# CforestFilePath = "variables_glm_clipped_ST_minus_NP/forest.tif"
#   # mixed forest --> band 313
# MixforestFilePath = "variables_glm_clipped_ST_minus_NP/landuse.tif"
# transitionalLandFilePath = "forGLM/rastersized_data/smallWoodyFeatures_EPSG25832.tif"
#   #pastures --> band 231
# pasturesFilePath = "variables_glm_clipped_ST_minus_NP/landuse.tif"
# grasslandsFilePath = "variables_glm_clipped_ST_minus_NP/grassland_cover_100m.tif"
# # !!!!!! sdmFilePath = "/forGLM/
# watercoursesFilePath = "variables_glm_clipped_ST_minus_NP/watercourses.tif"
# lakesFilePath = "variables_glm_clipped_ST_minus_NP/lakes.tif"
# roadTypeFilePath = "variables_glm_clipped_ST_minus_NP/RoadType.tif"
# humansFilePath = "variables_glm_clipped_ST_minus_NP/human_influence.tif"
# 
# 
# ST_borderFilePath = "GIS/border_southTyrol_withoutNP.shp"
# 
# # # # # Load all the files # # #
# 
#  # subds allows to load specific band
# roadKills <- rast(roadKillFilePath)
# roadKillsCount <- rast(roadKillCountsFilePath)
# 
# broadleave_unmasked <- rast(BLforestFilePath)
# mask_layer_broadleave <- broadleave_unmasked == 1
# broadleave <- mask(broadleave_unmasked, mask_layer_broadleave)
# # broadleave <- rast(BLforestFilePath, subds = 1)
# 
# coniferous_unmasked <- rast(CforestFilePath)
# mask_layer_coniferous <- coniferous_unmasked == 2
# coniferous <- mask(coniferous_unmasked, mask_layer_coniferous)
# # coniferous <- rast(CforestFilePath, subds = 2)
# 
# mixedForest_unmasked <- rast(MixforestFilePath)
# mask_layer_mixedForest <- mixedForest_unmasked == 25
# mixedForest <- mask(mixedForest_unmasked, mask_layer_mixedForest)
# # mixedForest <- rast(MixforestFilePath, subds = 25)
# 
# transitionalLand <- rast(transitionalLandFilePath)
# 
# pastures_unmasked <- rast(pasturesFilePath)
# mask_layer_pasture <- pastures_unmasked == 18
# pastures <- mask(pastures_unmasked, mask_layer_pasture)
# 
# grasslands <- rast(grasslandsFilePath)
# roadType <- rast(roadTypeFilePath)
# #roadNetwork <- rast(roadTypeFilePath)
# humanInfluence <- rast(humansFilePath)
# watercourses <- rast(watercoursesFilePath)
# lakes<- rast(lakesFilePath)
# #sdm <- rast(sdmFilePath)
# 
# # list of all variables             !!!!!!! add sdm !!!!!!!
# all_variables <- list(roadKills, roadKillsCount, broadleave, coniferous,
#                       mixedForest, transitionalLand, pastures, grasslands,
#                       roadType, humanInfluence, watercourses, lakes)
# names(all_variables) <- c("roadKills", "roadKillsCount", "broadleave","coniferous",
#                           "mixedForest", "transitionalLand", "pastures", "grasslands",
#                           "roadType", "humanInfluence", "watercourses", "lakes")
# 
# # get an overview over raster properties - function crs(raster) gives
# # coordinate reference system, function res(r) gives resolution
# print(grasslands)
# 
# # ----------------------------------------------------------------------------
# # BEFORE checking assumptions and running the GLM
# # ----------------------------------------------------------------------------
# # Make sure that...
# # a. all layers have the same CRS (can be made sure of in R) --> EPSG 25832
# # b. all layers have the same resolution (can be done in R) --> 100m by 100m
# # c. all layers have the same spatial extent (can be done in R, except if minimum
# #     extent is not given) --> border South Tyrol (not the one from NUTS!)
# 
# # ----------------------------------------------------------------------------
# #                               (A)  - CRS
# 
# # Define the target CRS
# target_crs <- "EPSG:25832"
# 
# 
# roadKills <- project(roadKills, target_crs)
# roadKillsCount <- project(roadKillsCount, target_crs)
# broadleave <- project(broadleave, target_crs)
# coniferous <- project(coniferous, target_crs)
# mixedForest <- project(mixedForest, target_crs)
# transitionalLand <- project(transitionalLand, target_crs)
# pastures <- project(pastures , target_crs)
# grasslands <- project(grasslands, target_crs)
# roadType <- project(roadType, target_crs)
# humanInfluence <- project(humanInfluence, target_crs)
# watercourses <- project(watercourses, target_crs)
# lakes <- project(lakes, target_crs)
# #sdm <- project(sdm, target_crs)
# 
# # ----------------------------------------------------------------------------
# #                               (B)  - Resolution
# 
# # Define the target resolution (100 meters)
# target_res <- 100  # Target resolution in meters
# template_raster <- rast(roadKillFilePath)  # Copy extent and CRS from the reference raster (roadKills)
# res(template_raster) <- target_res  # Ensure resolution is set to 100x100 meters
# 
# 
# roadKills <- resample(roadKills,template_raster, method = "near")
# roadKillsCount <- resample(roadKillsCount,template_raster, method = "near")
# broadleave <- resample(broadleave,template_raster, method = "near")
# coniferous <- resample(coniferous,template_raster, method = "near")
# mixedForest <- resample(mixedForest,template_raster, method = "near")
# transitionalLand <- resample(transitionalLand,template_raster, method = "near")
# pastures <- resample(pastures,template_raster, method = "near")
# grasslands <- resample(grasslands,template_raster, method = "near")
# roadType <- resample(roadType,template_raster, method = "near")
# humanInfluence <- resample(humanInfluence,template_raster, method = "near")
# watercourses <- resample(watercourses,template_raster, method = "near")
# lakes<- resample(lakes,template_raster, method = "near")
# #sdm <- resample(sdm,template_raster, method = "bilinear")
# 
# 
# # ----------------------------------------------------------------------------
# #                               (C)  - extent
# 
# #start with setting the vector layer of south tyrols border as extent we want to crop the other layers to
# ST_border <- vect(ST_borderFilePath)
# crs(ST_border) <- target_crs
# 
# roadKills_crop <- mask(crop (roadKills, ST_border), ST_border)
# roadKillsCount_crop <- mask(crop (roadKillsCount, ST_border), ST_border)
# broadleave_crop <- mask(crop (broadleave, ST_border), ST_border)
# coniferous_crop <- mask(crop (coniferous, ST_border), ST_border)
# mixedForest_crop <- mask(crop (mixedForest, ST_border), ST_border)
# transitionalLand_crop <- mask(crop (transitionalLand, ST_border), ST_border)
# pastures_crop <- mask(crop (pastures, ST_border), ST_border)
# grasslands_crop <- mask(crop (grasslands, ST_border), ST_border)
# roadType_crop <- mask(crop (roadType, ST_border), ST_border)
# humanInfluence_crop <- mask(crop (humanInfluence, ST_border), ST_border)
# watercourses_crop <- mask(crop (watercourses, ST_border), ST_border)
# lakes_crop <- mask(crop (lakes, ST_border), ST_border)
# #sdm_crop <- mask(crop (sdm_crop, ST_border), ST_border)
# 
# # ----------------------------------------------------------------------------
# # save all layers and in the future use them!
# # ----------------------------------------------------------------------------
# 
# writeRaster(roadKills_crop, file.path(rasterExportPath, "roadKills.tif"), overwrite = TRUE)
# writeRaster(roadKillsCount_crop, file.path(rasterExportPath, "roadKillsCount.tif"), overwrite = TRUE)
# writeRaster(broadleave_crop, file.path(rasterExportPath, "broadleave.tif"), overwrite = TRUE)
# writeRaster(coniferous_crop, file.path(rasterExportPath, "coniferous.tif"), overwrite = TRUE)
# writeRaster(mixedForest_crop, file.path(rasterExportPath, "mixedForest.tif"), overwrite = TRUE)
# writeRaster(transitionalLand_crop, file.path(rasterExportPath, "transitionalLand.tif"), overwrite = TRUE)
# writeRaster(pastures_crop, file.path(rasterExportPath, "pastures.tif"), overwrite = TRUE)
# writeRaster(grasslands_crop, file.path(rasterExportPath, "grasslands.tif"), overwrite = TRUE)
# writeRaster(roadType_crop, file.path(rasterExportPath, "roadType.tif"), overwrite = TRUE)
# writeRaster(humanInfluence_crop, file.path(rasterExportPath, "humanInfluence.tif"), overwrite = TRUE)
# writeRaster(watercourses_crop, file.path(rasterExportPath, "watercourses.tif"), overwrite = TRUE)
# writeRaster(lakes_crop, file.path(rasterExportPath, "lakes.tif"), overwrite = TRUE)
# #writeRaster(sdm_crop, file.path(rasterExportPath, "sdm.tif"), overwrite = TRUE)

######
#############       SKIPPING OF CODE - E N D                    ##############
######

# # # Load all the files # # #
roadKills <- rast(file.path(rasterExportPath, "roadKills.tif"))
roadKillsCount <- rast(file.path(rasterExportPath, "roadKillsCount.tif"))
broadleave <- rast(file.path(rasterExportPath, "broadleave.tif"))
coniferous <- rast(file.path(rasterExportPath, "coniferous.tif"))
mixedForest <- rast(file.path(rasterExportPath, "mixedForest.tif"))
transitionalLand <- rast(file.path(rasterExportPath, "transitionalLand.tif"))
pastures <- rast(file.path(rasterExportPath, "pastures.tif"))
grasslands <- rast(file.path(rasterExportPath, "grasslands.tif"))
roadType <- rast(file.path(rasterExportPath, "roadType.tif"))
humanInfluence <- rast(file.path(rasterExportPath, "roadType.tif"))
watercourses <- rast(file.path(rasterExportPath, "humanInfluence.tif"))
lakes <- rast(file.path(rasterExportPath, "lakes_crop.tif"))
#sdm <- rast(rasterExportPath, "sdm.tif")


# ----------------------------------------------------------------------------
# Visualization of the data
# ----------------------------------------------------------------------------
library(ggplot2) # needed for some graphs

##############################################################################
# Description of the data:
# most variables are binary, the response variable (road kill counts) is poisson
# the distribution of human influence is XXXXXX,
# the distribution of road type is XXXXX

# ----------------- (1) Extract values from rasters
# Extract values and raster names

roadKills_values <- values(roadKills)
roadKills_values <- as.factor(ifelse(roadKills_values == 1, "Yes", "No")) #binary data
roadKillsCount_values <- values(roadKillsCount)
broadleave_values <- values(broadleave)
broadleave_values <- as.factor(ifelse(broadleave_values == 1, "Yes", "No")) #binary data
coniferous_values <- values(coniferous)
coniferous_values <- as.factor(ifelse(coniferous_values == 1, "Yes", "No")) #binary data
mixedForest_values <- values(mixedForest)
mixedForest_values <- as.factor(ifelse(mixedForest_values == 1, "Yes", "No")) #binary data
transitionalLand_values <- values(transitionalLand)
transitionalLand_values <- as.factor(ifelse(transitionalLand_values == 1, "Yes", "No")) #binary data
pastures_values <- values(pastures)
pastures_values <- as.factor(ifelse(pastures_values == 1, "Yes", "No")) #binary data
grasslands_values <- values(grasslands)
grasslands_values <- as.factor(ifelse(grasslands_values == 1, "Yes", "No")) #binary data
roadType_values <- values(roadType)
humanInfluence_values <- values(humanInfluence)
watercourses_values <- values(watercourses)
watercourses_values <- as.factor(ifelse(watercourses_values == 1, "Yes", "No")) #binary data
lakes_values <- values(lakes)
lakes_values <- as.factor(ifelse(lakes_values == 1, "Yes", "No")) #binary data

# Create the factor with appropriate levels
roadType_factor <- factor(roadType_values, 
                          levels = c("0", "1", "2", "3", "4", "5"),
                          labels = c("No Road", "Highway", "State Road", 
                                     "Country Road", "Municipal Road", 
                                     "Adress Road"))

# Check the levels to confirm
levels(roadType_factor)

# Create the data frame
df <- data.frame(roadKills_values, roadKillsCount_values, broadleave_values,
                 coniferous_values, mixedForest_values, transitionalLand_values,
                 pastures_values, grasslands_values, roadType_factor, 
                 humanInfluence_values, watercourses_values, lakes_values)

# Assign column names
colnames(df) <- c("Road_kills", "Road_kills_Count", "Broadleave", "Coniferous", 
                  "Mixed_Forest", "Small_Woody_Features", "Pastures",
                  "Grasslands", "Road_Type", "Human_Influence", "Watercourses",
                  "Lakes")

# Inspect the first few rows
head(df)

# Remove rows where "Road kills Count" is NA
df_na <- subset(df, !is.na(df$Road_kills_Count))

# Inspect the filtered data
head(df_na)

# Remove rows where "Road kill Count" is 0
df_greater0 <- subset(df_na, df_na$Road_kills_Count != 0)

# Make simple boxplots
par(mfrow = c(2,2))
boxplot(Road_kills_Count ~ Road_Type,
        data = df_greater0,
        # change axes labels
        xlab = "Road type",
        ylab = "Road kill count")

boxplot(Road_kills_Count ~ Broadleave,
        data = df_greater0,
        # change axes labels
        xlab = "Broadleaved Forest",
        ylab = "Road kill count")

boxplot(Road_kills_Count ~ Coniferous,
        data = df_greater0,
        # change axes labels
        xlab = "Coniferous forest",
        ylab = "Road kill count")

boxplot(Road_kills_Count ~ Mixed_Forest,
        data = df_greater0,
        # change axes labels
        xlab = "Mixed forest",
        ylab = "Road kill count")

boxplot(Road_kills_Count ~ Small_Woody_Features,
        data = df_greater0,
        # change axes labels
        xlab = "Small Woody Features",
        ylab = "Road kill count")

boxplot(Road_kills_Count ~ Pastures,
        data = df_greater0,
        # change axes labels
        xlab = "Pastures",
        ylab = "Road kill count")

boxplot(Road_kills_Count ~ Grasslands,
        data = df_greater0,
        # change axes labels
        xlab = "Grasslands",
        ylab = "Road kill count")

boxplot(Road_kills_Count ~ Human_Influence,
        data = df_greater0,
        # change axes labels
        xlab = "Human Influence",
        ylab = "Road kill count")

boxplot(Road_kills_Count ~ Watercourses,
        data = df_greater0,
        # change axes labels
        xlab = "Watercourses",
        ylab = "Road kill count")

boxplot(Road_kills_Count ~ Lakes,
        data = df_greater0,
        # change axes labels
        xlab = "Lakes",
        ylab = "Road kill count")


# ----------------------------------------------------------------------------
#                                   glm
# ----------------------------------------------------------------------------
# Activate the libraries needed for this 
library(car) # needed for the Anova() function, contains function vif(), which offers an easy alternative way to check predictor independence.


