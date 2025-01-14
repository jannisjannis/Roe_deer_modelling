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

###############################################################################
# Filepaths for data
###############################################################################
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

# # Reproject all raster variables to the target CRS
# for (i in seq_along(all_variables)) {
#   obj <- all_variables[[i]]
# 
#   # Reproject to target CRS
#   all_variables[[i]] <- project(obj, target_crs)  # Reproject to target CRS
#   cat("Reprojected raster", i, "to the target CRS.\n")
# }


roadKills <- project(roadKills, target_crs)
roadKillsCount <- project(roadKillsCount, target_crs)
broadleave <- project(broadleave, target_crs)
coniferous <- project(coniferous, target_crs)
mixedForest <- project(mixedForest, target_crs)
transitionalLand <- project(transitionalLand, target_crs)
pastures <- project(pastures , target_crs)
grasslands <- project(grasslands, target_crs)
roadType <- project(roadType, target_crs)
humanInfluence <- project(humanInfluence, target_crs)
watercourses <- project(watercourses, target_crs)
lakes <- project(lakes, target_crs)
#sdm <- project(sdm, target_crs)

# ----------------------------------------------------------------------------
#                               (B)  - Resolution

# Define the target resolution (100 meters)
target_res <- 100  # Target resolution in meters
template_raster <- rast(roadKillFilePath)  # Copy extent and CRS from the reference raster (roadKills)
res(template_raster) <- target_res  # Ensure resolution is set to 100x100 meters

# # Loop through all rasters and resample to the target resolution
# for (i in seq_along(all_variables)) {
#   obj <- all_variables[[i]]
# 
#   # Resample to the target resolution (100 meters)
#   all_variables[[i]] <- resample(x=obj, y = template_raster, method = "near")  # 'bilinear' for continuous, 'near' for categorical
#   cat("Resampled raster", i, "to 100 meters resolution.\n")
# }

roadKills <- resample(roadKills,template_raster, method = "near")
roadKillsCount <- resample(roadKillsCount,template_raster, method = "near")
broadleave <- resample(broadleave,template_raster, method = "near")
coniferous <- resample(coniferous,template_raster, method = "near")
mixedForest <- resample(mixedForest,template_raster, method = "near")
transitionalLand <- resample(transitionalLand,template_raster, method = "near")
pastures <- resample(pastures,template_raster, method = "near")
grasslands <- resample(grasslands,template_raster, method = "near")
roadType <- resample(roadType,template_raster, method = "near")
humanInfluence <- resample(humanInfluence,template_raster, method = "near")
watercourses <- resample(watercourses,template_raster, method = "near")
lakes<- resample(lakes,template_raster, method = "near")
#sdm <- resample(sdm,template_raster, method = "bilinear")


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
# save all layers and in the future use them!
# ----------------------------------------------------------------------------

writeRaster(roadKills_crop, file.path(rasterExportPath, "roadKills.tif"), overwrite = TRUE)
writeRaster(roadKillsCount_crop, file.path(rasterExportPath, "roadKillsCount.tif"), overwrite = TRUE)
writeRaster(broadleave_crop, file.path(rasterExportPath, "broadleave.tif"), overwrite = TRUE)
writeRaster(coniferous_crop, file.path(rasterExportPath, "coniferous.tif"), overwrite = TRUE)
writeRaster(mixedForest_crop, file.path(rasterExportPath, "mixedForest.tif"), overwrite = TRUE)
writeRaster(transitionalLand_crop, file.path(rasterExportPath, "transitionalLand.tif"), overwrite = TRUE)
writeRaster(pastures_crop, file.path(rasterExportPath, "pastures.tif"), overwrite = TRUE)
writeRaster(grasslands_crop, file.path(rasterExportPath, "grasslands.tif"), overwrite = TRUE)
writeRaster(roadType_crop, file.path(rasterExportPath, "roadType.tif"), overwrite = TRUE)
writeRaster(humanInfluence_crop, file.path(rasterExportPath, "humanInfluence.tif"), overwrite = TRUE)
writeRaster(watercourses_crop, file.path(rasterExportPath, "watercourses.tif"), overwrite = TRUE)
writeRaster(lakes_crop, file.path(rasterExportPath, "lakes.tif"), overwrite = TRUE)
#writeRaster(sdm_crop, file.path(rasterExportPath, "sdm.tif"), overwrite = TRUE)

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
roadType <- rast(file.path(rasterExportPath, "mixedForest.tif"))
humanInfluence <- rast(file.path(rasterExportPath, "roadType.tif"))
watercourses <- rast(file.path(rasterExportPath, "watercourses.tif"))
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
roadKills_values <- ifelse(roadKills_values == 1, "Yes", "No") #binary data
roadKillsCount_values <- values(roadKillsCount)
broadleave_values <- values(broadleave)
broadleave_values <- ifelse(broadleave_values == 1, "Yes", "No") #binary data
coniferous_values <- values(coniferous)
coniferous_values <- ifelse(coniferous_values == 1, "Yes", "No") #binary data
mixedForest_values <- values(mixedForest)
mixedForest_values <- ifelse(mixedForest_values == 1, "Yes", "No") #binary data
transitionalLand_values <- values(transitionalLand)
transitionalLand_values <- ifelse(transitionalLand_values == 1, "Yes", "No") #binary data
pastures_values <- values(pastures)
pastures_values <- ifelse(pastures_values == 1, "Yes", "No") #binary data
grasslands_values <- values(grasslands)
grasslands_values <- ifelse(grasslands_values == 1, "Yes", "No") #binary data
roadType_values <- values(roadType)
humanInfluence_values <- values(humanInfluence)
watercourses_values <- values(watercourses)
watercourses_values <- ifelse(watercourses_values == 1, "Yes", "No") #binary data
lakes_values <- values(lakes)
lakes_values <- ifelse(lakes_values == 1, "Yes", "No") #binary data


roadType_factor <- factor(roadType_values, levels = c("No Road", "Highway", "State Road", 
                                                      "Country Road", "Municipal Road",
                                                      "Adress Road"))
levels(roadType_factor)
# Replace ALL numeric entries by road type
roadType_factor[roadType_factor == "0"] <- "No Road"
roadType_factor[roadType_factor == "1"] <- "Highway"
roadType_factor[roadType_factor == "2"] <- "State Road"
roadType_factor[roadType_factor == "3"] <- "Country Road"
roadType_factor[roadType_factor == "4"] <- "Municipal Road"
roadType_factor[roadType_factor == "5"] <- "Adress Road"
levels(roadType_factor)

df <- data.frame(roadKills_values , roadKillsCount_values, broadleave_values,
                 coniferous_values,mixedForest_values,transitionalLand_values,
                 transitionalLand_values,pastures_values, grasslands_values,
                 roadType_factor, humanInfluence_values, watercourses_values,
                 lakes_values)
head(df)



# Beispiel für einen Boxplot
boxplot(roadKillsCount_values, horizontal = TRUE)
boxplot(roadType_factor, horizontal = TRUE)
boxplot(humanInfluence_values, horizontal = TRUE)

# Make simple boxplots
boxplot(roadKillsCount_values ~ roadType_factor,
        # change axes labels
        xlab = "Road type",
        ylab = "Road kill count")


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


