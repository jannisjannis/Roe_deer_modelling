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
# 
# ###############################################################################
# # Filepaths for data
# ###############################################################################
# 
# # response variable
# roadKillFilePath = "GIS/roadKills_presence_res100.tif"     
# roadKillCountsFilePath = "GIS/roadKills_count_raster_res100.tif"
# 
# ## predictor variables
# BLforestFilePath = "variables_glm_clipped_ST_minus_NP/forest.tif"
# CforestFilePath = "variables_glm_clipped_ST_minus_NP/forest.tif"
# MixforestFilePath = "variables_glm_clipped_ST_minus_NP/landuse.tif"
# transitionalLandFilePath = "forGLM/rastersized_data/smallWoodyFeatures_EPSG25832.tif"
# pasturesFilePath = "variables_glm_clipped_ST_minus_NP/landuse.tif"
# grasslandsFilePath = "variables_glm_clipped_ST_minus_NP/grassland_cover_100m.tif"
# sdmFilePath = "variables_glm_clipped_ST_minus_NP/species_distribution_RF.tif"
# watercoursesFilePath = "variables_glm_clipped_ST_minus_NP/watercourses.tif"
# lakesFilePath = "variables_glm_clipped_ST_minus_NP/lakes.tif"
# roadTypeFilePath = "variables_glm_clipped_ST_minus_NP/RoadType.tif"
# roadNetworkFilePath = "GIS/relevant_roads_raster_10by10.tif"
# humansFilePath = "variables_glm_clipped_ST_minus_NP/human_influence.tif"
# # predictor variables - topography
# demFilePath = "variables_glm_clipped_ST_minus_NP/dem_100m.tif"
# aspectFilePath = "variables_glm_clipped_ST_minus_NP/aspect_100m.tif"
# slopeFilePath = "variables_glm_clipped_ST_minus_NP/slope_100m.tif"
# 
# # Extent file (South Tyrol - NP)
# ST_borderFilePath = "GIS/border_southTyrol_withoutNP.shp"
# 
# 
# # # # # Load all the files # # #
# 
# 
# roadKills <- rast(roadKillFilePath)
# roadKillsCount <- rast(roadKillCountsFilePath)
# 
# 
# broadleave_unmasked <- rast(BLforestFilePath)
# broadleave <- app(broadleave_unmasked, function(x) ifelse(x == 1, x, 0))
# 
# coniferous_unmasked <- rast(CforestFilePath)
# coniferous <- app(coniferous_unmasked, function(x) ifelse(x == 2, x, 0))
# 
# mixedForest_unmasked <- rast(MixforestFilePath)
# mixedForest <- app(mixedForest_unmasked, function(x) ifelse(x == 25, x, 0))
# 
# transitionalLand_unmasked<- rast(transitionalLandFilePath)
# transitionalLand <- app(transitionalLand_unmasked, function(x) ifelse(is.na(x), 0, x))
# 
# pastures_unmasked <- rast(pasturesFilePath)
# pastures <- app(pastures_unmasked, function(x) ifelse(x == 18, x, 0))
# 
# grasslands <- rast(grasslandsFilePath)
# roadType <- rast(roadTypeFilePath)
# 
# roadNetwork_unmasked <- rast(roadNetworkFilePath)
# roadNetwork <- app(roadNetwork_unmasked, function(x) ifelse(is.na(x), 0, x))
# 
# humanInfluence <- rast(humansFilePath)
# watercourses <- rast(watercoursesFilePath)
# lakes<- rast(lakesFilePath)
# sdm <- rast(sdmFilePath)
# 
# dem <- rast(demFilePath)
# aspect <- rast(aspectFilePath)
# slope <- rast(slopeFilePath)
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
# roadNetwork <- project(roadNetwork, target_crs)
# humanInfluence <- project(humanInfluence, target_crs)
# watercourses <- project(watercourses, target_crs)
# lakes <- project(lakes, target_crs)
# sdm <- project(sdm, target_crs)
# dem <- project(dem, target_crs)
# aspect <- project(aspect, target_crs)
# slope <- project(slope, target_crs)
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
# roadNetwork <- resample(roadNetwork,template_raster, method = "near")
# humanInfluence <- resample(humanInfluence,template_raster, method = "near")
# watercourses <- resample(watercourses,template_raster, method = "near")
# lakes<- resample(lakes,template_raster, method = "near")
# sdm <- resample(sdm,template_raster, method = "bilinear")
# dem <- resample(dem, template_raster, method = "bilinear")
# aspect <- resample(aspect, template_raster, method = "bilinear")
# slope <- resample(slope, template_raster, method = "bilinear")
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
# roadNetwork_crop <- mask(crop (roadNetwork, ST_border), ST_border)
# humanInfluence_crop <- mask(crop (humanInfluence, ST_border), ST_border)
# watercourses_crop <- mask(crop (watercourses, ST_border), ST_border)
# lakes_crop <- mask(crop (lakes, ST_border), ST_border)
# sdm_crop <- mask(crop (sdm, ST_border), ST_border)
# dem_crop <- mask(crop (dem, ST_border), ST_border)
# aspect_crop <- mask(crop (aspect, ST_border), ST_border)
# slope_crop <- mask(crop (slope, ST_border), ST_border)
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
# writeRaster(roadNetwork_crop, file.path(rasterExportPath, "roadNetwork.tif"), overwrite = TRUE)
# writeRaster(humanInfluence_crop, file.path(rasterExportPath, "humanInfluence.tif"), overwrite = TRUE)
# writeRaster(watercourses_crop, file.path(rasterExportPath, "watercourses.tif"), overwrite = TRUE)
# writeRaster(lakes_crop, file.path(rasterExportPath, "lakes.tif"), overwrite = TRUE)
# writeRaster(sdm_crop, file.path(rasterExportPath, "sdm.tif"), overwrite = TRUE)
# writeRaster(dem_crop, file.path(rasterExportPath, "dem.tif"), overwrite = TRUE)
# writeRaster(aspect_crop, file.path(rasterExportPath, "aspect.tif"), overwrite = TRUE)
# writeRaster(slope_crop, file.path(rasterExportPath, "slope.tif"), overwrite = TRUE)
# 
# ####### DISTANCE VARIABLES
# 
# # Human influence as distance from nearest town
# humanInfluence[humanInfluence==0] <- NA
# distance_humanInfluence <- distance(humanInfluence)
# writeRaster(distance_humanInfluence, file.path(rasterExportPath, "distance_humanInfluence.tif"), overwrite = TRUE)
# 
# broadleave[broadleave==0] <- NA
# distance_broadleave <- distance(broadleave)
# writeRaster(distance_broadleave, file.path(rasterExportPath, "distance_broadleave.tif"), overwrite = TRUE)
# 
# coniferous[coniferous==0] <- NA
# distance_coniferous <- distance(coniferous)
# writeRaster(distance_coniferous, file.path(rasterExportPath, "distance_coniferous.tif"), overwrite = TRUE)
# 
# mixedForest[mixedForest==0] <- NA
# distance_mixedForest <- distance(mixedForest)
# writeRaster(distance_mixedForest, file.path(rasterExportPath, "distance_mixedForest.tif"), overwrite = TRUE)
# 
# transitionalLand[transitionalLand==0] <- NA
# distance_transitionalLand <- distance(transitionalLand)
# writeRaster(distance_transitionalLand, file.path(rasterExportPath, "distance_transitionalLand.tif"), overwrite = TRUE)
# 
# pastures[pastures==0] <- NA
# distance_pastures <- distance(pastures)
# writeRaster(distance_pastures, file.path(rasterExportPath, "distance_pastures.tif"), overwrite = TRUE)
# 
# grasslands[grasslands==0] <- NA
# distance_grasslands <- distance(grasslands)
# writeRaster(distance_grasslands, file.path(rasterExportPath, "distance_grasslands.tif"), overwrite = TRUE)
# 
# watercourses[watercourses==0] <- NA
# distance_watercourses <- distance(watercourses)
# writeRaster(distance_watercourses, file.path(rasterExportPath, "distance_watercourses.tif"), overwrite = TRUE)
#
# lakes[lakes==0] <- NA
# distance_lakes <- distance(lakes)
# writeRaster(distance_lakes, file.path(rasterExportPath, "distance_lakes.tif"), overwrite = TRUE)
# 


######
#############       SKIPPING OF CODE - E N D                    ##############
######

# # # Load all the files # # #
roadKills_rast <- rast(file.path(rasterExportPath, "roadKills.tif"))
roadKillsCount_rast  <- rast(file.path(rasterExportPath, "roadKillsCount.tif"))
broadleave_rast  <- rast(file.path(rasterExportPath, "broadleave.tif"))
coniferous_rast  <- rast(file.path(rasterExportPath, "coniferous.tif"))
mixedForest_rast  <- rast(file.path(rasterExportPath, "mixedForest.tif"))
transitionalLand_rast  <- rast(file.path(rasterExportPath, "transitionalLand.tif"))
pastures_rast  <- rast(file.path(rasterExportPath, "pastures.tif"))
grasslands_rast  <- rast(file.path(rasterExportPath, "grasslands.tif"))
roadType_rast  <- rast(file.path(rasterExportPath, "roadType.tif"))
roadNetwork_rast  <- rast(file.path(rasterExportPath, "roadNetwork.tif"))
humanInfluence_rast  <- rast(file.path(rasterExportPath, "humanInfluence.tif"))
watercourses_rast  <- rast(file.path(rasterExportPath, "watercourses.tif"))
lakes_rast  <- rast(file.path(rasterExportPath, "lakes.tif"))
sdm_rast  <- rast(file.path(rasterExportPath,"sdm.tif"))
dem_rast  <- rast(file.path(rasterExportPath, "dem.tif"))
aspect_rast  <- rast(file.path(rasterExportPath, "aspect.tif"))
slope_rast  <- rast(file.path(rasterExportPath, "slope.tif"))

# all DISTANCE files
broadleave_dist <- rast(file.path(rasterExportPath, "distance_broadleave.tif"))
coniferous_dist <- rast(file.path(rasterExportPath, "distance_coniferous.tif"))
mixedForest_dist <- rast(file.path(rasterExportPath, "distance_mixedForest.tif"))
transitionalLand_dist <- rast(file.path(rasterExportPath, "distance_transitionalLand.tif"))
pastures_dist <- rast(file.path(rasterExportPath, "distance_pastures.tif"))
grasslands_dist <- rast(file.path(rasterExportPath, "distance_grasslands.tif"))
watercourses_dist <- rast(file.path(rasterExportPath, "distance_watercourses.tif"))
lakes_dist <- rast(file.path(rasterExportPath, "distance_lakes.tif"))
humanInfluence_dist <- rast(file.path(rasterExportPath, "distance_humanInfluence.tif"))


# ----------------------------------------------------------------------------
# Visualization of the data
# ----------------------------------------------------------------------------
library(ggplot2) # needed for some graphs

##############################################################################
# Description of the data:
# most variables are binary or continuous, the response variable (road kill counts) is poisson
# the distribution of human influence is XXXXXX,
# the distribution of road type is XXXXX

# ----------------- (1) Extract values from rasters
# Extract values and raster names

roadKills_values <- values(roadKills_rast)
roadKills_values <- as.factor(ifelse(roadKills_values == 1, "Yes", "No")) #binary data
roadKillsCount_values <- as.integer(values(roadKillsCount_rast ))
broadleave_values <- values(broadleave_rast)
broadleave_values <- as.factor(ifelse(broadleave_values == 1, "Yes", "No")) #binary data
coniferous_values <- values(coniferous_rast)
coniferous_values <- as.factor(ifelse(coniferous_values == 2, "Yes", "No")) #binary data
mixedForest_values <- values(mixedForest_rast)
mixedForest_values <- as.factor(ifelse(mixedForest_values == 25, "Yes", "No")) #binary data
transitionalLand_values <- values(transitionalLand_rast)
transitionalLand_values <- as.factor(ifelse(transitionalLand_values == 1, "Yes", "No")) #binary data
pastures_values <- values(pastures_rast)
pastures_values <- as.factor(ifelse(pastures_values == 18, "Yes", "No")) #binary data
grasslands_values <- values(grasslands_rast)
grasslands_values <- as.factor(ifelse(grasslands_values == 1, "Yes", "No")) #binary data
roadType_values <- values(roadType_rast)
roadNetwork_values <- values(roadNetwork_rast)
roadNetwork_values <- as.factor(ifelse(roadNetwork_values == 1, "Yes", "No")) #binary data
humanInfluence_values <- values(humanInfluence_dist)
watercourses_values <- values(watercourses_rast)
watercourses_values <- as.factor(ifelse(watercourses_values == 1, "Yes", "No")) #binary data
lakes_values <- values(lakes_rast)
lakes_values <- as.factor(ifelse(lakes_values == 1, "Yes", "No")) #binary data
sdm_values <- values(sdm_rast)
dem_values <- values(dem_rast)
aspect_values <- values(aspect_rast)
slope_values <- values(slope_rast)

####### EXTRACT DISTANCE VALUES
broadleave_dist_values <- values(broadleave_dist)
coniferous_dist_values <- values(coniferous_dist)
mixedForest_dist_values <- values(mixedForest_dist)
transitionalLand_dist_values <- values(transitionalLand_dist)
pastures_dist_values <- values(pastures_dist)
grasslands_dist_values <- values(grasslands_dist)
watercourses_dist_values <- values(watercourses_dist)
lakes_dist_values <- values(lakes_dist)


# Create the factor with appropriate levels
roadType_factor <- factor(roadType_values, 
                          levels = c("0", "1", "2", "3", "4", "5"),
                          labels = c("No Road", "Highway", "State Road", 
                                     "Country Road", "Municipal Road", 
                                     "Adress Road"))

# Check the levels to confirm
levels(roadType_factor)

# Create the data frame (without the binary variable roadKills (roadKills_values))
df <- data.frame(roadKillsCount_values, broadleave_values,
                 coniferous_values, mixedForest_values, transitionalLand_values,
                 pastures_values, grasslands_values, roadType_factor, 
                 roadNetwork_values,humanInfluence_values, watercourses_values, 
                 lakes_values, sdm_values, dem_values, aspect_values, slope_values)

# Assign column names
colnames(df) <- c("Road_kills_Count", "Broadleave", "Coniferous", 
                  "Mixed_Forest", "Small_Woody_Features", "Pastures",
                  "Grasslands", "Road_Type", "Road_network", "Human_Influence", "Watercourses",
                  "Lakes", "SDM", "DEM", "aspect", "slope")

# Create data frame with continuous distance values
df_dist <- data.frame(roadKillsCount_values, broadleave_dist_values,
                      coniferous_dist_values, mixedForest_dist_values, 
                      transitionalLand_dist_values, pastures_dist_values, 
                      grasslands_dist_values, roadType_factor, roadNetwork_values,
                      humanInfluence_values, watercourses_dist_values, 
                      lakes_dist_values, sdm_values, dem_values, aspect_values,
                      slope_values)
# Assign column names
colnames(df_dist) <- c("Road_kills_Count", "Broadleave", "Coniferous", 
                  "Mixed_Forest", "Small_Woody_Features", "Pastures",
                  "Grasslands", "Road_Type", "Road_network", "Human_Influence", "Watercourses",
                  "Lakes", "SDM","DEM", "aspect", "slope")



# Inspect the first few rows
head(df)

# Remove rows where "Road kills Count" is NA
df_na <- subset(df, !is.na(df$Road_kills_Count))
df_dist_na <- subset(df_dist, !is.na(df$Road_kills_Count))

# Remove rows where "Road kill Count" is 0
df_greater0 <- subset(df_na, df_na$Road_kills_Count != 0)

# ------------------------- Make simple plots ----------------------------------
####                Boxplots for binary/factorial variables                 ####

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

# all factorial: Broadleave + Coniferous + Mixed_Forest + Small_Woody_Features + Pastures + Grasslands + Watercourses + Lakes


####                  Histogram for numerical variable                     ####

hist(df_na$Human_Influence,
     xlab = "Distance to human settlements [m]",
     col = "lightblue")

######## histograms of all my numeric variables
lapply(names(df_dist_na)[sapply(df_dist_na, is.numeric)], function(var) {
  hist(df_dist_na[[var]], main = paste("Histogram of", var), xlab = var, col = "lightblue", breaks = 30)
})

########          Scatterplots
# List of numeric variables (excluding the response variable)
numeric_vars <- names(df_dist_na)[sapply(df_dist_na, is.numeric)]

# Create scatterplots for each numeric predictor against the response variable
lapply(numeric_vars, function(var) {
  plot(df_dist_na[[var]], df_dist_na$Road_kills_Count, main = paste("Scatterplot of", var, "vs Response"), 
       xlab = var, ylab = "Response Variable", col = "blue", pch = 16)
})



# ----------------------------------------------------------------------------
#                              Check assumptions
# ----------------------------------------------------------------------------
# Activate the libraries needed for this 
library(car) # needed for the Anova() function, contains function vif(), which offers an easy alternative way to check predictor independence.

# ----------------- (A) - Independence of predictors ---------------------------
# To find out: check by pairwise "correlations" among predictor values
library(corrplot) 

df_dist_na_numeric <- df_dist_na[sapply(df_dist_na, is.numeric)]
cor_matrix <- cor(df_dist_na_numeric, method = "spearman", use = "pairwise.complete.obs") #spearman correltion to test for correlation
colors <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA")) #define colorpalette

par(mfrow=c(1,1))
#show correlation plot
corrplot(cor_matrix, 
         method = "color",                     # Use color tiles for visualization
         col = colors(200),                    # Apply your custom color palette
         type = "upper",                       # Show only the upper triangle
         order = "hclust",                     # Cluster variables based on correlations
         addCoef.col = "black",                # Add correlation coefficients in black
         number.cex = 0.7,                     # Adjust size of correlation numbers
         tl.cex = 0.3,                         # Adjust size of text labels
         title = "Spearman Correlation Matrix", # Add a title
         mar = c(0, 0, 2, 0),                  # Adjust margins for the title
         diag = FALSE                          # Do not show the diagonal
)
# Find variable pairs with |correlation| > 0.8
high_corr_pairs <- which(abs(cor_matrix) > 0.8, arr.ind = TRUE)
# Exclude self-correlations
high_corr_pairs <- high_corr_pairs[high_corr_pairs[, 1] != high_corr_pairs[, 2], ]
# Display highly correlated pairs
print(high_corr_pairs) # yaaay! No highly correlated pairs!


# ----------------- (B) - Dispersion parameter ---------------------------
# To obtain model validation plots, We first need to calculate the glm
# ... WITHOUT checking its statistical results!

# First model option without interaction and binary predictors
mod1.counts <- glm(Road_kills_Count ~  Broadleave + Coniferous + Mixed_Forest + 
                  Small_Woody_Features + Pastures + Grasslands + Watercourses +
                  Lakes + Human_Influence + Road_Type + SDM + DEM + aspect + slope, 
                data = df_na, 
                family = poisson(link = "log"))
mod2.counts <- glm(Road_kills_Count ~  Broadleave + Coniferous + Mixed_Forest + 
                  Small_Woody_Features + Pastures + Grasslands + Watercourses +
                  Lakes + Human_Influence + Road_Type + SDM + DEM + aspect + slope, 
                data = df_na, 
                family = quasipoisson(link = "log"))

mod3.distance <- glm(Road_kills_Count ~  Broadleave + Coniferous + Mixed_Forest + 
                       Small_Woody_Features + Pastures + Grasslands + Watercourses +
                       Lakes + Human_Influence + Road_Type + SDM + DEM + aspect + slope, 
                     data = df_dist_na, 
                     family = poisson(link = "log"))

# the variance inflation factor VIF
# ... offers an easy alternative way to check for predictor independence.
# Without further explanation: Kick out predictors that have a value > 5 in the following list:
# ... They are too heavily correlated with other predictors.
vif(mod1.counts) # all looks good --> Predictors are independent

# model validation: residual plots only partially informative 
# ...many road kill counts are 0 or close to 0
#par(mfrow = c(2,2))
#plot(mod1.glm) 

# --> MODEL VALIDATION with dispersion parameter
# What the model believes:
# mean value = data variance --> obs. / exp. dispersion = 1
# True dispersion in the data: Deviance / df
summary(mod1.counts)
dispersion.parameter <- mod1.counts$deviance / mod1.counts$df.residual
dispersion.parameter

summary(mod2.counts)
dispersion.parameter <- mod2.counts$deviance / mod2.counts$df.residual
dispersion.parameter

summary(mod3.distance)
dispersion.parameter <- mod3.distance$deviance / mod3.distance$df.residual
dispersion.parameter

# dispersion parameters a way too small --> UNDERDISPERSION
# -->  zero-inflated models


# ----------------- (C) - check for near-zero variance -------------------------
# is needed for zero inflated model, else it can cause numerical instability
# check for near-zero variance:
library(caret)
nzv <- nearZeroVar(df_dist_na, saveMetrics = TRUE)
print(nzv)
# --> Road_network has near-zero variance

# ----------------- (D) - check NaNs  ---------------------------------
# Proportion of NaNs in each variable
prop_slope_na <- sum(is.na(df_dist_na$slope)) / nrow(df_dist_na)
prop_aspect_na <- sum(is.na(df_dist_na$aspect)) / nrow(df_dist_na)
prop_dem_na <- sum(is.na(df_dist_na$DEM)) / nrow(df_dist_na)

# Print the proportions
prop_slope_na
prop_aspect_na
prop_dem_na

# biggest proportion is 0.0019 (1333 pixels)
# Solution: remove them

df_dist_na_clean <- df_dist_na[complete.cases(df_dist_na$slope, df_dist_na$aspect, df_dist_na$DEM), ]

# ----------------- (D) - check Multicollinearity ------------------------------
# see below with the vif() function


# ----------------------------------------------------------------------------
#                                 Model selection
# ----------------------------------------------------------------------------
# Load required packages
library(pscl)
library(caret)
library(MASS)
# in general a zero inlated model with poisson distribution was chosen
# poisson due to counts, zero inflated due to many cells with 0


# part behind | --> Explains structural zeros 
#...(e.g., locations where road kills cannot occur due to unsuitable conditions, 
#...such as areas with no roads)
mod.zip <- zeroinfl(Road_kills_Count ~  Broadleave + Coniferous + Mixed_Forest + 
                      Small_Woody_Features + Pastures + Grasslands + Watercourses +
                      Lakes  + Human_Influence + Road_Type + SDM + DEM + aspect + slope       
                    | SDM + DEM + aspect + slope                                                        
                    data = df_dist_na_clean, 
                    dist = "poisson")
vif(mod.zip)
summary(mod.zip)

mod.zip2 <- glmmTMB(Road_kills_Count ~ Road_Type+Broadleave + Coniferous + Mixed_Forest+ 
                      Small_Woody_Features + Pastures + Grasslands + Watercourses+
                      Lakes + Human_Influence+ DEM + aspect + slope,
                    zi = ~Road_Type + DEM + aspect + slope,
                    family = "poisson",
                    data = df_dist_na_clean)
summary(mod.zip2)


# model with binary variables
mod.zip.binary <- zeroinfl(Road_kills_Count ~  Broadleave + Coniferous + Mixed_Forest + 
                     Small_Woody_Features + Pastures + Grasslands + Watercourses +
                     Lakes  + Human_Influence + Road_Type #+ DEM + aspect + slope       # AM ENDE SDM hinzufügen!!
                   | Road_Type,                                                        # AM ENDE SDM hinzufügen!!
                   data = df_na,
                   dist = "poisson")

vif(mod.zip.binary)
summary(mod.zip.binary)

# compare AIC values (model fit)
AIC(mod1.counts, mod.zip) # mod.zip has the much better fit

# stepwise selection
step_model <- stepAIC(mod.zip, direction ="both", trace = TRUE)
summary(step_model)

library(MuMIn)
# selection by calculation of AIC for all combinations of predictors
# Generate all model subsets
model_set <- dredge(mod.zip, rank = "AIC")

# View best models
head(model_set)

# Get the best model
best_model <- get.models(model_set, subset = 1)[[1]]

summary(best_model)


# ----------------------------------------------------------------------------
#                                   glm - significance
# ----------------------------------------------------------------------------
library(emmeans)
# Anova() allows estimating the "overall" contribution of each predictor.
#  We use - and report - results of the Likelihood-ratio chi-square test (type III):
# LR: Likelihood Ratio Test compares models by their likelihood
# Anova(mod.zip, type = "III",test ="LR") 

# ----------------------------------------------------------------------------
#                   generate raster for resulting roadkill risk map
# ----------------------------------------------------------------------------
