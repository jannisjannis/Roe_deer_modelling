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
# sdmFilePath = "variables_glm_clipped_ST_minus_NP/species_distribution_RF.tif"
# watercoursesFilePath = "variables_glm_clipped_ST_minus_NP/watercourses.tif"
# lakesFilePath = "variables_glm_clipped_ST_minus_NP/lakes.tif"
# roadTypeFilePath = "variables_glm_clipped_ST_minus_NP/RoadType.tif"
# roadNetworkFilePath = "GIS/relevant_roads_raster_10by10.tif"
# humansFilePath = "variables_glm_clipped_ST_minus_NP/human_influence.tif"
# 
# 
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
# 
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
# distance_watercourses <- distance(grasslands)
# writeRaster(distance_watercourses, file.path(rasterExportPath, "distance_watercourses.tif"), overwrite = TRUE)
# 
# lakes[lakes==0] <- NA 
# distance_lakes <- distance(lakes)
# writeRaster(distance_lakes, file.path(rasterExportPath, "distance_lakes.tif"), overwrite = TRUE)



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
roadNetwork <- rast(file.path(rasterExportPath, "roadNetwork.tif"))
humanInfluence <- rast(file.path(rasterExportPath, "humanInfluence.tif"))
watercourses <- rast(file.path(rasterExportPath, "watercourses.tif"))
lakes <- rast(file.path(rasterExportPath, "lakes.tif"))
sdm <- rast("/Users/carlabehringer/iCloud Drive (Archive)/Documents/Documents – Carlas MacBook Air/dokumente/Uni/Master/second_year/year2_sem1/Com_EnvMan/Projectstudy/data/Export_R/raster/sdm.tif")

# all DISTANCE files
broadleave_dist <- rast(file.path(rasterExportPath, "distance_broadleave.tif"))
coniferous_dist <- rast(file.path(rasterExportPath, "distance_coniferous.tif"))
mixedForest_dist <- rast(file.path(rasterExportPath, "distance_mixedForest.tif"))
transitionalLand_dist <- rast(file.path(rasterExportPath, "distance_transitionalLand.tif"))
pastures_dist <- rast(file.path(rasterExportPath, "distance_pastures.tif"))
grasslands_dist <- rast(file.path(rasterExportPath, "distance_grasslands.tif"))
watercourses_dist <- rast(file.path(rasterExportPath, "distance_watercourses.tif"))
lakes_dist <- rast(file.path(rasterExportPath, "distance_lakes.tif"))


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

roadKills_values <- values(roadKills)
roadKills_values <- as.factor(ifelse(roadKills_values == 1, "Yes", "No")) #binary data
roadKillsCount_values <- as.integer(values(roadKillsCount))
broadleave_values <- values(broadleave)
broadleave_values <- as.factor(ifelse(broadleave_values == 1, "Yes", "No")) #binary data
coniferous_values <- values(coniferous)
coniferous_values <- as.factor(ifelse(coniferous_values == 2, "Yes", "No")) #binary data
mixedForest_values <- values(mixedForest)
mixedForest_values <- as.factor(ifelse(mixedForest_values == 25, "Yes", "No")) #binary data
transitionalLand_values <- values(transitionalLand)
transitionalLand_values <- as.factor(ifelse(transitionalLand_values == 1, "Yes", "No")) #binary data
pastures_values <- values(pastures)
pastures_values <- as.factor(ifelse(pastures_values == 18, "Yes", "No")) #binary data
grasslands_values <- values(grasslands)
grasslands_values <- as.factor(ifelse(grasslands_values == 1, "Yes", "No")) #binary data
roadType_values <- values(roadType)
roadNetwork_values <- values(roadNetwork)
roadNetwork_values <- as.factor(ifelse(roadNetwork_values == 1, "Yes", "No")) #binary data
humanInfluence_values <- values(distance_humanInfluence)
watercourses_values <- values(watercourses)
watercourses_values <- as.factor(ifelse(watercourses_values == 1, "Yes", "No")) #binary data
lakes_values <- values(lakes)
lakes_values <- as.factor(ifelse(lakes_values == 1, "Yes", "No")) #binary data
sdm_values <- values(sdm)

####### EXTRACT DISTANCE VALUES
humanInfluence_values <- values(distance_humanInfluence)
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
                 pastures_values, grasslands_values, roadType_factor, roadNetwork_values,
                 humanInfluence_values, watercourses_values, lakes_values, sdm_values)

# Assign column names
colnames(df) <- c("Road_kills_Count", "Broadleave", "Coniferous", 
                  "Mixed_Forest", "Small_Woody_Features", "Pastures",
                  "Grasslands", "Road_Type", "Road_network", "Human_Influence", "Watercourses",
                  "Lakes", "SDM")

# Create data frame with continuous distance values
df_dist <- data.frame(roadKillsCount_values, broadleave_dist_values,
                      coniferous_dist_values, mixedForest_dist_values, transitionalLand_dist_values,
                      pastures_dist_values, grasslands_dist_values, roadType_factor, roadNetwork_values,
                      humanInfluence_values, watercourses_dist_values, lakes_dist_values, sdm_values)
# Assign column names
colnames(df_dist) <- c("Road_kills_Count", "Broadleave", "Coniferous", 
                  "Mixed_Forest", "Small_Woody_Features", "Pastures",
                  "Grasslands", "Road_Type", "Road_network", "Human_Influence", "Watercourses",
                  "Lakes", "SDM")



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

# COMMENT: quite a lot of data logically at 0, therefore filtered out. 
# Plus an outlier: Bolzano with extremely high population density --> also filtered out

######## histograms of all my numeric variables
lapply(names(df_dist_na)[sapply(df_dist_na, is.numeric)], function(var) {
  hist(df_dist_na[[var]], main = paste("Histogram of", var), xlab = var, col = "lightblue", breaks = 30)
})

########          Scatterplots
# List of numeric variables (excluding the response variable)
numeric_vars <- names(df_dist_na)[sapply(df_dist_na, is.numeric)]

# Create scatterplots for each numeric predictor against the response variable
lapply(numeric_vars, function(var) {
  plot(df[[var]], df_dist_na$Road_kills_Count, main = paste("Scatterplot of", var, "vs Response"), 
       xlab = var, ylab = "Response Variable", col = "blue", pch = 16)
})



# ----------------------------------------------------------------------------
#                              Check assumptions
# ----------------------------------------------------------------------------
# Activate the libraries needed for this 
library(car) # needed for the Anova() function, contains function vif(), which offers an easy alternative way to check predictor independence.

# ----------------- (A) - Independence of predictors ---------------------------
# To find out: check by pairwise "correlations" among predictor values

# To obtain model validation plots, We first need to calculate the glm
# ... WITHOUT checking its statistical results!

# First model option without interaction and binary predictors
mod1.counts <- glm(Road_kills_Count ~  Broadleave + Coniferous + Mixed_Forest + 
                  Small_Woody_Features + Pastures + Grasslands + Watercourses +
                  Lakes + Human_Influence + Road_Type + SDM, 
                data = df_na, 
                family = poisson(link = "log"))
mod2.counts <- glm(Road_kills_Count ~  Broadleave + Coniferous + Mixed_Forest + 
                  Small_Woody_Features + Pastures + Grasslands + Watercourses +
                  Lakes + Human_Influence + Road_Type + SDM, 
                data = df_na, 
                family = quasipoisson(link = "log"))

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

# Let's try out zero-inflated models
# Load required package
library(pscl)
library(caret)
# check for near-zero variance:
nzv <- nearZeroVar(df_na, saveMetrics = TRUE)
print(nzv)

# part behind | --> Explains structural zeros 
#...(e.g., locations where road kills cannot occur due to unsuitable conditions, 
#...such as areas with no roads)
mod.zip <- zeroinfl(Road_kills_Count ~  Broadleave + Coniferous + Mixed_Forest + 
                      Small_Woody_Features + Pastures + Grasslands + Watercourses +
                      Lakes  + Human_Influence + Road_Type + SDM | Road_network + SDM, 
                    data = df_na, 
                    dist = "poisson")
summary(mod.zip)

# compare AIC values (model fit)
AIC(mod1.counts, mod.zip)


# ----------------------------------------------------------------------------
#                                   glm - significance
# ----------------------------------------------------------------------------
library(emmeans)
# Anova() allows estimating the "overall" contribution of each predictor.
#  We use - and report - results of the Likelihood-ratio chi-square test (type III):
# LR: Likelihood Ratio Test compares models by their likelihood
Anova(mod.zip, type = "III",test ="LR") 

# ----------------------------------------------------------------------------
#                   generate raster for resulting roadkill risk map
# ----------------------------------------------------------------------------
