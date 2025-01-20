# #-------------------------------------------------------------------#
# # Course:
# # Topic: Road kill risk Roe deer South Tyrol
# # Script for glm
# # Author: Carla Behringer
# # Date: 20.01.2024
# #-------------------------------------------------------------------#
# 
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
rasterExportPath <- "Export_R/raster_dist_buf"

######
#############       SKIPPING OF CODE - S T A R T                    ##############
######
#
###############################################################################
# Filepaths for data
###############################################################################

# response variable
roadKillCountsFilePath = "GIS/roadKills_count_raster_res100.tif"

## predictor variables
BLforestFilePath = "forGLM/landcover/europe/clipped_30km/FTY_clip.tif"
CforestFilePath = "forGLM/landcover/europe/clipped_30km/FTY_clip.tif"
MixforestFilePath = "forGLM/landcover/europe/clipped_30km/FTY_clip.tif"
transitionalLandFilePath = "forGLM/landcover/europe/clipped_30km/SWF_clip.tif"
pasturesFilePath = "forGLM/landcover/europe/clipped_30km/CLC_clip.tif"
grasslandsFilePath = "forGLM/landcover/europe/clipped_30km/GRA_clip.tif"
sdmFilePath = "variables_glm_clipped_ST_minus_NP/species_distribution_RF.tif"
waterFilePath = "forGLM/landcover/europe/clipped_30km/WAW_clip.tif"
roadTypeFilePath = "variables_glm_clipped_ST_minus_NP/RoadType.tif"
roadNetworkFilePath = "GIS/relevant_roads_raster_10by10.tif"
humansFilePath = "forGLM/landcover/europe/clipped_30km/CLC_clip.tif"
# predictor variables - topography
demFilePath = "variables_glm_clipped_ST_minus_NP/dem_100m.tif"
aspectFilePath = "variables_glm_clipped_ST_minus_NP/aspect_100m.tif"
slopeFilePath = "variables_glm_clipped_ST_minus_NP/slope_100m.tif"

# Extent files (South Tyrol - NP, road buffer)
ST_borderFilePath = "GIS/border_southTyrol_withoutNP.shp"
roads_buffer_vector_FilePath ="GIS/100mBuffer_roads_clip.shp"


# # # # Load all the files # # #

roadKillsCount <- rast(roadKillCountsFilePath)


broadleave_unmasked <- rast(BLforestFilePath)
broadleave <- app(broadleave_unmasked, function(x) ifelse(x == 1, x, 0))

coniferous_unmasked <- rast(CforestFilePath)
coniferous <- app(coniferous_unmasked, function(x) ifelse(x == 2, x, 0))

mixedForest_unmasked <- rast(MixforestFilePath)
mixedForest <- app(mixedForest_unmasked, function(x) ifelse(x == 3, x, 0))

transitionalLand_unmasked<- rast(transitionalLandFilePath)
transitionalLand <- app(transitionalLand_unmasked, function(x) ifelse(is.na(x), 0, x))

pastures_unmasked <- rast(pasturesFilePath)
pastures <- app(pastures_unmasked, function(x) ifelse(x == 18, x, 0))

grasslands <- rast(grasslandsFilePath)
roadType <- rast(roadTypeFilePath)

humanInfluence_unmasked <- rast(humansFilePath)
humanInfluence <- app(humanInfluence_unmasked, function(x) ifelse(x %in% c(1, 2, 3, 4), x, 0))

water_unmasked <- rast(waterFilePath)
water <- app(water_unmasked, function(x) ifelse(x == 1, x, 0))

sdm <- rast(sdmFilePath)

dem <- rast(demFilePath)
aspect <- rast(aspectFilePath)
slope <- rast(slopeFilePath)

roads_buf_vect <- vect(roads_buffer_vector_FilePath)

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

roadKillsCount <- project(roadKillsCount, target_crs)
broadleave <- project(broadleave, target_crs)
coniferous <- project(coniferous, target_crs)
mixedForest <- project(mixedForest, target_crs)
transitionalLand <- project(transitionalLand, target_crs)
pastures <- project(pastures , target_crs)
grasslands <- project(grasslands, target_crs)
roadType <- project(roadType, target_crs)
humanInfluence <- project(humanInfluence, target_crs)
water <- project(water, target_crs)
sdm <- project(sdm, target_crs)
dem <- project(dem, target_crs)
aspect <- project(aspect, target_crs)
slope <- project(slope, target_crs)

# ----------------------------------------------------------------------------
#                               (B)  - Resolution

# Define the target resolution (100 meters)
target_res <- 100  # Target resolution in meters
template_raster <- rast(roadKillCountsFilePath)  # Copy extent and CRS from the reference raster (roadKills)
res(template_raster) <- target_res  # Ensure resolution is set to 100x100 meters

roadKillsCount <- resample(roadKillsCount,template_raster, method = "near")
broadleave <- resample(broadleave,template_raster, method = "near")
coniferous <- resample(coniferous,template_raster, method = "near")
mixedForest <- resample(mixedForest,template_raster, method = "near")
transitionalLand <- resample(transitionalLand,template_raster, method = "near")
pastures <- resample(pastures,template_raster, method = "near")
grasslands <- resample(grasslands,template_raster, method = "near")
roadType <- resample(roadType,template_raster, method = "near")
humanInfluence <- resample(humanInfluence,template_raster, method = "near")
water <- resample(water,template_raster, method = "near")
sdm <- resample(sdm,template_raster, method = "bilinear")
dem <- resample(dem, template_raster, method = "bilinear")
aspect <- resample(aspect, template_raster, method = "bilinear")
slope <- resample(slope, template_raster, method = "bilinear")

####### DISTANCE VARIABLES

# Human influence as distance from urban fabric
humanInfluence[humanInfluence==0] <- NA
distance_humanInfluence <- distance(humanInfluence)

broadleave[broadleave==0] <- NA
distance_broadleave <- distance(broadleave)

coniferous[coniferous==0] <- NA
distance_coniferous <- distance(coniferous)

mixedForest[mixedForest==0] <- NA
distance_mixedForest <- distance(mixedForest)

pastures[pastures==0] <- NA
distance_pastures <- distance(pastures)

grasslands[grasslands==0] <- NA
distance_grasslands <- distance(grasslands)

water[water==0] <- NA
distance_water <- distance(water)
#
#
# # # ----------------------------------------------------------------------------
# # #                               (C)  - extent
# #
# #start with setting the vector layer of the roads buffer as extent we want to crop the other layers to
# roads_buf_vect <- vect(roads_buffer_vector_FilePath)
# crs(roads_buf_vect) <- target_crs
# 
# roadKillsCount_crop <- mask(crop (roadKillsCount, roads_buf_vect), roads_buf_vect)
# roadType_crop <- mask(crop (roadType, roads_buf_vect), roads_buf_vect)
# sdm_crop <- mask(crop (sdm, roads_buf_vect), roads_buf_vect)
# dem_crop <- mask(crop (dem, roads_buf_vect), roads_buf_vect)
# aspect_crop <- mask(crop (aspect, roads_buf_vect), roads_buf_vect)
# slope_crop <- mask(crop (slope, roads_buf_vect), roads_buf_vect)
# transitionalLand_crop <- mask(crop(transitionalLand, roads_buf_vect), roads_buf_vect)
# 
# distance_humanInfluence_crop<- mask(crop (distance_humanInfluence, roads_buf_vect), roads_buf_vect)
# distance_broadleave_crop <- mask(crop (distance_broadleave, roads_buf_vect), roads_buf_vect)
# distance_coniferous_crop <- mask(crop (distance_coniferous, roads_buf_vect), roads_buf_vect)
# distance_mixedForest_crop <- mask(crop (distance_mixedForest, roads_buf_vect), roads_buf_vect)
# distance_pastures_crop <- mask(crop (distance_pastures, roads_buf_vect), roads_buf_vect)
# distance_grasslands_crop <- mask(crop (distance_grasslands, roads_buf_vect), roads_buf_vect)
# distance_water_crop <- mask(crop (distance_water, roads_buf_vect), roads_buf_vect)
# 
# # ----------------------------------------------------------------------------
# # save all layers and in the future use them!
# # ----------------------------------------------------------------------------
# 
# writeRaster(roadKillsCount_crop, file.path(rasterExportPath, "roadKillsCount.tif"), overwrite = TRUE)
# writeRaster(roadType_crop, file.path(rasterExportPath, "roadType.tif"), overwrite = TRUE)
# writeRaster(sdm_crop, file.path(rasterExportPath, "sdm.tif"), overwrite = TRUE)
# writeRaster(dem_crop, file.path(rasterExportPath, "dem.tif"), overwrite = TRUE)
# writeRaster(aspect_crop, file.path(rasterExportPath, "aspect.tif"), overwrite = TRUE)
# writeRaster(slope_crop, file.path(rasterExportPath, "slope.tif"), overwrite = TRUE)
# writeRaster(transitionalLand_crop, file.path(rasterExportPath, "transitionalLand.tif"), overwrite = TRUE)
# 
# # write distance rasters
# writeRaster(distance_humanInfluence_crop, file.path(rasterExportPath, "distance_humanInfluence.tif"), overwrite = TRUE)
# writeRaster(distance_broadleave_crop, file.path(rasterExportPath, "distance_broadleave.tif"), overwrite = TRUE)
# writeRaster(distance_coniferous_crop, file.path(rasterExportPath, "distance_coniferous.tif"), overwrite = TRUE)
# writeRaster(distance_mixedForest_crop, file.path(rasterExportPath, "distance_mixedForest.tif"), overwrite = TRUE)
# writeRaster(distance_pastures_crop, file.path(rasterExportPath, "distance_pastures.tif"), overwrite = TRUE)
# writeRaster(distance_grasslands_crop, file.path(rasterExportPath, "distance_grasslands.tif"), overwrite = TRUE)
# writeRaster(distance_water_crop, file.path(rasterExportPath, "distance_water.tif"), overwrite = TRUE)

######
#############       SKIPPING OF CODE - E N D                    ##############
######

# # # Load all the files # # #
roadKillsCount_rast <- rast(file.path(rasterExportPath, "roadKillsCount.tif"))
roadType_rast <- rast(file.path(rasterExportPath, "roadType.tif"))
sdm_rast <- rast(file.path(rasterExportPath,"sdm.tif"))
dem_rast <- rast(file.path(rasterExportPath, "dem.tif"))
aspect_rast <- rast(file.path(rasterExportPath, "aspect.tif"))
slope_rast <- rast(file.path(rasterExportPath, "slope.tif"))
transitionalLand <- rast(file.path(rasterExportPath, "transitionalLand.tif"))

# all DISTANCE files
broadleave_dist <- rast(file.path(rasterExportPath, "distance_broadleave.tif"))
coniferous_dist <- rast(file.path(rasterExportPath, "distance_coniferous.tif"))
mixedForest_dist <- rast(file.path(rasterExportPath, "distance_mixedForest.tif"))
pastures_dist <- rast(file.path(rasterExportPath, "distance_pastures.tif"))
grasslands_dist <- rast(file.path(rasterExportPath, "distance_grasslands.tif"))
water_dist <- rast(file.path(rasterExportPath, "distance_water.tif"))
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

roadKillsCount_values <- as.integer(values(roadKillsCount_rast))
roadType_values <- values(roadType_rast)
sdm_values <- values(sdm_rast)
dem_values <- values(dem_rast)
aspect_values <- values(aspect_rast)
slope_values <- values(slope_rast)
transitionalLand_values <- values(transitionalLand)

####### EXTRACT DISTANCE VALUES

broadleave_dist_values <- values(broadleave_dist)
coniferous_dist_values <- values(coniferous_dist)
mixedForest_dist_values <- values(mixedForest_dist)

pastures_dist_values <- values(pastures_dist)
grasslands_dist_values <- values(grasslands_dist)
water_dist_values <- values(water_dist)
humanInfluence_dist_values <- values(humanInfluence_dist)

# Create the factor with appropriate levels
roadType_factor <- factor(roadType_values, 
                          levels = c("0", "1", "2", "3", "4", "5"),
                          labels = c("No Road", "Highway", "State Road", 
                                     "Country Road", "Municipal Road", 
                                     "Adress Road"))

# Check the levels to confirm
levels(roadType_factor)


# Create data frame with continuous distance values
df_dist <- data.frame(roadKillsCount_values, broadleave_dist_values,
                      coniferous_dist_values, mixedForest_dist_values, 
                      transitionalLand_values, pastures_dist_values, 
                      grasslands_dist_values, roadType_factor,
                      humanInfluence_dist_values, water_dist_values, sdm_values, 
                      dem_values, aspect_values, slope_values)
# Assign column names
colnames(df_dist) <- c("Road_kills_Count", "Broadleave", "Coniferous", 
                  "Mixed_Forest", "Small_Woody_Features", "Pastures",
                  "Grasslands", "Road_Type", "Human_Influence", "Water", "SDM",
                  "DEM", "aspect", "slope")

# Inspect the first few rows
head(df_dist)

# Remove rows where "Road kills Count" is NA or other variable is 
df_dist_na <- na.omit(df_dist)
head(df_dist_na)

# ------------------------- Make simple plots ----------------------------------
####                Boxplots for binary/factorial variables                 ####

par(mfrow = c(1,1))
boxplot(Road_kills_Count ~ Road_Type,
        data = df_dist_na,
        # change axes labels
        xlab = "Road type",
        ylab = "Road kill count")

# Count the frequency of each level in Road_Type
road_type_counts <- table(df_dist_na$Road_Type)

# Create a bar plot of the counts
barplot(road_type_counts, main = "Road Type Distribution", col = "skyblue", 
        xlab = "Road Type", ylab = "Count")

####                  Histogram for all numerical variables                     ####
par(mfrow = c(2,2))

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
mod1.poi <- glm(Road_kills_Count ~  Broadleave + Coniferous + Mixed_Forest + 
                  Small_Woody_Features + Pastures + Grasslands + Water + 
                  Human_Influence + Road_Type + DEM + aspect + slope,             # SDM HINZUFÜGEN
                data = df_dist_na, 
                family = poisson(link = "log"))


# the variance inflation factor VIF
# ... offers an easy alternative way to check for predictor independence.
# Without further explanation: Kick out predictors that have a value > 5 in the following list:
# ... They are too heavily correlated with other predictors.
vif(mod1.poi) # all looks good --> Predictors are independent

# model validation: residual plots only partially informative 
# ...many road kill counts are 0 or close to 0
#par(mfrow = c(2,2))
#plot(mod1.glm) 

# --> MODEL VALIDATION with dispersion parameter
# What the model believes:
# mean value = data variance --> obs. / exp. dispersion = 1
# True dispersion in the data: Deviance / df
summary(mod1.poi)
dispersion.parameter <- mod1.poi$deviance / mod1.poi$df.residual
dispersion.parameter # UNDERDISPERSION

# Let's try out zero-inflated models

# ----------------- (C) - check for near-zero variance -------------------------
# is needed for zero inflated model, else it can cause numerical instability
# check for near-zero variance:
library(caret)
nzv <- nearZeroVar(df_dist_na, saveMetrics = TRUE)
print(nzv)
# --> PREDICTOR has near-zero variance

# ----------------- (D) - check Multicollinearity ------------------------------
# see above with the vif() function

# ----------------------------------------------------------------------------
#                                 Model selection
# ----------------------------------------------------------------------------
# Load required packages
library(pscl)
library(caret)
library(MASS)
library(glmmTMB)
# in general a zero inflated model with poisson distribution was chosen
# poisson due to counts, zero inflated due to many cells with 0


# part behind | --> Explains structural zeros 
#...(e.g., locations where road kills cannot occur due to unsuitable conditions, 
#...such as areas with no roads)
mod.zip <- with(
  df_dist_na, 
  zeroinfl(
    Road_kills_Count ~ Broadleave + Coniferous + Mixed_Forest + 
      Small_Woody_Features + Pastures + Grasslands + Water+
      Human_Influence + Road_Type + SDM + DEM + aspect + slope       
    | SDM + DEM + aspect + slope,                                                         
    dist = "poisson"
  )
)
vif(mod.zip)
summary(mod.zip)

mod.zip2 <- glmmTMB(Road_kills_Count ~ Road_Type+Broadleave + Coniferous + Mixed_Forest+ 
                      Small_Woody_Features + Pastures + Grasslands + Water+
                      Human_Influence+ DEM+ aspect + slope,
                    zi = ~ DEM+ aspect + slope,
                    family = "poisson",
                    data = df_dist_na,
                    na.action = na.fail
                    )
diagnose(mod.zip2)
summary(mod.zip2)

# stepwise selection
step_model <- stepAIC(mod.zip2, direction ="both", trace = TRUE)
summary(step_model)
formula(step_model) # final model

library(MuMIn)
# selection by calculation of AIC for all combinations of predictors
# Generate all model subsets
model_set <- dredge(mod.zip2, rank = "AIC")

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
