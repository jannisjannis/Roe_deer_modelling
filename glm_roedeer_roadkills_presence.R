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
library(terra)
library(dplyr)

###############################################################################
# Filepaths to export and use data
###############################################################################
rasterExportPath <- "Export_R/final_rasters"

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
sdmFilePath = "forGLM/new_230125/species_distribution_sre_5000_pa1_run1_final.tif"
waterFilePath = "forGLM/landcover/europe/clipped_30km/WAW_clip.tif"
roadTypeFilePath = "variables_glm_clipped_ST_minus_NP/RoadType.tif"
humansFilePath = "forGLM/landcover/europe/clipped_30km/CLC_clip.tif"
# predictor variables - topography
demFilePath = "variables_glm_clipped_ST_minus_NP/dem_100m.tif"
aspectFilePath = "variables_glm_clipped_ST_minus_NP/aspect_100m.tif"
slopeFilePath = "variables_glm_clipped_ST_minus_NP/slope_100m.tif"
# new variables
RD_densityFilePath = "forGLM/new_230125/density_roedeer_r1km.tif"
forestFilePath = "forGLM/landcover/europe/clipped_30km/FTY_clip.tif"
roadNetworkFilePath = "forGLM/new_230125/road_network_total_rast.tif"
road_densityFilePath = "forGLM/new_230125/road_density_multi1000.tif"



# Extent files (South Tyrol - NP, road buffer)
ST_borderFilePath = "GIS/border_southTyrol_withoutNP.shp"
#roads_buffer_vector_FilePath ="GIS/100mBuffer_roads_clip.shp"

# Rasterize the buffered relevant roads
# Load your vector data
roads_buf_vect <- vect("forGLM/new_230125/relevant_roads_buf20m.shp")

# Define a raster template
template <- rast(ext(roads_buf_vect), resolution = 100, crs = "EPSG:25832")

# Rasterize the vector
roads_buf_rast <- rasterize(roads_buf_vect, template, touches = TRUE)
roads_buf_rast[is.na(roads_buf_rast)] <- 0

plot(roads_buf_rast)


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

#roads_buf_vect <- vect(roads_buffer_vector_FilePath)

RD_density <- rast(RD_densityFilePath)

forest_unmasked <- rast(forestFilePath)
forest <- app(forest_unmasked, function(x) ifelse(x %in% c(1, 2, 3), x, 0))

roadNetwork_unmasked <- rast(roadNetworkFilePath)
roadNetwork <- app(roadNetwork_unmasked, function(x) ifelse(is.na(x), 0, x))

roadDensity <- rast(road_densityFilePath)

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
RD_density <- project(RD_density, target_crs)
forest <- project(forest, target_crs)
roadNetwork <- project(roadNetwork, target_crs)
roadDensity <- project(roadDensity, target_crs)
roadNetwork_buf <- project(roads_buf_rast, target_crs)

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
RD_density <- resample(RD_density, template_raster, method = "bilinear")
forest <- resample(forest,template_raster, method = "near")
roadNetwork <- resample(roadNetwork,template_raster, method = "near")
roadDensity <- resample(roadDensity,template_raster, method = "bilinear")
roadNetwork_buf <-resample(roadNetwork_buf,template_raster, method = "near")

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

forest[forest==0] <- NA
distance_forest <- distance(forest)

# # ----------------------------------------------------------------------------
# #                               (C)  - extent
#
#start with setting the vector layer of the roads buffer as extent we want to crop the other layers to
# roads_buf_vect <- vect(roads_buffer_vector_FilePath)
# crs(roads_buf_vect) <- target_crs
ST_border <- vect(ST_borderFilePath)
crs(ST_border) <- target_crs

roadKillsCount_crop <- mask(crop (roadKillsCount, ST_border), ST_border)
roadType_crop <- mask(crop (roadType, ST_border), ST_border)
sdm_crop <- mask(crop (sdm, ST_border), ST_border)
dem_crop <- mask(crop (dem, ST_border), ST_border)
aspect_crop <- mask(crop (aspect, ST_border), ST_border)
slope_crop <- mask(crop (slope, ST_border), ST_border)
transitionalLand_crop <- mask(crop(transitionalLand, ST_border), ST_border)
RD_density_crop <- mask(crop(RD_density, ST_border), ST_border)
roadNetwork_crop <- mask(crop(roadNetwork, ST_border), ST_border)
roadDensity_crop <- mask(crop(roadDensity, ST_border), ST_border)
roadNetwork_buf_crop <- mask(crop(roadNetwork_buf, ST_border), ST_border)

distance_humanInfluence_crop<- mask(crop (distance_humanInfluence, ST_border), ST_border)
distance_broadleave_crop <- mask(crop (distance_broadleave, ST_border), ST_border)
distance_coniferous_crop <- mask(crop (distance_coniferous, ST_border), ST_border)
distance_mixedForest_crop <- mask(crop (distance_mixedForest, ST_border), ST_border)
distance_pastures_crop <- mask(crop (distance_pastures, ST_border), ST_border)
distance_grasslands_crop <- mask(crop (distance_grasslands, ST_border), ST_border)
distance_water_crop <- mask(crop (distance_water, ST_border), ST_border)
distance_forest_crop <- mask(crop (distance_forest, ST_border), ST_border)

# ----------------------------------------------------------------------------
# save all layers and in the future use them!
# ----------------------------------------------------------------------------

writeRaster(roadKillsCount_crop, file.path(rasterExportPath, "roadKillsCount.tif"), overwrite = TRUE)
writeRaster(roadType_crop, file.path(rasterExportPath, "roadType.tif"), overwrite = TRUE)
writeRaster(sdm_crop, file.path(rasterExportPath, "sdm.tif"), overwrite = TRUE)
writeRaster(dem_crop, file.path(rasterExportPath, "dem.tif"), overwrite = TRUE)
writeRaster(aspect_crop, file.path(rasterExportPath, "aspect.tif"), overwrite = TRUE)
writeRaster(slope_crop, file.path(rasterExportPath, "slope.tif"), overwrite = TRUE)
writeRaster(transitionalLand_crop, file.path(rasterExportPath, "transitionalLand.tif"), overwrite = TRUE)
writeRaster(RD_density_crop, file.path(rasterExportPath, "RD_density.tif"), overwrite = TRUE)
writeRaster(roadNetwork_crop, file.path(rasterExportPath, "roadNetwork.tif"), overwrite = TRUE)
writeRaster(roadDensity_crop, file.path(rasterExportPath, "roadDensity.tif"), overwrite = TRUE)
writeRaster(roadNetwork_buf_crop, file.path(rasterExportPath, "roadNetworkBuffered.tif"), overwrite = TRUE)

# write distance rasters
writeRaster(distance_humanInfluence_crop, file.path(rasterExportPath, "distance_humanInfluence.tif"), overwrite = TRUE)
writeRaster(distance_broadleave_crop, file.path(rasterExportPath, "distance_broadleave.tif"), overwrite = TRUE)
writeRaster(distance_coniferous_crop, file.path(rasterExportPath, "distance_coniferous.tif"), overwrite = TRUE)
writeRaster(distance_mixedForest_crop, file.path(rasterExportPath, "distance_mixedForest.tif"), overwrite = TRUE)
writeRaster(distance_pastures_crop, file.path(rasterExportPath, "distance_pastures.tif"), overwrite = TRUE)
writeRaster(distance_grasslands_crop, file.path(rasterExportPath, "distance_grasslands.tif"), overwrite = TRUE)
writeRaster(distance_water_crop, file.path(rasterExportPath, "distance_water.tif"), overwrite = TRUE)
writeRaster(distance_forest_crop, file.path(rasterExportPath, "distance_forest.tif"), overwrite = TRUE)

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
transitionalLand_rast <- rast(file.path(rasterExportPath, "transitionalLand.tif"))
RD_density_rast <- rast(file.path(rasterExportPath, "RD_density.tif"))
roadNetwork_rast <- rast(file.path(rasterExportPath, "roadNetwork.tif"))
roadDensity_rast <- rast(file.path(rasterExportPath, "roadDensity.tif"))
roadNetwork_buf_rast <- rast(file.path(rasterExportPath, "roadNetworkBuffered.tif"))

# all DISTANCE files
broadleave_dist <- rast(file.path(rasterExportPath, "distance_broadleave.tif"))
coniferous_dist <- rast(file.path(rasterExportPath, "distance_coniferous.tif"))
mixedForest_dist <- rast(file.path(rasterExportPath, "distance_mixedForest.tif"))
pastures_dist <- rast(file.path(rasterExportPath, "distance_pastures.tif"))
grasslands_dist <- rast(file.path(rasterExportPath, "distance_grasslands.tif"))
water_dist <- rast(file.path(rasterExportPath, "distance_water.tif"))
humanInfluence_dist <- rast(file.path(rasterExportPath, "distance_humanInfluence.tif"))
forest_dist <- rast(file.path(rasterExportPath, "distance_forest.tif"))

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
transitionalLand_values <- values(transitionalLand_rast)
RD_density_values<- values(RD_density_rast)
roadNetwork_values <- values(roadNetwork_rast)
roadDensity_values <- values(roadDensity_rast)
roadNetwork_buf_values <- values(roadNetwork_buf_rast)

####### EXTRACT DISTANCE VALUES

broadleave_dist_values <- values(broadleave_dist)
coniferous_dist_values <- values(coniferous_dist)
mixedForest_dist_values <- values(mixedForest_dist)

pastures_dist_values <- values(pastures_dist)
grasslands_dist_values <- values(grasslands_dist)
water_dist_values <- values(water_dist)
humanInfluence_dist_values <- values(humanInfluence_dist)
forest_dist_values <- values(forest_dist)

# Create the factor with appropriate levels
roadType_factor <- factor(roadType_values, 
                          levels = c("0", "1", "2", "3", "4", "5"),
                          labels = c("Other", "Highway", "State Road", 
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
                      dem_values, aspect_values, slope_values, RD_density_values,
                      roadNetwork_values,forest_dist_values, roadDensity_values)
# Assign column names
colnames(df_dist) <- c("Road_kills_Count", "Broadleave", "Coniferous", 
                  "Mixed_Forest", "Small_Woody_Features", "Pastures",
                  "Grasslands", "Road_Type", "Human_Influence", "Water", "SDM",
                  "DEM", "aspect", "slope","Roe_deer_density", "Road_Network", 
                  "Forest", "Road_density")

# Inspect the first few rows
head(df_dist)

# Remove rows where "Road kills Count" is NA or other variable is 
df_dist_na <- na.omit(df_dist)
head(df_dist_na)
#df_roads <- subset(df_dist_na, Road_Type != "No Road")
#df_roads$Road_Type <- droplevels(df_roads$Road_Type)
df_roads <- subset(df_dist_na, df_dist_na$Road_kills_Count != 0)
head(df_roads)

df_roads_yes <- subset(df_roads, df_roads$Road_Network != 0) # only all the kills on Roads (exclude unclear coordinated)
df_roads_clean <- subset(df_roads_yes, select = - Road_Network)

# transform variables
# List of variables to transform
variables_to_tf <- c("Broadleave", "Coniferous", 
                      "Mixed_Forest", "Small_Woody_Features", "Pastures",
                      "Grasslands", "Human_Influence", "Road_density", "SDM", "slope", 
                      "Forest", "Water", "Roe_deer_density")

# Check which variables in `variables_to_log` have values < 1
variables_with_values_below_1 <- sapply(variables_to_tf, function(var) {
  if (var %in% names(df_roads_clean)) { # Ensure the variable exists in the dataframe
    any(df_roads_clean[[var]] < 1, na.rm = TRUE)
  } else {
    FALSE # If the variable is not in the dataframe
  }
})

# Output variables with values < 1
variables_to_tf[variables_with_values_below_1]

# Transformed data frame
df_transformed <- df_roads_clean %>%
  mutate(across(all_of(variables_to_tf), ~ log(. + 3), .names = "log_{.col}"))

# Identify the variables not in variables_to_log
variables_not_transformed <- setdiff(names(df_roads_clean), 
                                     c(variables_to_tf, paste0("log_", variables_to_tf)))

# Select the untransformed variables
df_untransformed <- df_roads_clean %>%
  select(all_of(variables_not_transformed))

# Combine the untransformed variables with the transformed variables
df_tf <- bind_cols(df_transformed, df_untransformed)

# Check the structure of the final dataframe
str(df_tf)


# ------------------------- Make simple plots ----------------------------------
####                Boxplots for binary/factorial variables                 ####

par(mfrow = c(1,1))
boxplot(Road_kills_Count ~ Road_Type,
        data = df_roads_clean,
        # change axes labels
        xlab = "Road type",
        ylab = "Road kill count")

# Count the frequency of each level in Road_Type
road_type_counts <- table(df_roads_clean$Road_Type)

# Create a bar plot of the counts
barplot(road_type_counts, main = "Road Type Distribution", col = "skyblue", 
        xlab = "Road Type", ylab = "Count")

            ####                  Histogram for all numerical variables                     ####
par(mfrow = c(2,2))

lapply(names(df_roads_clean)[sapply(df_roads_clean, is.numeric)], function(var) {
  hist(df_roads_clean[[var]], main = paste("Histogram of", var), xlab = var, col = "lightblue", breaks = 30)
})

# transformed df
lapply(names(df_tf)[sapply(df_tf, is.numeric)], function(var) {
  hist(df_tf[[var]], main = paste("Histogram of", var), xlab = var, col = "lightblue", breaks = 30)
})

df_tf$sdm_tf_squared <- df_roads_clean$SDM^2
hist(df_tf$sdm_tf_squared)

########          Scatterplots
# List of numeric variables (excluding the response variable)
numeric_vars <- names(df_roads_clean)[sapply(df_roads_clean, is.numeric)]

# Create scatterplots for each numeric predictor against the response variable
lapply(numeric_vars, function(var) {
  plot(df_roads_clean[[var]], df_roads_clean$Road_kills_Count, main = paste("Scatterplot of", var, "vs Response"), 
       xlab = var, ylab = "Response Variable", col = "blue", pch = 16)
})

df_final <- bind_cols (df_roads_clean %>% select(-Roe_deer_density, -Road_density), df_tf %>% select (log_Roe_deer_density, log_Road_density)) 
str(df_final)
# ----------------------------------------------------------------------------
#                              Check assumptions
# ----------------------------------------------------------------------------
# Activate the libraries needed for this 
library(car) # needed for the Anova() function, contains function vif(), which offers an easy alternative way to check predictor independence.

# ----------------- (A) - Independence of predictors ---------------------------
# To find out: check by pairwise "correlations" among predictor values
library(corrplot) 

df_final_numeric <- df_final[sapply(df_final, is.numeric)]
# df_roads_numeric <- subset(df_roads_numeric, select = - Road_density)
cor_matrix <- cor(df_final_numeric, method = "spearman", use = "pairwise.complete.obs") #spearman correltion to test for correlation
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
         tl.cex = 0.7,                         # Adjust size of text labels
         title = "Spearman Correlation Matrix", # Add a title
         mar = c(0, 0, 2, 0),                  # Adjust margins for the title
         diag = FALSE                          # Do not show the diagonal
)
# Find variable pairs with |correlation| > 0.5
high_corr_pairs <- which(abs(cor_matrix) > 0.5, arr.ind = TRUE)
# Exclude self-correlations
high_corr_pairs <- high_corr_pairs[high_corr_pairs[, 1] != high_corr_pairs[, 2], ]
# Display highly correlated pairs
print(high_corr_pairs) 
# Roe deer density and SDM correlate with 0.6, remove roe deer density

# ----------------- (B) - Dispersion parameter ---------------------------
# To obtain model validation plots, We first need to calculate the glm
# ... WITHOUT checking its statistical results!

# Model WITHOUT roe deer density (Correlation)
mod1.poi <- glm(Road_kills_Count ~  Broadleave + Coniferous + Mixed_Forest + 
                  Small_Woody_Features + Pastures + Grasslands + Water + 
                  Human_Influence + Road_Type + DEM + aspect + slope+ 
                  Forest +SDM + Road_density, 
                data = df_roads_clean, 
                family = poisson(link = "log"),
                na.action = na.fail)

# the variance inflation factor VIF
# ... offers an easy alternative way to check for predictor independence.
# Without further explanation: Kick out predictors that have a value > 5 in the following list:
# ... They are too heavily correlated with other predictors.
vif(mod1.poi) # all looks good --> Predictors are independent

# model validation: residual plots only partially informative 
par(mfrow = c(2,2))
plot(mod1.poi) 

# --> MODEL VALIDATION with dispersion parameter
# What the model believes:
# mean value = data variance --> obs. / exp. dispersion = 1
# True dispersion in the data: Deviance / df
summary(mod1.poi)
dispersion.parameter <- mod1.poi$deviance / mod1.poi$df.residual
dispersion.parameter # UNDERDISPERSION

mod.poi.tf <- glm(Road_kills_Count ~  Broadleave + Coniferous + Mixed_Forest + 
                    Small_Woody_Features + Pastures + Grasslands + Water + 
                    Human_Influence + Road_Type + DEM + aspect + slope+ 
                    Forest +SDM + log_Road_density, 
                  data = df_final, 
                  family = poisson(link = "log"),
                  na.action = na.fail)

vif(mod.poi.tf) # all looks good --> Predictors are independent

# model validation: residual plots only partially informative 
# ...many road kill counts are 0 or close to 0
par(mfrow = c(2,2))
plot(mod.poi.tf) 

summary(mod.poi.tf)
dispersion.parameter.tf <- mod.poi.tf$deviance / mod.poi.tf$df.residual
dispersion.parameter.tf # UNDERDISPERSION

Anova(mod.poi.tf, type = "III")

# Let's try out quasipoisson

mod2.quasipoi <- glm(Road_kills_Count ~  Broadleave + Coniferous + Mixed_Forest + 
                  Small_Woody_Features + Pastures + Grasslands + Water + 
                  Human_Influence + Road_Type + DEM + aspect + slope + 
                  Roe_deer_density + Forest +SDM + Road_density,             # SDM HINZUFÜGEN
                data = df_roads_clean, 
                family = quasipoisson(link = "log"))

summary(mod2.quasipoi)


#### MODEL SELECTION
# Extract p-values
p_values <- summary(mod1.poi)$coefficients[, "Pr(>|z|)"]

# Select predictors with p < 0.1 (excluding the intercept)
selected_predictors <- names(p_values)[p_values < 0.1 & names(p_values) != "(Intercept)"]
# Ensure 'Road_Type' is a factor and not split into dummy variables
selected_predictors <- setdiff(selected_predictors, grep("Road_Type", selected_predictors, value = TRUE))
selected_predictors <- c(selected_predictors, "Road_Type")  # Add the factor variable itself

# Print selected predictors
print(selected_predictors)

# Construct the formula for the reduced model
reduced_formula <- as.formula(paste("Road_kills_Count ~", paste(selected_predictors, collapse = " + ")))

mod3.poi.reduced <- glm(reduced_formula,
                                 data = df_roads_clean,
                                 family = poisson(link = "log"),
                                 na.action = na.fail)


summary(mod3.poi.reduced)
dispersion.parameter.red <- mod3.poi.reduced$deviance / mod3.poi.reduced$df.residual
dispersion.parameter.red # UNDERDISPERSION

## Conclusio: Stick with poisson distribution to select a model (estimates are
# more conservative, that should be ok)


# As previously, we CAN calculate MacFadden's Pseudo-R2.
# !! NOVELTY: #####
# Pseudo-R2 CANNOT anymore be interpreted as "proportion of explained variance"
# ... because this requires residual normality (which we don't have in non-gaussian models)
# Still, the value indicates how "close" you get from a null model towards a model that explains ALL variation in the data.
# Hence, pseudo-R2 = 0.19 means a "more informative" model than pseudo-R2 = 0.12
macFadden.pseudoR2.mod1 <- 1 - mod1.poi$deviance / mod1.poi$null.deviance
macFadden.pseudoR2.mod1

macFadden.pseudoR2.mod1tf <- 1 - mod.poi.tf$deviance / mod.poi.tf$null.deviance
macFadden.pseudoR2.mod1tf

macFadden.pseudoR2.mod3 <- 1 - mod3.poi.reduced$deviance / mod3.poi.reduced$null.deviance
macFadden.pseudoR2.mod3

# ----------------- (D) - check Multicollinearity ------------------------------
# see above with the vif() function

# # ----------------------------------------------------------------------------
# #                                 Model selection
# # ----------------------------------------------------------------------------
library(MASS)
# stepwise selection
step_model <- stepAIC(mod1.poi, direction ="both", trace = TRUE)
summary(step_model)

dispersion.parameter.step <- step_model$deviance / step_model$df.residual
dispersion.parameter.step # UNDERDISPERSION

macFadden.pseudoR2.step <- 1 - step_model$deviance / step_model$null.deviance
macFadden.pseudoR2.step

Anova(step_model, type = "III")

library(MuMIn)
# selection by calculation of AIC for all combinations of predictors
# Generate all model subsets
model_set <- dredge(mod1.poi, rank = "AIC")

# View best models
head(model_set)

# Get the best model
best_model <- get.models(model_set, subset = 1)[[1]]

summary(best_model)

Anova(best_model, type = "III") 


# Zero truncated model

# Install pscl package if you haven't already
# install.packages("pscl")

# Load the pscl package
library(pscl)
library(VGAM) 

# Fit the Zero-Truncated Poisson model
mod1.ztp <-vglm(formula = Road_kills_Count ~ Broadleave + Coniferous + Mixed_Forest + 
                        Small_Woody_Features + Pastures + Grasslands + Water + 
                        Human_Influence + Road_Type + DEM + aspect + slope + 
                        Forest + SDM + Road_density, 
                      data = df_roads_clean, 
                family = pospoisson())

# Check the model summary
summary(mod1.ztp)

# Extract deviance and degrees of freedom
deviance <- deviance(mod1.ztp)
df <- df.residual(mod1.ztp)

# Calculate the dispersion parameter
dispersion <- deviance / df
dispersion
