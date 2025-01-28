# #-------------------------------------------------------------------#
# # Course:
# # Topic: Road kill risk Roe deer South Tyrol
# # Script for glm
# # Author: Carla Behringer
# # Date: 24.01.2024
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
rasterExportPath <- "Export_R/rasters_binomial_model"

######
#############       SKIPPING OF CODE - S T A R T                    ##############
######
#
###############################################################################
# Filepaths for data
###############################################################################

# response variable
roadKillBinomialFilePath = "GIS/pseudoAbsence/presence_absence_5000_random_rast.tif"

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

# roadKills_vect <- vect(roadKillBinomialFilePath)
# # Extract coordinates using geom()
# coords <- geom(roadKills_vect)
# summary(coords)
# 
# # Create an empty raster with desired resolution and extent
# r <- rast(ext = ext(roadKills_vect), res = 100, crs="EPSG:25832")
# 
# roadKills <- rasterize(roadKills_vect, r, field="presence", fun = max)

roadKills <- rast(roadKillBinomialFilePath)

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

roadKills <- project(roadKills, target_crs)
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
ST_border <- vect(ST_borderFilePath)
template_raster <-rast(ext = ST_border, res = 100, crs = "EPSG:25832")  # Copy extent and CRS from the reference raster (roadKills)
res(template_raster) <- target_res  # Ensure resolution is set to 100x100 meters

roadKills <- resample(roadKills,template_raster, method = "near")
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

roadKills_crop <- mask(crop (roadKills, ST_border), ST_border)
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

writeRaster(roadKills_crop, file.path(rasterExportPath, "roadKills.tif"), overwrite = TRUE)
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
roadKills_rast <- rast(file.path(rasterExportPath, "roadKills.tif"))
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

roadKills_values <- as.integer(values(roadKills_rast))
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
                                     "Address Road"))

# Check the levels to confirm
levels(roadType_factor)


# Create data frame with continuous distance values
df_dist <- data.frame(roadKills_values, broadleave_dist_values,
                      coniferous_dist_values, mixedForest_dist_values, 
                      transitionalLand_values, pastures_dist_values, 
                      grasslands_dist_values, roadType_factor,
                      humanInfluence_dist_values, water_dist_values, sdm_values, 
                      dem_values, slope_values, RD_density_values,
                      roadNetwork_values,forest_dist_values, roadDensity_values)
# Assign column names
colnames(df_dist) <- c("Road_kills", "Broadleave", "Coniferous", 
                  "Mixed_Forest", "Small_Woody_Features", "Pastures",
                  "Grasslands", "Road_Type", "Human_Influence", "Water", "SDM",
                  "DEM", "slope","Roe_deer_density", "Road_Network", 
                  "Forest", "Road_density")

# Inspect the first few rows
head(df_dist)

# Remove rows where "Road kills Count" is NA or other variable is 
#df_dist_na <- na.omit(df_dist)
df_dist_na <- df_dist[!is.na(df_dist$Road_kill), ]

head(df_dist_na)

df_roads_yes <- subset(df_dist_na, df_dist_na$Road_Network != 0) # only all the kills on Roads (exclude unclear coordinated)
df_clean <- subset(df_roads_yes, select = - Road_Network)
sum(df_clean$Road_kill == 1)
sum(df_clean$Road_kill == 0)

# ------------------------- Make simple plots to check assumptions-------------

###### Assumption A: values are more or less leveled across predictor range

# List of continuous predictors
continuous_predictors <- c("Pastures", "Grasslands", "Water", "DEM", 
                            "slope", "Forest", "SDM", 
                           "Road_density", "Broadleave", "Coniferous", 
                           "Mixed_Forest", "Small_Woody_Features")

# Loop through each continuous predictor to create a boxplot
for (predictor in continuous_predictors) {
  boxplot(as.formula(paste(predictor, "~ Road_kills")), 
          data = df_clean, 
          main = paste("Boxplot of", predictor, "by Road_kills"), 
          xlab = "Road_kills", 
          ylab = predictor)
}
# ----------------------------------------------------------------------------
#                              Check assumptions
# ----------------------------------------------------------------------------
# Activate the libraries needed for this 
library(car) # needed for the Anova() function, contains function vif(), which offers an easy alternative way to check predictor independence.

# ----------------- (A) - Independence of predictors ---------------------------
# To find out: check by pairwise "correlations" among predictor values
library(corrplot) 

df_clean_numeric <- df_clean[sapply(df_clean, is.numeric)]
# df_roads_numeric <- subset(df_roads_numeric, select = - Road_density)
cor_matrix <- cor(df_clean_numeric, method = "spearman", use = "pairwise.complete.obs") #spearman correltion to test for correlation
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
# Find variable pairs with |correlation| > 0.6
high_corr_pairs <- which(abs(cor_matrix) > 0.6, arr.ind = TRUE)
# Exclude self-correlations
high_corr_pairs <- high_corr_pairs[high_corr_pairs[, 1] != high_corr_pairs[, 2], ]
# Display highly correlated pairs
print(high_corr_pairs) 
# Pairs with Spearman correlation r > 0.6:
# Distance to Forest – Distance to Coniferous Forest (0.65) -> removed distance to Forest
# Roe deer density – Species Distribution Model (0.67) -> removed roe deer density
# Road density – distance to Human Influence (-0.60) -> removed  Road density


# ----------------- (B) - Dispersion parameter ---------------------------
# To obtain model validation plots, We first need to calculate the glm
# ... WITHOUT checking its statistical results!

# Model WITHOUT roe deer density and coniferous forest (Correlation)
mod1.binom <- glm(Road_kills ~   Pastures + Grasslands + Water + 
                  Human_Influence + Road_Type + DEM + slope +SDM + 
                  Broadleave + Forest + Mixed_Forest + Small_Woody_Features
                  + Road_density,
                data = df_clean, 
                family = binomial(link = "logit"),
                na.action = na.fail)

# the variance inflation factor VIF
# ... offers an easy alternative way to check for predictor independence.
# Without further explanation: Kick out predictors that have a value > 5 in the following list:
# ... They are too heavily correlated with other predictors.
vif(mod1.binom) # all looks good --> Predictors are independent

# model validation: residual plots only partially informative 
par(mfrow = c(2,2))
plot(mod1.binom) 

macFadden.pseudoR2 <- 1 - mod1.binom$deviance / mod1.binom$null.deviance
macFadden.pseudoR2

# --> MODEL VALIDATION with dispersion parameter
# What the model believes:
# mean value = data variance --> obs. / exp. dispersion = 1
# True dispersion in the data: Deviance / df
library(MASS)

summary(mod1.binom)
dispersion.parameter <- mod1.binom$deviance / mod1.binom$df.residual
dispersion.parameter # 1.18

step_model <- stepAIC(mod1.binom, direction ="both", trace = TRUE)
summary(step_model)

dispersion.parameter.step <- step_model$deviance / step_model$df.residual
dispersion.parameter.step # UNDERDISPERSION

macFadden.pseudoR2.step <- 1 - step_model$deviance / step_model$null.deviance
macFadden.pseudoR2.step

Anova(step_model, type = "III")

mod3.forest <- glm(Road_kills ~   Pastures + Grasslands +
                    Human_Influence + Road_Type + DEM + slope +SDM + 
                    Forest+ Broadleave + Mixed_Forest+ Road_density,
                  data = df_clean, 
                  family = binomial(link = "logit"),
                  na.action = na.fail)
summary(mod3.forest)
macFadden.pseudoR2.forest <- 1 - mod3.forest$deviance / mod3.forest$null.deviance
macFadden.pseudoR2.forest

# -----------------------------------------------------------------------------#
#                                 GRAPHS
# -----------------------------------------------------------------------------#
# Required libraries
library(ggplot2)
library(ggeffects)
library(gridExtra)



effect_data <- ggpredict(step_model, terms = c("Grasslands [all]"))  # Adjust for your predictors
# Plot predicted effects
plot1 <- ggplot(effect_data, aes(x = x, y = predicted)) +
  geom_line(color = "orange") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "orange") +
  geom_point(data = df_clean, aes(x = Grasslands, y = Road_kills), 
             color = "grey", alpha = 0.5, size = 2) +  # Add actual data points
  labs(x = "Distance to Grasslands [m]", 
       y = "Predicted Roadkill Risk") +
  #title = "Effect of Distance to Grasslands on Road Kill Risk") +
  theme_bw() +  # Apply theme first
  theme(
    axis.title = element_text(size = 16),  # Increase axis labels font size
    axis.text = element_text(size = 14)    # Increase axis ticks font size
  )+
  scale_y_continuous(limits = c(0, 1))  # Then apply scale


# effect_data <- ggpredict(step_model, terms = c("Water [all]"))  # Adjust for your predictors
# # Plot predicted effects
# plot2 <-ggplot(effect_data, aes(x = x, y = predicted)) +
#   geom_line(color = "blue") +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill ="blue") +
#   geom_point(data = df_clean, aes(x = Water, y = Road_kills), 
#              color = "grey", alpha = 0.5, size = 2) +  # Add actual data points
#   labs(x = "Distance to Water [m]", 
#        y = "Predicted Roadkill Risk") +
#   #title = "Effect of Distance to Water on Road Kill Risk") +
#   theme_bw() +  # Apply theme first
#   theme(
#     axis.title = element_text(size = 16),  # Increase axis labels font size
#     axis.text = element_text(size = 14)    # Increase axis ticks font size
#   )+
#   scale_y_continuous(limits = c(0, 1))  # Then apply scale



effect_data <- ggpredict(step_model, terms = c("Human_Influence [all]"))  # Adjust for your predictors
# Plot predicted effects
plot3 <-ggplot(effect_data, aes(x = x, y = predicted)) +
  geom_line(color = "red") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill ="orange") +
  geom_point(data = df_clean, aes(x = Human_Influence, y = Road_kills), 
             color = "grey", alpha = 0.5, size = 2) +  # Add actual data points
  labs(x = "Distance to Human Settlements [m]", 
       y = "Predicted Roadkill Risk") +
  #title = "Effect of Distance to Human Settlements on Road Kill Risk") +
  theme_bw() +  # Apply theme first
  theme(
    axis.title = element_text(size = 16),  # Increase axis labels font size
    axis.text = element_text(size = 14)    # Increase axis ticks font size
  )+
  scale_y_continuous(limits = c(0, 1))  # Then apply scale


effect_data <- ggpredict(step_model, terms = c("DEM [all]"))  # Adjust for your predictors
# Plot predicted effects
plot4 <-ggplot(effect_data, aes(x = x, y = predicted)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill ="gray25") +
  geom_point(data = df_clean, aes(x = DEM, y = Road_kills), 
             color = "grey", alpha = 0.5, size = 2) +  # Add actual data points
  labs(x = "Elevation [m]", 
       y = "Predicted Roadkill Risk") +
  #title = "Effect of Elevation on Road Kill Risk") +
  theme_bw() +  # Apply theme first
  theme(
    axis.title = element_text(size = 16),  # Increase axis labels font size
    axis.text = element_text(size = 14)    # Increase axis ticks font size
  )+
  scale_y_continuous(limits = c(0, 1))  # Then apply scale


# effect_data <- ggpredict(step_model, terms = c("aspect [all]"))  # Adjust for your predictors
# # Plot predicted effects
# plot5 <-ggplot(effect_data, aes(x = x, y = predicted)) +
#   geom_line(color = "black") +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill ="gray25") +
#   geom_point(data = df_clean, aes(x = aspect, y = Road_kills), 
#              color = "grey", alpha = 0.5, size = 2) +  # Add actual data points
#   labs(x = "Aspect [°]", 
#        y = "Predicted Roadkill Risk") +
#   #title = "Effect of Aspect on Road Kill Risk") +
#   theme_bw() +  # Apply theme first
#   theme(
#     axis.title = element_text(size = 16),  # Increase axis labels font size
#     axis.text = element_text(size = 14)    # Increase axis ticks font size
#   )+
#   scale_y_continuous(limits = c(0, 1))  # Then apply scale


effect_data <- ggpredict(step_model, terms = c("slope [all]"))  # Adjust for your predictors
# Plot predicted effects
plot6 <-ggplot(effect_data, aes(x = x, y = predicted)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill ="gray25") +
  geom_point(data = df_clean, aes(x = slope, y = Road_kills), 
             color = "grey", alpha = 0.5, size = 2) +  # Add actual data points
  labs(x = "Slope [°]", 
       y = "Predicted Roadkill Risk") +
  #title = "Effect of Slope on Road Kill Risk") +
  theme_bw() +  # Apply theme first
  theme(
    axis.title = element_text(size = 16),  # Increase axis labels font size
    axis.text = element_text(size = 14)    # Increase axis ticks font size
  )+
  scale_y_continuous(limits = c(0, 1))  # Then apply scale


effect_data <- ggpredict(step_model, terms = c("SDM [all]"))  # Adjust for your predictors
# Plot predicted effects
plot7 <-ggplot(effect_data, aes(x = x, y = predicted)) +
  geom_line(color = "chocolate4") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill ="chocolate") +
  geom_point(data = df_clean, aes(x = SDM, y = Road_kills), 
             color = "grey", alpha = 0.5, size = 2) +  # Add actual data points
  labs(x = "Habitat suitability", 
       y = "Predicted Roadkill Risk") +
  #title = "Effect of Habitat suitability on Road Kill Risk") +
  theme_bw() +  # Apply theme first
  theme(
    axis.title = element_text(size = 16),  # Increase axis labels font size
    axis.text = element_text(size = 14)    # Increase axis ticks font size
  )+
  scale_y_continuous(limits = c(0, 1))  # Then apply scale

# effect_data <- ggpredict(step_model, terms = c("Broadleave [all]"))  # Adjust for your predictors
# # Plot predicted effects
# plot8 <-ggplot(effect_data, aes(x = x, y = predicted)) +
#   geom_line(color = "limegreen") +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill ="lawngreen") +
#   geom_point(data = df_clean, aes(x = Broadleave, y = Road_kills), 
#              color = "grey", alpha = 0.5, size = 2) +  # Add actual data points
#   labs(x = "Distance to Deciduous Forest [m]", 
#        y = "Predicted Roadkill Risk") +
#   #title = "Effect of Distance to Deciduous Forest on Road Kill Risk") +
#   theme_bw() +  # Apply theme first
#   theme(
#     axis.title = element_text(size = 16),  # Increase axis labels font size
#     axis.text = element_text(size = 14)    # Increase axis ticks font size
#   )+
#   scale_y_continuous(limits = c(0, 1))  # Then apply scale


effect_data <- ggpredict(step_model, terms = c("Forest [all]"))  # Adjust for your predictors
# Plot predicted effects
plot9 <-ggplot(effect_data, aes(x = x, y = predicted)) +
  geom_line(color = "forestgreen") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill ="yellowgreen") +
  geom_point(data = df_clean, aes(x = Coniferous, y = Road_kills), 
             color = "grey", alpha = 0.5, size = 2) +  # Add actual data points
  labs(x = "Distance to Forest [m]", 
       y = "Predicted Roadkill Risk")+
       #title = "Effect of Distance to Forest on Road Kill Risk") +
  theme_bw() +  # Apply theme first
  theme(
    axis.title = element_text(size = 16),  # Increase axis labels font size
    axis.text = element_text(size = 14)    # Increase axis ticks font size
  )+
  scale_y_continuous(limits = c(0, 1))  # Then apply scale


effect_data <- ggpredict(step_model, terms = c("Mixed_Forest [all]"))  # Adjust for your predictors
# Plot predicted effects
plot10 <-ggplot(effect_data, aes(x = x, y = predicted)) +
  geom_line(color = "olivedrab4") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill ="olivedrab4") +
  geom_point(data = df_clean, aes(x = Mixed_Forest, y = Road_kills), 
             color = "grey", alpha = 0.5, size = 2) +  # Add actual data points
  labs(x = "Distance to Mixed Forest [m]", 
       y = "Predicted Roadkill Risk") +
  #title = "Effect of Distance to Mixed Forest on Road Kill Risk") +
  theme_bw() +  # Apply theme first
  theme(
    axis.title = element_text(size = 16),  # Increase axis labels font size
    axis.text = element_text(size = 14)    # Increase axis ticks font size
  )+
  scale_y_continuous(limits = c(0, 1))  # Then apply scale

effect_data <- ggpredict(step_model, terms = c("Road_density [all]"))  # Adjust for your predictors
# Plot predicted effects
plot11 <-ggplot(effect_data, aes(x = x, y = predicted)) +
  geom_line(color = "deeppink4") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill ="deeppink3") +
  geom_point(data = df_clean, aes(x = Road_density, y = Road_kills), 
             color = "grey", alpha = 0.5, size = 2) +  # Add actual data points
  labs(x = expression("Road density [km/km"^2*"]"), 
       y = "Predicted Roadkill Risk") +
       #title = "Effect of Road density on Road Kill Risk") +
  theme_bw() +  # Apply theme first
  theme(
    axis.title = element_text(size = 16),  # Increase axis labels font size
    axis.text = element_text(size = 14)    # Increase axis ticks font size
  )+
  scale_y_continuous(limits = c(0, 1))  # Then apply scale

# Generate predicted values for Road_Type
effect_data <- ggpredict(step_model, terms = "Road_Type")
effect_data$x <- factor(effect_data$x, levels = c("Highway", "State Road","Country Road",  "Municipal Road", "Address Road", "Other"))

# Dot plot of predictions
dot_plot <- ggplot(effect_data, aes(x = x, y = predicted)) +
  geom_point(size = 3, color = "chocolate4") +  # Predicted values as dots
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, color = "chocolate4") +  # Confidence intervals
  labs(x = "Road Type", 
       y = "Predicted Roadkill Risk", 
       title = "Predicted Roadkill Risk by Road Type") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels
  )

violin_plot <- ggplot(df_clean, aes(x = Road_Type, y = Road_kills)) +
  geom_violin(fill = "darkblue", alpha = 0.5) +  # Violin for observed data
  geom_point(data = effect_data, aes(x = x, y = predicted), 
             color = "orange", size = 3, shape = 17, inherit.aes = FALSE) +  # Predicted probabilities
  geom_errorbar(data = effect_data, aes(x = x, ymin = conf.low, ymax = conf.high), 
                width = 0.2, color = "orange", inherit.aes = FALSE) +  # Confidence intervals
  labs(x = "Road Type", 
       y = "Roadkill Risk") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels
  )



all <- grid.arrange(plot1, plot3, plot4, plot6, plot7, plot9, plot10, plot11, 
             ncol = 3, nrow = 3)

vegetation <- grid.arrange(plot1, plot10, plot9,
             ncol = 3, nrow = 1)

topography <- grid.arrange(plot4, plot6,
             ncol = 2, nrow = 1)

habitatandhumans <- grid.arrange(plot7, plot3, plot11, 
             ncol = 3, nrow = 1)

anthro <- grid.arrange(plot3, plot11, dot_plot,
                       ncol = 3, nrow = 1)

# Save the grid to a file
export_pngPath <- "Export_R/figures/"

ggsave(paste0(export_pngPath, "allPredictors.png"), all, width = 15, height = 15, dpi = 300)
ggsave(paste0(export_pngPath, "vegetation.png"), vegetation, width = 15, height = 5, dpi = 300)
ggsave(paste0(export_pngPath, "topography.png"), topography, width = 10, height = 5, dpi = 300)
ggsave(paste0(export_pngPath, "habitatandhumans.png"), habitatandhumans, width = 15, height = 5, dpi = 300)
ggsave(paste0(export_pngPath, "sdm.png"), plot7, width = 7, height = 6, dpi = 300)
ggsave(paste0(export_pngPath, "anthro.png"), anthro, width = 15, height = 5, dpi = 300)
ggsave(paste0(export_pngPath, "roadtype.png"), dot_plot, width = 9, height = 6, dpi = 300)
ggsave(paste0(export_pngPath, "roadtype2.png"), violin_plot, width = 6, height = 6, dpi = 300)




