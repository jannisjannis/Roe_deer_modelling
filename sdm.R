#install.packages("biomod2")
#install.packages("readr")
#install.packages("terra")
#install.packages("spThin")
#install.packages("tidyterra")
#install.packages("ggtext")
#install.packages("rio")
#install.packages("corrplot")
#install.packages("randomForest")
#install.packages("xgboost")
#install.packages("dismo")
#install.packages("nnet")
#install.packages("mgcv")
#install.packages("gbm")
#install.packages("gam")


setwd("C:/Users/Jannis/OneDrive - Scientific Network South Tyrol/Documents/Master - EMMA/3. Semester/Southtyrol-hunting data")

library(biomod2)
library(readr)
library(terra)
library(dplyr)
library(spThin)
library(ggplot2)
library(tidyterra)
library(ggtext)
library(rio) 
library(corrplot)
library(randomForest)
library(xgboost)
library(dismo)
library(nnet)  
library(mgcv)  
library(gbm)
library(gam)
#skip to line 151 from here

#### Worlclim raster layer creation ####
#create precipitation sum per year raster for southtyrol
rastersprecip <- rast(paste0("precip_raster/precip_clip_", 1:12, ".tif"))
annual_precip <- sum(rastersprecip)
#writeRaster(annual_precip, "annual_precip_southtyrol.tif", overwrite = TRUE)

#create precipitation sum WINTER raster for southtyrol
rastersprecip_WINTER <- rast(paste0("precip_raster/precip_clip_", c(12, 1, 2), ".tif"))
winter_precip <- sum(rastersprecip_WINTER)
#writeRaster(winter_precip, "precip_raster/winter_precip_southtyrol.tif", overwrite = TRUE)

#create mean annual temp raster for southtyrol
rasters <- rast(paste0("climate_raster/temp_clip_mean_", 1:12, ".tif"))
mean_annual_temp <- sum(rasters) / 12
#writeRaster(mean_annual_temp, "mean_annual_temp_southtyrol.tif", overwrite = TRUE)

#create mean WINTER (Dez, jan, feb) temp raster for southtyrol
rasters_WINTER <- rast(paste0("climate_raster/temp_clip_mean_", c(12, 1, 2), ".tif"))
mean_WINTER_temp <- sum(rasters_WINTER) / 3
#writeRaster(mean_WINTER_temp, "climate_raster/mean_winter_temp_southtyrol.tif", overwrite = TRUE)

#create Tmax raster for southtyrol
rasters_tmax <- rast(paste0("Tmax_raster/Tmax_mean_clipped_", 1:12, ".tif"))
Tmax_annual_temp <- sum(rasters_tmax) / 12
#writeRaster(Tmax_annual_temp, "Tmax_raster/Tmax_annual_temp_southtyrol.tif", overwrite = TRUE)

#create Tmin raster for southtyrol
rasters_tmin <- rast(paste0("Tmin_raster/Tmin_mean_clipped_", 1:12, ".tif"))
Tmin_annual_temp <- sum(rasters_tmin) / 12
#writeRaster(Tmin_annual_temp, "Tmin_raster/Tmin_annual_temp_southtyrol.tif", overwrite = TRUE)

#### clc raster creation ####
clc_raster <- rast("Layer/CLC_clip.tif")
# Extract cells with values 1, 2, 3, and 4
# Create a mask where these values are retained
desired_values <- c(1, 2, 3, 4)
clc_1234 <- ifel(clc_raster %in% desired_values, 1, 0)#writeRaster(clc_1234, "Layer/clc_1234.tif", overwrite = TRUE)
#plot(clc_1234)

#### set CRS ####
dem <- rast("dem_alps_100m.tif")
aspect <- rast("Layer/aspect.tif")
slope <- rast("Layer/slope.tif")
bio1 <-rast("Layer/bio1.tif")
bio2 <-rast("Layer/bio2.tif")
bio11 <-rast("Layer/bio11.tif")
bio12 <-rast("Layer/bio12.tif")
bio19 <-rast("Layer/bio19.tif")
forest_cover <- rast("Layer/FTY_2018_010m_03035_V1_0.tif")
grassland_cover <- rast("Layer/GRA_2018_010m_03035_V1_0.tif")
heat_map <- rast("Layer/heatmap_filled.tif")

#set EPSG:25832 as target crs into which we want to reproject all other layers
target_crs <- crs(dem)
#dem, aspect and slope are already in EPSG:25832
bio1 <- project(bio1, target_crs)
bio2 <- project(bio2, target_crs)
bio11 <- project(bio11, target_crs)
bio12 <- project(bio12, target_crs)
bio19 <- project(bio19, target_crs)
forest_cover <- project(forest_cover, target_crs)
grassland_cover <- project(grassland_cover, target_crs)
clc_cover <- project(clc_1234, target_crs)
crs(heat_map)

#### Match resolutions of the layers ####
#we start with downscaling the climate variables from 1km2 to 100m2
target_res <- 100
template_raster <- rast(dem)  # Copy extent and CRS from the reference raster (dem)
res(template_raster) <- target_res  # Ensure resolution is set to 100x100 meters

# Step 3: Resample bio1 to match the template raster
bio12_mm_per_year <- rast("Layer/bio12_mm_per_year.tif")
bio1_100m <- resample(bio1, template_raster, method = "bilinear")
bio2_100m <- resample(bio2, template_raster, method = "bilinear")
bio11_100m <- resample(bio11, template_raster, method = "bilinear")
bio12_100m <- resample(bio12, template_raster, method = "bilinear")
bio19_100m <- resample(bio19, template_raster, method = "bilinear")
heat_map_100m <- resample(heat_map, template_raster, method = "bilinear")
clc_100m <- resample(clc_cover, template_raster, method = "bilinear")

forest_cover_100m <- aggregate(forest_cover, fact = 10, fun = mean) #factor 10 due to grid size before is 10m
forest_cover_100m_aligned <- resample(forest_cover_100m, template_raster, method = "bilinear")

grassland_cover_100m <- aggregate(grassland_cover, fact = 10, fun = mean) #factor 10 due to grid size before is 10m
grassland_cover_100m_aligned <- resample(grassland_cover_100m, template_raster, method = "bilinear")

#### Cropping all layers to the same extent ####
#start with setting the vector layer of south tyrols border as extent we want to crop the other layers to
border_southtyrol <- vect("Layer/border_southTyrol_withoutNP.shp")

dem_crop <- mask(crop (dem, border_southtyrol), border_southtyrol)
aspect_crop <- mask(crop(aspect, border_southtyrol), border_southtyrol)
slope_crop <- mask(crop(slope, border_southtyrol), border_southtyrol)
bio1_crop <- mask(crop(bio1_100m, border_southtyrol), border_southtyrol)
bio2_crop <- mask(crop(bio2_100m, border_southtyrol), border_southtyrol)
bio11_crop <- mask(crop(bio11_100m, border_southtyrol), border_southtyrol)
bio12_crop <- mask(crop(bio12_100m, border_southtyrol), border_southtyrol)
bio19_crop <- mask(crop(bio19_100m, border_southtyrol), border_southtyrol)
forest_cover_crop <- mask(crop(forest_cover_100m_aligned, border_southtyrol), border_southtyrol)
grassland_cover_crop <- mask(crop(grassland_cover_100m_aligned, border_southtyrol), border_southtyrol)
heat_map_crop <- mask(crop(heat_map_100m, border_southtyrol), border_southtyrol)
clc_crop <- mask(crop(clc_100m, border_southtyrol), border_southtyrol)

bio12_mm_per_year_crop <- mask(crop(bio12_mm_per_year, border_southtyrol), border_southtyrol)
# due to the realignment the raster it not binary anymore. we make it binary again with values < 0.3 are 0
clc_crop_binary <- ifel(clc_crop < 0.5, 0, 1)

#### Save aligned rasters ####
writeRaster(dem_crop, "aligned_rasters/dem_100m.tif", overwrite = TRUE)
writeRaster(aspect_crop, "aligned_rasters/aspect_100m.tif", overwrite = TRUE)
writeRaster(slope_crop, "aligned_rasters/slope_100m.tif", overwrite = TRUE)
writeRaster(bio1_crop, "aligned_rasters/bio1_100m.tif", overwrite = TRUE)
writeRaster(bio2_crop, "aligned_rasters/bio2_100m.tif", overwrite = TRUE)
writeRaster(bio11_crop, "aligned_rasters/bio11_100m.tif", overwrite = TRUE)
writeRaster(bio12_crop, "aligned_rasters/bio12_100m.tif", overwrite = TRUE)
writeRaster(bio19_crop, "aligned_rasters/bio19_100m.tif", overwrite = TRUE)
writeRaster(forest_cover_crop, "aligned_rasters/forest_cover_100m.tif", overwrite = TRUE)
writeRaster(grassland_cover_crop, "aligned_rasters/grassland_cover_100m.tif", overwrite = TRUE)
writeRaster(heat_map_crop, "aligned_rasters/heat_map_100m.tif", overwrite = TRUE)
writeRaster(clc_crop_binary, "aligned_rasters/clc_binary_100m.tif", overwrite = TRUE)

writeRaster(bio12_mm_per_year_crop, "aligned_rasters/bio12_mm_per_year_crop.tif", overwrite = TRUE)

#### Read in Aligned rasters ####
dem <- rast("aligned_rasters/dem_100m.tif")
aspect <- rast("aligned_rasters/aspect_100m.tif")
slope <- rast("aligned_rasters/slope_100m.tif")
bio1 <-rast("aligned_rasters/bio1_100m.tif")
#bio2 <-rast("aligned_rasters/bio2_100m.tif")
bio11 <-rast("aligned_rasters/bio11_100m.tif")
bio12 <-rast("aligned_rasters/bio12_100m.tif")
bio19 <-rast("aligned_rasters/bio19_100m.tif")
forest_cover <- rast("aligned_rasters/forest_cover_100m.tif")
grassland_cover <- rast("aligned_rasters/grassland_cover_100m.tif")
heat_map <- rast("aligned_rasters/heat_map_100m.tif")
human_settlement <- rast("aligned_rasters/clc_binary_100m.tif")

#### Test for correlation ####
env_stack <- c(dem, aspect, slope, bio1, bio11, bio12, bio19, forest_cover, grassland_cover, human_settlement)

env_values <- as.data.frame(env_stack, na.rm = TRUE) #extract values from rasters
colors <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA")) #define colorpalette
colnames(env_values) <- c(
  "DEM", "Aspect", "Slope", "Bio1", "Bio11", "Bio12", "Bio19", 
  "ForestCover", "GrasslandCover", "HumanSettlement")

#Col names with correlated parameters removed
colnames(env_values) <- c(
  "Aspect", "Slope", "Bio11", "Bio12",
  "ForestCover", "GrasslandCover", "HeatMap", "HumanSettlement")

cor_matrix <- cor(env_values, method = "spearman")

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

# Find variable pairs with |correlation| > 0.8
high_corr_pairs <- which(abs(cor_matrix) > 0.6, arr.ind = TRUE)
# Exclude self-correlations
high_corr_pairs <- high_corr_pairs[high_corr_pairs[, 1] != high_corr_pairs[, 2], ]
# Display highly correlated pairs
print(high_corr_pairs)
# -> therefore exclude dem_alps_100m and bio1
env_stack <- env_stack[[names(env_stack) != "dem_alps_100m"]]
env_stack <- env_stack[[names(env_stack) != "bio1"]]
env_stack <- env_stack[[names(env_stack) != "bio19"]]
names(env_stack)

#### Read in presence data ####
presence_data <- import("Rehwilddaten_punktverortet.csv", header = TRUE)

colnames(presence_data)[colnames(presence_data) == "GPS (lat)"] <- "lat"
colnames(presence_data)[colnames(presence_data) == "GPS (lng)"] <- "lon"
colnames(presence_data)

#change from character to numeric the coordinates
presence_data$lat <- as.numeric(gsub(",", ".", presence_data$lat))
presence_data$lon <- as.numeric(gsub(",", ".", presence_data$lon))

#removing Nas
presence_data <- presence_data %>%
  filter(!is.na(lat), !is.na(lon))

#Change from old WGS 84 to EPSG:35832
#create  Vector from the old dataframe
presence_vect <- vect(presence_data, geom = c("lon", "lat"), crs = "EPSG:4326")  
presence_vect_25832 <- project(presence_vect, "EPSG:25832") # Transform to EPSG:25832
coords <- crds(presence_vect_25832) # extract transformed coordinates
# Add the transformed coordinates back to the data frame
presence_data$X_25832 <- coords[, 1]  # X coordinate
presence_data$y_25832 <- coords[, 2]  # Y coordinate

#Remove all the points that are outside our southtyrol boundary
border_southtyrol <- vect("Layer/border_southTyrol_withoutNP.shp")
presence_vect_25832 <- vect(presence_data, geom = c("X_25832", "y_25832"), crs = "EPSG:25832")  
presence_25832_within_southtyrol <- presence_vect_25832[border_southtyrol, ]
#crs(presence_25832_within_southtyrol)
clipped_data <- as.data.frame(presence_25832_within_southtyrol)
#add new coords back to the dataframe
# Ensure the new coordinates are included in the data frame
clipped_data$X_25832 <- crds(presence_25832_within_southtyrol)[, 1]  # X coordinate (EPSG:25832)
clipped_data$y_25832 <- crds(presence_25832_within_southtyrol)[, 2]  # Y coordinate (EPSG:25832)

#create table for biomod2 with just x,y coordinates, the presence indicator and the species name
presence_data_biomod <- clipped_data[, c("y_25832", "X_25832")]
presence_data_biomod$presence <- 1
presence_data_biomod$species <- "Roe_deer"
colnames(presence_data_biomod)[colnames(presence_data_biomod) == "y_25832"] <- "y"
colnames(presence_data_biomod)[colnames(presence_data_biomod) == "X_25832"] <- "x"
#export(presence_data_biomod, "Roedeer_within_ST_25832.csv")

#erstelle random subset von 5.000 punkten
set.seed(245) 
presence_data_biomod_subset_5000 <- presence_data_biomod[sample(nrow(presence_data_biomod), 5000), ]
export(presence_data_biomod_subset_5000, "presence_data_biomod_subset_5000.csv", overwrite = TRUE)

#random subset with 10.000 points
set.seed(513) 
presence_data_biomod_subset_10000 <- presence_data_biomod[sample(nrow(presence_data_biomod), 10000), ]
export(presence_data_biomod_subset_10000, "presence_data_biomod_subset_10000.csv", overwrite = TRUE)


#### Format data and generate pseudo-absences SRE 5000 points ####
myBiomodData_r_sre_5000 <- BIOMOD_FormatingData(
  resp.var = presence_data_biomod_subset_5000$presence,           # Presence data
  expl.var = env_stack,                        # Environmental variables
  resp.xy = presence_data_biomod_subset_5000[, c("x", "y")],      # Coordinates of presences
  resp.name = "Roe_deer",                      # Name of the species
  PA.nb.rep = 1,                               # Number of pseudo-absence replicates
  PA.nb.absences = 5000,                       # Number of pseudo-absences
  PA.strategy = "sre", # Strategy for generating pseudo-absences
  filter.raster = TRUE # Enable automatic filtering
)

#look at generated pseudo-absences
# Extract the pseudo-absence table
pseudo_absence_table_sre_5000 <- myBiomodData_r_sre_5000@PA.table
# Extract the coordinates (x, y) of pseudo-absences
coordinates_sre_5000 <- myBiomodData_r_sre_5000@coord
# Combine the pseudo-absence data with coordinates
pseudo_absence_data_sre_5000 <- cbind(coordinates_sre_5000, pseudo_absence_table_sre_5000)
write.csv(pseudo_absence_data_sre_5000, "pseudo_absences_sre_5000.csv", row.names = FALSE)
plot(myBiomodData_r_sre_5000)

#### Format data and generate pseudo-absences SRE 10.000 points ####
myBiomodData_r_sre_10000 <- BIOMOD_FormatingData(
  resp.var = presence_data_biomod_subset_10000$presence,           # Presence data
  expl.var = env_stack,                        # Environmental variables
  resp.xy = presence_data_biomod_subset_10000[, c("x", "y")],      # Coordinates of presences
  resp.name = "Roe_deer",                      # Name of the species
  PA.nb.rep = 1,                               # Number of pseudo-absence replicates
  PA.nb.absences = 10000,                       # Number of pseudo-absences
  PA.strategy = "sre", # Strategy for generating pseudo-absences
  filter.raster = TRUE # Enable automatic filtering
)

#look at generated pseudo-absences
# Extract the pseudo-absence table
pseudo_absence_table_sre_10000 <- myBiomodData_r_sre_10000@PA.table
# Extract the coordinates (x, y) of pseudo-absences
coordinates_sre_10000 <- myBiomodData_r_sre_10000@coord
# Combine the pseudo-absence data with coordinates
pseudo_absence_data_sre_10000 <- cbind(coordinates_sre_10000, pseudo_absence_table_sre_10000)
write.csv(pseudo_absence_data_sre_10000, "pseudo_absences_sre_10000.csv", row.names = FALSE)
plot(myBiomodData_r__10000)

#### Format data and generate pseudo-absences RANDOM 5000 points####
myBiomodData_r_random_5000 <- BIOMOD_FormatingData(
  resp.var = presence_data_biomod_subset_5000$presence,           # Presence data
  expl.var = env_stack,                        # Environmental variables
  resp.xy = presence_data_biomod_subset_5000[, c("x", "y")],      # Coordinates of presences
  resp.name = "Roe_deer",                      # Name of the species
  PA.nb.rep = 1,                               # Number of pseudo-absence replicates
  PA.nb.absences = 5000,                       # Number of pseudo-absences
  PA.strategy = "random", # Strategy for generating pseudo-absences
  filter.raster = TRUE # Enable automatic filtering
)

#look at generated pseudo-absences
# Extract the pseudo-absence table
pseudo_absence_table_random_5000 <- myBiomodData_r_random_5000@PA.table
# Extract the coordinates (x, y) of pseudo-absences
coordinates_random_5000 <- myBiomodData_r_random_5000@coord
# Combine the pseudo-absence data with coordinates
pseudo_absence_data_random_5000 <- cbind(coordinates_random_5000, pseudo_absence_table_random_5000)
write.csv(pseudo_absence_data_random_5000, "pseudo_absences_random_5000.csv", row.names = FALSE)
plot(myBiomodData_r_random)

#### Format data and generate pseudo-absences RANDOM 10.000 points####
myBiomodData_r_random_10000 <- BIOMOD_FormatingData(
  resp.var = presence_data_biomod_subset_10000$presence,           # Presence data
  expl.var = env_stack,                        # Environmental variables
  resp.xy = presence_data_biomod_subset_10000[, c("x", "y")],      # Coordinates of presences
  resp.name = "Roe_deer",                      # Name of the species
  PA.nb.rep = 1,                               # Number of pseudo-absence replicates
  PA.nb.absences = 10000,                       # Number of pseudo-absences
  PA.strategy = "random", # Strategy for generating pseudo-absences
  filter.raster = TRUE # Enable automatic filtering
)

#look at generated pseudo-absences
# Extract the pseudo-absence table
pseudo_absence_table_random_10000 <- myBiomodData_r_random_10000@PA.table
# Extract the coordinates (x, y) of pseudo-absences
coordinates_random_10000 <- myBiomodData_r_random_10000@coord
# Combine the pseudo-absence data with coordinates
pseudo_absence_data_random_10000 <- cbind(coordinates_random_10000, pseudo_absence_table_random_10000)
write.csv(pseudo_absence_data_random_10000, "pseudo_absences_random_10000.csv", row.names = FALSE)
plot(myBiomodData_r_random)

#### Thin Presence Dataset down from the 52.900 occurences -> not needed right now ####
#r <- rast(ext(border_southtyrol), resolution = 1000, crs = "EPSG:25832") #Set the extent and resolution for the raster grid
#r_points <- rasterize(presence_25832_within_southtyrol, r, fun = "first", background = NA) #Rasterize the points (assign each point to a grid cell)
#unique_points <- as.points(r_points, na.rm = TRUE) #Extract unique points based on the raster cells
#writeVector(unique_points, "Layer/thinned_occurence_data_1km.shp", overwrite = TRUE)
#thinned_occurence_data <- as.data.frame(crds(unique_points))
#thinned_occurence_data$presence <- 1
#thinned_occurence_data$species <- "Roe_deer"

#occur_thin <- thin(
  #loc.data = presence_data_biomod_subset,
  #lat.col = "x",
  #long.col = "y",
  #spec.col = "species",
  #thin.par = 1000,--> not 
  #reps = 5,
  #write.files = TRUE,
  #out.dir = "Layer/",
  # = "thinned_data_1_home"
#)

#### Modelling ####
#SRE_5000 modelling
biomod_model_sre5000 <- BIOMOD_Modeling(
  bm.format = myBiomodData_r_sre_5000,
  modeling.id = "Example",
  models = c('RF'),
  #models = c('RF', 'GLM', 'ANN', 'XGBOOST', 'SRE'),
  CV.strategy = 'random',
  CV.nb.rep = 3,
  CV.perc = 0.8,
  OPT.strategy = 'bigboss',
  metric.eval = c('TSS','ROC'),
  var.import = 2,
  seed.val = 42)

biomod_model_scores_sre5000 <- get_evaluations(biomod_model_sre5000)
dim(biomod_model_scores_sre5000)
dimnames(biomod_model_scores_sre5000) 
biomod_model_variable_importance_sre5000 <- get_variables_importance(biomod_model_sre5000)
bm_PlotEvalBoxplot(bm.out = biomod_model_sre5000, group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = biomod_model_sre5000, group.by = c('expl.var', 'algo', 'run'))

#SRE_10000 modelling
biomod_model_sre10000 <- BIOMOD_Modeling(
  bm.format = myBiomodData_r_sre_10000,
  modeling.id = "Example",
  models = c('RF'),
  CV.strategy = 'random',
  CV.nb.rep = 3,
  CV.perc = 0.8,
  OPT.strategy = 'bigboss',
  metric.eval = c('TSS','ROC'),
  var.import = 2,
  seed.val = 42)

biomod_model_scores_sre10000 <- get_evaluations(biomod_model_sre10000)
dim(biomod_model_scores_sre10000)
dimnames(biomod_model_scores_sre10000) 
biomod_model_variable_importance_sre10000 <- get_variables_importance(biomod_model_sre10000)
bm_PlotEvalBoxplot(bm.out = biomod_model_sre10000, group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = biomod_model_sre10000, group.by = c('expl.var', 'algo', 'run'))

#random_5000 modelling
biomod_model_random_5000 <- BIOMOD_Modeling(
  bm.format = myBiomodData_r_random_5000 ,
  modeling.id = "Example",
  models = c('RF'),
  CV.strategy = 'random',
  CV.nb.rep = 3,
  CV.perc = 0.8,
  OPT.strategy = 'bigboss',
  metric.eval = c('TSS','ROC'),
  var.import = 2,
  seed.val = 42)

biomod_model_scores_random_5000  <- get_evaluations(biomod_model_random_5000 )
dim(biomod_model_scores_random_5000 )
dimnames(biomod_model_scores_random_5000 ) 
biomod_model_variable_importance_random_5000  <- get_variables_importance(biomod_model_random_5000 )
bm_PlotEvalBoxplot(bm.out = biomod_model_random_5000 , group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = biomod_model_random_5000 , group.by = c('expl.var', 'algo', 'run'))

#random_10000 modelling
biomod_model_random_10000 <- BIOMOD_Modeling(
  bm.format = myBiomodData_r_random_10000 ,
  modeling.id = "Example",
  models = c('RF'),
  CV.strategy = 'random',
  CV.nb.rep = 3,
  CV.perc = 0.8,
  OPT.strategy = 'bigboss',
  metric.eval = c('TSS','ROC'),
  var.import = 2,
  seed.val = 42)

biomod_model_scores_random_10000  <- get_evaluations(biomod_model_random_10000)
dim(biomod_model_scores_random_10000)
dimnames(biomod_model_scores_random_10000) 
biomod_model_variable_importance_random_10000  <- get_variables_importance(biomod_model_random_10000)
bm_PlotEvalBoxplot(bm.out = biomod_model_random_10000 , group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = biomod_model_random_10000 , group.by = c('expl.var', 'algo', 'run'))


biomod_projection <- BIOMOD_Projection(
  bm.mod = biomod_model_sre5000,
  new.env = env_stack,        # Your environmental variable stack
  proj.name = "species_projection",
  selected.models = c("Roe.deer_PA1_Run3_RF"),    # You can specify models, e.g., RF or GLM
  binary.meth = "TSS",        # Use a binary method based on TSS or ROC
  compress = FALSE,
  clamping.mask = FALSE,
  output.format = ".img"      # You can choose GeoTIFF or other formats
)
    
proj_files <- get_predictions(biomod_projection)
names(proj_files)
rf_projection <- proj_files[["Roe.deer_PA1_RUN3_RF"]]
writeRaster(rf_projection, filename = "species_distribution_sre_5000_pa1_run3_final.tif", overwrite = TRUE)

# Corellation Heat_map and final model
corr_den_sdm <- as.data.frame(c(rf_projection, heat_map), na.rm = TRUE) #extract values from rasters
cor_matrix_sdm <- cor(corr_den_sdm, method = "spearman") #spearman correltion to test for correlation

corrplot(cor_matrix_sdm, 
         method = "color",                     # Use color tiles for visualization
         col = colors(200),                    # Apply your custom color palette
         type = "upper",                       # Show only the upper triangle
         order = "hclust",                     # Cluster variables based on correlations
         addCoef.col = "black",                # Add correlation coefficients in black
         number.cex = 0.7,                     # Adjust size of correlation numbers
         tl.cex = 0.5,                         # Adjust size of text labels
         title = "Spearman Correlation Matrix", # Add a title
         mar = c(0, 0, 2, 0),                  # Adjust margins for the title
         diag = FALSE                          # Do not show the diagonal
)


r <- rasterize(bio19)  # For terra (use raster() if using the raster package)

# Multiply by the factor
factor <- 3600 * 24 * 91.3 * 1000
bio19_new <- bio19 * factor

# Save the modified raster
writeRaster(bio19_new, "Layer/bio19_per3months.tif", overwrite=TRUE)
