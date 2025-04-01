library(adehabitatHR)
library(sf)
library(raster)
library(future.apply)
library(parallel)
library(ggplot2)
library(plyr)
library(dplyr)
library(pbapply)
library(compiler)
library(drc)

#### 1. LOAD FUNCTIONS #####

source("./Codes/1_predicting_minimum_sample_sizes/functions.R")

#### 2. LOADING IN GLS FILES ####

in_direx <- "./Data_inputs/gls_data/"
f <- dir(in_direx)
f <- f[-grep("apboth", f)]
f <- f[-grep("bsboth", f)]
f <- f[-grep("waall", f)]
f <- f[-grep("wcpboth", f)]

in_df  <- list()
for (i in 1:length(f)) {
  in_df[[i]] <- read.csv(paste0(in_direx, f[i]))
  in_df[[i]]$id <- as.factor(as.character(in_df[[i]]$id))
}

# set projections
proj_utm <- "+proj=laea +lat_0=-90 +lon_0=-38.05 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj_dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"   


#### 3. RESAMPLING AND PASTING OUT DATASETS - DO NOT PASTE OUT FOR NOW ####

# output directory
out_direx <- "./Data_outputs/resampled_area_outputs/gls_data/"

for (i in 1:length(f)) {
  print(paste(in_df[[i]]$species[1], in_df[[i]]$year[1], sep = "_"))
  # resampling 50% UD
  res_50 <- resampling_UD_area(x = in_df[[i]], n_iter = 1000, 
                                iso_level = 50, cellsize = 50000, h = 186000, proj_in = proj_dec,
                                proj_out = proj_utm, n_core = 6)
  # pasting out 50% UD iterations
  out_name <- paste0(out_direx, in_df[[i]]$species[1], substr(in_df[[i]]$year[1], 3, 4), "_", "50.csv")
  #write.csv(res_50, out_name)
  # resampling 95% UD
  res_95 <- resampling_UD_area(x = in_df[[i]], n_iter = 1000, 
                                iso_level = 95, cellsize = 50000, h = 186000, proj_in = proj_dec,
                                proj_out = proj_utm, n_core = 6)
  # pasting out 95% UD iterations
  out_name <- paste0(out_direx, in_df[[i]]$species[1], substr(in_df[[i]]$year[1], 3, 4), "_", "95.csv")
  #write.csv(res_95, out_name)
}


#### 4. RUN ASYMPTOTIC MODELS, COMPARE AIC AND PREDICT MINIMUM SAMPLE SIZES FOR A RANGE OF POP SIZES ####

in_direx <- "./Data_outputs/resampled_area_outputs/gls_data/"
out_direx <- "./Data_outputs/"

pops <- c(100, 1000, 10000, 100000, 1000000)


##### 4A. 50% UD #####

# load in resampled datasets 
f_50 <- dir(in_direx, "_50")

out_aics_50 <- list()
mod_mm2_50 <- list()
mod_mm3_50 <- list()
mod_ar2_50 <- list()
mod_ar3_50 <- list()
out_preds_50 <- list()
for (i in 1:length(f_50)) {
  print(gsub("_50.csv", "", f_50[i])) 
  dat <- read.csv(paste0(in_direx, f_50[i]))
  # running asymptotic model code
  print("running asymptotic models...")
  dat$x <- dat$sample
  dat$area <- dat$tot_area
  mod_mm2_50[[i]] <- run_asymptotic_model(dat = dat, type = "MM2")
  mod_mm3_50[[i]] <- run_asymptotic_model(dat = dat, type = "MM3")
  mod_ar2_50[[i]] <- run_asymptotic_model(dat = dat, type = "AR2")
  mod_ar3_50[[i]] <- run_asymptotic_model(dat = dat, type = "AR3")
  # compare AICs
  out_aics_50[[i]] <- AIC(mod_mm2_50[[i]], mod_mm3_50[[i]], mod_ar2_50[[i]], mod_ar3_50[[i]])
  out_aics_50[[i]]$model <- c("MM2", "MM3", "AR2", "AR3")
  out_aics_50[[i]] <- out_aics_50[[i]][order(out_aics_50[[i]]$AIC),]
  out_aics_50[[i]]$dAIC <- out_aics_50[[i]]$AIC-out_aics_50[[i]]$AIC[1]
  out_aics_50[[i]]$dataset <- gsub("_50.csv", "", f_50[i])
  # predicting minimum sample sizes
  print("predicting minimum sample sizes...")
  pred_mm2 <- sample_prediction(model = mod_mm2_50[[i]], pop_size = pops, percent_cover = 95)
  pred_mm2$model <- "MM2"
  pred_mm3 <- sample_prediction(model = mod_mm3_50[[i]], pop_size = pops, percent_cover = 95)
  pred_mm3$model <- "MM3"
  pred_ar2 <- sample_prediction(model = mod_ar2_50[[i]], pop_size = pops, percent_cover = 95)
  pred_ar2$model <- "AR2"
  pred_ar3 <- sample_prediction(model = mod_ar3_50[[i]], pop_size = pops, percent_cover = 95)
  pred_ar3$model <- "AR3"
  out_preds_50[[i]] <- rbind(pred_mm2, pred_mm3, pred_ar2, pred_ar3)
  out_preds_50[[i]]$dataset <- gsub("_50.csv", "", f_50[i])
}
all_aics_50 <- do.call(rbind, out_aics_50)
all_preds_50 <- do.call(rbind, out_preds_50)

# paste out AICs and minimum sample sizes - DON'T PASTE OUT FOR NOW
out_name <- paste0(out_direx, "area_occupied_50_aics.csv")
#write.csv(all_aics_50, out_name)
out_name <- paste0(out_direx, "area_occupied_50_predictions.csv")
#write.csv(all_preds_50, out_name)

# identifying best model for each dataset
aic_best <- all_aics_50 %>%
  group_by(dataset) %>%
  slice_min(AIC)
aic_best
nrow(aic_best[aic_best$model == "MM3",])
# 13 of 21 have MM3 as best model

# pasting out MM3 models for each dataset - - DON'T PASTE OUT FOR NOW
out_mod <- "./Data_outputs/resampled_area_outputs/gls_models/"
for (i in 1:length(f_50)) {
  out_name <- paste0(out_mod, gsub(".csv", "", f_50[i]), ".RData")
  mod <- mod_mm3_50[[i]]
  #save(mod, file = out_name)
}


##### 4B. 95% UD #####

# load in resampled datasets 
f_95 <- dir(in_direx, "_95")

out_aics_95 <- list()
mod_mm2_95 <- list()
mod_mm3_95 <- list()
mod_ar2_95 <- list()
mod_ar3_95 <- list()
out_preds_95 <- list()
for (i in 1:length(f_95)) {
  print( gsub("_95.csv", "", f_95[i])) 
  dat <- read.csv(paste0(in_direx, f_95[i]))
  # running asymptotic model code
  print("running asymptotic models...")
  dat$x <- dat$sample
  dat$area <- dat$tot_area
  mod_mm2_95[[i]] <- run_asymptotic_model(dat = dat, type = "MM2")
  mod_mm3_95[[i]] <- run_asymptotic_model(dat = dat, type = "MM3")
  mod_ar2_95[[i]] <- run_asymptotic_model(dat = dat, type = "AR2")
  mod_ar3_95[[i]] <- run_asymptotic_model(dat = dat, type = "AR3")
  # compare AICs
  out_aics_95[[i]] <- AIC(mod_mm2_95[[i]], mod_mm3_95[[i]], mod_ar2_95[[i]], mod_ar3_95[[i]])
  out_aics_95[[i]]$model <- c("MM2", "MM3", "AR2", "AR3")
  out_aics_95[[i]] <- out_aics_95[[i]][order(out_aics_95[[i]]$AIC),]
  out_aics_95[[i]]$dAIC <- out_aics_95[[i]]$AIC-out_aics_95[[i]]$AIC[1]
  out_aics_95[[i]]$dataset <- gsub("_95.csv", "", f_95[i])
  # predicting minimum sample sizes
  print("predicting minimum sample sizes...")
  pred_mm2 <- sample_prediction(model = mod_mm2_95[[i]], pop_size = pops, percent_cover = 95)
  pred_mm2$model <- "MM2"
  pred_mm3 <- sample_prediction(model = mod_mm3_95[[i]], pop_size = pops, percent_cover = 95)
  pred_mm3$model <- "MM3"
  pred_ar2 <- sample_prediction(model = mod_ar2_95[[i]], pop_size = pops, percent_cover = 95)
  pred_ar2$model <- "AR2"
  pred_ar3 <- sample_prediction(model = mod_ar3_95[[i]], pop_size = pops, percent_cover = 95)
  pred_ar3$model <- "AR3"
  out_preds_95[[i]] <- rbind(pred_mm2, pred_mm3, pred_ar2, pred_ar3)
  out_preds_95[[i]]$dataset <- gsub("_95.csv", "", f_95[i])
}
all_aics_95 <- do.call(rbind, out_aics_95)
all_preds_95 <- do.call(rbind, out_preds_95)

# paste out AICs and minimum sample sizes - DON'T PASTE OUT FOR NOW
out_name <- paste0(out_direx, "area_occupied_95_aics.csv")
#write.csv(all_aics_95, out_name)
out_name <- paste0(out_direx, "area_occupied_95_predictions.csv")
#write.csv(all_preds_95, out_name)

# identifying best model for each dataset
aic_best <- all_aics_95 %>%
  group_by(dataset) %>%
  slice_min(AIC)
aic_best
nrow(aic_best[aic_best$model == "MM3",])
# 18 of 21 have MM3 as best model
