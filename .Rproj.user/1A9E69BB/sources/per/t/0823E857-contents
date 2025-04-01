library(sf)
library(terra)
library(dplyr)
library(adehabitatHR)
library(raster)

#### 1. LOAD FUNCTIONS #####

source("./Codes/1_predicting_minimum_sample_sizes/functions.R")

#### 2. LOAD IN GLS DATASETS ####

in_direx <- "./Data_inputs/gls_data/"
f <- dir(in_direx)

in_df  <- list()
for (i in 1:length(f)) {
  in_df[[i]] <- read.csv(paste0(in_direx, f[i]))
  in_df[[i]]$id <- as.factor(as.character(in_df[[i]]$id))
}

# set projections
proj_utm <- "+proj=laea +lat_0=-90 +lon_0=-38.05 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj_dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"   


#### 3. CALCULATE OVERLAP SCORES ####

out_df <- data.frame(dataset = gsub(".csv", "", f), species = NA, year = NA, ba = NA)
for (i in 1:length(f)) {
  print(gsub(".csv", "", f[i]))
  ba <- pairwise_BA_overlap(x = in_df[[i]], h = 186000, cellsize = 50000, iso_level = 50, proj_in = proj_dec, proj_out = proj_utm)
  out_df$ba[i] <- mean(ba) 
  out_df$species[i] <- in_df[[i]]$species[1]
  if (grepl("both", f[i])) {
    out_df$year[i] <- "BOTH" 
  } else if (grepl("all", f[i])) {
    out_df$year[i] <- "ALL" 
  } else {out_df$year[i] <- in_df[[i]]$year[1] }
}

# paste out
out_name <- "./Data_outputs/ba_overlap_scores.csv"
write.csv(out_df, out_name)
