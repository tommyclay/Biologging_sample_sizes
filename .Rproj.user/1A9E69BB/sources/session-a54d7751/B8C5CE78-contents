library(sf)
library(terra)
library(dplyr)
library(adehabitatHR)
library(raster)

#### 1. LOAD FUNCTIONS #####

source("./Codes/1_predicting_minimum_sample_sizes/functions.R")

#### 2. LOAD IN EXAMPLE DATASETS ####

in_direx <- "./Data_inputs/gls_data/"

# grey-headed albatrosses tracked in 1999-2000
gha <- read.csv(paste0(in_direx, "gha99.csv"))

# white-chinned petrels tracked in 2003
wcp <- read.csv(paste0(in_direx, "wcp03.csv"))

#### 3. CALCULATE AMONG-INDIVIDUAL OVERLAP ####

# set projections
# for gha, lambert azimuthal equal area projection centred around south pole
proj_gha <- "+proj=laea +lat_0=-90 +lon_0=-38.05 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# for gha, lambert azimuthal equal area projection centred around median lon and lat
proj_wcp <- paste0("+proj=laea +lat_0=", median(wcp$latitude), " +lon_0=", median(wcp$longitude), " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
# WGS84 projection
proj_dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"   

# GHA
gha_ba <- pairwise_BA_overlap(x = gha, h = 186000, cellsize = 50000, iso_level = 50, proj_in = proj_dec, proj_out = proj_gha)
gha_ba_avg <- mean(gha_ba) 
gha_ba_avg
# 0.2145811

# WCP
wcp_ba <- pairwise_BA_overlap(x = wcp, h = 186000, cellsize = 50000, iso_level = 50, proj_in = proj_dec, proj_out = proj_wcp)
wcp_ba_avg <- mean(wcp_ba) 
wcp_ba_avg
# 0.2262891

#### 4. PREDICT MINIMUM SAMPLE SIZE BASED ON RELATIONSHIP WITH BA ####

# GHA
prediction_from_ba(gha_ba_avg) 
#        ba prediction lower upper
# 0.2145811         17    15    18
# 17 individuals is slightly lower than that predicted using area occupied method = 19 (from 1a)

# WCP
prediction_from_ba(wcp_ba_avg) 
#       ba prediction lower upper
# 0.2262891         16    14    17
# 16 individuals is slightly lower than that predicted using area occupied method = 17 (from 1a)
