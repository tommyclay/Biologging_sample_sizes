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

#### 2. LOADING IN EXAMPLE DATASETS ####

in_direx <- "./Data_inputs/gls_data/"

# grey-headed albatrosses tracked in 1999-2000
gha <- read.csv(paste0(in_direx, "gha99.csv"))

# white-chinned petrels tracked in 2003
wcp <- read.csv(paste0(in_direx, "wcp03.csv"))


#### 3. RESAMPLING AND PASTING OUT DATASETS ####

# set projections
# for gha, lambert azimuthal equal area projection centred around south pole
proj_gha <- "+proj=laea +lat_0=-90 +lon_0=-38.05 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# for gha, lambert azimuthal equal area projection centred around median lon and lat
proj_wcp <- paste0("+proj=laea +lat_0=", median(wcp$latitude), " +lon_0=", median(wcp$longitude), " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
# WGS84 projection
proj_dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"   

# resampling for 50% UD - for 500 iterations as an example (recommend using 1,000+ for final run)
# 50 km grid cell size and 186 km smoothing factor for GLS data (change based on type of biologging data)

set.seed(123) # for reproducibility 

# GHA
gha_res <- resampling_UD_area(x = gha, n_iter = 500, 
                             iso_level = 50, cellsize = 50000, h = 186000, proj_in = proj_dec,
                             proj_out = proj_gha, n_core = 6)
ggplot(gha_res, aes(sample, tot_area)) + geom_line(aes(group = factor(iteration)))


# WCP 
wcp_res <- resampling_UD_area(x = wcp, n_iter = 500, 
                                 iso_level = 50, cellsize = 50000, h = 186000, proj_in = proj_dec,
                                 proj_out = proj_wcp, n_core = 6)
ggplot(wcp_res, aes(sample, tot_area)) + geom_line(aes(group = factor(iteration)))


#### 4. RUN ASYMPTOTIC MODELS AND COMPARE AIC ####

# GHA 
gha_res$x <- gha_res$sample
gha_res$area <- gha_res$tot_area
gha_mm2 <- run_asymptotic_model(dat = gha_res, type = "MM2")
gha_mm3 <- run_asymptotic_model(dat = gha_res, type = "MM3")
gha_ar2 <- run_asymptotic_model(dat = gha_res, type = "AR2")
gha_ar3 <- run_asymptotic_model(dat = gha_res, type = "AR3")
# compare AICs
gha_aics <- AIC(gha_mm2, gha_mm3, gha_ar2, gha_ar3)
gha_aics$model <- c("MM2", "MM3", "AR2", "AR3")
gha_aics <- gha_aics[order(gha_aics$AIC),]
gha_aics$dAIC <- gha_aics$AIC-gha_aics$AIC[1]
gha_aics
#         df      AIC model      dAIC
# gha_mm3  4 330233.7   MM3   0.00000
# gha_mm2  3 330246.7   MM2  13.05636
# gha_ar3  4 330322.4   AR3  88.74593
# gha_ar2  3 330692.8   AR2 459.10652

# MM3 model is best

plot(gha_mm3)
# looks good

# WCP 
wcp_res$x <- wcp_res$sample
wcp_res$area <- wcp_res$tot_area
wcp_mm2 <- run_asymptotic_model(dat = wcp_res, type = "MM2")
wcp_mm3 <- run_asymptotic_model(dat = wcp_res, type = "MM3")
wcp_ar2 <- run_asymptotic_model(dat = wcp_res, type = "AR2")
wcp_ar3 <- run_asymptotic_model(dat = wcp_res, type = "AR3")
# compare AICs
wcp_aics <- AIC(wcp_mm2, wcp_mm3, wcp_ar2, wcp_ar3)
wcp_aics$model <- c("MM2", "MM3", "AR2", "AR3")
wcp_aics <- wcp_aics[order(wcp_aics$AIC),]
wcp_aics$dAIC <- wcp_aics$AIC-wcp_aics$AIC[1]
wcp_aics
#         df    AIC   model    dAIC
# wcp_mm3  4 138061.3   MM3  0.000000
# wcp_ar3  4 138066.3   AR3  5.025178
# wcp_ar2  3 138089.2   AR2 27.899487
# wcp_mm2  3 138098.0   MM2 36.695000

# MM3 model best

plot(wcp_mm3)
# looks good



#### 5. PREDICT MINIMUM SAMPLE SIZES FOR A GIVEN POPULATION SIZE ####

# just predicting to a population size of 100 for now

# GHA 
gha_pred <- sample_prediction(model = gha_mm3, pop_size = 100, percent_cover = 95)
gha_pred
#    prediction lower upper pop_size
#          19    18    20      100

# WCP
wcp_pred <- sample_prediction(model = wcp_mm3, pop_size = 100, percent_cover = 95)
wcp_pred
#    prediction lower upper pop_size
#          17    16    18      100