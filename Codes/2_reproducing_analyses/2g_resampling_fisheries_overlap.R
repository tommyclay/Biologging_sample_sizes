library(plyr)
library(raster)
library(terra)
library(progress)
library(future.apply)
library(raster)
library(sf)
library(adehabitatHR)
library(future.apply)
library(ggplot2)

#### 1. LOAD FUNCTIONS #####

source("./Codes/1_predicting_minimum_sample_sizes/functions.R")

#### 2. LOAD IN FISHERIES DATA ####

fish_df <- read.csv("./Data_inputs/fishing_effort/gfw_pll_effort_2_degrees.csv")

#### 3. LOAD IN TRACKING DATA ####

in_gls <- "./Data_inputs/gls_data/"
bba <- read.csv(paste0(in_gls, "bba02.csv"))
gha <- read.csv(paste0(in_gls, "gha99.csv"))
wa <- read.csv(paste0(in_gls, "waall.csv"))
wcp <- read.csv(paste0(in_gls, "wcpboth.csv"))
bba$id <- as.factor(as.character(bba$id))
gha$id <- as.factor(as.character(gha$id))
wa$id <- as.factor(as.character(wa$id))
wcp$id <- as.factor(as.character(wcp$id))

# utm centred around lon of colony and lat of south pole
proj_utm <- paste0("+proj=laea +lat_0=", -90, " +lon_0=", -38.05, " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj_dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"   

#### 4. RUNNING RESAMPLING PROCESS AND PASTING OUT ####

out_direx <- "./Data_outputs/resampled_fisheries_overlap_outputs/"

##### 4A. BBA ######

bba_run <- resampling_fisheries_overlap(x = bba, n_iter = 1000, UD_cell = 50000, UD_h = 186000, 
                                                   threat_df = fish_df, grid_res = 2, n_cores = 4, 
                                                   proj_utm = proj_utm, proj_dec = proj_dec) 
# 7.8 mins for 1000 iterations/4 cores

# plotting output
ggplot(bba_run, aes(sample, overlap)) + geom_line(aes(group = iteration))
ggplot(bba_run, aes(sample, o_change_p)) + geom_line(aes(group = iteration))

# pasting out
out_name <- paste0(out_direx, "bba02_50.csv")
write.csv(bba_run, out.name)


##### 4B. GHA ######

gha_run <- resampling_fisheries_overlap(x = gha, n_iter = 1000, UD_cell = 50000, UD_h = 186000, 
                                         threat_df = fish_df, grid_res = 2, n_cores = 4, 
                                         proj_utm = proj_utm, proj_dec = proj_dec) 

# plotting output
ggplot(gha_run, aes(sample, overlap)) + geom_line(aes(group = iteration))
ggplot(gha_run, aes(sample, o_change_p)) + geom_line(aes(group = iteration))

# pasting out
out_name <- paste0(out_direx, "gha99_50.csv")
write.csv(gha_run, out.name)


##### 4C. WA ######

wa_run <- resampling_fisheries_overlap(x = wa, n_iter = 1000, UD_cell = 50000, UD_h = 186000, 
                                        threat_df = fish_df, grid_res = 2, n_cores = 4, 
                                        proj_utm = proj_utm, proj_dec = proj_dec) 

# plotting output
ggplot(wa_run, aes(sample, overlap)) + geom_line(aes(group = iteration))
ggplot(wa_run, aes(sample, o_change_p)) + geom_line(aes(group = iteration))

# pasting out
out_name <- paste0(out_direx, "waall_50.csv")
write.csv(wa_run, out.name)


##### 4D. WCP ######

wcp_run <- resampling_fisheries_overlap(x = wcp, n_iter = 1000, UD_cell = 50000, UD_h = 186000, 
                                       threat_df = fish_df, grid_res = 2, n_cores = 4, 
                                       proj_utm = proj_utm, proj_dec = proj_dec) 

# plotting output
ggplot(wcp_run, aes(sample, overlap)) + geom_line(aes(group = iteration))
ggplot(wcp_run, aes(sample, o_change_p)) + geom_line(aes(group = iteration))

# pasting out
out_name <- paste0(out_direx, "wcpboth_50.csv")
write.csv(wcp_run, out.name)

