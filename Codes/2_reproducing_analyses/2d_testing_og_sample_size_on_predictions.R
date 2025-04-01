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
library(stringr)
library(lme4)
library(lmtest)


#### 1. LOAD FUNCTIONS #####

source("./Codes/1_predicting_minimum_sample_sizes/functions.R")

#### 2. LOAD IN DATA ####

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


#### 3. RANDOMLY SAMPLE DATASETS AT DIFFERENT SAMPLE SIZE INTERVALS AND RUN RESAMPLING CODE  - DON'T PASTE OUT FOR NOW ####

# set seed
set.seed(512)

# output directory
out_direx <- "./Data_outputs/sensitivity_analysis_outputs/"

# this takes a really long time!

for (i in 1:length(f)) {
  # set number of levels to resample at
  n_sample <- seq(5, length(unique(in_df[[i]]$id)), 5)
  # if the max n_sample is the same as the max number of individuals, then remove
  if(max(n_sample) == length(unique(in_df[[i]]$id))) { 
    n_sample = n_sample[-length(n_sample)]
  } else {}
  # iterate through levels of 5
  for (a in 1:length(n_sample)) {
    print(paste(paste(in_df[[i]]$species[1], in_df[[i]]$year[1], sep = "_"), n_sample[a], sep = "..."))
    # subsetting random individuals
    random <- sample(levels(in_df[[i]]$id), n_sample[a], replace = FALSE)
    dat_sub <- subset(in_df[[i]], id %in% random,)
    dat_sub$id <- as.factor(as.character(dat_sub$id))
    # resampling 
    res_out <- resampling_UD_area(x = dat_sub, n_iter = 1000, 
                                   iso_level = 50, cellsize = 50000, h = 186000, proj_in = proj_dec,
                                   proj_out = proj_utm, n_core = 6)
    # pasting out
    out_name <- paste0(out_direx, in_df[[i]]$species[1], substr(in_df[[i]]$year[1], 3, 4), "_", n_sample[a], "ind.csv")
    #write.csv(res_out, out_name)
  }
}


#### 4. RUN ASYMPTOTIC MODELS AND GENERATE PREDICTIONS BASED ON MM3 MODEL ####


in_direx <- "./Data_outputs/sensitivity_analysis_outputs/"
f_df <- data.frame(file = dir(in_direx))
f_df$dataset <- sub("_.*", "", f_df$file)
f_df$n_ind <- sub("ind.csv*", "", f_df$file)
f_df$n_ind <- sub("^[^_]*_", "", f_df$n_ind)

out_preds <- list()
for (i in 1:nrow(f_df)) {
  print(paste0(f_df$dataset[i], " ", f_df$n_ind[i], " individuals"))
  # load in file
  dat <- read.csv(paste0(in_direx, f_df$file[i]))
  # running asymptotic model code
  print("running asymptotic models...")
  dat$x <- dat$sample
  dat$area <- dat$tot_area
  mod_mm3 <- run_asymptotic_model(dat = dat, type = "MM3")
  # predicting minimum sample sizes
  print("predicting minimum sample sizes...")
  out_preds[[i]] <- sample_prediction(model = mod_mm3, pop_size = 100, percent_cover = 95)
  out_preds[[i]]$dataset <- f_df$dataset[i]
  out_preds[[i]]$n_ind <-  f_df$n_ind[i]
}
# combine all preds
all_preds <- do.call(rbind, out_preds)

# paste out
out_name <- "./Data_outputs/sensitivity_analysis_predictions.csv"
write.csv(all_preds, out_name)


#### 5. GLMM TO EXAMINE EFFECT OF ORIGINAL SAMPLE SIZE ON PREDICTIONS ####


# loading in sensitivity analysis predictions
pred_red <- read.csv("./Data_outputs/sensitivity_analysis_predictions.csv")
pred_og <- read.csv("./Data_outputs/area_occupied_50_predictions.csv")
pred_og <- subset(pred_og, model == "MM3",)
pred_og <- subset(pred_og, pop_size == 100,)
# pasting in original sample sizes
pred_og$n_ind <- c(9, 6, 14, 25, 24, 9, 6, 19, 25, 22, 25, 11, 25, 29, 15, 37, 33, 78, 10, 9, 19)
# combine into one dataframe 
preds <- data.frame(prediction = c(pred_og$prediction, pred_red$prediction), dataset = c(pred_og$dataset, pred_red$dataset), 
                    sample_og = c(pred_og$n_ind, pred_red$n_ind))
preds$sp <- gsub("[0-9]", "", preds$dataset)
preds$yr <-gsub("[a-zA-Z]", "", preds$dataset)
preds$sp_yr <- as.factor(as.character(paste(preds$sp, preds$yr, sep = "_")))
# remove combined years of data
preds <- subset(preds, dataset != "APBOTH",)
preds <- subset(preds, dataset != "WCPBOTH",)
preds <- subset(preds, dataset != "WAALL",)
preds <- subset(preds, dataset != "BSBOTH",)

preds$sp <- as.factor(as.character(preds$sp))
preds$yr <- as.factor(as.character(preds$yr))
# scaling variable to help convergence
preds$sample_og <- scale(preds$sample_og)  

# running model
mod1 <- glmer(prediction ~ sample_og + sp + (1|yr), family = "poisson",
                   data = preds)
summary(mod1)

# removing effect of original sample
mod2 <- update(mod1,~.- sample_og)
lrtest(mod1, mod2)
# Model 1: prediction ~ sample_og + sp + (1 | yr)
# Model 2: prediction ~ sp + (1 | yr)
# #Df  LogLik Df  Chisq Pr(>Chisq)  
# 1  12 -163.10                       
# 2  11 -165.02 -1 3.8335    0.05024 

# removing effect of species
mod3 <- update(mod1,~.- sp)
lrtest(mod1, mod3)
# Model 1: prediction ~ sample_og + sp + (1 | yr)
# Model 2: prediction ~ sample_og + (1 | yr)
# #Df  LogLik Df  Chisq Pr(>Chisq)    
# 1  12 -163.10                         
# 2   3 -215.11 -9 104.02  < 2.2e-16 ***

plot(mod3)
