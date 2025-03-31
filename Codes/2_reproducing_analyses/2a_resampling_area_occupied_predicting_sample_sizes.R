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

#### 1. SPECIFYING FUNCTIONS - IN FUTURE THIS WILL BE IN SEPARATE R FILE WE CAN LINK TO #####

# function to project from dec to utm
proj_fun <- function(df_in, proj_in, proj_out) {
  # Convert to sf object
  df_sf <- st_as_sf(df_in, coords = c("longitude", "latitude"), crs = proj_in)
  # Transform projection
  df_out <- st_transform(df_sf, crs = proj_out)
  # Extract coordinates and rename to Longitude & Latitude
  coords <- st_coordinates(df_out)
  colnames(coords) <- c("longitude", "latitude")
  # Convert back to data frame (if needed, with coordinates extracted)
  df_out <- cbind(st_drop_geometry(df_out), coords)
  return(df_out)
}

# 1. FUNCTION TO RESAMPLE HOME RANGE AREA WITH SAMPLE SIZE

# if utm, put smoothing factor in metres not km
resampling_UD_area <- function(x, n_iter, iso_level, cellsize, h, proj_in, proj_out, n_core, verbose = TRUE) {
  # re-projecting to utm to calculate kernel, if needed
  if (proj_in != proj_out) {
    x <- proj_fun(df_in = x, proj_in = proj_in, proj_out = proj_out) 
  } else { proj_out <- proj_in }
  
  if (verbose) start_time <- Sys.time()
  
  # 1. SETTING UP NULL GRID
  buffer <- h * 2
  long_seq <- seq(min(x$longitude) - buffer, max(x$longitude) + buffer, by = cellsize)
  lat_seq <- seq(min(x$latitude) - buffer, max(x$latitude) + buffer, by = cellsize)
  null_grid <- expand.grid(x = long_seq, y = lat_seq)
  coordinates(null_grid) <- ~x + y
  gridded(null_grid) <- TRUE
  
  # 2. CALCULATE KERNELS & MATRICES PER INDIVIDUAL
  id_levels <- levels(as.factor(as.character(x$id)))
  n_levels <- length(id_levels)
  kern_in <- vector("list", n_levels)
  matrix_in <- vector("list", n_levels)
  names(kern_in) <- names(matrix_in) <- id_levels
  
  for (j in seq_along(id_levels)) {
    id_sub <- x[x$id == id_levels[j], ]
    lox <- id_sub[, c("longitude", "latitude")]
    coordinates(lox) <- ~ longitude + latitude
    kern_in[[j]] <- kernelUD(lox, grid = null_grid, h = h)
    matrix_in[[j]] <- as.matrix(raster(kern_in[[j]]))
  }
  
  # 3. ITERATION FUNCTION (optimized inner loop)
  iteration_func <- function(i) {
    #set.seed(i)  # For reproducibility
    random_ids <- sample(id_levels, n_levels, replace = FALSE)
    # Get indices once:
    idxs <- match(random_ids, names(matrix_in))
    tot_v <- numeric(n_levels)
    
    running_sum <- matrix_in[[ idxs[1] ]]
    kern_base <- kern_in[[ idxs[1] ]]
    iso <- getverticeshr(kern_base, percent = iso_level)
    tot_v[1] <- as.numeric(st_area(st_as_sf(iso)))
    
    for (a in 2:n_levels) {
      running_sum <- running_sum + matrix_in[[ idxs[a] ]]
      avg_matrix <- running_sum / a
      
      # Copy and update the kernel object (S4 copying overhead is minimized here)
      kern_new <- kern_base
      kern_new@data$ud <- as.vector(avg_matrix)
      iso_new <- getverticeshr(kern_new, percent = iso_level)
      tot_v[a] <- as.numeric(st_area(st_as_sf(iso_new)))
    }
    
    data.frame(
      sample = seq_len(n_levels),
      tot_area = tot_v,
      isopleth = iso_level,
      iteration = i
    )
  }
  
  # Compile the iteration function for a speed boost
  iteration_func <- cmpfun(iteration_func)
  
  if (verbose) message("Starting resampling iterations...")
  
  # 4. SET UP CLUSTER AND EXPORT NECESSARY OBJECTS & PACKAGES
  available_cores <- detectCores() - 1
  n_core <- min(n_core, available_cores)
  
  if (n_core > 1) {
    cl <- makeCluster(n_core)
    on.exit(stopCluster(cl))
    clusterEvalQ(cl, { library(adehabitatHR); library(raster); library(sf) })
    clusterExport(cl, varlist = c("id_levels", "matrix_in", "kern_in", "iso_level", "getverticeshr"), envir = environment())
    result_list <- parLapply(cl, seq_len(n_iter), iteration_func)
  } else {
    result_list <- lapply(seq_len(n_iter), iteration_func)  # Run sequentially
  }
  
  result <- do.call(rbind, result_list)
  result$tot_area <- result$tot_area / 1e6
  
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message(paste0("Resampling completed in ", round(total_time, 2), " ", attr(total_time, "units")))
  }
  
  return(result)
}

# 2. FUNCTION TO RUN ASMYPTOTIC MODEL

# sample size should be "x", y should be "area". 
run_asymptotic_model <- function(dat, type) {
  # Define the function form
  fun_form <- switch(type,
                 "MM2" = MM.2(),
                 "MM3" = MM.3(),
                 "AR2" = AR.2(),
                 "AR3" = AR.3(),
                 stop("Invalid model type specified"))
  
  # Try fitting the model, handle errors
  tryCatch({
    # Extract the max sample area
    a <- dat$area[dat$x == max(dat$x)][1]
    # Calculate median areas
    dat_med <- ddply(dat, .(x), summarize, med_area = median(area, na.rm = TRUE))
    dat_med <- na.omit(dat_med)
    # Find closest x value to half of max area
    b <- dat$x[which.min(abs(dat_med$med_area - a / 2))]
    # Define starting values
    start <- if (type %in% c("MM2", "AR2")) c(a, b) else c(0, a, b)
    # Retry fitting model with starting values
    mod <- drm(area ~ x, data = dat, fct = fun_form, start = start)
  }, error = function(e2) {
    message("Model fitting failed: ", e2$message)
    return(NULL)
    })
return(mod)
}

# 3. FUNCTION TO PREDICT MINIMUM SAMPLE SIZE
sample_prediction <- function(model, pop_size, percent_cover) {
  
  pred.ci <- lapply(seq_along(pop_size), function(i) {
    # Create data frame for sampling
    pred.df <- data.frame(sample = seq_len(pop_size[i]), groupings = 1)
    # Get prediction results with confidence intervals
    pred.df2 <- as.data.frame(predict(model, newdata = pred.df, se.fit = TRUE))
    # Extract predictions and standard errors
    pred.df$Prediction <- pred.df2$Prediction
    pred.df$SE  <- pred.df2$SE
    # Compute confidence intervals manually (assuming normal distribution) as we get a warning message otherwise
    tquan <- qt(0.975, df = model$df.residual)  # 95% CI quantile
    pred.df$Upper <- pred.df$Prediction + tquan * pred.df$SE
    pred.df$Lower <- pred.df$Prediction - tquan * pred.df$SE
    # Directly get max values
    max.value.mean  <- max(na.omit(pred.df$Prediction))
    max.value.upper <- max(na.omit(pred.df$Upper))
    max.value.lower <- max(na.omit(pred.df$Lower))
    # Calculate threshold values related to percent cover
    value.mean  <- max.value.mean  * percent_cover/100
    value.upper <- max.value.upper * percent_cover/100
    value.lower <- max.value.lower * percent_cover/100
    # Find the first value greater than or equal to the threshold using logical indexing
    number.mean  <- which(pred.df$Prediction >= value.mean)[1]
    number.upper <- which(pred.df$Upper >= value.upper)[1]
    number.lower <- which(pred.df$Lower >= value.lower)[1]
    # If thresholds are not met, use the maximum sample size
    number.mean  <- ifelse(is.na(number.mean),  pop_size[i], number.mean)
    number.upper <- ifelse(is.na(number.upper), pop_size[i], number.upper)
    number.lower <- ifelse(is.na(number.lower), pop_size[i], number.lower)
    # Return the results as a data frame
    data.frame(
      prediction = number.mean,
      lower      = number.lower,
      upper      = number.upper,
      pop_size = pop_size[i]
    )
  })
  # Combine all results into one data frame
  predict.df <- do.call(rbind, pred.ci)
  return(predict.df)
}


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
  res_50 <- resampling_UD_area(x = in_df[[i]], n_iter = 100, 
                                iso_level = 50, cellsize = 50000, h = 186000, proj_in = proj_dec,
                                proj_out = proj_utm, n_core = 6)
  # pasting out 50% UD iterations
  out_name <- paste0(out_direx, in_df[[i]]$species[1], substr(in_df[[i]]$year[1], 3, 4), "_", "50.csv")
  #write.csv(res_50, out_name)
  # resampling 95% UD
  res_95 <- resampling_UD_area(x = in_df[[i]], n_iter = 100, 
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

# paste out AICs and minimum sample sizes
out_name <- paste0(out_direx, "area_occupied_50_aics.csv")
write.csv(all_aics_50, out_name)
out_name <- paste0(out_direx, "area_occupied_50_predictions.csv")
write.csv(all_preds_50, out_name)

# identifying best model for each dataset
aic_best <- all_aics_50 %>%
  group_by(dataset) %>%
  slice_min(AIC)
aic_best
nrow(aic_best[aic_best$model == "MM3",])
# 13 of 21 have MM3 as best model

# pasting out MM3 models for each dataset
out_mod <- "./Data_outputs/resampled_area_outputs/gls_models/"
for (i in 1:length(f_50)) {
  out_name <- paste0(out_mod, gsub(".csv", "", f_50[i]), ".RData")
  mod <- mod_mm3_50[[i]]
  save(mod, file = out_name)
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

# paste out AICs and minimum sample sizes
out_name <- paste0(out_direx, "area_occupied_95_aics.csv")
write.csv(all_aics_95, out_name)
out_name <- paste0(out_direx, "area_occupied_95_predictions.csv")
write.csv(all_preds_95, out_name)

# identifying best model for each dataset
aic_best <- all_aics_95 %>%
  group_by(dataset) %>%
  slice_min(AIC)
aic_best
nrow(aic_best[aic_best$model == "MM3",])
# 18 of 21 have MM3 as best model
