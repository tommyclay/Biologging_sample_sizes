library(sf)
library(terra)
library(parallel)
library(adehabitatHR)
library(compiler)

library(drc)

#### 1. SPECIFYING FUNCTIONS - IN FUTURE THIS WILL BE IN SEPARATE R FILE WE CAN LINK TO #####

# 1. FUNCTION TO RESAMPLE HOME RANGE AREA OCCUPIED 

# function to project from dec to utm
proj_fun <- function(df_in, proj_in, proj_out) {
  # Convert to sf object
  df_sf <- st_as_sf(df_in, coords = c("Longitude", "Latitude"), crs = proj_in)
  # Transform projection
  df_out <- st_transform(df_sf, crs = proj_out)
  # Extract coordinates and rename to Longitude & Latitude
  coords <- st_coordinates(df_out)
  colnames(coords) <- c("Longitude", "Latitude")
  # Convert back to data frame (if needed, with coordinates extracted)
  df_out <- cbind(st_drop_geometry(df_out), coords)
  return(df_out)
}

# if utm, put smoothing factor in metres not km
resampling_UD_area <- function(x, n_iter, iso_level, cellsize, h, proj_in, proj_out, n_core, verbose = TRUE) {
  # re-projecting to utm to calculate kernel, if needed
  if (proj_in != proj_out) {
    x <- proj_fun(df_in = x, proj_in = proj_in, proj_out = proj_out) 
  } else { proj_out <- proj_in }
  
  if (verbose) start_time <- Sys.time()
  
  # pre-calculate constants
  h_m <- h            # smoothing factor in meters
  buffer <- h_m * 2          # buffer distance
  cellsize_m <- cellsize
  
  # 1. SETTING UP NULL GRID (terra + sf)
  bbox <- st_bbox(c(
    xmin = min(x$longitude) - buffer,
    xmax = max(x$longitude) + buffer,
    ymin = min(x$latitude)  - buffer,
    ymax = max(x$latitude)  + buffer
  ), crs = st_crs(proj_out))
  # Create empty raster grid
  grid_raster <- raster(xmn = bbox["xmin"], xmx = bbox["xmax"],
    ymn = bbox["ymin"], ymx = bbox["ymax"],
    resolution = cellsize_m,crs = proj_out)
  
  # 2. CALCULATE KERNELS & MATRICES PER INDIVIDUAL
  id_levels <- levels(as.factor(as.character(x$Id)))
  n_levels <- length(id_levels)
  kern_in <- vector("list", n_levels)
  matrix_in <- vector("list", n_levels)
  names(kern_in) <- names(matrix_in) <- id_levels
  
  for (j in seq_along(id_levels)) {
    id_sub <- x[x$Id == id_levels[j], ]
    lox <- id_sub[, c("Longitude", "Latitude")]
    coordinates(lox) <- ~ longitude + latitude
    kern_in[[j]] <- kernelUD(lox, grid = null.grid, h = h_m)
    matrix_in[[j]] <- as.matrix(raster(kern))
    kern_in[[j]] <- kern
  }
  
  # 3. ITERATION FUNCTION (optimized inner loop)
  iteration_func <- function(i) {
    set.seed(i)  # For reproducibility
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
      Sample = seq_len(n_levels),
      Tot_area = tot_v,
      Isopleth = iso_level,
      Iteration = i
    )
  }
  
  # Compile the iteration function for a speed boost
  iteration_func <- cmpfun(iteration_func)
  
  if (verbose) message("Starting resampling iterations...")
  
  # 4. SET UP CLUSTER AND EXPORT NECESSARY OBJECTS & PACKAGES
  cl <- makeCluster(n_core)
  clusterEvalQ(cl, {
    library(adehabitatHR)
    library(raster)
    library(sf)
  })
  # Export required variables and functions to each worker:
  clusterExport(cl, varlist = c("id_levels", "matrix_in", "kern_in", "iso_level", 
                                "gArea", "getverticeshr"), envir = environment())
  
  # 5. RUN PARALLEL ITERATIONS
  result_list <- parLapply(cl, seq_len(n_iter), iteration_func)
  stopCluster(cl)
  result <- do.call(rbind, result_list)
  
  # 6. POST-PROCESSING: Convert areas from m² to km²
  result$Tot_area <- result$Tot_area / 1e6
  
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message(paste0("Resampling completed in ", round(total_time, 2), " ", attr(total_time, "units")))
  }
  
  return(result)
}

# 2. FUNCTION TO RUN ASYMPTOTIC MODELS 
# sample size should be "x", y should be "area". 
run_asymptotic_model <- function(x, type, start) {
  if (type == "MM2") { form = "MM.2()"
  } else if (type == "MM3") { form = "MM.3()"
  } else if (type == "AR2") { form = "AR.2()"
  } else if (type == "AR3") { form = "AR.3()" }
  mod <- drm(area ~ x, data = x, fct = type)
}

# NEED TO MODIFY THIS CODE!!!
run_asymptotic_model2 <- function(x, type, start = NULL) {
  
  # Determine the model function and default start values based on type
  if (type == "MM2") {
    model_func <- MM.2()
    if (is.null(start)) start <- c(a = max(x$area), b = 1)
  } else if (type == "MM3") {
    model_func <- MM.3()
    if (is.null(start)) start <- c(a = max(x$area), b = 1, c = 0)
  } else if (type == "AR2") {
    model_func <- AR.2()
    if (is.null(start)) start <- c(a = max(x$area), b = 1)
  } else if (type == "AR3") {
    model_func <- AR.3()
    if (is.null(start)) start <- c(a = max(x$area), b = 1, c = 0)
  } else {
    stop("Invalid model type specified.")
  }
  
  # Fit model
  mod <- drm(area ~ x, data = x, fct = model_func, start = start)
  
  return(mod)
}


# 3. FUNCTION TO PREDICT MINIMUM SAMPLE SIZE

# can specify multiple pop sizes
sample_prediction <- function(model, pop_size, percent_cover) {
  
  pred.ci <- lapply(seq_along(pop_size), function(i) {
    # Create data frame for sampling
    pred.df <- data.frame(sample = seq_len(pop_size[i]), groupings = 1)
    # Get prediction results with confidence intervals
    pred.df2 <- as.data.frame(predict(model, newdata = pred.df, se.fit = FALSE, interval = "confidence", level = 0.95))
    # Directly get max values
    max.value.mean  <- max(pred.df2$Prediction)
    max.value.upper <- max(pred.df2$Upper)
    max.value.lower <- max(pred.df2$Lower)
    # Calculate threshold values related to percent cover
    value.mean  <- max.value.mean  * percent_cover
    value.upper <- max.value.upper * percent_cover
    value.lower <- max.value.lower * percent_cover
    # Find the first value greater than or equal to the threshold using logical indexing
    number.mean  <- which(pred.df2$Prediction >= value.mean)[1]
    number.upper <- which(pred.df2$Upper >= value.upper)[1]
    number.lower <- which(pred.df2$Lower >= value.lower)[1]
    # If thresholds are not met, use the maximum sample size
    number.mean  <- ifelse(is.na(number.mean),  pop_size[i], number.mean)
    number.upper <- ifelse(is.na(number.upper), pop_size[i], number.upper)
    number.lower <- ifelse(is.na(number.lower), pop_size[i], number.lower)
    # Return the results as a data frame
    data.frame(
      Prediction = number.mean,
      Upper      = number.upper,
      Lower      = number.lower,
      Predict.To = pop_size[i]
    )
  })
  # Combine all results into one data frame
  predict.df <- do.call(rbind, pred.ci)
  return(predict.df)
}





#### 2. RUNNING ON EXAMPLE CODE ####

# 



