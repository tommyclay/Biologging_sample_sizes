
#### FUNCTION TO REPROJECT (USUALLY FROM WGS84 TO UTM TO CREATE UDS) ####

proj_fun <- function(x, proj_in, proj_out) {
  # Convert to sf object
  df_sf <- st_as_sf(x, coords = c("longitude", "latitude"), crs = proj_in)
  # Transform projection
  df_out <- st_transform(df_sf, crs = proj_out)
  # Extract coordinates and rename to Longitude & Latitude
  coords <- st_coordinates(df_out)
  colnames(coords) <- c("longitude", "latitude")
  # Convert back to data frame (if needed, with coordinates extracted)
  df_out <- cbind(st_drop_geometry(df_out), coords)
  return(df_out)
}

# x = input dataframe
# proj_in = input projection (usually WGS84)
# proj_out = what to project to (usually utm coords)


### FUNCTION TO RESAMPLE HOME RANGE AREA WITH SAMPLE SIZE ####

resampling_UD_area <- function(x, n_iter, iso_level, cellsize, h, proj_in, proj_out, n_core, verbose = TRUE) {
  # re-projecting to utm to calculate kernel, if needed
  if (proj_in != proj_out) {
    x <- proj_fun(x = x, proj_in = proj_in, proj_out = proj_out) 
  } else { proj_out <- proj_in }
  
  if (verbose) start_time <- Sys.time()
  
  # 1. SETTING UP NULL GRID
  buffer <- h * 3
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

# x = input dataframe (id = individual identity, longitude and latitude)
# n_iter = number of resampled iterations
# iso_level = UD isopleth level
# cellsize = the grid cell size to calculate UDs (if utm then put in meters, i.e. X 1000)
# h = smoothing factor (if utm then put in meters, i.e. X 1000)
# proj_in = input projection (usually WGS84)
# proj_out = what to project to (usually utm coords)
# n_core = the number of cores to use
# verbose = if true, then plots time taken to run code


#### FUNCTION TO RUN ASMYPTOTIC MODELS ####

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
    # Remove invalid observations:
    dat2 <- subset(dat, is.finite(area) & is.finite(x) & x >= 0)
    # specifying starting values
    c_start <- mean(dat2$area[dat2$x == min(dat2$x)], na.rm = TRUE)
    a_start <- max(dat2$area, na.rm = TRUE)
    half <- c_start + 0.5 * (a_start - c_start)
    b_start <- dat2$x[which.min(abs(dat2$area - half))]
    # Define starting values
    start <- if (type %in% c("MM2", "AR2")) c(a_start, b_start) else c(c_start, a_start, b_start)
    # Retry fitting model with starting values
    mod <- drm(area ~ x, data = dat, fct = fun_form, start = start)
  }, error = function(e2) {
    message("Model fitting failed: ", e2$message)
    return(NULL)
  })
  return(mod)
}

# dat = input dataframe (x = sample size, area)
# type = type of asymptotic model ("MM2" = Michaelis Menten 2 parameter, "MM3" = MM 3-parameter, 
# "AR2" = asymptotic exponential 2 parameter, "AR3" = AR 3 parameter)
# function outputs model object

#### FUNCTION TO PREDICT MINIMUM SAMPLE SIZE FROM ASYMPTOTIC MODEL ####

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

# model = model output from run_asymptotic_model code
# pop_size = size of population to project to, can be a vector of multiple population sizes. 
# percent_cover = percentage of predicted asymptote to use as threshold (usually 95).

#### FUNCTION TO CALCULATE BA OVERLAP FOR ALL PAIRWISE COMBINATIONS OF INDIVIDUALS ####

pairwise_BA_overlap <- function(x, h, cellsize, iso_level, proj_in, proj_out) {
  if (proj_in != proj_out) {
    x <- proj_fun(x = x, proj_in = proj_in, proj_out = proj_out) 
  } else { proj_out <- proj_in }
  # calculating null grid
  coordinates(x) <- ~ longitude+latitude
  null.grid <- expand.grid(x=seq(min(coordinates(x)[,1])-(h*2), max(coordinates(x)[,1])+(h*2), by = cellsize),
                           y=seq(min(coordinates(x)[,2])-(h*2), max(coordinates(x)[,2])+(h*2), by = cellsize))
  coordinates(null.grid) <- ~x+y
  gridded(null.grid) <- TRUE
  # calculating all individual kernels
  lox <- data.frame(x= x$longitude, y = x$latitude, id = x$id)
  coordinates(lox) <-  ~x + y
  kerns <- kernelUD(lox, grid = null.grid, h = h, same4all = FALSE)
  BA <- kerneloverlaphr(kerns, method = "BA", conditional = TRUE, percent = iso_level)
  # clean matrix (remove diagonal, duplicates, NAs)
  BA[lower.tri(BA, diag = TRUE)] <- NA
  BA_v <- na.omit(as.vector(BA))
  # scale to 1
  BA_v <- as.numeric(BA_v / (iso_level / 100))
  return(BA_v)
}

# x = input dataframe (id = individual identity, longitude and latitude)
# iso_level = UD isopleth level
# cellsize = the grid cell size to calculate UDs (if utm then put in meters, i.e. X 1000)
# h = smoothing factor (if utm then put in meters, i.e. X 1000)
# proj_in = input projection (usually WGS84)
# proj_out = what to project to (usually utm coords)

#### FUNCTION TO PREDICT MINIMUM SAMPLE SIZES BASED ON AVERAGE BA VALUES ####

prediction_from_ba <- function(ba) {
  # Check if vector or dataframe
  if (is.data.frame(ba)) {
    x <- ba$ba
  } else if (is.vector(ba)) {
    x <- ba 
  }
  # Compute predictions
  y <- -3.9926/log(1-x)
  # Compute the gradient manually: derivative w.r.t k
  grad_matrix <- matrix(1 / log(1 - x), ncol = 1)  # Ensure it's a matrix
  # Compute standard error of predictions from vcov of model
  se_fit <- sqrt(grad_matrix^2 * 0.05120739)
  # Compute 95% confidence intervals
  alpha <- 0.05
  z_score <- qnorm(1 - alpha / 2)  # 1.96 for 95% CI
  upper <- y + z_score * se_fit
  lower <- y - z_score * se_fit
  # Output dataframe as whole numbers
  out_df <- data.frame(ba = x, prediction = round(y), lower = round(lower), upper = round(upper))
  return(out_df)
}

# either a vector or dataframe with column ba with average BA values (between 0 and 1) to predict from 


#### FUNCTION TO RESAMPLE OVERLAP WITH THREATS (E.G. FISHERIES) WITH SAMPLE SIZE ####

# function to resample home range area with sample size
resampling_fisheries_overlap <- function(x, n_iter, UD_cell, UD_h, threat_df, grid_res, n_cores, proj_dec, proj_utm) {
  
  # Start timing
  start_time <- Sys.time()
  # 1. DATA PROCESSING
  
  # converting fishing effort data to raster grid
  coordinates(threat_df) <- ~longitude + latitude
  proj4string(threat_df) <- proj_dec
  ext <- extent(-180, 180, -90, 90) 
  threat_r <- raster(ext, res = grid_res)
  threat_r <- rasterize(threat_df, threat_r, field = "value")
  threat_r[is.na(threat_r)] <- 0
  
  # create empty grid to convert to kernels to 
  world.grid <- expand.grid(longitude = seq(-179, 179, grid_res), latitude = seq(-89, 89, grid_res))
  world.grid2 <- world.grid 
  coordinates(world.grid2) <- ~longitude + latitude
  gridded(world.grid2) <- TRUE
  world_r <- raster(world.grid2)
  proj4string(world_r) <- proj_dec
  
  # just taking the ID and longitude and latitude
  x <- x[, c("longitude", "latitude", "id")]
  # re-projecting points using sf
  x_sf <- st_as_sf(x, coords = c("longitude", "latitude"))
  st_crs(x_sf) <- proj_dec
  if (proj_dec != proj_utm) {
    x_sf <- st_transform(x_sf, crs = proj_utm)
  }
  coords <- st_coordinates(x_sf)
  x_df <- cbind(st_drop_geometry(x_sf), coords)
  names(x_df) <- c("id", "longitude", "latitude")
  
  # 1. SETTING UP NULL GRID - using buffer of 2 times the smoothing factor around min and max points 
  null.grid <- expand.grid(x = seq(min(x_df$longitude) - (UD_h * 2), max(x_df$longitude) + (UD_h * 2), by = UD_cell),
                           y = seq(min(x_df$latitude) - (UD_h * 2), max(x_df$latitude) + (UD_h * 2), by = UD_cell))
  coordinates(null.grid) <- ~x + y
  gridded(null.grid) <- TRUE
  
  # 2. CALCULATE KERNELS FOR EACH INDIVIDUAL AND TURN TO MATRICES (TO AVERAGE ACROSS LATER)
  kern_in <- list()
  matrix_in <- list()
  n_ind <- length(unique(x_df$id))
  for (j in 1:n_ind) {
    id.sub <- subset(x_df, id == unique(x_df$id)[j])   
    id.sub$id <- as.factor(as.character(id.sub$id))
    lox <- id.sub[, c("longitude", "latitude")]
    coordinates(lox) <- ~longitude + latitude
    kern_in[[j]] <- adehabitatHR::kernelUD(lox, grid = null.grid, h = UD_h)
    matrix_in[[j]] <- as.matrix(raster(kern_in[[j]]))
  }
  names(kern_in) <- unique(x_df$id)
  names(matrix_in) <- unique(x_df$id)
  
  # Setting up parallel backend
  plan(multisession, workers = n_cores)
  on.exit(plan(sequential), add = TRUE)  # Ensures cleanup after function exits
  
  # Parallelized iterations
  iter.summary <- future_lapply(
    seq_len(n_iter),
    function(i) { 
      # Initialize lists and vectors for each iteration
      kern_out <- list()
      matrix_out <- list()
      overlap_v <- numeric(n_ind)
      o_change_v <- numeric(n_ind)
      
      # Randomize order of IDs
      random.ids <- sample(unique(x_df$id), n_ind, replace = FALSE)
      
      # First level of sample size (N = 1)
      matrix_out[[1]] <- matrix_in[[random.ids[1]]]
      kern_out[[1]] <- kern_in[[random.ids[1]]]
      
      # Resample to desired resolution
      r_out <- raster(kern_out[[1]])
      proj4string(r_out) <- proj_utm
      r_2 <- projectRaster(from = r_out, to = world_r)
      c.factor <- 1 / sum(na.omit(values(r_2)))
      values(r_2) <- values(r_2) * c.factor
      
      # Calculate overlap
      threat_overlap <- overlay(r_2, threat_r, fun = "*")
      overlap_df <- as.data.frame(threat_overlap, xy = TRUE)
      names(overlap_df) <- c("longitude", "latitude", "overlap")
      overlap_v[1] <- sum(na.omit(overlap_df$overlap))
      o_change_v[1] <- NA
      
      # Subsequent levels of sample size (N = 2+)
      for (a in 2:n_ind) { 
        matrix_out[[a]] <- matrix_in[[random.ids[a]]]
        matrix_sum <- Reduce('+', matrix_out) / a
        
        kern_out[[a]] <- kern_out[[1]]
        kern_out[[a]]@data$ud <- as.vector(matrix_sum)
        
        r_out <- raster(kern_out[[a]])
        proj4string(r_out) <- proj_utm
        r_2 <- projectRaster(from = r_out, to = world_r)
        c.factor <- 1 / sum(na.omit(values(r_2)))
        values(r_2) <- values(r_2) * c.factor
        
        threat_overlap <- overlay(r_2, threat_r, fun = "*")
        overlap_df <- as.data.frame(threat_overlap, xy = TRUE)
        names(overlap_df) <- c("longitude", "latitude", "overlap")
        overlap_v[a] <- sum(na.omit(overlap_df$overlap))
        o_change_v[a] <- (overlap_v[a] - overlap_v[a - 1]) / overlap_v[a - 1] * 100
      }
      
      data.frame(sample = seq_len(n_ind), overlap = overlap_v, o_change_p = o_change_v, iteration = i)
    },
    future.seed = TRUE,
    future.packages = c("adehabitatHR", "raster", "sp")  # <-- Pre-load packages
  )
  
  # Combine all iteration summaries
  all.iter <- do.call(rbind, iter.summary)
  
  # End timing
  end_time <- Sys.time()
  total_time <- end_time - start_time
  message(paste0("Resampling completed in ", round(total_time, 2), " ", attr(total_time, "units")))
  
  return(all.iter)
}

# x = input dataframe (id = individual identity, longitude and latitude)
# UD_cell = the grid cell size to calculate UDs (if utm then put in meters, i.e. X 1000)
# UD_h = smoothing factor for UD (if utm then put in meters, i.e. X 1000)
# threat_df = dataframe of threat grid (longitude, latitude, value = e.g number of hooks)
# grid_res = resolution of threat grid that we can to reproject UD to (in WGS84, e.g. 2 degrees) 
# n_cores = the number of cores to use
# proj_in = input projection (usually WGS84)
# proj_out = what to project to (usually utm coords)

#### FUNCTIONS TO CREATE INDIVIDUAL AND POPULATION KERNELS FOR PLOTTING ####

# calculate individual kernels
ind_kernels_poly <- function(x, h, cellsize, iso_level) {
  # calculating null grid
  coordinates(x) <- ~ longitude+latitude
  null.grid <- expand.grid(x=seq(min(coordinates(x)[,1])-(h*1000*2), max(coordinates(x)[,1])+(h*1000*2), by = cellsize*1000),
                           y=seq(min(coordinates(x)[,2])-(h*1000*2), max(coordinates(x)[,2])+(h*1000*2), by = cellsize*1000))
  coordinates(null.grid) <- ~x+y
  gridded(null.grid) <- TRUE
  # calculating all individual kernels
  lox <- data.frame(x= x$longitude, y = x$latitude, id = x$id)
  coordinates(lox) <-  ~x + y
  kerns <- kernelUD(lox, grid = null.grid, h = h*1000, same4all = FALSE)
  # getting isopleths
  isos <- getverticeshr(x = kerns, percent = iso_level)
  return(isos)
}

# calculate population kernels
pop_kernel_poly <- function(x, h, cellsize, iso_level) {
  
  # 1. SETTING UP NULL GRID - using buffer of 2 times the smoothing factor around min and max points 
  null.grid <- expand.grid(x=seq(min(x$longitude)-(h*1000*2), max(x$longitude)+(h*1000*2), by = cellsize*1000),
                           y=seq(min(x$latitude)-(h*1000*2), max(x$latitude)+(h*1000*2), by = cellsize*1000))
  coordinates(null.grid) <- ~x+y
  gridded(null.grid) <- TRUE
  
  # 2. CALCULATE KERNELS FOR EACH INDIVIDUAL AND TURN TO MATRICES (TO AVERAGE ACROSS LATER)
  kern_in <- list()
  matrix_in <- list()
  x$Id <- as.factor(as.character(x$id))
  for (j in 1:nlevels(x$id)) {
    id.sub <- subset(x, id==levels(x$id)[j])   
    id.sub$Id <- as.factor(as.character(id.sub$id))
    lox <- id.sub[,c("longitude", "latitude")]
    coordinates(lox) <-  ~longitude + latitude
    kern_in[[j]] <- adehabitatHR::kernelUD(lox, grid = null.grid, h = h*1000)
    matrix_in[[j]] <- as.matrix(raster(kern_in[[j]]))
  }
  
  # 3. MERGE AND PASTE OUT
  kern_out <- kern_in[[1]]
  matrix.summed <- Reduce('+', matrix_in)
  matrix.summed2 <- matrix.summed/length(kern_in)
  kern_out@data$ud <- as.vector(matrix.summed2)
  iso_out <- getverticeshr(kern_out, percent = iso_level)
  return(iso_out)
}

# x = input dataframe (id = individual identity, longitude and latitude)
# cellsize = the grid cell size to calculate UDs (if utm then put in meters, i.e. X 1000)
# h = smoothing factor for UD (if utm then put in meters, i.e. X 1000)
# iso_level = isopleth level (e.g. 50%)

