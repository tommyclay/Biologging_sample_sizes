
library(plyr)
library(raster)
library(terra)
library(ks)
library(progress)
library(future.apply)
library(raster)
library(sf)
library(adehabitatHR)
library(future.apply)

# function to resample home range area with sample size
resampling_fisheries_overlap_og_changed <- function(x, n_iter, UD_cell, UD_h, threat_df, grid_res, n_cores, proj_dec, proj_utm) {
  
  # Start timing
  start_time <- Sys.time()
  
  # 1. DATA PROCESSING
  
  # converting fishing effort data to raster grid
  coordinates(threat_df)<- ~longitude+latitude
  proj4string(threat_df) <- proj_dec
  ext <- extent(-180, 180, -90, 90) 
  threat_r <- raster(ext, res = grid_res)
  threat_r <- rasterize(threat_df, threat_r, field = "value")
  threat_r[is.na(threat_r)] <- 0
  
  # create empty grid to convert to kernels to 
  world.grid <- expand.grid(longitude = seq(-179, 179, grid_res), latitude = seq(-89, 89, grid_res))
  world.grid2 <- world.grid 
  coordinates(world.grid2) <- ~longitude+latitude
  gridded(world.grid2) <- TRUE
  world_r <- raster(world.grid2)
  proj4string(world_r) <- proj_dec
  
  # just taking the ID and longitude and latitude
  x <- x[,c("longitude", "latitude", "Id")]
  # re-projecting points using sf
  x_sf <- st_as_sf(x, coords = c("longitude", "latitude"))
  st_crs(x_sf) <- proj_dec
  if (proj_in != proj_out) {
    x_sf <- st_transform(x_sf, crs = proj_utm)
  }
  coords <- st_coordinates(x_sf)
  x_df <- cbind(st_drop_geometry(x_sf), coords)
  names(x_df) <- c("Id", "longitude", "latitude")
  
  # 1. SETTING UP NULL GRID - using buffer of 2 times the smoothing factor around min and max points 
  null.grid <- expand.grid(x=seq(min(x_df$longitude)-(UD_h*2), max(x_df$longitude)+(UD_h*2), by = UD_cell),
                           y=seq(min(x_df$latitude)-(UD_h*2), max(x_df$latitude)+(UD_h*2), by = UD_cell))
  coordinates(null.grid) <- ~x+y
  gridded(null.grid) <- TRUE
  
  # 2. CALCULATE KERNELS FOR EACH INDIVIDUAL AND TURN TO MATRICES (TO AVERAGE ACROSS LATER)
  kern_in <- list()
  matrix_in <- list()
  n_ind <- length(unique(x_df$Id))
  for (j in 1:n_ind) {
    id.sub <- subset(x_df, Id==unique(x_df$Id)[j])   
    id.sub$Id <- as.factor(as.character(id.sub$Id))
    lox <- id.sub[,c("longitude", "latitude")]
    coordinates(lox) <-  ~longitude + latitude
    kern_in[[j]] <- adehabitatHR::kernelUD(lox, grid = null.grid, h = UD_h)
    matrix_in[[j]] <- as.matrix(raster(kern_in[[j]]))
  }
  names(kern_in) <- unique(x_df$Id)
  names(matrix_in) <- unique(x_df$Id)
  
  # set up progress bar
  pb <- progress_bar$new(format = "[:bar] :current/:total (Time remaining::eta)", total = n_iter)
  
  iter.summary <- list() # output summary dataframe for each iteration
  for (i in 1:n_iter) {
   
    # for progress bar
    pb$tick()
    Sys.sleep(1 / n_iter)
    # set up lists and output dataframes per level of sample size
    kern_out <- list()
    matrix_out <- list()
    overlap_v <- vector(length = n_ind)
    o_change_v <- vector(length = n_ind)
    # randomize order of ids
    random.ids <- sample(unique(x_df$Id),n_ind, replace = FALSE, prob = NULL) # random ids without replacement
    
    #  FOR FIRST LEVEL OF SAMPLE SIZE (IE, N = 1)
    matrix_out[[1]] <- subset(matrix_in, names(matrix_in)==random.ids[1])[[1]]
    kern_out[[1]] <- subset(kern_in, names(kern_in)==random.ids[1])[[1]]
    # resample to 2 degrees
    r_out <- raster::raster(kern_out[[1]])
    proj4string(r_out) <- proj_utm
    r_2 <- projectRaster(from = r_out, to = world_r)
    c.factor <- 1/sum(na.omit(raster::values(r_2)))
    raster::values(r_2) <- raster::values(r_2)*c.factor
    
    # calculate overlap
    threat_overlap <- raster::overlay(r_2, threat_r, fun ="*")
    overlap_df <- raster::as.data.frame(threat_overlap, xy = TRUE)
    names(overlap_df) <- c("longitude", "Lat", "overlap")
    overlap_v[1] <- sum(na.omit(overlap_df$overlap))
    o_change_v[1] <- NA
    
    #  FOR OTHER LEVEL OF SAMPLE SIZE (IE, N = 2+)
    for (a in 2:n_ind) { 
      # averaging values across matrices
      matrix_out[[a]] <- subset(matrix_in, names(matrix_in)==random.ids[a])[[1]]
      matrix_sum <- Reduce('+', matrix_out)
      matrix_sum <- matrix_sum/a
      # create kernel output, paste in merged values and calculate isopleth 
      kern_out[[a]] <- kern_out[[1]]
      kern_out[[a]]@data$ud <- as.vector(matrix_sum)
      r_out <- raster::raster(kern_out[[a]])
      proj4string(r_out) <- proj_utm
      r_2 <- projectRaster(from = r_out, to = world_r)
      c.factor <- 1/sum(na.omit(raster::values(r_2)))
      raster::values(r_2) <- raster::values(r_2)*c.factor
      # calculate overlap
      threat_overlap <- raster::overlay(r_2, threat_r, fun ="*")
      overlap_df <- raster::as.data.frame(threat_overlap, xy = TRUE)
      names(overlap_df) <- c("longitude", "latitude", "overlap")
      overlap_v[a] <- sum(na.omit(overlap_df$overlap))
      o_change <- overlap_v[a]-overlap_v[a-1]
      o_change_v[a] <- o_change/overlap_v[a-1]*100
    }
    iter.summary[[i]] <- data.frame(sample = seq(1, n_ind, 1), overlap = overlap_v, 
                                    o_change_p = o_change_v, iteration = i)
  }
  all.iter <- do.call(rbind, iter.summary)
  
  # End timing
  end_time <- Sys.time()
  total_time <- end_time - start_time
  message(paste0("Resampling completed in ", round(total_time, 2), " ", attr(total_time, "units")))
  
  return(all.iter)
}

resampling_fisheries_overlap_og_faster <- function(x, n_iter, UD_cell, UD_h, threat_df, grid_res, n_cores, proj_dec, proj_utm) {
  
  # Start timing
  start_time <- Sys.time()
  # 1. DATA PROCESSING
  
  # converting fishing effort data to raster grid
  coordinates(threat_df)<- ~longitude+latitude
  proj4string(threat_df) <- proj_dec
  ext <- extent(-180, 180, -90, 90) 
  threat_r <- raster(ext, res = grid_res)
  threat_r <- rasterize(threat_df, threat_r, field = "value")
  threat_r[is.na(threat_r)] <- 0
  
  # create empty grid to convert to kernels to 
  world.grid <- expand.grid(longitude = seq(-179, 179, grid_res), latitude = seq(-89, 89, grid_res))
  world.grid2 <- world.grid 
  coordinates(world.grid2) <- ~longitude+latitude
  gridded(world.grid2) <- TRUE
  world_r <- raster(world.grid2)
  proj4string(world_r) <- proj_dec
  
  # just taking the ID and longitude and latitude
  x <- x[,c("longitude", "latitude", "Id")]
  # re-projecting points using sf
  x_sf <- st_as_sf(x, coords = c("longitude", "latitude"))
  st_crs(x_sf) <- proj_dec
  if (proj_in != proj_out) {
    x_sf <- st_transform(x_sf, crs = proj_utm)
  }
  coords <- st_coordinates(x_sf)
  x_df <- cbind(st_drop_geometry(x_sf), coords)
  names(x_df) <- c("Id", "longitude", "latitude")
  
  # 1. SETTING UP NULL GRID - using buffer of 2 times the smoothing factor around min and max points 
  null.grid <- expand.grid(x=seq(min(x_df$longitude)-(UD_h*2), max(x_df$longitude)+(UD_h*2), by = UD_cell),
                           y=seq(min(x_df$latitude)-(UD_h*2), max(x_df$latitude)+(UD_h*2), by = UD_cell))
  coordinates(null.grid) <- ~x+y
  gridded(null.grid) <- TRUE
  
  # 2. CALCULATE KERNELS FOR EACH INDIVIDUAL AND TURN TO MATRICES (TO AVERAGE ACROSS LATER)
  kern_in <- list()
  matrix_in <- list()
  n_ind <- length(unique(x_df$Id))
  for (j in 1:n_ind) {
    id.sub <- subset(x_df, Id==unique(x_df$Id)[j])   
    id.sub$Id <- as.factor(as.character(id.sub$Id))
    lox <- id.sub[,c("longitude", "latitude")]
    coordinates(lox) <-  ~longitude + latitude
    kern_in[[j]] <- adehabitatHR::kernelUD(lox, grid = null.grid, h = UD_h)
    matrix_in[[j]] <- as.matrix(raster(kern_in[[j]]))
  }
  names(kern_in) <- unique(x_df$Id)
  names(matrix_in) <- unique(x_df$Id)
  
  # Setting up parallel backend
  plan(multisession, workers = n_cores)
  future::clusterEvalQ(future::getCluster(), {
    suppressPackageStartupMessages({
      library(adehabitatHR)
      library(raster)
      library(sp)
    })
  })
  
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
    random.ids <- sample(unique_ids, n_ind, replace = FALSE)
    
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
    
    data.frame(sample = seq_len(n_ind), overlap = overlap_v, o_change_p = o_change_v, Iteration = i)
    },
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

resampling_fisheries_overlap_og_faster2 <- function(x, n_iter, UD_cell, UD_h, threat_df, grid_res, n_cores, proj_dec, proj_utm) {
  
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
  x <- x[, c("longitude", "latitude", "Id")]
  # re-projecting points using sf
  x_sf <- st_as_sf(x, coords = c("longitude", "latitude"))
  st_crs(x_sf) <- proj_dec
  if (proj_in != proj_out) {
    x_sf <- st_transform(x_sf, crs = proj_utm)
  }
  coords <- st_coordinates(x_sf)
  x_df <- cbind(st_drop_geometry(x_sf), coords)
  names(x_df) <- c("Id", "longitude", "latitude")
  
  # 1. SETTING UP NULL GRID - using buffer of 2 times the smoothing factor around min and max points 
  null.grid <- expand.grid(x = seq(min(x_df$longitude) - (UD_h * 2), max(x_df$longitude) + (UD_h * 2), by = UD_cell),
                           y = seq(min(x_df$latitude) - (UD_h * 2), max(x_df$latitude) + (UD_h * 2), by = UD_cell))
  coordinates(null.grid) <- ~x + y
  gridded(null.grid) <- TRUE
  
  # 2. CALCULATE KERNELS FOR EACH INDIVIDUAL AND TURN TO MATRICES (TO AVERAGE ACROSS LATER)
  kern_in <- list()
  matrix_in <- list()
  n_ind <- length(unique(x_df$Id))
  for (j in 1:n_ind) {
    id.sub <- subset(x_df, Id == unique(x_df$Id)[j])   
    id.sub$Id <- as.factor(as.character(id.sub$Id))
    lox <- id.sub[, c("longitude", "latitude")]
    coordinates(lox) <- ~longitude + latitude
    kern_in[[j]] <- adehabitatHR::kernelUD(lox, grid = null.grid, h = UD_h)
    matrix_in[[j]] <- as.matrix(raster(kern_in[[j]]))
  }
  names(kern_in) <- unique(x_df$Id)
  names(matrix_in) <- unique(x_df$Id)
  
  # Setting up parallel backend
  plan(multisession, workers = n_cores)
  
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
      random.ids <- sample(unique_ids, n_ind, replace = FALSE)
      
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

resampling_fisheries_overlap_mod <- function(x, n_iter, UD_cell, UD_h, threat_df, grid_res, n_cores, proj_in, proj_out) {
  
  # 1. DATA PROCESSING
  
  # converting fishing effort data to raster grid
  coordinates(threat_df)<- ~longitude+latitude
  proj4string(threat_df) <- proj_out
  ext <- extent(-180, 180, -90, 90) 
  threat_r <- raster(ext, res = grid_res)
  threat_r <- rasterize(threat_df, threat_r, field = "value")
  threat_r[is.na(threat_r)] <- 0
  
  # create empty grid to convert to kernels to 
  world.grid <- expand.grid(longitude = seq(-179, 179, grid_res), latitude = seq(-89, 89, grid_res))
  world.grid2 <- world.grid 
  coordinates(world.grid2) <- ~longitude+latitude
  gridded(world.grid2) <- TRUE
  world_r <- raster(world.grid2)
  proj4string(world_r) <- proj_out
  
  # just taking the ID and longitude and latitude
  x <- x[,c("longitude", "latitude", "Id")]
  # re-projecting points using sf
  x_sf <- st_as_sf(x, coords = c("longitude", "latitude"))
  st_crs(x_sf) <- proj_in
  if (proj_in != proj_out) {
    x_sf <- st_transform(x_sf, crs = proj_out)
  }
  
  # 1. SETTING UP NULL GRID - using buffer of 2 times the smoothing factor around min and max points 
  null.grid <- expand.grid(x = seq(min(x$longitude) - (UD_h * 2), 
                                   max(x$longitude) + (UD_h * 2), 
                                   by = UD_cell),
                           y = seq(min(x$latitude) - (UD_h * 2), 
                                   max(x$latitude) + (UD_h * 2), 
                                   by = UD_cell))
  
  coordinates(null.grid) <- ~x + y
  gridded(null.grid) <- TRUE
  
  # 2. CALCULATE KERNELS FOR EACH INDIVIDUAL AND TURN TO MATRICES (TO AVERAGE ACROSS LATER)
  kern_in <- list()
  matrix_in <- list()
  x$Id <- as.factor(as.character(x$Id))
  
  # Parallelize kernel calculations using future_lapply
  plan(multisession, workers = n_cores) # Setup parallel plan using multisession
  kern_in <- future_lapply(1:nlevels(x$Id), function(j) {
    id.sub <- subset(x, Id == levels(x$Id)[j])   
    lox <- id.sub[, c("longitude", "latitude")]
    coordinates(lox) <- ~longitude + latitude
    kern <- kernelUD(lox, grid = null.grid, h = UD_h)
    matrix <- as.matrix(raster(kern))
    list(kern = kern, matrix = matrix)
  })
  
  # Extract kernels and matrices into separate lists
  kern_in <- lapply(kern_in, `[[`, 1)
  matrix_in <- lapply(kern_in, `[[`, 2)
  
  names(kern_in) <- levels(x$Id)
  names(matrix_in) <- levels(x$Id)
  
  # Set up progress bar
  pb <- progress_bar$new(format = "[:bar] :current/:total (Time remaining::eta)", total = n_iter)
  
  iter.summary <- list() # Output summary dataframe for each iteration
  
  # Parallelize the iteration loop
  for (i in 1:n_iter) {
    
    # For progress bar
    pb$tick()
    Sys.sleep(1 / n_iter)
    
    # Set up lists and output dataframes per level of sample size
    kern_out <- list()
    matrix_out <- list()
    overlap_v <- numeric(nlevels(x$Id))
    o_change_v <- numeric(nlevels(x$Id))
    
    # Randomize order of ids
    random.ids <- sample(levels(x$Id), nlevels(x$Id), replace = FALSE)
    
    # FOR FIRST LEVEL OF SAMPLE SIZE (IE, N = 1)
    matrix_out[[1]] <- matrix_in[[random.ids[1]]]
    kern_out[[1]] <- kern_in[[random.ids[1]]]
    
    # Resample to grid resolution
    r_out <- raster(kern_out[[1]])
    proj4string(r_out) <- proj_in
    r_2 <- projectRaster(from = r_out, to = world_r)
    c.factor <- 1 / sum(na.omit(values(r_2)))
    values(r_2) <- values(r_2) * c.factor
    
    # Calculate overlap
    threat_overlap <- overlay(r_2, threat_r, fun = "*")
    overlap_df <- as.data.frame(threat_overlap, xy = TRUE)
    names(overlap_df) <- c("longitude", "latitude", "overlap")
    overlap_v[1] <- sum(na.omit(overlap_df$overlap))
    o_change_v[1] <- NA
    
    # FOR OTHER LEVEL OF SAMPLE SIZE (IE, N = 2+)
    for (a in 2:nlevels(x$Id)) { 
      
      # Averaging values across matrices
      matrix_out[[a]] <- matrix_in[[random.ids[a]]]
      matrix_sum <- Reduce('+', matrix_out) / a
      
      # Create kernel output, paste in merged values and calculate isopleth 
      kern_out[[a]] <- kern_out[[1]]
      kern_out[[a]]@data$ud <- as.vector(matrix_sum)
      
      r_out <- raster(kern_out[[a]])
      proj4string(r_out) <- proj.utm
      r_2 <- projectRaster(from = r_out, to = world_r)
      c.factor <- 1 / sum(na.omit(values(r_2)))
      values(r_2) <- values(r_2) * c.factor
      
      # Calculate overlap
      threat_overlap <- overlay(r_2, threat_r, fun = "*")
      overlap_df <- as.data.frame(threat_overlap, xy = TRUE)
      names(overlap_df) <- c("longitude", "latitude", "overlap")
      overlap_v[a] <- sum(na.omit(overlap_df$overlap))
      
      o_change_v[a] <- (overlap_v[a] - overlap_v[a - 1]) / overlap_v[a - 1] * 100
    }
    
    iter.summary[[i]] <- data.frame(Sample = 1:nlevels(x$Id), 
                                    Overlap = overlap_v, 
                                    O_change_p = o_change_v, 
                                    Iteration = i)
  }
  
  all.iter <- do.call(rbind, iter.summary)
  return(all.iter)
}


#### LOAD IN FISHERIES DATA ####


all.ll <- read.csv("C:/Users/44753/Documents/Kalinka WCP study/Data_outputs/GFW/GFW_v2_resample_2x2_grids/GFW_Pel_LL_resample5x5_CHN_rm.csv")
all.ll$flag <- as.factor(as.character(all.ll$flag))
# SUM FOR ALL FLEETS #
all.ll.sum <- ddply(all.ll,.(month, lon, lat), summarize, Hours = sum(Hours))
# AVERAGE ACROSS MONTHS
all.ll.avg <- ddply(all.ll,.(lon, lat), summarize, Hours = mean(Hours))
fish_df <- all.ll.avg
names(fish_df) <- c("longitude", "latitude", "value")


#### LOAD IN EXAMPLE TRACKING DATA ####

gha.99 <- read.csv("./Data_inputs/gls_data/gha.1999.csv")
gha.99$Id <- as.factor(as.character(gha.99$Id))
gha.99$longitude <- gha.99$Lon
gha.99$latitude <- gha.99$Lat

proj.utm <- paste0("+proj=laea +lat_0=", -90, " +lon_0=", -38.05, " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj.dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"   

test <- resampling_fisheries_overlap_mod(x = gha.99, n_iter = 5, UD_cell = 50000, UD_h = 186000, 
                                         threat_df = fish_df, grid_res = 2, n_cores = 4, 
                                         proj_in = proj.utm, proj_out = proj.dec) 
  
test <- resampling_fisheries_overlap_og_changed(x = gha.99, n_iter = 100, UD_cell = 50000, UD_h = 186000, 
                                         threat_df = fish_df, grid_res = 2, n_cores = 4, 
                                         proj_utm = proj.utm, proj_dec = proj.dec) 
# 18 seconds with progress bar for 10 iterations/4 cores
# 2.5 minutes for 100 iterations/4 cores

test <- resampling_fisheries_overlap_og_faster(x = gha.99, n_iter = 10, UD_cell = 50000, UD_h = 186000, 
                                                threat_df = fish_df, grid_res = 2, n_cores = 4, 
                                                proj_utm = proj.utm, proj_dec = proj.dec) 
# 24 seconds for 10 iterations/4 cores
# 1.04 for 100 iterations/4 cores

test <- resampling_fisheries_overlap_og_faster2(x = gha.99, n_iter = 500, UD_cell = 50000, UD_h = 186000, 
                                               threat_df = fish_df, grid_res = 2, n_cores = 4, 
                                               proj_utm = proj.utm, proj_dec = proj.dec) 
# 24 seconds for 10 iterations/4 cores
# 1.12 mins for 100 iterations/4 cores
# 3.88 mins for 500 iterations/4 cores

# need to remove future apply warning messages

ggplot(test, aes(sample, overlap)) + geom_line(aes(group = Iteration))
ggplot(test, aes(sample, o_change_p)) + geom_line(aes(group = Iteration))
