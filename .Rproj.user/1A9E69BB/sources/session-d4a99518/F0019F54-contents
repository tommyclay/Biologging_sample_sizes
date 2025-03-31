library(sf)
library(terra)
library(dplyr)
library(adehabitatHR)
library(raster)

#### 1. FUNCTIONS TO CALCULATE INDIVIDUAL OVERLAP ####

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

pairwise_BA_overlap <- function(x, h, cellsize, iso, proj_in, proj_out) {
  if (proj_in != proj_out) {
    x <- proj_fun(df_in = x, proj_in = proj_in, proj_out = proj_out) 
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
  BA <- kerneloverlaphr(kerns, method = "BA", conditional = TRUE, percent = iso)
  # clean matrix (remove diagonal, duplicates, NAs)
  BA[lower.tri(BA, diag = TRUE)] <- NA
  BA_v <- na.omit(as.vector(BA))
  # scale to 1
  BA_v <- as.numeric(BA_v / (iso / 100))
  return(BA_v)
}

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
  ba <- pairwise_BA_overlap(x = in_df[[i]], h = 186000, cellsize = 50000, iso = 50, proj_in = proj_dec, proj_out = proj_utm)
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
