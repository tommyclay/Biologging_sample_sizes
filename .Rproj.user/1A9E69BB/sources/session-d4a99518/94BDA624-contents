library(sf)
library(terra)
library(dplyr)
#library(ks)        # For KDE
#library(proxy)     # For pairwise comparisons
library(adehabitatHR)
library(raster)
#library(dplyr)

# projection formats - utm coordinates centered colony location
proj.dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"   

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

ap10 <- read.csv("./Data_inputs/gls_data/ap10.csv")
ap10$id <- as.factor(as.character(ap.10$id))
# setting utm projection 
proj.utm <- paste0("+proj=laea +lat_0=", -90, " +lon_0=", -38.05, " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
ap10_BA <- pairwise_BA_overlap(x = ap.10, h = 186000, cellsize = 50000, iso = 50, proj_in = proj.dec, proj_out = proj.utm)
median(ap10_BA) # 0.4371286
> 
