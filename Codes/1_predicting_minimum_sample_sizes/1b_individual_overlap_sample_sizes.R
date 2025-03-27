library(sf)
library(terra)
library(dplyr)
library(ks)        # For KDE
library(proxy)     # For pairwise comparisons
library(adehabitatHR)
library(raster)
library(dplyr)

# projection formats - utm coordinates centered colony location
proj.dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"   

pairwise_BA_og <- function(x, h, cellsize, iso_level) {
  # calculating null grid
  coordinates(x) <- ~ Longitude+Latitude
  null.grid <- expand.grid(x=seq(min(coordinates(x)[,1])-(h*1000*2), max(coordinates(x)[,1])+(h*1000*2), by = cellsize*1000),
                           y=seq(min(coordinates(x)[,2])-(h*1000*2), max(coordinates(x)[,2])+(h*1000*2), by = cellsize*1000))
  coordinates(null.grid) <- ~x+y
  gridded(null.grid) <- TRUE
  # calculating all individual kernels
  lox <- data.frame(x= x$Longitude, y = x$Latitude, id = x$Id)
  coordinates(lox) <-  ~x + y
  kerns <- kernelUD(lox, grid = null.grid, h = h*1000, same4all = FALSE)
  BA <- kerneloverlaphr(kerns, method = "BA", conditional = TRUE, percent = iso_level)
  # making diagonal elements NA and delete repeated values and nas
  diag(BA) <- NA
  BA_v <- as.vector(BA)
  BA.df <- data.frame(BA = BA_v, Dup = duplicated(BA_v))
  BA.df <- subset(BA.df, Dup != TRUE,)
  BA.df <- BA.df[!is.na(BA.df$BA),]
  # converting to out of 1 rather than 0.5
  BA.df$BA <- BA.df$BA/(iso_level/100)*1 
  return(BA.df$BA)
}
  
pairwise_BA_ade <- function(x, h, cellsize, iso, proj_in, proj_out, buffer_dist = h*2) {
  # just taking the ID and longitude and latitude
  x <- x[,c("longitude", "latitude", "Id")]
  # re-projecting points using sf
  x_sf <- st_as_sf(x, coords = c("longitude", "latitude"))
  st_crs(x_sf) <- proj_in
  if (proj_in != proj_out) {
    x_sf <- st_transform(x_sf, crs = proj_out)
  }
  # get bbox as numeric
  bbox <- st_bbox(x_sf)
  xmin <- as.numeric(bbox["xmin"]) - buffer_dist
  xmax <- as.numeric(bbox["xmax"]) + buffer_dist
  ymin <- as.numeric(bbox["ymin"]) - buffer_dist
  ymax <- as.numeric(bbox["ymax"]) + buffer_dist
  # create raster grid
  grid_raster <- raster(
    xmn = xmin, xmx = xmax,
    ymn = ymin, ymx = ymax,
    res = cellsize, crs = proj_out
  )
  grid_sp <- as(grid_raster, "SpatialPixelsDataFrame")
  proj4string(grid_sp) <- proj_out
  # convert sf to SpatialPointsDataFrame with ID
  x_sp <- as(x_sf, "Spatial")
  x_sp@data$Id <- as.factor(x$Id)  # Ensure ID is a factor
  proj4string(x_sp) <- proj_out
  # calculate kernel
  kerns <- kernelUD(x_sp, grid = grid_sp, h = h, same4all = FALSE)
  # calculate BA overlap
  BA <- kerneloverlaphr(kerns, method = "BA", conditional = TRUE, percent = iso)
  # clean matrix (remove diagonal, duplicates, NAs)
  BA[lower.tri(BA, diag = TRUE)] <- NA
  BA_v <- na.omit(as.vector(BA))
  # scale BA to 1
  BA_v <- as.numeric(BA_v / (iso / 100))
  return(BA_v)
}
# produces different values to function above

pairwise_BA_wo_ade <- pairwise_BA_ade <- function(x, h, cellsize, iso, proj_in, proj_out, buffer_dist = h*2) {
  # just taking the ID and longitude and latitude
  x <- x[,c("longitude", "latitude", "Id")]
  # re-projecting points using sf
  x_sf <- st_as_sf(x, coords = c("longitude", "latitude"))
  st_crs(x_sf) <- proj_in
  if (proj_in != proj_out) {
    x_sf <- st_transform(x_sf, crs = proj_out)
  }
  # get bbox as numeric
  bbox <- st_bbox(x_sf)
  xmin <- as.numeric(bbox["xmin"]) - buffer_dist
  xmax <- as.numeric(bbox["xmax"]) + buffer_dist
  ymin <- as.numeric(bbox["ymin"]) - buffer_dist
  ymax <- as.numeric(bbox["ymax"]) + buffer_dist
  # create raster grid
  grid_raster <- rast(
    xmin = xmin, xmax = xmax,
    ymin = ymin, ymax = ymax,
    resolution = cellsize, crs = proj_out
  )
  # convert sf to SpatVector with ID
  x_vect <- vect(x_sf)
  x_vect$Id <- as.factor(x$Id)  # Ensure ID is a factor
  crs(x_vect) <- proj_out
  
  # Manually calculate kernel density estimation
  kde <- function(points, grid, h) {
    density <- terra::density(points, grid, bandwidth = h)
    return(density)
  }
  
  # Calculate kernel density for each ID
  ids <- unique(x_vect$Id)
  kerns <- lapply(ids, function(id) {
    points <- x_vect[x_vect$Id == id, ]
    kde(points, grid_raster, h)
  })
  
  # Calculate BA overlap manually
  kernel_overlap <- function(kerns, iso) {
    n <- length(kerns)
    BA <- matrix(NA, n, n)
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        overlap <- sum(pmin(kerns[[i]], kerns[[j]])) / sum(kerns[[i]])
        BA[i, j] <- overlap * (iso / 100)
      }
    }
    return(BA)
  }
  
  BA <- kernel_overlap(kerns, iso)
  # clean matrix (remove diagonal, duplicates, NAs)
  BA[lower.tri(BA, diag = TRUE)] <- NA
  BA_v <- na.omit(as.vector(BA))
  return(BA_v)
}
# doesnt work!


ap.10 <- read.csv("./Data_inputs/gls_data/ap.2010.csv")
ap.10$Id <- as.factor(as.character(ap.10$Id))
ap.10$longitude <- ap.10$Lon
ap.10$latitude <- ap.10$Lat
# setting utm projection 
proj.utm <- paste0("+proj=laea +lat_0=", median(ap.10$latitude), " +lon_0=", median(ap.10$longitude), " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
# trying process without adehabitat
ap10_ade <- pairwise_BA_ade(x = ap.10, h = 186000, cellsize = 50000, iso = 50,
                               proj_in = proj.dec, proj_out = proj.utm)
median(ap10_ade)

ap10_wo_ade <- function(x, h, cellsize, iso, proj_in, proj_out, buffer_dist = h*2) {
  # just taking the ID and longitude and latitude
  x <- x[,c("longitude", "latitude", "Id")]
  # re-projecting points using sf
  x_sf <- st_as_sf(x, coords = c("longitude", "latitude"))
  st_crs(x_sf) <- proj_in
  if (proj_in != proj_out) {
    x_sf <- st_transform(x_sf, crs = proj_out)
  }
  # get bbox as numeric
  bbox <- st_bbox(x_sf)
  xmin <- as.numeric(bbox["xmin"]) - buffer_dist
  xmax <- as.numeric(bbox["xmax"]) + buffer_dist
  ymin <- as.numeric(bbox["ymin"]) - buffer_dist
  ymax <- as.numeric(bbox["ymax"]) + buffer_dist
  # create raster grid
  grid_raster <- rast(
    xmin = xmin, xmax = xmax,
    ymin = ymin, ymax = ymax,
    resolution = cellsize, crs = proj_out
  )
  values(grid_raster) <- 0  # Initialize raster with zeros
  
  # convert sf to SpatVector with ID
  x_vect <- vect(x_sf)
  x_vect$Id <- as.factor(x$Id)  # Ensure ID is a factor
  crs(x_vect) <- proj_out
  
  # Manually calculate kernel density estimation
  kde <- function(points, grid, h) {
    density <- terra::density(points, grid, bandwidth = h)
    return(density)
  }
  
  # Calculate kernel density for each ID
  ids <- unique(x_vect$Id)
  kerns <- lapply(ids, function(id) {
    points <- x_vect[x_vect$Id == id, ]
    coords <- as.matrix(geom(points)[, c("x", "y")])
    kde(coords, h, grid_raster)
  })
  
  # Calculate BA overlap manually
  kernel_overlap <- function(kerns, iso) {
    n <- length(kerns)
    BA <- matrix(NA, n, n)
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        overlap <- sum(pmin(values(kerns[[i]]), values(kerns[[j]]))) / sum(values(kerns[[i]]))
        BA[i, j] <- overlap * (iso / 100)
      }
    }
    return(BA)
  }
  
  BA <- kernel_overlap(kerns, iso)
  # clean matrix (remove diagonal, duplicates, NAs)
  BA[lower.tri(BA, diag = TRUE)] <- NA
  BA_v <- na.omit(as.vector(BA))
  return(BA_v)
}

median(ap10_wo_ade)



# 0.4371286