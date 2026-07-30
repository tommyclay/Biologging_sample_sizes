library(track2KBA)
library(adehabitatHR)
library(sf)

#### 1. SET UP DIRECTORY AND PROJECTIONS ####

in_direx <- "./Data_inputs/gls_data/"
f <- dir(in_direx)

proj_utm <- "+proj=laea +lat_0=-90 +lon_0=-38.05 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj_dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"   


#### 2. FOR EACH DATASET (JUST GLS DATA FOR NOW) = LOAD IN DATA, CONVERT TO KERNEL AND RUN REPRESENTATIVENESS CODE ####


iter = 1000
iso = 50
h = 186000
cellsize = 50000

out_l <- list()
for (i in 1:length(f)) {
  print(gsub(".csv", "", f[i]))
  # read in file
  dat <- read.csv(paste0(in_direx, f[i]))
  dat$id <- as.factor(as.character(dat$id))
  # project to utm
  lox <- data.frame(x = dat$longitude, y = dat$latitude, ID = dat$id)
  coordinates(lox) <-  ~x + y
  proj4string(lox) <- proj_dec
  lox <- spTransform(lox, CRS(proj_utm))
  # create null grid
  null.grid <- expand.grid(x=seq(min(coordinates(lox)[,1])-(h*2), max(coordinates(lox)[,1])+(h*2), by = cellsize),
                           y=seq(min(coordinates(lox)[,2])-(h*2), max(coordinates(lox)[,2])+(h*2), by = cellsize))
  coordinates(null.grid) <- ~x+y
  gridded(null.grid) <- TRUE
  # create individual kernel using adehabitatHR
  UD_in <- kernelUD(lox, grid = null.grid, h = h, same4all = FALSE)
  # running representativeness
  out_l[[i]] <- repAssess(tracks = lox, KDE = UD_in, iteration = iter, 
                            levelUD = iso, avgMethod = "mean", nCores = 4, bootTable = FALSE)
  # output
  out_l[[i]]$Species <- dat$species[1]
  out_l[[i]]$Year <- dat$year[1]
  out_l[[i]]$Dataset <- gsub(".csv", "", f[i])
}
all_out <- do.call(rbind, out_l)


#### 3. PASTE OUT - DON'T PASTE OUT FOR NOW ####

out_name <- "./Data_outputs/birdlife_predictions.csv"
# write.csv(all_out, out_name)
