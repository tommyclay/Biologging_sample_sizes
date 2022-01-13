library(adehabitatHR)
library(raster)
library(maptools)
library(rgeos)
library(plyr)
library(ggplot2)
library(sf)
library(progress)


######### 1. CODE TO RESAMPLE BOTH NEW AND TOTAL HOME-RANGE AREA WITH INCREASING SAMPLE SIZE ##################



# projection formats - utm coordinates centered around South Pole
colony <- data.frame(Longitude = -38.05, Latitude = -90)
proj.utm <- paste0("+proj=laea +lat_0=", colony$Latitude, " +lon_0=", colony$Longitude, " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj.dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"   

# function to project from dec to utm
proj_fun <- function(df_in, proj_in, proj_out) {
  coordinates(df_in) <- ~Longitude+Latitude
  proj4string(df_in) <- proj_in
  df_out <- spTransform(df_in, CRS(proj_out))
  df_out <- data.frame(df_out)
  return(df_out)
}


#### PARALLELIZE WITHIN LOOP ####


# function to resample home range area with sample size
resampling_UD_area <- function(x, n_iter, iso_level, cellsize, h) {
  
  # 1. SETTING UP NULL GRID - using buffer of 2 times the smoothing factor around min and max points 
  null.grid <- expand.grid(x=seq(min(x$Longitude)-(h*1000*2), max(x$Longitude)+(h*1000*2), by = cellsize*1000),
                           y=seq(min(x$Latitude)-(h*1000*2), max(x$Latitude)+(h*1000*2), by = cellsize*1000))
  coordinates(null.grid) <- ~x+y
  gridded(null.grid) <- TRUE

  # 2. CALCULATE KERNELS FOR EACH INDIVIDUAL AND TURN TO MATRICES (TO AVERAGE ACROSS LATER)
  kern_in <- list()
  matrix_in <- list()
  x$Id <- as.factor(as.character(x$Id))
  for (j in 1:nlevels(x$Id)) {
    id.sub <- subset(x, Id==levels(x$Id)[j])   
    id.sub$Id <- as.factor(as.character(id.sub$Id))
    lox <- id.sub[,c("Longitude", "Latitude")]
    coordinates(lox) <-  ~Longitude + Latitude
    kern_in[[j]] <- adehabitatHR::kernelUD(lox, grid = null.grid, h = h*1000)
    matrix_in[[j]] <- as.matrix(raster(kern_in[[j]]))
  }
  names(kern_in) <- levels(x$Id) 
  names(matrix_in) <- levels(x$Id)
  
  # 3. RESAMPLING HR AREA THROUGH LEVELS OF SAMPLE SIZE
  # randomly sample an individual without replacement and consecutively add another individual and merge together
  
  # set up progress bar
  pb <- progress_bar$new(format = "[:bar] :current/:total (Time remaining::eta)", total = n_iter)
  
  iter.summary <- list() # output summary dataframe for each iteration
  
  # run through iterations...
  for (i in 1:n_iter) {
    # for progress bar
    pb$tick()
    Sys.sleep(1 / n_iter)
    # set up lists and output dataframes per level of sample size
    kern_out <- list()
    iso_out <- list()
    matrix_out <- list()
    iter.summary[[i]] <- data.frame(Sample = seq(1, nlevels(x$Id), 1),Tot_area = NA, New_area = NA, Perc_new = NA, 
                                Iteration = rep(i), Isopleth = rep(iso_level))
    # randomize order of ids
    random.ids <- sample(levels(x$Id),nlevels(x$Id), replace = FALSE, prob = NULL) # random ids without replacement
 
    
    #  FOR FIRST LEVEL OF SAMPLE SIZE (IE, N = 1)
    matrix_out[[1]] <- subset(matrix_in, names(matrix_in)==random.ids[1])[[1]]
    kern_out[[1]] <- subset(kern_in, names(kern_in)==random.ids[1])[[1]]
    iso_out[[1]] <- adehabitatHR::getverticeshr(kern_out[[1]], percent = iso_level)
    iter.summary[[i]]$Tot_area[1] <- st_area(st_as_sf(iso_out[[1]]))
    iter.summary[[i]]$New_area[1] <- st_area(st_as_sf(iso_out[[1]]))
    iter.summary[[i]]$Perc_new[1] <- iter.summary[[i]]$New_area[1]/iter.summary[[i]]$Tot_area[1]*100
    
    #  FOR OTHER LEVEL OF SAMPLE SIZE (IE, N = 2+)
    for (a in 2:nlevels(x$Id)) { 
      # averaging values across matrices
      matrix_out[[a]] <- subset(matrix_in, names(matrix_in)==random.ids[a])[[1]]
      matrix_sum <- Reduce('+', matrix_out)
      matrix_sum <- matrix_sum/a
      # create kernel output, paste in merged values and calculate isopleth 
      kern_out[[a]] <- kern_out[[1]]
      kern_out[[a]]@data$ud <- as.vector(matrix_sum)
      iso_out[[a]] <- adehabitatHR::getverticeshr(kern_out[[a]], percent = iso_level)
      # intersect with previous polygon to get new area occupied
      iter.summary[[i]]$Tot_area[a] <- st_area(st_as_sf(iso_out[[a]]))
      polydiff <- rgeos::gDifference(iso_out[[a]], iso_out[[a-1]])
      if (length(polydiff) == 0) { polydiff <- gDifference(iso_out[[a-1]], iso_out[[a]]) # if new area is contained within previous polygon, then value is minus
          iter.summary[[i]]$New_area[a] <- st_area(st_as_sf(polydiff))*-1
      } else { iter.summary[[i]]$New_area[a] <- st_area(st_as_sf(polydiff)) }
      # paste out percentage of new area as a percentage of total area
      iter.summary[[i]]$Perc_new[a] <- iter.summary[[i]]$New_area[a]/iter.summary[[i]]$Tot_area[a]*100
    }
  }
  all.iter <- do.call(rbind, iter.summary)
  # convert area measurements back to km2 from m2
  all.iter$Tot_area <- all.iter$Tot_area/1000000
  all.iter$New_area <- all.iter$New_area/1000000
  return(all.iter)
}


# x = data.frame with 3 columns, "Longitude" and "Latitude" and individual "Id"
# n_iterations =  number of times you want to resample data - usually 1,000
# iso_level = isopleth level you want to cut at, e.g. 50. If you want at multiple levels then specify
# e.g. "c(50, 95)
# cellsize = cellsize of grid you want to input into kernels in km
# h = smoothing factor used in analysis in km (for GLS usually around 200 km)


############# LOADING DATA IN ####################


##### AP 2010 ######

ap.10 <- read.csv("./Data/ap.2010.csv")
ap.10$Id <- as.factor(as.character(ap.10$Id))
ap.10$Longitude <- ap.10$Lon
ap.10$Latitude <- ap.10$Lat
# reproject
ap.10.df <- proj_fun(df_in = ap.10, proj_in = proj.dec, proj_out = proj.utm) 

ap10_out <- resampling_UD_area(x = ap.10.df, n_iter = 1000, 
                                         iso_level = 50, cellsize = 50, h = 186)
# takes 4-5 minutes to run for 1000 iterations

# plotting
ap10_out$Tot_area <- ap10_out$Tot_area/1000000
ap10_out$New_area <- ap10_out$New_area/1000000
ggplot(ap10_out, aes(Sample, Tot_area)) + geom_line(aes(group = as.factor(Iteration)))+
  stat_smooth() + ylab("Home range area (106 km2)") + xlab("Sample size")
ggplot(ap10_out, aes(Sample, New_area)) + geom_line(aes(group = as.factor(Iteration)))+
  stat_smooth() + ylab("New home range area (106 km2)") + xlab("Sample size")


# setwd("C:/Users/tclay/Documents/Methods paper/Boots")
# write.csv(resample.out, "AP2010_Resample_50_1k.csv")
# 
# resample.out <- resampling.sample.size(x = ap.10.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
# setwd("C:/Users/tclay/Documents/Methods paper/Boots")
# write.csv(resample.out, "AP2010_Resample_95_1k.csv")
# 




##### AP 2011 ######


ap.11 <- read.csv("./Data/ap.2011.csv")
ap.11$Id <- as.factor(as.character(ap.11$Id))
ap.11$Longitude <- ap.11$Lon
ap.11$Latitude <- ap.11$Lat
# reproject
ap.11.df <- proj_fun(df_in = ap.11, proj_in = proj.dec, proj_out = proj.utm) 

ap11_out <- resampling_UD_area(x = ap.11.df, n_iter = 1000, 
                               iso_level = 50, cellsize = 50, h = 186)
# takes 3-4 minutes to run for 1000 iterations

# plotting
ap11_out$Tot_area <- ap11_out$Tot_area/1000000
ap11_out$New_area <- ap11_out$New_area/1000000
ggplot(ap11_out, aes(Sample, Tot_area)) + geom_line(aes(group = as.factor(Iteration)))+
  stat_smooth() + ylab("Home range area (106 km2)") + xlab("Sample size")
ggplot(ap11_out, aes(Sample, New_area)) + geom_line(aes(group = as.factor(Iteration)))+
  stat_smooth() + ylab("New home range area (106 km2)") + xlab("Sample size")


# 
# resample.out <- resampling.sample.size(x = ap.11.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
# setwd("C:/Users/tclay/Documents/Methods paper/Boots")
# write.csv(resample.out, "AP2011_Resample_95_1k.csv")
# 



######## BBA 2002 #######

bba.02 <- read.csv("./Data/bba.2002.csv")
bba.02$Id <- as.factor(as.character(bba.02$Id))
bba.02$Longitude <- bba.02$Lon
bba.02$Latitude <- bba.02$Lat
bba.02.df <- proj_fun(df_in = bba.02, proj_in = proj.dec, proj_out = proj.utm) 

bba02_out <- resampling_UD_area(x = bba.02.df, n_iter = 1000, 
                                         iso_level = 50, cellsize =50, h = 186)
# takes 45-55 mins for 1000 iterations

# plotting
bba02_out$Tot_area <- bba02_out$Tot_area/1000000
bba02_out$New_area <- bba02_out$New_area/1000000

p1 <- ggplot(bba02_out, aes(Sample, Tot_area)) + geom_line(aes(group = as.factor(Iteration)),
                                                           col = "grey70")+
  stat_smooth() + ylab("Home range area (106 km 2)") + xlab("Sample size") + theme_bw()+
  scale_y_continuous(limits = c(0, 3.2))
p1
p2 <- ggplot(bba02_out, aes(Sample, New_area)) + geom_line(aes(group = as.factor(Iteration)), 
                                                          col = "grey70")+
  stat_smooth() + ylab("New home range area (106 km 2)") + xlab("Sample size") + theme_bw()
p2
p3 <- ggplot(bba02_out, aes(Sample, Perc_new)) + geom_line(aes(group = as.factor(Iteration)), 
                                                           col = "grey70")+
  stat_smooth() + ylab("New home range area (% of total)") + xlab("Sample size") + theme_bw()+
  scale_y_continuous(limits = c(-30, 100), breaks = seq(-30, 100, 10)) +
  geom_hline(yintercept=5, linetype="dashed")
p3




# resample.out <- resampling.sample.size(x = bba.02.3, niter = 1000, iso = 50, cellsize =50000, h = 200000)
# setwd("C:/Users/tclay/Documents/Methods paper/Boots")
# write.csv(resample.out, "BBA2002_Resample_50_1k.csv")
# 
# resample.out <- resampling.sample.size(x = bba.02.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
# setwd("C:/Users/tclay/Documents/Methods paper/Boots")
# write.csv(resample.out, "BBA2002_Resample_95_1k.csv")




######## BBA 2003 #######


bba.03 <- read.csv("./Data/bba.2003.csv")
bba.03$Id <- as.factor(as.character(bba.03$Id))
bba.03$Longitude <- bba.03$Lon
bba.03$Latitude <- bba.03$Lat
bba.03.df <- proj_fun(df_in = bba.03, proj_in = proj.dec, proj_out = proj.utm) 

bba03_out <- resampling_UD_area(x = bba.03.df, n_iter = 1000, 
                                iso_level = 50, cellsize = 50, h = 186)
# takes 45-50 mins for 1000 iterations

# plotting
bba03_out$Tot_area <- bba03_out$Tot_area/1000000
bba03_out$New_area <- bba03_out$New_area/1000000

p1 <- ggplot(bba03_out, aes(Sample, Tot_area)) + geom_line(aes(group = as.factor(Iteration)),
                                                           col = "grey70")+
  stat_smooth() + ylab("Home range area (106 km 2)") + xlab("Sample size") + theme_bw()+
  scale_y_continuous(limits = c(0, 3.2))
p1
p2 <- ggplot(bba03_out, aes(Sample, New_area)) + geom_line(aes(group = as.factor(Iteration)), 
                                                           col = "grey70")+
  stat_smooth() + ylab("New home range area (106 km 2)") + xlab("Sample size") + theme_bw()
p2
p3 <- ggplot(bba03_out, aes(Sample, Perc_new)) + geom_line(aes(group = as.factor(Iteration)), 
                                                           col = "grey70")+
  stat_smooth() + ylab("New home range area (% of total)") + xlab("Sample size") + theme_bw()+
  scale_y_continuous(limits = c(-30, 100), breaks = seq(-30, 100, 10)) +
  geom_hline(yintercept=5, linetype="dashed")
p3

# setwd("C:/Users/tclay/Documents/Methods paper/Boots")
# write.csv(resample.out, "BBA2003_Resample_50_1k.csv")
# 
# resample.out <- resampling.sample.size(x = bba.03.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
# setwd("C:/Users/tclay/Documents/Methods paper/Boots")
# write.csv(resample.out, "BBA2003_Resample_95_1k.csv")



##### BP 2011 ######


bp.11 <- read.csv("./Data/bp.2011.csv")
bp.11$Id <- as.factor(as.character(bp.11$Id))
bp.11$Longitude <- bp.11$Lon
bp.11$Latitude <- bp.11$Lat
bp.11.df <- proj_fun(df_in = bp.11, proj_in = proj.dec, proj_out = proj.utm) 

bp11_out <- resampling_UD_area(x = bp.11.df, n_iter = 1000, 
                                iso_level = 50, cellsize = 50, h = 186)
# takes 11-13 mins for 1000 iterations

# plotting
bp11_out$Tot_area <- bp11_out$Tot_area/1000000
bp11_out$New_area <- bp11_out$New_area/1000000

p1 <- ggplot(bp11_out, aes(Sample, Tot_area)) + geom_line(aes(group = as.factor(Iteration)),
                                                           col = "grey70")+
  stat_smooth() + ylab("Home range area (106 km 2)") + xlab("Sample size") + theme_bw()+
  scale_y_continuous(limits = c(0, 3.2))
p1
p2 <- ggplot(bp11_out, aes(Sample, New_area)) + geom_line(aes(group = as.factor(Iteration)), 
                                                           col = "grey70")+
  stat_smooth() + ylab("New home range area (106 km 2)") + xlab("Sample size") + theme_bw()
p2
p3 <- ggplot(bp11_out, aes(Sample, Perc_new)) + geom_line(aes(group = as.factor(Iteration)), 
                                                           col = "grey70")+
  stat_smooth() + ylab("New home range area (% of total)") + xlab("Sample size") + theme_bw()+
  scale_y_continuous(limits = c(-30, 100), breaks = seq(-30, 100, 10)) +
  geom_hline(yintercept=5, linetype="dashed")
p3

# setwd("C:/Users/tclay/Documents/Methods paper/Boots")
# write.csv(resample.out, "BP2011_Resample_50_1k.csv")
# 
# resample.out <- resampling.sample.size(x = bp.11.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
# setwd("C:/Users/tclay/Documents/Methods paper/Boots")
# write.csv(resample.out, "BP2011_Resample_95_1k.csv")




##### BS 2002 ######


bs.02 <- read.csv("./Data/bs.2002.csv")
bs.02$Id <- as.factor(as.character(bs.02$Id))
bs.02$Longitude <- bs.02$Lon
bs.02$Latitude <- bs.02$Lat
bs.02.df <- proj_fun(df_in = bs.02, proj_in = proj.dec, proj_out = proj.utm) 

bs02_out <- resampling_UD_area(x = bs.02.df, n_iter = 1000, 
                               iso_level = 50, cellsize = 50, h = 186)
# takes 3 mins for 1000 iterations

# plotting
bs02_out$Tot_area <- bs02_out$Tot_area/1000000
bs02_out$New_area <- bs02_out$New_area/1000000

p1 <- ggplot(bs02_out, aes(Sample, Tot_area)) + geom_line(aes(group = as.factor(Iteration)),
                                                          col = "grey70")+
  stat_smooth() + ylab("Home range area (106 km 2)") + xlab("Sample size") + theme_bw()+
  scale_y_continuous(limits = c(0, 3.2))
p1
p2 <- ggplot(bs02_out, aes(Sample, New_area)) + geom_line(aes(group = as.factor(Iteration)), 
                                                          col = "grey70")+
  stat_smooth() + ylab("New home range area (106 km 2)") + xlab("Sample size") + theme_bw()
p2
p3 <- ggplot(bs02_out, aes(Sample, Perc_new)) + geom_line(aes(group = as.factor(Iteration)), 
                                                          col = "grey70")+
  stat_smooth() + ylab("New home range area (% of total)") + xlab("Sample size") + theme_bw()+
  scale_y_continuous(limits = c(-30, 100), breaks = seq(-30, 100, 10)) +
  geom_hline(yintercept=5, linetype="dashed")
p3

# setwd("C:/Users/tclay/Documents/Methods paper/Boots")
# write.csv(resample.out, "BS2002_Resample_50_1k.csv")
# 
# resample.out <- resampling.sample.size(x = bs.02.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
# setwd("C:/Users/tclay/Documents/Methods paper/Boots")
# write.csv(resample.out, "BS2002_Resample_95_1k.csv")





##### BS 2012 ######


bs.12 <- read.csv("./Data/bs.2012.csv")
bs.12$Id <- as.factor(as.character(bs.12$Id))
bs.12$Longitude <- bs.12$Lon
bs.12$Latitude <- bs.12$Lat
bs.12.df <- proj_fun(df_in = bs.12, proj_in = proj.dec, proj_out = proj.utm) 

bs12_out <- resampling_UD_area(x = bs.12.df, n_iter = 1000, 
                               iso_level = 50, cellsize = 50, h = 186)
# takes 17-20 mins for 1000 iterations

# plotting
bs12_out$Tot_area <- bs12_out$Tot_area/1000000
bs12_out$New_area <- bs12_out$New_area/1000000

p1 <- ggplot(bs12_out, aes(Sample, Tot_area)) + geom_line(aes(group = as.factor(Iteration)),
                                                          col = "grey70")+
  stat_smooth() + ylab("Home range area (106 km 2)") + xlab("Sample size") + theme_bw()+
  scale_y_continuous(limits = c(0, 3.2))
p1
p2 <- ggplot(bs12_out, aes(Sample, New_area)) + geom_line(aes(group = as.factor(Iteration)), 
                                                          col = "grey70")+
  stat_smooth() + ylab("New home range area (106 km 2)") + xlab("Sample size") + theme_bw()
p2
p3 <- ggplot(bs12_out, aes(Sample, Perc_new)) + geom_line(aes(group = as.factor(Iteration)), 
                                                          col = "grey70")+
  stat_smooth() + ylab("New home range area (% of total)") + xlab("Sample size") + theme_bw()+
  scale_y_continuous(limits = c(-30, 100), breaks = seq(-30, 100, 10)) +
  geom_hline(yintercept=5, linetype="dashed")
p3



##### GHA 1999 ######

gha.99 <- read.csv("./Data/gha.1999.csv")
gha.99$Id <- as.factor(as.character(gha.99$Id))
gha.99$Longitude <- gha.99$Lon
gha.99$Latitude <- gha.99$Lat
gha.99.df <- proj_fun(df_in = gha.99, proj_in = proj.dec, proj_out = proj.utm) 

gha99_out <- resampling_UD_area(x = gha.99.df, n_iter = 1000, 
                               iso_level = 50, cellsize = 50, h = 186)
# takes over 1 hour for 1000 iterations

# plotting
bs12_out$Tot_area <- bs12_out$Tot_area/1000000
bs12_out$New_area <- bs12_out$New_area/1000000

p1 <- ggplot(bs12_out, aes(Sample, Tot_area)) + geom_line(aes(group = as.factor(Iteration)),
                                                          col = "grey70")+
  stat_smooth() + ylab("Home range area (106 km 2)") + xlab("Sample size") + theme_bw()+
  scale_y_continuous(limits = c(0, 3.2))
p1
p2 <- ggplot(bs12_out, aes(Sample, New_area)) + geom_line(aes(group = as.factor(Iteration)), 
                                                          col = "grey70")+
  stat_smooth() + ylab("New home range area (106 km 2)") + xlab("Sample size") + theme_bw()
p2
p3 <- ggplot(bs12_out, aes(Sample, Perc_new)) + geom_line(aes(group = as.factor(Iteration)), 
                                                          col = "grey70")+
  stat_smooth() + ylab("New home range area (% of total)") + xlab("Sample size") + theme_bw()+
  scale_y_continuous(limits = c(-30, 100), breaks = seq(-30, 100, 10)) +
  geom_hline(yintercept=5, linetype="dashed")
p3


### PLOTTING FOR 95% KERNEL ###

# median for each sample size
gha95_med <- ddply(gha95_out,.(Sample), summarize, Tot_area_med = median(Tot_area), New_area_med = median(New_area),
                   Perc_new_med = median(Perc_new))

p1 <- ggplot() + geom_line(data = gha95_out, aes(x=Sample, y=Tot_area, group = as.factor(Iteration)),
                           col = "grey70")+
  stat_smooth(data = gha95_out, aes(x=Sample, y=Tot_area)) + 
  geom_point(data = gha95_med, aes(x = Sample, y = Tot_area_med), pch = 21) + 
  ylab("Home range area (millions km 2)") + xlab("Sample size") + theme_bw() +
  scale_y_continuous(limits = c(0, 40))
p1

p2 <- ggplot() + geom_line(data = gha95_out, aes(x = Sample, y = New_area, group = as.factor(Iteration)), 
                           col = "grey70")+
  stat_smooth(data = gha95_out, aes(x=Sample, y=New_area)) + 
  geom_point(data = gha95_med, aes(x = Sample, y = New_area_med), pch = 21) + 
  ylab("New home range area (millions km 2)") + xlab("Sample size") + theme_bw()
#geom_hline(yintercept=0)
p2

p3 <- ggplot() + geom_line(data = gha95_out, aes(x = Sample, y = Perc_new, group = as.factor(Iteration)), 
                           col = "grey70")+
  stat_smooth(data = gha95_out, aes(x=Sample, y=Perc_new)) + 
  geom_point(data = gha95_med, aes(x = Sample, y = Perc_new_med), pch = 21) + 
  ylab("New home range area (% of total)") + xlab("Sample size") + theme_bw()+
  scale_y_continuous(limits = c(-5, 100), breaks = seq(-10, 100, 10)) +
  geom_hline(yintercept=5, linetype="dashed")+geom_hline(yintercept=0)
p3


# ls.mod <- nls(Perc_new ~ Sample, data = gha99_out)
# nls(Perc_new~Sample*exp(-b*x),start=list(a=a_start,b=b_start))


setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "GHA1999_Resample_95_1k.csv")




######## LMA 2004 #######


setwd("C:/Users/tclay/Documents/Methods paper/Data")
lma.04 <- read.csv("lmsa.2004.csv")
lma.04$Id <- as.factor(as.character(lma.04$Id))

lma.04.2 <- lma.04
coordinates(lma.04.2) <- ~Lon+Lat
proj4string(lma.04.2) <- projection.dec
lma.04.2 <- spTransform(lma.04.2, CRS(projection))
lma.04.3 <- data.frame(lma.04.2)
lma.04.3$Longitude <- lma.04.3$Lon
lma.04.3$Latitude <- lma.04.3$Lat

lma.04.df <- as.data.frame(lma.04.3)
lma.04.df$Id <- as.factor(as.character(lma.04.df$Id))

resample.out <- resampling.sample.size(x = lma.04.3, niter = 1000, iso = 50, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "LMA2004_Resample_50_1k.csv")

resample.out <- resampling.sample.size(x = lma.04.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "LMA2004_Resample_95_1k.csv")



######## LMA 2005 #######


setwd("C:/Users/tclay/Documents/Methods paper/Data")
lma.05 <- read.csv("lmsa.2005.csv")
lma.05$Id <- as.factor(as.character(lma.05$Id))

lma.05.2 <- lma.05
coordinates(lma.05.2) <- ~Lon+Lat
proj4string(lma.05.2) <- projection.dec
lma.05.2 <- spTransform(lma.05.2, CRS(projection))
lma.05.3 <- data.frame(lma.05.2)
lma.05.3$Longitude <- lma.05.3$Lon
lma.05.3$Latitude <- lma.05.3$Lat

lma.05.df <- as.data.frame(lma.05.3)
lma.05.df$Id <- as.factor(as.character(lma.05.df$Id))

resample.out <- resampling.sample.size(x = lma.05.3, niter = 1000, iso = 50, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "LMA2005_Resample_50_1k.csv")

resample.out <- resampling.sample.size(x = lma.05.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "LMA2005_Resample_95_1k.csv")





######## NGP 2000 #######


setwd("C:/Users/tclay/Documents/Methods paper/Data")
ngp.00 <- read.csv("ngp.2000.csv")
ngp.00$Id <- as.factor(as.character(ngp.00$Id))

ngp.00.2 <- ngp.00
coordinates(ngp.00.2) <- ~Lon+Lat
proj4string(ngp.00.2) <- projection.dec
ngp.00.2 <- spTransform(ngp.00.2, CRS(projection))
ngp.00.3 <- data.frame(ngp.00.2)
ngp.00.3$Longitude <- ngp.00.3$Lon
ngp.00.3$Latitude <- ngp.00.3$Lat

ngp.00.df <- as.data.frame(ngp.00.3)
ngp.00.df$Id <- as.factor(as.character(ngp.00.df$Id))

resample.out <- resampling.sample.size(x = ngp.00.3, niter = 1000, iso = 50, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "NGP2000_Resample_50_1k.csv")

resample.out <- resampling.sample.size(x = ngp.00.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "NGP2000_Resample_95_1k.csv")





######## SGP 2000 #######


setwd("C:/Users/tclay/Documents/Methods paper/Data")
sgp.00 <- read.csv("sgp.2000.csv")
sgp.00$Id <- as.factor(as.character(sgp.00$Id))

sgp.00.2 <- sgp.00
coordinates(sgp.00.2) <- ~Lon+Lat
proj4string(sgp.00.2) <- projection.dec
sgp.00.2 <- spTransform(sgp.00.2, CRS(projection))
sgp.00.3 <- data.frame(sgp.00.2)
sgp.00.3$Longitude <- sgp.00.3$Lon
sgp.00.3$Latitude <- sgp.00.3$Lat

sgp.00.df <- as.data.frame(sgp.00.3)
sgp.00.df$Id <- as.factor(as.character(sgp.00.df$Id))

resample.out <- resampling.sample.size(x = sgp.00.3, niter = 1000, iso = 50, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "SGP2000_Resample_50_1k.csv")

resample.out <- resampling.sample.size(x = sgp.00.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "SGP2000_Resample_95_1k.csv")






######## WA 2004 #######


setwd("C:/Users/tclay/Documents/Methods paper/Data")
wa.04 <- read.csv("wa.2004.csv")
wa.04$Id <- as.factor(as.character(wa.04$Id))

wa.04.2 <- wa.04
coordinates(wa.04.2) <- ~Lon+Lat
proj4string(wa.04.2) <- projection.dec
wa.04.2 <- spTransform(wa.04.2, CRS(projection))
wa.04.3 <- data.frame(wa.04.2)
wa.04.3$Longitude <- wa.04.3$Lon
wa.04.3$Latitude <- wa.04.3$Lat

wa.04.df <- as.data.frame(wa.04.3)
wa.04.df$Id <- as.factor(as.character(wa.04.df$Id))

resample.out <- resampling.sample.size(x = wa.04.3, niter = 1000, iso = 50, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "WA2004_Resample_50_1k.csv")

resample.out <- resampling.sample.size(x = wa.04.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "WA2004_Resample_95_1k.csv")


######## WA 2008 #######


setwd("C:/Users/tclay/Documents/Methods paper/Data")
wa.08 <- read.csv("wa.2008.csv")
wa.08$Id <- as.factor(as.character(wa.08$Id))

wa.08.2 <- wa.08
coordinates(wa.08.2) <- ~Lon+Lat
proj4string(wa.08.2) <- projection.dec
wa.08.2 <- spTransform(wa.08.2, CRS(projection))
wa.08.3 <- data.frame(wa.08.2)
wa.08.3$Longitude <- wa.08.3$Lon
wa.08.3$Latitude <- wa.08.3$Lat

wa.08.df <- as.data.frame(wa.08.3)
wa.08.df$Id <- as.factor(as.character(wa.08.df$Id))

resample.out <- resampling.sample.size(x = wa.08.3, niter = 1000, iso = 50, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "WA2008_Resample_50_1k.csv")

resample.out <- resampling.sample.size(x = wa.08.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "WA2008_Resample_95_1k.csv")



######## WA 2009 #######


setwd("C:/Users/tclay/Documents/Methods paper/Data")
wa.09 <- read.csv("wa.2009.csv")
wa.09$Id <- as.factor(as.character(wa.09$Id))

wa.09.2 <- wa.09
coordinates(wa.09.2) <- ~Lon+Lat
proj4string(wa.09.2) <- projection.dec
wa.09.2 <- spTransform(wa.09.2, CRS(projection))
wa.09.3 <- data.frame(wa.09.2)
wa.09.3$Longitude <- wa.09.3$Lon
wa.09.3$Latitude <- wa.09.3$Lat

wa.09.df <- as.data.frame(wa.09.3)
wa.09.df$Id <- as.factor(as.character(wa.09.df$Id))

resample.out <- resampling.sample.size(x = wa.09.3, niter = 1000, iso = 50, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "WA2009_Resample_50_1k.csv")

resample.out <- resampling.sample.size(x = wa.09.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "WA2009_Resample_95_1k.csv")





######## WCP 2003 ########



setwd("C:/Users/tclay/Documents/Methods paper/Data")
wcp.03 <- read.csv("wcp.2003.csv")
wcp.03$Id <- as.factor(as.character(wcp.03$Id))
wcp.03$Date.time <- as.POSIXct(wcp.03$Date.time, tz = "GMT")
wcp.03$Date <- as.Date(wcp.03$Date, tz = "GMT")


wcp.03.2 <- wcp.03
coordinates(wcp.03.2) <- ~Lon+Lat
proj4string(wcp.03.2) <- projection.dec
wcp.03.2 <- spTransform(wcp.03.2, CRS(projection))
wcp.03.3 <- data.frame(wcp.03.2)
wcp.03.3$Longitude <- wcp.03.3$Lon
wcp.03.3$Latitude <- wcp.03.3$Lat

wcp.03.df <- as.data.frame(wcp.03.3)
wcp.03.df$Id <- as.factor(as.character(wcp.03.df$Id))

resample.out <- resampling.sample.size(x = wcp.03.3, niter = 1000, iso = 50, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "WCP2003_Resample_50_1k.csv")

resample.out <- resampling.sample.size(x = wcp.03.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "WCP2003_Resample_95_1k.csv")




######## WCP 2013 ########

setwd("C:/Users/tclay/Documents/Methods paper/Data")
wcp.13 <- read.csv("wcp.2013.csv")
wcp.13$Year <- 2013
wcp.13$Id <- as.factor(as.character(wcp.13$Id))
wcp.13$Date.time <- as.POSIXct(wcp.13$Date.time, tz = "GMT")
wcp.13$Date <- as.Date(wcp.13$Date, tz = "GMT")
ddply(wcp.13,.(Id), summarize, Start = Date[1], End = Date[length(Date)])
# remove c358 as only lasts a month
wcp.13 <- subset(wcp.13, Id != "C538",)
wcp.13$Id <- as.factor(as.character(wcp.13$Id))
nlevels(wcp.13$Id)
# 9

wcp.13.2 <- wcp.13
coordinates(wcp.13.2) <- ~Lon+Lat
proj4string(wcp.13.2) <- projection.dec
wcp.13.2 <- spTransform(wcp.13.2, CRS(projection))
wcp.13.3 <- data.frame(wcp.13.2)
wcp.13.3$Longitude <- wcp.13.3$Lon
wcp.13.3$Latitude <- wcp.13.3$Lat

wcp.13.df <- as.data.frame(wcp.13.3)
wcp.13.df$Id <- as.factor(as.character(wcp.13.df$Id))

resample.out <- resampling.sample.size(x = wcp.13.3, niter = 1000, iso = 50, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "WCP2013_Resample_50_1k.csv")

resample.out <- resampling.sample.size(x = wcp.13.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "WCP2013_Resample_95_1k.csv")



##### AP ALL YRS #####

setwd("C:/Users/tclay/Documents/Methods paper/Data")
ap.10 <- read.csv("ap.2010.csv")
ap.10$Id <- as.factor(as.character(ap.10$Id))
ap.11 <- read.csv("ap.2011.csv")
ap.11$Id <- as.factor(as.character(ap.11$Id))
levels(ap.10$Id)
# "DK33302" "DK33303" "DK33304" "DK33305" "DK33306" "DK33308" "DK33309" "DK33310" "DK33314"
levels(ap.11$Id)
# "DK33309" "DK33386" "DK33389" "DK33395" "DK33398" "DK33409"
# there is one individual DK33309 tracked over 2 years - remove from 2010 data
#ggplot(ap.10, aes(Lon, Lat)) + geom_point() + facet_wrap(~Id)
#ggplot(ap.11, aes(Lon, Lat)) + geom_point() + facet_wrap(~Id)
ap.10.2 <- subset(ap.10, Id != "DK33309",)
ap.10.2$Id <- as.factor(as.character(ap.10.2$Id))

ap.both <- rbind(ap.10.2, ap.11)
ap.both$Id  <- as.factor(as.character(ap.both$Id))
ap.both.2 <- ap.both
coordinates(ap.both.2) <- ~Lon+Lat
proj4string(ap.both.2) <- projection.dec
ap.both.2 <- spTransform(ap.both.2, CRS(projection))
ap.both.3 <- data.frame(ap.both.2)
ap.both.3$Longitude <- ap.both.3$Lon
ap.both.3$Latitude <- ap.both.3$Lat

ap.both.df <- as.data.frame(ap.both.3)
ap.both.df$Id <- as.factor(as.character(ap.both.df$Id))

resample.out <- resampling.sample.size(x = ap.both.3, niter = 1000, iso = 50, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "APBOTH_Resample_50_1k.csv")

resample.out <- resampling.sample.size(x = ap.both.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "APBOTH_Resample_95_1k.csv")



##### BS ALL YRS ######



setwd("C:/Users/tclay/Documents/Methods paper/Data")
bs.02 <- read.csv("bs.2002.csv")
bs.02$Id <- as.factor(as.character(bs.02$Id))
bs.12 <- read.csv("bs.2012.csv")
bs.12$Id <- as.factor(as.character(bs.12$Id))
levels(bs.02$Id)
# "27706" "67760" "67768" "67799" "67912" "67916"
levels(bs.12$Id)
# "RA00" "RA01" "RA02" "RA03" "RA05" "RA06" "RA07" "RA10" "RA11" "RA14" "RA17" "RA20" "RA36" "RA37" "RA38" "RA39" "RA44"
# "RA45" "RA46"

bs.both <- rbind(bs.02, bs.12)
bs.both$Id  <- as.factor(as.character(bs.both$Id))

bs.both.2 <- bs.both
coordinates(bs.both.2) <- ~Lon+Lat
proj4string(bs.both.2) <- projection.dec
bs.both.2 <- spTransform(bs.both.2, CRS(projection))
bs.both.3 <- data.frame(bs.both.2)
bs.both.3$Longitude <- bs.both.3$Lon
bs.both.3$Latitude <- bs.both.3$Lat

bs.both.df <- as.data.frame(bs.both.3)
bs.both.df$Id <- as.factor(as.character(bs.both.df$Id))

resample.out <- resampling.sample.size(x = bs.both.3, niter = 1000, iso = 50, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "BSBOTH_Resample_50_1k.csv")

resample.out <- resampling.sample.size(x = bs.both.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "BSBOTH_Resample_95_1k.csv")





##### WA ALL YRS #####

# set seed
#set.seed(1)
setwd("C:/Users/tclay/Documents/Methods paper/Data")
wa.04 <- read.csv("wa.2004.csv")
wa.04$Id <- as.factor(as.character(wa.04$Id))
wa.08 <- read.csv("wa.2008.csv")
wa.08$Id <- as.factor(as.character(wa.08$Id))
wa.09 <- read.csv("wa.2009.csv")
wa.09$Id <- as.factor(as.character(wa.09$Id))
levels(wa.04$Id)
#"5072082" "5077825" "5098231" "5098618" "5116345" "5116408" "5127247" "5127268" "5144775" "5144948" "5146578"
# "5147364" "5147491" "5202916" "5222657"
levels(wa.08$Id)
# "5077825" "5097085" "5098718" "5109020" "5109491" "5116188" "5116345" "5116408" "5116432" "5116480" "5121065"
# "5127232" "5127310" "5127446" "5127710" "5142556" "5144070" "5144241" "5146515" "5147364" "5148220" "5156003"
# "5156044" "5156109" "5156225" "5156335" "5156396" "5156460" "5164487" "5183024" "5183919" "5187824" "5202916"
# "5217268" "5222556" "5222647" "5222685
levels(wa.09$Id)
# "5054713" "5072450" "5077863" "5097256" "5098231" "5108659" "5109292" "5116011" "5116184" "5117985" "5127247"
# "5127610" "5127614" "5132339" "5144003" "5144582" "5145276" "5145808" "5146315" "5147871" "5147960" "5155793"
# "5155827" "5156372" "5164488" "5183481" "5184820" "5184821" "5184834" "5184976" "5202320" "5202427" "5216987"

wa4in8 <- subset(wa.04, Id %in% levels(wa.08$Id),)
wa4in8$Id <- as.factor(as.character(wa4in8$Id))
levels(wa4in8$Id)
# "5077825" "5116345" "5116408" "5147364" "5202916" - 5 individuals
wa4in9 <- subset(wa.04, Id %in% levels(wa.09$Id),)
wa4in9$Id <- as.factor(as.character(wa4in9$Id))
levels(wa4in9$Id)
# "5098231" "5127247"

# randomly remove these
wa.all <- rbind(wa.04, wa.08, wa.09)
wa.all$Id_yr <- as.factor(as.character(paste(wa.all$Id, wa.all$Year, sep = "_")))
levels(wa.all$Id_yr)
# summarize id years
wa.sum <- ddply(wa.all,.(Id, Year), summarize, Nlox = length(Year))
wa.sum2 <- ddply(wa.sum,.(Id), summarize, Nyrs = length(Id))
wa.sum3 <- subset(wa.sum2, Nyrs > 1,)
wa.repeat <- subset(wa.all, Id %in% wa.sum3$Id,)
wa.repeat$Id <- as.factor(as.character(wa.repeat$Id))
wa.repeat$Id_yr <- as.factor(as.character(wa.repeat$Id_yr))
wa.levs <- vector(length = nrow(wa.sum3))
for (i in 1:nlevels(wa.repeat$Id)) {
  wa.sub <- subset(wa.repeat, Id == levels(wa.repeat$Id)[i],)
  wa.levs[i] <- sample(as.character(unique(wa.sub$Id_yr)), 1, replace = FALSE)
}
# not repeat
wa.norepeat <- wa.all[!(wa.all$Id %in% wa.sum3$Id),]
wa.norepeat$Id <- as.factor(as.character(wa.norepeat$Id))
wa.repeat <- subset(wa.repeat, Id_yr %in% wa.levs,)

# merge together
wa.all2 <- rbind(wa.repeat, wa.norepeat)
wa.all2$Id <- as.factor(as.character(wa.all2$Id))
wa.all2$Id_yr <- as.factor(as.character(wa.all2$Id_yr))
nlevels(wa.all2$Id) # 78
levels(wa.all2$Id_yr)
# "5054713_2009" "5072082_2004" "5072450_2009" "5077825_2004" "5077863_2009" "5097085_2008" "5097256_2009"
#  "5098231_2009" "5098618_2004" "5098718_2008" "5108659_2009" "5109020_2008" "5109292_2009" "5109491_2008"
#  "5116011_2009" "5116184_2009" "5116188_2008" "5116345_2008" "5116408_2004" "5116432_2008" "5116480_2008"
#  "5117985_2009" "5121065_2008" "5127232_2008" "5127247_2009" "5127268_2004" "5127310_2008" "5127446_2008"
#  "5127610_2009" "5127614_2009" "5127710_2008" "5132339_2009" "5142556_2008" "5144003_2009" "5144070_2008"
#  "5144241_2008" "5144582_2009" "5144775_2004" "5144948_2004" "5145276_2009" "5145808_2009" "5146315_2009"
#  "5146515_2008" "5146578_2004" "5147364_2004" "5147491_2004" "5147871_2009" "5147960_2009" "5148220_2008"
#  "5155793_2009" "5155827_2009" "5156003_2008" "5156044_2008" "5156109_2008" "5156225_2008" "5156335_2008"
#  "5156372_2009" "5156396_2008" "5156460_2008" "5164487_2008" "5164488_2009" "5183024_2008" "5183481_2009"
#  "5183919_2008" "5184820_2009" "5184821_2009" "5184834_2009" "5184976_2009" "5187824_2008" "5202320_2009"
#  "5202427_2009" "5202916_2004" "5216987_2009" "5217268_2008" "5222556_2008" "5222647_2008" "5222657_2004"
#  "5222685_2008"

wa.all3 <- wa.all2
coordinates(wa.all3) <- ~Lon+Lat
proj4string(wa.all3) <- projection.dec
wa.all3 <- spTransform(wa.all3, CRS(projection))
wa.all4 <- data.frame(wa.all3)
wa.all4$Longitude <- wa.all4$Lon
wa.all4$Latitude <- wa.all4$Lat

wa.all.df <- as.data.frame(wa.all4)
wa.all.df$Id <- as.factor(as.character(wa.all.df$Id))

resample.out <- resampling.sample.size(x = wa.all4, niter = 1000, iso = 50, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "WAALL_Resample_50_1k.csv")

resample.out <- resampling.sample.size(x = wa.all4, niter = 1000, iso = 95, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "WAALL_Resample_95_1k.csv")
# 
# # need to record the individuals you used so you can use again!
# wa.all4$Id_yr <- as.factor(as.character(wa.all4$Id_yr))
# levels(wa.all4$Id_yr)
# # "5054713_2009" "5072082_2004" "5072450_2009" "5077825_2004" "5077863_2009" "5097085_2008" "5097256_2009"
# #  "5098231_2009" "5098618_2004" "5098718_2008" "5108659_2009" "5109020_2008" "5109292_2009" "5109491_2008"
# #  "5116011_2009" "5116184_2009" "5116188_2008" "5116345_2008" "5116408_2004" "5116432_2008" "5116480_2008"
# #  "5117985_2009" "5121065_2008" "5127232_2008" "5127247_2009" "5127268_2004" "5127310_2008" "5127446_2008"
# #  "5127610_2009" "5127614_2009" "5127710_2008" "5132339_2009" "5142556_2008" "5144003_2009" "5144070_2008"
# #  "5144241_2008" "5144582_2009" "5144775_2004" "5144948_2004" "5145276_2009" "5145808_2009" "5146315_2009"
# #  "5146515_2008" "5146578_2004" "5147364_2004" "5147491_2004" "5147871_2009" "5147960_2009" "5148220_2008"
# #  "5155793_2009" "5155827_2009" "5156003_2008" "5156044_2008" "5156109_2008" "5156225_2008" "5156335_2008"
# #  "5156372_2009" "5156396_2008" "5156460_2008" "5164487_2008" "5164488_2009" "5183024_2008" "5183481_2009"
# #  "5183919_2008" "5184820_2009" "5184821_2009" "5184834_2009" "5184976_2009" "5187824_2008" "5202320_2009"
# #  "5202427_2009" "5202916_2004" "5216987_2009" "5217268_2008" "5222556_2008" "5222647_2008" "5222657_2004"
# #  "5222685_2008"
# 

#####  WCP ALL #####


setwd("C:/Users/tclay/Documents/Methods paper/Data")
wcp.13 <- read.csv("wcp.2013.csv")
wcp.13$Year <- 2013
wcp.13$Id <- as.factor(as.character(wcp.13$Id))
wcp.13$Date.time <- as.POSIXct(wcp.13$Date.time, tz = "GMT")
wcp.13$Date <- as.Date(wcp.13$Date, tz = "GMT")
ddply(wcp.13,.(Id), summarize, Start = Date[1], End = Date[length(Date)])
# remove c358 as only lasts a month
wcp.13 <- subset(wcp.13, Id != "C538",)
wcp.13$Id <- as.factor(as.character(wcp.13$Id))
nlevels(wcp.13$Id)
# 9
wcp.03 <- read.csv("wcp.2003.csv")
wcp.03$Id <- as.factor(as.character(wcp.03$Id))
wcp.03$Date.time <- as.POSIXct(wcp.03$Date.time, tz = "GMT")
wcp.03$Date <- as.Date(wcp.03$Date, tz = "GMT")

levels(wcp.03$Id)
# "4178" "4216" "6016" "6751" "6756" "6758" "6759" "6766" "6769" "6773"
levels(wcp.13$Id)
# "C535" "C541" "C542" "C546" "C548" "C552" "C553" "C554" "C555"

wcp.03 <- wcp.03[,-9]
wcp.03 <- wcp.03[,-9]
wcp.both <- rbind(wcp.03, wcp.13)

wcp.both.2 <- wcp.both
coordinates(wcp.both.2) <- ~Lon+Lat
proj4string(wcp.both.2) <- projection.dec
wcp.both.2 <- spTransform(wcp.both.2, CRS(projection))
wcp.both.3 <- data.frame(wcp.both.2)
wcp.both.3$Longitude <- wcp.both.3$Lon
wcp.both.3$Latitude <- wcp.both.3$Lat

wcp.both.df <- as.data.frame(wcp.both.3)
wcp.both.df$Id <- as.factor(as.character(wcp.both.df$Id))

resample.out <- resampling.sample.size(x = wcp.both.3, niter = 1000, iso = 50, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "WCPBOTH_Resample_50_1k.csv")

resample.out <- resampling.sample.size(x = wcp.both.3, niter = 1000, iso = 95, cellsize =50000, h = 200000)
setwd("C:/Users/tclay/Documents/Methods paper/Boots")
write.csv(resample.out, "WCPBOTH_Resample_95_1k.csv")









########## CALCULATING SAMPLE SIZE PREDICTIONS - 95% COVER #################


# from the model
# c is c
# d is a
# e is b


library(drc)

# function to use asymptotic model to predict beyond range of data to 50 individuals, 100 individuals,
# and to the size of the population. 
pred.number.coverage <- function(model, pop.size, sample.size, percent.cover) {
  forecast <- c(50, 100, pop.size)
  pred.ci <- list()
  for (i in 1:length(forecast)) {
    pred.df <- data.frame(sample = seq(1, forecast[i], length.out = forecast[i]), groupings = 1) 
    pred.df2 <- data.frame(predict(model, newdata = pred.df, se.fit = FALSE, interval = "confidence", level = 0.95) )
    max.no.mean <- which.max(pred.df2$Prediction)
    max.no.upper <- which.max(pred.df2$Upper)
    max.no.lower <- which.max(pred.df2$Lower)
    max.value.mean <- pred.df2$Prediction[max.no.mean]
    max.value.upper <- pred.df2$Upper[max.no.upper]
    max.value.lower <- pred.df2$Lower[max.no.lower]
    value.mean <- max.value.mean*percent.cover
    value.upper <- max.value.upper*percent.cover
    value.lower <- max.value.lower*percent.cover
    # choosing next value above your threshold value
    sample.subset.mean <- subset(pred.df2, Prediction < value.mean)
    sample.subset.upper <- subset(pred.df2, Upper < value.upper)
    sample.subset.lower <- subset(pred.df2, Lower < value.lower)
    number.mean <- nrow(sample.subset.mean)+1
    number.upper <- nrow(sample.subset.upper)+1
    number.lower <- nrow(sample.subset.lower)+1
    pred.ci[[i]] <- data.frame(Prediction = number.mean, Upper = number.upper, Lower = number.lower,
                               Predict.To = forecast[i])
  }
  predict.df <- do.call(rbind, pred.ci)
  return(predict.df)
}

# model = asymptotic model run on resampled sample size area relationships (above).
# pop.size = the size of the population you want to predict to.
# sample.size = the original tracking sample size.
# percent.cover = the percentage cut-off of the predicted population area you think you need to represent. 



##### AP 2010 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("AP2010_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.044987e+08/1000000000 = 0.1044987
# b =   3.459769e-01 = 0.3459769
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -7.308799e+07/100000000 = -0.7308799
# a = 1.038420e+08/100000000 = 1.03842
# b = 1.770828e-01  = 0.1770828

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#           6     6     6         50
#           6     6     6        100
#           6     7     6       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#           6     6     6         50
#           7     7     6        100
#           7     7     7       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("AP2010_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 5.806558e+08/1000000000 = 0.5806558
# b = 6.486085e-01   = 0.6486085
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 7.232077e+07/1000000000 = 0.07232077
# a = 5.857572e+08/1000000000 = 0.5857572
# b = 8.187522e-01    = 0.8187522

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         11    12    10         50
# 2         12    13    12        100
# 3         14    15    13       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# #Prediction Upper Lower Predict.To
# 1         10    10    10         50
# 2         11    12    11        100
# 3         13    13    12       1000




##### AP 2011 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("AP2011_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.312365e+08/1000000000 = 0.1312365
# b =   1.024148
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -8.858028e+07/100000000 = -0.8858028
# a =  1.244097e+08/100000000 = 1.244097
# b = 4.052868e-01   = 0.4052868

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 6, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         11    12    10         50
# 2         12    14    11        100
# 3         14    15    12       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 6, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         14    15    14         50
# 2         17    17    16        100
# 3         20    20    19       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("AP2011_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 6.418790e+08/1000000000 = 0.641879
# b = 7.523508e-01 = 0.7523508
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 1.500464e+08/1000000000 =  0.1500464
# a =  6.698964e+08/1000000000 = 0.6698964
# b = 1.345309

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 6, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#         14    17    11         50
#         17    20    13        100
#         20    25    14       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 6, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          11    12    11         50
#          13    14    12        100
#          15    16    14       1000



##### BBA 2002 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("BBA2002_Resample_50_1k.csv")
#mm2 <- drm(Area ~ Sample, data = x1, fct = MM.2())
mm2<- drm(Area ~ Sample, data = x1, start = c(1, 1),fct = MM.2())
#mm2<- drm(Area ~ Sample, data = x1, start = c(1.425, 0.448),fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.448768e+08/1000000000 = 0.1448768
# b =  4.595700e-01 =  0.45957
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 4.387249e+07/100000000 = 0.4387249
# a = 1.457196e+08/100000000 = 1.457196
# b = 7.557581e-01   = 0.7557581

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 25, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#           9     9     8         50
#           9    10     9        100
#          10    10    10       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 25, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1          8     8     8         50
# 2          8     9     8        100
# 3          9     9     9       1000

# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("BBA2002_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 2.239177e+09/1000000000 = 2.239177
# b = 2.156264
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 3.826470e+08/1000000000 = 0.382647
# a = 2.315505e+09/1000000000 = 2.315505
# b = 3.402400

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 25, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          25    26    25         50
#          34    35    33        100
#          51    53    49       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 25, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          22    23    22         50
#          29    30    29        100
#          40    41    39       1000



##### BBA 2003 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("BBA2003_Resample_50_1k.csv")
#mm2 <- drm(Area ~ Sample, data = x1, fct = MM.2())
mm2<- drm(Area ~ Sample, data = x1, start = c(1, 1),fct = MM.2())
#mm2 <- drm(Area ~ Sample, data = x1, start = c(1.425, 1), fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.516382e+08/1000000000 = 0.1516382
# b =  3.834835e-01 = 0.3834835
# mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
mm3<- drm(Area ~ Sample, data = x1, start = c(2, 4, 2),fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -4.974825e+08 /1000000000 = -0.4974825
# a = 1.509715e+08/100000000 =  1.509715
# b = 7.404177e-02 = 0.07404177

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 24, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1          6     6     6         50
# 2          6     6     6        100
# 3          7     7     7       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 24, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1          7     7     7         50
# 2          7     7     7        100
# 3          8     8     8       1000

# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("BBA2003_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.896184e+09/1000000000 = 1.896184
# b = 1.577947
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 2.476874e+08/1000000000 = 0.2476874
# a = 1.925351e+09/1000000000 = 1.925351
# b = 2.120891

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 24, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          20    21    20         50
#          26    27    25        100
#          34    36    32       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 24, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          19    19    19         50
#          23    24    23        100
#          30    30    29       1000




##### BP 2011 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("BP2011_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
#mm2<- drm(Area ~ Sample, data = x1, start = c(0.1, 2),fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.703845e+08/1000000000 = 0.1703845
# b =   1.653245e-01  = 0.1653245
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
mm3<- drm(Area ~ Sample, data = x1, start = c(0.1, 2, 2),fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 1.441681e+08/100000000 =  1.441681
# a =  2.009868e+08/100000000 = 2.009868
# b =  9.704486

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         20    32     6         50
# 2         28    52    10        100
# 3         43   115   910       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#           3     4     3         50
#           4     4     3        100
#           4     4     3       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("BP2011_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.910972e+09/1000000000 = 1.910972
# b = 2.035725
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 1.241700e+08 /1000000000 =  0.12417
# a =   1.969770e+09/1000000000 = 1.96977
# b = 2.512087

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          23    25    22         50
#          31    33    29        100
#          43    47    39       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          22    22    21         50
#          28    29    27        100
#          38    39    36       1000



##### BS 2002 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("BS2002_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.404924e+08/1000000000 = 0.1404924
# b =   8.757407e-01 = 0.8757407
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -4.286115e+07/100000000 = -0.4286115
# a =  1.367040e+08/100000000 =  1.36704
# b =   5.375201e-01 = 0.5375201

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          11    12    10         50
#          12    13    11        100
#          14    15    13       100
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          13    13    13         50
#          15    15    14        100
#          17    17    16       1000

# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("BS2002_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 5.260252e+08/1000000000 = 0.5260252
# b = 7.797855e-01 = 0.7797855
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 2.575677e+07/1000000000 =  0.02575677
# a =    5.292744e+08/1000000000 = 0.5292744
# b = 8.612683e-01 = 0.8612683

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          12    13    11         50
#          14    15    13        100
#          16    17    14       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          12    12    12         50
#          13    14    13        100
#          15    15    15       1000



##### BS 2012 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("BS2012_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 2.426024e+08/1000000000 = 0.2426024
# b =    1.152997
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -1.677930e+08 /100000000 = -1.67793
# a =  2.370692e+08 /100000000 = 2.370692
# b =  5.002910e-01  =  0.500291

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          13    13    13         50
#          15    15    14        100
#          17    17    16       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          15    16    15         50
#          18    19    18        100
#          22    22    21       1000

# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("BS2012_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.298218e+09 /1000000000 = 1.298218
# b = 2.006139
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c =  9.045549e+07/1000000000 =  0.09045549
# a =    1.315640e+09/1000000000 = 1.31564
# b = 2.384635

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          23    23    22         50
#          30    30    29        100
#          41    42    39       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          22    22    21         50
#          28    28    27        100
#          37    38    37       1000
# 


##### GHA 1999 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("GHA1999_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a =  6.512144e+08/1000000000 =  0.6512144
# b =    1.278042
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -5.625786e+07 /100000000 = -0.5625786
# a =  6.476343e+08 /100000000 =  6.476343
# b = 1.099108

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          16    16    16         50
#          19    19    19        100
#          23    23    22       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          17    17    16         50
#          20    20    20        100
#          24    24    24       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("GHA1999_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 3.548325e+09 /1000000000 =  3.548325
# b = 1.432193
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c =  5.683460e+08/1000000000 =  0.568346
# a =    3.615449e+09/1000000000 = 3.615449
# b = 2.066562

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          20    20    19         50
#          25    25    24        100
#          32    33    31       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          18    18    18         50
#          22    22    21        100
#          27    27    27       1000





##### LMA 2004 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("LMA2004_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 3.685211e+08/1000000000 = 0.3685211
# b =   5.469613e-01 = 0.5469613
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -5.195824e+08/100000000 = -5.195824
# a = 3.656137e+08/100000000 =  3.656137
# b = 1.805479e-01 = 0.1805479

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#           8     8     8         50
#           8     8     8        100
#           9     9     9       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#           9     9     9         50
#          10    10    10        100
#          11    11    11       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("LMA2004_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a =  2.126798e+09/1000000000 = 2.126798
# b = 8.070287e-01 = 0.8070287
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 4.074771e+08/1000000000 = 0.4074771
# a = 2.167567e+09/1000000000 = 2.167567
# b = 1.215267

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          14    14    13         50
#          16    17    15        100
#          19    20    17       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#          12    12    12         50
#          14    14    13        100
#          16    16    15       1000





##### LMA 2005 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("LMA2005_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 2.724118e+08/1000000000 =  0.2724118
# b =   3.702634e-01 = 0.3702634
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -3.201521e+07/100000000 = -0.3201521
# a = 2.719334e+08/100000000 =  2.719334
# b = 3.187142e-01  = 0.3187142

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#           6     7     6         50
#           7     7     7        100
#           7     8     7       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
#           7     7     7         50
#           7     7     7        100
#           7     8     7       1000

# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("LMA2005_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a =  2.724118e+08/1000000000 = 0.2724118
# b = 3.702634e-01= 0.3702634
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -3.201521e+07/1000000000 = -0.03201521
# a = 2.719334e+08/1000000000 = 0.2719334
# b =  3.187142e-01 = 0.3187142

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1          6     7     6         50
# 2          7     7     7        100
# 3          7     8     7       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1          7     7     7         50
# 2          7     7     7        100
# 3          7     8     7       1000





##### NGP 2000 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("NGP2000_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.734073e+08/1000000000 =  0.1734073
# b =   4.893666e-01  =  0.4893666
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -3.071376e+07/100000000 = -0.3071376
# a = 1.730703e+08/100000000 =   1.730703
# b =  3.925292e-01  =  0.3925292

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1          8     8     8         50
# 2          9     9     9        100
# 3          9     9     9       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
Prediction Upper Lower Predict.To
# 1          8     8     8         50
# 2          9     9     9        100
# 3         10    10    10       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("NGP2000_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a =  8.504598e+08/1000000000 = 0.8504598
# b = 4.624149e-01  = 0.4624149
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 3.825515e+08/1000000000 = 0.3825515
# a = 8.618198e+08/1000000000 =  0.8618198
# b = 1.186985

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         10    10    10         50
# 2         11    11    11        100
# 3         12    13    12       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1          8     8     8         50
# 2          9     9     9        100
# 3          9     9     9       1000




##### SGP 2000 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("SGP2000_Resample_50_1k.csv")
#mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
mm2<- drm(Area ~ Sample, data = x1, start = c(1, 1),fct = MM.2())
# mm2<- drm(Area ~ Sample, data = x1, start = c(1.425, 0.448),fct = MM.2())
plot(mm2) 
mm2$coefficients
# a =  1.435301e+08/1000000000 = 0.1435301
# b = 1.996108e-01 = 0.1996108
#mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
mm3<- drm(Area ~ Sample, data = x1, start = c(2, 4, 2),fct = MM.3())
# mm3<- drm(Area ~ Sample, data = x1, start = c(2, 2, 1),fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -1.076698e+08/1000000000 = -0.1076698
# a = 1.434225e+08/1000000000 = 0.1434225
# b = 1.071805e-01 = 0.1071805

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1          4     4     4         50
# 2          4     4     4        100
# 3          4     4     4       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1          4     4     4         50
# 2          4     4     4        100
# 3          4     4     4       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("SGP2000_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a =  1.473385e+09/1000000000 = 1.473385
# b =  1.081275
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 3.169101e+08/1000000000 = 0.3169101
# a = 1.493825e+09/1000000000 =  1.493825
# b = 1.695110

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         17    17    16         50
# 2         20    21    20        100
# 3         25    26    24       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         15    15    15         50
# 2         17    18    17        100
# 3         21    21    20       1000






##### WA 2004 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("WA2004_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 5.489386e+09/1000000000 = 5.489386
# b =  1.958169
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 2.533478e+08/100000000 = 0.2533478
# a = 5.545731e+09/100000000 = 5.545731
# b = 2.20015

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         22    22    22         50
# 2         29    29    28        100
# 3         39    40    38       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         21    22    21         50
# 2         27    27    27        100
# 3         36    37    36       1000
# # In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("WA2004_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 5.458988e+09/1000000000 = 5.458988
# b = 1.883517e+00  = 1.883517
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 1.342372e+08/1000000000 = 0.1342372
# a = 5.485105e+09/1000000000 = 5.485105
# b = 1.997712e+00   = 1.997712

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         21    22    21         50
# 2         27    28    27        100
# 3         36    37    35       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         21    21    21         50
# 2         26    27    26        100
# 3         35    35    35       1000




##### WA 2008 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("WA2008_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 6.994537e+08/1000000000 = 0.6994537
# b =  1.575760
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -7.408580e+07/100000000 = -0.740858
# a = 6.958049e+08/100000000 = 6.958049
# b =  1.318036

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
wa.pred.mm3 <- pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
wa.pred.mm3
# Prediction Upper Lower Predict.To
# 1         18    18    18         50
# 2         22    22    22        100
# 3         28    28    27       1000
wa.pred.mm2 <- pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
wa.pred.mm2
# Prediction Upper Lower Predict.To
# 1         19    19    19         50
# 2         23    23    23        100
# 3         30    30    29       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("WA2008_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 5.568490e+09/1000000000 = 5.56849
# b = 2.230625
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 2.633474e+08/1000000000 = 0.2633474
# a = 5.594948e+09/1000000000 = 5.594948
# b = 2.466244

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         23    24    23         50
# 2         31    31    31        100
# 3         43    43    43       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         23    23    23         50
# 2         30    30    30        100
# 3         41    41    41       1000


##### WA 2009 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("WA2009_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 5.558802e+08/1000000000 =  0.5558802
# b =  1.066079
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -7.041043e+07/100000000 = -0.7041043
# a = 5.537291e+08/100000000 = 5.537291
# b =  8.798918e-01 = 0.8798918

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         14    14    14         50
# 2         16    17    16        100
# 3         19    20    19       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         15    15    15         50
# 2         17    17    17        100
# 3         20    21    20       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("WA2009_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 4.955721e+09/1000000000 = 4.955721
# b = 2.523067
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 7.399283e+08/1000000000 =  0.7399283
# a = 5.081190e+09/1000000000 = 5.08119
# b = 3.691168

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         27    27    26         50
# 2         37    37    36        100
# 3         56    58    55       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         24    24    24         50
# 2         32    33    32        100
# 3         46    47    46       1000





#### WCP 2003 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("WCP2003_Resample_50_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.889647e+08/1000000000 =  0.1889647
# b = 1.384131
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -8.165973e+07/1000000000 =  -0.08165973
# a = 1.810094e+08/1000000000 = 0.1810094
# b = 6.876653e-01 =  0.6876653


# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         14    15    14         50
# 2         17    17    16        100
# 3         19    21    18       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         17    18    17         50
# 2         21    22    21        100
# 3         26    27    25       1000
# # In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("WCP2003_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 8.178225e+08/1000000000 = 0.8178225
# b =  1.279069
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 7.446207e+07/1000000000 =  0.07446207
# a = 8.327305e+08/1000000000 = 0.8327305
# b = 1.597956

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         18    18    17         50
# 2         22    23    21        100
# 3         27    29    26       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         17    17    16         50
# 2         20    20    20        100
# 3         24    25    24       1000




#### WCP 2013 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("WCP2013_Resample_50_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2)
mm2$coefficients
# a = 1.474499e+08 /1000000000 =  0.1474499
# b = 9.888987e-01 = 0.9888987
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -1.495247e+07/1000000000 =  -0.01495247
# a = 1.459660e+08/1000000000 = 0.145966
# b = 8.216576e-01 = 0.8216576

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         13    14    12         50
# 2         15    16    14        100
# 3         17    19    16       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         14    14    14         50
# 2         16    17    16        100
# 3         19    20    18       1000

# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("WCP2013_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 6.568766e+08 /1000000000 = 0.6568766
# b = 8.291787e-01  = 0.8291787
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 5.631578e+07/1000000000 =  0.05631578
# a = 6.626760e+08/1000000000 =  0.662676
# b = 9.849061e-01  = 0.9849061

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         13    14    12         50
# 2         15    16    14        100
# 3         17    18    16       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         12    13    12         50
# 2         14    14    14        100
# 3         16    16    16       1000
# 





##### AP ALL #####



# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("APBOTH_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.209011e+08/1000000000 = 0.1209011
# b =   6.388617e-01 = 0.6388617
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -1.150595e+08/100000000 = -1.150595
# a = 1.191758e+08/100000000 = 1.191758
# b = 2.543918e-01   = 0.2543918

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 14, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1          9     9     8         50
# 2          9    10     9        100
# 3         10    11    10       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 14, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         10    10    10         50
# 2         11    11    11        100
# 3         12    13    12       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("APBOTH_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 6.642996e+08/1000000000 = 0.6642996
# b = 1.030524
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 2.258254e+08/1000000000 = 0.2258254
# a = 7.079439e+08/1000000000 = 0.7079439
# b = 2.694532

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 14, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         20    21    19         50
# 2         25    27    23        100
# 3         33    36    30       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 14, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         14    15    14         50
# 2         17    17    16        100
# 3         20    20    19       1000






##### BS ALL #####





# 50% KERNEL #


setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("BSBOTH_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
#mm2<- drm(Area ~ Sample, data = x1, start = c(1, 1),fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 2.251949e+08 /1000000000 = 0.2251949
# b =   1.184961
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -1.204011e+08/100000000 = -1.204011
# a = 2.217293e+08/100000000 =2.217293
# b = 6.051414e-01  = 0.6051414

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 25, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         14    14    14         50
# 2         16    16    16        100
# 3         18    19    18       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 25, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         16    16    16         50
# 2         19    19    18        100
# 3         22    23    22       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("BSBOTH_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
#mm2<- drm(Area ~ Sample, data = x2, start = c(1, 1),fct = MM.2())
plot(mm2) 
mm2$coefficients
# a =  1.249974e+09/1000000000 = 1.249974
# b = 2.339168
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 1.377153e+08/1000000000 = 0.1377153
# a = 1.276833e+09/1000000000 = 1.276833
# b = 3.105543

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 25, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         25    26    25         50
# 2         34    35    33        100
# 3         50    52    48       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 25, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         23    24    23         50
# 2         31    31    30        100
# 3         43    44    42       1000






#### WA ALL ####




# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("WAALL_Resample_50_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 6.632030e+08/1000000000 =  0.663203
# b = 1.531291
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -7.491080e+07/1000000000 = -0.0749108
# a = 6.614442e+08/1000000000 = 0.6614442
# b =  1.287230


# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 78, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         18    18    18         50
# 2         22    22    22        100
# 3         27    27    27       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 78, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         19    19    19         50
# 2         23    23    23        100
# 3         29    29    29       1000
# # In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("WAALL_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 5.425202e+09 /1000000000 = 5.425202
# b =  2.498601
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 4.999381e+08/1000000000 =   0.4999381
# a = 5.454190e+09/1000000000 = 5.45419
# b = 3.004000

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 78, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         25    25    25         50
# 2         34    34    34        100
# 3         49    50    49       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 78, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         24    24    24         50
# 2         32    32    32        100
# 3         46    46    46       1000






####  WCP ALL ######



# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("WCPBOTH_Resample_50_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 2.035840e+08/1000000000 =  0.203584
# b = 1.676345
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -1.305616e+07/1000000000 = -0.01305616
# a = 2.021146e+08/1000000000 = 0.2021146
# b = 1.473966


# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 19, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         19    19    19         50
# 2         23    24    23        100
# 3         30    30    29       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 19, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         20    20    19         50
# 2         24    25    24        100
# 3         31    32    31       1000
# # In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("WCPBOTH_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 8.146544e+08/1000000000 = 0.8146544
# b =  1.209986
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -5.423566e+07/1000000000 =  -0.05423566
# a = 8.108681e+08/1000000000 = 0.8108681
# b = 1.074571

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 19, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         15    16    15         50
# 2         18    19    18        100
# 3         22    22    21       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 19, percent.cover = 0.95)
# Prediction Upper Lower Predict.To
# 1         16    16    16         50
# 2         19    19    19        100
# 3         23    23    23       1000



########## CALCULATING SAMPLE SIZE PREDICTIONS - 70% COVER #################


# from the model
# c is c
# d is a
# e is b


library(drc)

# function to use asymptotic model to predict beyond range of data to 50 individuals, 100 individuals,
# and to the size of the population. 
pred.number.coverage <- function(model, pop.size, sample.size, percent.cover) {
  forecast <- c(50, 100, pop.size)
  pred.ci <- list()
  for (i in 1:length(forecast)) {
    pred.df <- data.frame(sample = seq(1, forecast[i], length.out = forecast[i]), groupings = 1) 
    pred.df2 <- data.frame(predict(model, newdata = pred.df, se.fit = FALSE, interval = "confidence", level = 0.95) )
    max.no.mean <- which.max(pred.df2$Prediction)
    max.no.upper <- which.max(pred.df2$Upper)
    max.no.lower <- which.max(pred.df2$Lower)
    max.value.mean <- pred.df2$Prediction[max.no.mean]
    max.value.upper <- pred.df2$Upper[max.no.upper]
    max.value.lower <- pred.df2$Lower[max.no.lower]
    value.mean <- max.value.mean*percent.cover
    value.upper <- max.value.upper*percent.cover
    value.lower <- max.value.lower*percent.cover
    # choosing next value above your threshold value
    sample.subset.mean <- subset(pred.df2, Prediction < value.mean)
    sample.subset.upper <- subset(pred.df2, Upper < value.upper)
    sample.subset.lower <- subset(pred.df2, Lower < value.lower)
    number.mean <- nrow(sample.subset.mean)+1
    number.upper <- nrow(sample.subset.upper)+1
    number.lower <- nrow(sample.subset.lower)+1
    pred.ci[[i]] <- data.frame(Prediction = number.mean, Upper = number.upper, Lower = number.lower,
                               Predict.To = forecast[i])
  }
  predict.df <- do.call(rbind, pred.ci)
  return(predict.df)
}

# model = asymptotic model run on resampled sample size area relationships (above).
# pop.size = the size of the population you want to predict to.
# sample.size = the original tracking sample size.
# percent.cover = the percentage cut-off of the predicted population area you think you need to represent. 



##### AP 2010 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("AP2010_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.044987e+08/1000000000 = 0.1044987
# b =   3.459769e-01 = 0.3459769
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -7.308799e+07/100000000 = -0.7308799
# a = 1.038420e+08/100000000 = 1.03842
# b = 1.770828e-01  = 0.1770828

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          1     1     1         50
# 2          1     1     1        100
# 3          1     1     1       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          1     1     1         50
# 2          1     1     1        100
# 3          1     1     1       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("AP2010_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 5.806558e+08/1000000000 = 0.5806558
# b = 6.486085e-01   = 0.6486085
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 7.232077e+07/1000000000 = 0.07232077
# a = 5.857572e+08/1000000000 = 0.5857572
# b = 8.187522e-01    = 0.8187522

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     2     2         50
# 2          2     2     2        100
# 3          2     2     2       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     2     2         50
# 2          2     2     2        100
# 3          2     2     2       1000




##### AP 2011 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("AP2011_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.312365e+08/1000000000 = 0.1312365
# b =   1.024148
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -8.858028e+07/100000000 = -0.8858028
# a =  1.244097e+08/100000000 = 1.244097
# b = 4.052868e-01   = 0.4052868

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 6, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     2     2         50
# 2          2     2     2        100
# 3          2     2     2       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 6, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          3     3     3         50
# 2          3     3     3        100
# 3          3     3     3       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("AP2011_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 6.418790e+08/1000000000 = 0.641879
# b = 7.523508e-01 = 0.7523508
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 1.500464e+08/1000000000 =  0.1500464
# a =  6.698964e+08/1000000000 = 0.6698964
# b = 1.345309

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 6, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     3     2         50
# 2          3     3     2        100
# 3          3     3     2       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 6, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     2     2         50
# 2          2     2     2        100
# 3          2     2     2       1000



##### BBA 2002 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("BBA2002_Resample_50_1k.csv")
#mm2 <- drm(Area ~ Sample, data = x1, fct = MM.2())
mm2<- drm(Area ~ Sample, data = x1, start = c(1, 1),fct = MM.2())
#mm2<- drm(Area ~ Sample, data = x1, start = c(1.425, 0.448),fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.448768e+08/1000000000 = 0.1448768
# b =  4.595700e-01 =  0.45957
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 4.387249e+07/100000000 = 0.4387249
# a = 1.457196e+08/100000000 = 1.457196
# b = 7.557581e-01   = 0.7557581

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 25, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          1     1     2         50
# 2          1     1     2        100
# 3          2     1     2       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 25, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     2     2         50
# 2          2     2     2        100
# 3          2     2     2       1000

# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("BBA2002_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 2.239177e+09/1000000000 = 2.239177
# b = 2.156264
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 3.826470e+08/1000000000 = 0.382647
# a = 2.315505e+09/1000000000 = 2.315505
# b = 3.402400

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 25, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          6     6     6         50
# 2          6     6     6        100
# 3          7     7     6       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 25, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          5     5     5         50
# 2          5     5     5        100
# 3          5     5     5       1000



##### BBA 2003 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("BBA2003_Resample_50_1k.csv")
#mm2 <- drm(Area ~ Sample, data = x1, fct = MM.2())
mm2<- drm(Area ~ Sample, data = x1, start = c(1, 1),fct = MM.2())
#mm2 <- drm(Area ~ Sample, data = x1, start = c(1.425, 1), fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.516382e+08/1000000000 = 0.1516382
# b =  3.834835e-01 = 0.3834835
# mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
mm3<- drm(Area ~ Sample, data = x1, start = c(2, 4, 2),fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -4.974825e+08 /1000000000 = -0.4974825
# a = 1.509715e+08/100000000 =  1.509715
# b = 7.404177e-02 = 0.07404177

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 24, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          1     1     1         50
# 2          1     1     2        100
# 3          1     1     2       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 24, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          1     1     1         50
# 2          1     1     1        100
# 3          1     1     1       1000

# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("BBA2003_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.896184e+09/1000000000 = 1.896184
# b = 1.577947
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 2.476874e+08/1000000000 = 0.2476874
# a = 1.925351e+09/1000000000 = 1.925351
# b = 2.120891

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 24, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          4     4     4         50
# 2          4     4     4        100
# 3          5     5     5       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 24, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          4     4     4         50
# 2          4     4     4        100
# 3          4     4     4       1000




##### BP 2011 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("BP2011_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.703845e+08/1000000000 = 0.1703845
# b =   1.653245e-01  = 0.1653245
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
mm3<- drm(Area ~ Sample, data = x1, start = c(0.1, 2, 2),fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 1.441681e+08/100000000 =  1.441681
# a =  2.009868e+08/100000000 = 2.009868
# b =  9.704486

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          1     1     1         50
# 2          1     4     1        100
# 3          1     8     1       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          1     1     1         50
# 2          1     1     1        100
# 3          1     1     1       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("BP2011_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.910972e+09/1000000000 = 1.910972
# b = 2.035725
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 1.241700e+08 /1000000000 =  0.12417
# a =   1.969770e+09/1000000000 = 1.96977
# b = 2.512087

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          5     5     5         50
# 2          5     6     5        100
# 3          6     6     5       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          5     5     5         50
# 2          5     5     5        100
# 3          5     5     5       1000



##### BS 2002 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("BS2002_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.404924e+08/1000000000 = 0.1404924
# b =   8.757407e-01 = 0.8757407
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -4.286115e+07/100000000 = -0.4286115
# a =  1.367040e+08/100000000 =  1.36704
# b =   5.375201e-01 = 0.5375201

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     2     2         50
# 2          2     2     2        100
# 3          2     2     2       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     2     2         50
# 2          2     3     2        100
# 3          3     3     3       1000

# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("BS2002_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 5.260252e+08/1000000000 = 0.5260252
# b = 7.797855e-01 = 0.7797855
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 2.575677e+07/1000000000 =  0.02575677
# a =    5.292744e+08/1000000000 = 0.5292744
# b = 8.612683e-01 = 0.8612683

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     2     2         50
# 2          2     2     2        100
# 3          2     2     2       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     2     2         50
# 2          2     2     2        100
# 3          2     2     2       1000



##### BS 2012 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("BS2012_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 2.426024e+08/1000000000 = 0.2426024
# b =    1.152997
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -1.677930e+08 /100000000 = -1.67793
# a =  2.370692e+08 /100000000 = 2.370692
# b =  5.002910e-01  =  0.500291

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          3     3     3         50
# 2          3     3     3        100
# 3          3     3     3       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          3     3     3         50
# 2          3     3     3        100
# 3          3     3     3       1000

# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("BS2012_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.298218e+09 /1000000000 = 1.298218
# b = 2.006139
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c =  9.045549e+07/1000000000 =  0.09045549
# a =    1.315640e+09/1000000000 = 1.31564
# b = 2.384635

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          5     5     5         50
# 2          5     5     5        100
# 3          5     6     5       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          5     5     5         50
# 2          5     5     5        100
# 3          5     5     5       1000
# 


##### GHA 1999 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("GHA1999_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a =  6.512144e+08/1000000000 =  0.6512144
# b =    1.278042
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -5.625786e+07 /100000000 = -0.5625786
# a =  6.476343e+08 /100000000 =  6.476343
# b = 1.099108

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          3     3     3         50
# 2          3     3     3        100
# 3          3     3     3       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          3     3     3         50
# 2          3     3     3        100
# 3          3     3     3       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("GHA1999_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 3.548325e+09 /1000000000 =  3.548325
# b = 1.432193
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c =  5.683460e+08/1000000000 =  0.568346
# a =    3.615449e+09/1000000000 = 3.615449
# b = 2.066562

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          4     4     4         50
# 2          4     4     4        100
# 3          4     4     4       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          4     4     4         50
# 2          4     4     4        100
# 3          4     4     4       1000





##### LMA 2004 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("LMA2004_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 3.685211e+08/1000000000 = 0.3685211
# b =   5.469613e-01 = 0.5469613
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -5.195824e+08/100000000 = -5.195824
# a = 3.656137e+08/100000000 =  3.656137
# b = 1.805479e-01 = 0.1805479

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     2     2         50
# 2          2     2     2        100
# 3          2     2     2       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     2     2         50
# 2          2     2     2        100
# 3          2     2     2       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("LMA2004_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a =  2.126798e+09/1000000000 = 2.126798
# b = 8.070287e-01 = 0.8070287
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 4.074771e+08/1000000000 = 0.4074771
# a = 2.167567e+09/1000000000 = 2.167567
# b = 1.215267

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     2     2         50
# 2          3     3     2        100
# 3          3     3     3       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     2     2         50
# 2          2     2     2        100
# 3          2     2     2       1000





##### LMA 2005 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("LMA2005_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 2.724118e+08/1000000000 =  0.2724118
# b =   3.702634e-01 = 0.3702634
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -3.201521e+07/100000000 = -0.3201521
# a = 2.719334e+08/100000000 =  2.719334
# b = 3.187142e-01  = 0.3187142

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          1     1     1         50
# 2          1     1     1        100
# 3          1     1     1       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          1     1     1         50
# 2          1     1     1        100
# 3          1     1     1       1000

# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("LMA2005_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a =  2.724118e+08/1000000000 = 0.2724118
# b = 3.702634e-01= 0.3702634
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -3.201521e+07/1000000000 = -0.03201521
# a = 2.719334e+08/1000000000 = 0.2719334
# b =  3.187142e-01 = 0.3187142

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          1     1     1         50
# 2          1     1     1        100
# 3          1     1     1       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          1     1     1         50
# 2          1     1     1        100
# 3          1     1     1       1000





##### NGP 2000 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("NGP2000_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.734073e+08/1000000000 =  0.1734073
# b =   4.893666e-01  =  0.4893666
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -3.071376e+07/100000000 = -0.3071376
# a = 1.730703e+08/100000000 =   1.730703
# b =  3.925292e-01  =  0.3925292

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     2     2         50
# 2          2     2     2        100
# 3          2     2     2       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     2     2         50
# 2          2     2     2        100
# 3          2     2     2       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("NGP2000_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a =  8.504598e+08/1000000000 = 0.8504598
# b = 4.624149e-01  = 0.4624149
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 3.825515e+08/1000000000 = 0.3825515
# a = 8.618198e+08/1000000000 =  0.8618198
# b = 1.186985

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          1     1     1         50
# 2          1     1     1        100
# 3          2     1     2       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          2     2     2         50
# 2          2     2     2        100
# 3          2     2     2       1000




##### SGP 2000 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("SGP2000_Resample_50_1k.csv")
#mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
mm2<- drm(Area ~ Sample, data = x1, start = c(1, 1),fct = MM.2())
# mm2<- drm(Area ~ Sample, data = x1, start = c(1.425, 0.448),fct = MM.2())
plot(mm2) 
mm2$coefficients
# a =  1.435301e+08/1000000000 = 0.1435301
# b = 1.996108e-01 = 0.1996108
#mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
mm3<- drm(Area ~ Sample, data = x1, start = c(2, 4, 2),fct = MM.3())
# mm3<- drm(Area ~ Sample, data = x1, start = c(2, 2, 1),fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -1.076698e+08/1000000000 = -0.1076698
# a = 1.434225e+08/1000000000 = 0.1434225
# b = 1.071805e-01 = 0.1071805

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          1     1     1         50
# 2          1     1     1        100
# 3          1     1     1       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          1     1     1         50
# 2          1     1     1        100
# 3          1     1     1       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("SGP2000_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a =  1.473385e+09/1000000000 = 1.473385
# b =  1.081275
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 3.169101e+08/1000000000 = 0.3169101
# a = 1.493825e+09/1000000000 =  1.493825
# b = 1.695110

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          3     3     3         50
# 2          3     3     3        100
# 3          3     3     3       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          3     3     3         50
# 2          3     3     3        100
# 3          3     3     3       1000






##### WA 2004 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("WA2004_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 5.489386e+09/1000000000 = 5.489386
# b =  1.958169
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 2.533478e+08/100000000 = 0.2533478
# a = 5.545731e+09/100000000 = 5.545731
# b = 2.20015

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          5     5     5         50
# 2          5     5     5        100
# 3          5     5     5       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          5     5     5         50
#2          5     5     5        100
#3          5     5     5       1000
# # In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("WA2004_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 5.458988e+09/1000000000 = 5.458988
# b = 1.883517e+00  = 1.883517
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 1.342372e+08/1000000000 = 0.1342372
# a = 5.485105e+09/1000000000 = 5.485105
# b = 1.997712e+00   = 1.997712

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          4     4     4         50
#2          5     5     5        100
#3          5     5     5       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          4     4     4         50
#2          5     5     5        100
#3          5     5     5       1000




##### WA 2008 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("WA2008_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 6.994537e+08/1000000000 = 0.6994537
# b =  1.575760
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -7.408580e+07/100000000 = -0.740858
# a = 6.958049e+08/100000000 = 6.958049
# b =  1.318036

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
wa.pred.mm3 <- pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
wa.pred.mm3
# Prediction Upper Lower Predict.To
#1          4     4     4         50
#2          4     4     4        100
#3          4     4     4       1000
wa.pred.mm2 <- pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
wa.pred.mm2
# Prediction Upper Lower Predict.To
#1          4     4     4         50
#2          4     4     4        100
#3          4     4     4       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #
# Prediction Upper Lower Predict.To
x2 <- read.csv("WA2008_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 5.568490e+09/1000000000 = 5.56849
# b = 2.230625
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 2.633474e+08/1000000000 = 0.2633474
# a = 5.594948e+09/1000000000 = 5.594948
# b = 2.466244

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
# 1          5     5     5         50
# 2          5     5     5        100
# 3          6     6     6       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
##1          5     5     5         50
#2          5     5     5        100
#3          6     6     6       1000


##### WA 2009 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("WA2009_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 5.558802e+08/1000000000 =  0.5558802
# b =  1.066079
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -7.041043e+07/100000000 = -0.7041043
# a = 5.537291e+08/100000000 = 5.537291
# b =  8.798918e-01 = 0.8798918

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          3     3     3         50
#2          3     3     3        100
#3          3     3     3       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          3     3     3         50
#2          3     3     3        100
#3          3     3     3       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("WA2009_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 4.955721e+09/1000000000 = 4.955721
# b = 2.523067
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 7.399283e+08/1000000000 =  0.7399283
# a = 5.081190e+09/1000000000 = 5.08119
# b = 3.691168

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          6     6     6         50
#2          7     7     7        100
#3          7     7     7       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          6     6     6         50
#2          6     6     6        100
#3          6     6     6       1000





#### WCP 2003 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("WCP2003_Resample_50_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.889647e+08/1000000000 =  0.1889647
# b = 1.384131
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -8.165973e+07/1000000000 =  -0.08165973
# a = 1.810094e+08/1000000000 = 0.1810094
# b = 6.876653e-01 =  0.6876653


# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          3     3     3         50
#2          3     3     3        100
#3          3     3     3       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          3     3     3         50
#2          4     4     4        100
#3          4     4     4       1000
# # In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("WCP2003_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 8.178225e+08/1000000000 = 0.8178225
# b =  1.279069
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 7.446207e+07/1000000000 =  0.07446207
# a = 8.327305e+08/1000000000 = 0.8327305
# b = 1.597956

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          4     4     3         50
#2          4     4     4        100
#3          4     4     4       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          3     3     3         50
#2          3     3     3        100
#3          3     4     3       1000




#### WCP 2013 ####


# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("WCP2013_Resample_50_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2)
mm2$coefficients
# a = 1.474499e+08 /1000000000 =  0.1474499
# b = 9.888987e-01 = 0.9888987
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -1.495247e+07/1000000000 =  -0.01495247
# a = 1.459660e+08/1000000000 = 0.145966
# b = 8.216576e-01 = 0.8216576

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          3     3     3         50
#2          3     3     3        100
#3          3     3     3       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          3     3     3         50
#2          3     3     3        100
#3          3     3     3       1000

# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("WCP2013_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 6.568766e+08 /1000000000 = 0.6568766
# b = 8.291787e-01  = 0.8291787
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 5.631578e+07/1000000000 =  0.05631578
# a = 6.626760e+08/1000000000 =  0.662676
# b = 9.849061e-01  = 0.9849061

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          2     2     2         50
#2          2     2     2        100
#3          3     3     2       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 9, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          2     2     2         50
#2          2     2     2        100
#3          2     2     2       1000
# 





##### AP ALL #####



# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("APBOTH_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 1.209011e+08/1000000000 = 0.1209011
# b =   6.388617e-01 = 0.6388617
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -1.150595e+08/100000000 = -1.150595
# a = 1.191758e+08/100000000 = 1.191758
# b = 2.543918e-01   = 0.2543918

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 14, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          2     2     2         50
#2          2     2     2        100
#3          2     2     2       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 14, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          2     2     2         50
#2          2     2     2        100
#3          2     2     2       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("APBOTH_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 6.642996e+08/1000000000 = 0.6642996
# b = 1.030524
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 2.258254e+08/1000000000 = 0.2258254
# a = 7.079439e+08/1000000000 = 0.7079439
# b = 2.694532

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 14, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          3     4     3         50
#2          4     4     4        100
#3          4     4     4       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 14, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          3     3     3         50
#2          3     3     3        100
#3          3     3     3       1000






##### BS ALL #####





# 50% KERNEL #


setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("BSBOTH_Resample_50_1k.csv")
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
#mm2<- drm(Area ~ Sample, data = x1, start = c(1, 1),fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 2.251949e+08 /1000000000 = 0.2251949
# b =   1.184961
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -1.204011e+08/100000000 = -1.204011
# a = 2.217293e+08/100000000 =2.217293
# b = 6.051414e-01  = 0.6051414

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 25, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          3     3     3         50
#2          3     3     3        100
#3          3     3     3       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 25, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          3     3     3         50
#2          3     3     3        100
#3          3     3     3       1000
# In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("BSBOTH_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
#mm2<- drm(Area ~ Sample, data = x2, start = c(1, 1),fct = MM.2())
plot(mm2) 
mm2$coefficients
# a =  1.249974e+09/1000000000 = 1.249974
# b = 2.339168
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 1.377153e+08/1000000000 = 0.1377153
# a = 1.276833e+09/1000000000 = 1.276833
# b = 3.105543

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 25, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          6     6     6         50
#2          6     6     6        100
#3          7     7     7       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 25, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          5     5     5         50
#2          6     6     6        100
#3          6     6     6       1000






#### WA ALL ####




# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("WAALL_Resample_50_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 6.632030e+08/1000000000 =  0.663203
# b = 1.531291
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -7.491080e+07/1000000000 = -0.0749108
# a = 6.614442e+08/1000000000 = 0.6614442
# b =  1.287230


# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 78, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          4     4     4         50
#2          4     4     4        100
#3          4     4     4       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 78, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          4     4     4         50
#2          4     4     4        100
#3          4     4     4       1000
# # In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("WAALL_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 5.425202e+09 /1000000000 = 5.425202
# b =  2.498601
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = 4.999381e+08/1000000000 =   0.4999381
# a = 5.454190e+09/1000000000 = 5.45419
# b = 3.004000

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 78, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          6     6     6         50
#2          6     6     6        100
#3          7     7     7       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 78, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          5     5     6         50
#2          6     6     6        100
#3          6     6     6       1000






####  WCP ALL ######



# 50% KERNEL #

setwd("C:/Users/tclay/Documents/Methods paper/Boots") 
x1 <- read.csv("WCPBOTH_Resample_50_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x1, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 2.035840e+08/1000000000 =  0.203584
# b = 1.676345
mm3<- drm(Area ~ Sample, data = x1, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -1.305616e+07/1000000000 = -0.01305616
# a = 2.021146e+08/1000000000 = 0.2021146
# b = 1.473966


# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 19, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          4     4     4         50
#2          4     4     4        100
#3          4     4     4       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 19, percent.cover = 0.7)
# Prediction Upper Lower Predict.To
#1          4     4     4         50
#2          4     4     4        100
#3          4     4     4       1000
# # In this example, prediction made for Michaelis-Menten 3-parameter model based on a sample size of 9 individuals, 
# predicting to 95% of the total area occupied by a projected "population" size of 50, 100 and 1000 individuals.


# 95% KERNEL #

x2 <- read.csv("WCPBOTH_Resample_95_1k.csv")
#x <- subset(x, Isopleth == 50,)
# x$Area <- x$Area/1000000000000
mm2<- drm(Area ~ Sample, data = x2, fct = MM.2())
plot(mm2) 
mm2$coefficients
# a = 8.146544e+08/1000000000 = 0.8146544
# b =  1.209986
mm3<- drm(Area ~ Sample, data = x2, fct = MM.3())
plot(mm3) 
mm3$coefficients
# c = -5.423566e+07/1000000000 =  -0.05423566
# a = 8.108681e+08/1000000000 = 0.8108681
# b = 1.074571

# to see if model has converged or not. If it hasn't then you need to give the model starting values, 
# e.g. mm3<- drm(Area ~ Sample, data = x, start = c(1, 2, 1),fct = MM.3())
pred.number.coverage(model = mm3, pop.size = 1000, sample.size = 19, percent.cover = 0.7)
#  Prediction Upper Lower Predict.To
#1          3     3     3         50
#2          3     3     3        100
#3          3     3     3       1000
pred.number.coverage(model = mm2, pop.size = 1000, sample.size = 19, percent.cover = 0.7)
#  Prediction Upper Lower Predict.To
#1          3     3     3         50
#2          3     3     3        100
#3          3     3     3       1000
