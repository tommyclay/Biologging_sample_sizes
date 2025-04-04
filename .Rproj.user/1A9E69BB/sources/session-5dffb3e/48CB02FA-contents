library(adehabitatHR)
library(raster)
library(sf)
library(sfheaders)
library(ggplot2)
library(graticule)

#### 1. LOAD FUNCTIONS #####

source("./Codes/1_predicting_minimum_sample_sizes/functions.R")

#### 2. CALCULATING INDIVIDUAL AND POPULATION KERNELS FOR EACH DATASET ####

# setting projections
proj_utm <- paste0("+proj=laea +lat_0=", -90, " +lon_0=", -38.05, " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj_dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"   

# setting up directory
in_direx <- "./Data_inputs/gls_data/"
f <- dir(in_direx)
# remove multi-year datasets
f <- f[-grep("apboth.csv", f)]
f <- f[-grep("bsboth.csv", f)]
f <- f[-grep("waall.csv", f)]
f <- f[-grep("wcpboth.csv", f)]

# iterating for each dataset
ind_df <- list()
pop_df <- list()
for (i in 1:length(f)) {
  print(gsub(".csv", "", f[i]))
  dat <- read.csv(paste0(in_direx, f[i]))
  dat$id <- as.factor(as.character(dat$id))
  # reproject
  dat_proj <- proj_fun(x = dat, proj_in = proj_dec, proj_out = proj_utm) 
  # calculate population kernels
  dat_pop <- pop_kernel_poly(x = dat_proj, h = 186, cellsize = 50, iso_level = 50)
  # calculate individual kernels
  dat_ind <- ind_kernels_poly(x = dat_proj, h = 186, cellsize = 50, iso_level = 50)
  # project back to WGS84 if not WA
  if (dat$species[1] == "WA") {
  } else {
    proj4string(dat_ind) <- proj_utm
    dat_ind <- spTransform(dat_ind, CRS(proj_dec))
    proj4string(dat_pop) <- proj_utm
    dat_pop <- spTransform(dat_pop, CRS(proj_dec))
  }
  # convert to sf then to dataframes for plotting
  dat_ind_sf <- st_as_sf(dat_ind)
  ind_df[[i]] <- sf_to_df(dat_ind_sf, fill = TRUE, unlist = NULL)
  dat_pop_sf <- st_as_sf(dat_pop)
  pop_df[[i]] <- sf_to_df(dat_pop_sf, fill = TRUE, unlist = NULL)
  # create a merged individual polygon id
  ind_df[[i]]$id_polygon <- as.factor(as.character(paste(ind_df[[i]]$id, ind_df[[i]]$polygon_id, sep = "_")))
}

#### 3. SETTING UP PLOTTING MATERIALS ####

# colony location in WGS84 and utm coords
col_dec <- data.frame(x = -38.05, y = -54.00)
col_proj <- col_dec
coordinates(col_proj) <- ~x+y
proj4string(col_proj) <- proj_dec
col_proj <- spTransform(col_proj, CRS(proj_utm))
col_proj <- as.data.frame(col_proj)

# world file
country <- rgdal::readOGR("./Data_inputs/plotting/country.shp")
# wgs84
proj4string(country) <- proj_dec
country_dec <- st_as_sf(country)
country_dec <- sf_to_df(country_dec, fill = TRUE, unlist = NULL)
# utm
country_proj <- spTransform(country, CRS(proj_utm))
country_proj <- st_as_sf(country_proj)
country_proj <- sf_to_df(country_proj, fill = TRUE, unlist = NULL)

#### 4. PLOTTING FOR EACH SPECIES ####

plot_direx <- "./Plots/"

##### 4A. AP 2010  #####

# setting up graticule
ap_grat <- graticule(lons = seq(-120, 0, 10), lats = seq(-60, -20,10), xlim = c(-120, 00), ylim = c(-80, 0), proj = proj_dec)
ap_grat <- sf_to_df(st_as_sf(ap_grat), fill = TRUE, unlist = NULL)

p_ap10 <- ggplot() + geom_path(data = ap_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = country_dec, aes(x = x , y = y, group = polygon_id), fill = "grey80", colour = "grey30", size = 0.2) +
  geom_polygon(data = ind_df[[1]], aes(x=x, y=y, group = id_polygon), fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
  geom_polygon(data = pop_df[[1]], aes(x = x, y=y), fill = NA, col = "black", linetype = "dashed", size = 0.7)+
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(), legend.position = "none", axis.ticks=element_blank())  +
  coord_fixed(ratio = 1, xlim = c(-101, -37), ylim = c(-57, -33))+
  geom_point(data = col_dec, aes(x = x, y = y),  pch = 24, size = 3, fill = "yellow")  
p_ap10

out.name <- paste0(plot_direx, "4_ap10.png")
ggsave(out.name, dpi = 600, width = 3000, height = 3000, units = "px")
dev.off()

##### 4B. AP 2011  #####

p_ap11 <- ggplot() + geom_path(data = ap_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = country_dec, aes(x = x , y = y, group = polygon_id), fill = "grey80", colour = "grey30", size = 0.2) +
  geom_polygon(data = ind_df[[2]], aes(x=x, y=y, group = id_polygon), fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
  geom_polygon(data = pop_df[[3]], aes(x = x, y=y), fill = NA, col = "black", linetype = "dashed", size = 0.7)+
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),legend.position = "none", axis.ticks=element_blank())  +
  coord_fixed(ratio = 1, xlim = c(-101, -37), ylim = c(-57, -33))+
  geom_point(data = col_dec, aes(x = x, y = y),  pch = 24, size = 3, fill = "yellow") 
p_ap11

out.name <- paste0(plot_direx, "4_ap11.png")
ggsave(out.name, dpi = 600, width = 3000, height = 3000, units = "px")
dev.off()

##### 4C. BBA 2002  #####

bba_grat <- graticule(lons = seq(-70, 150, 10), lats = seq(-60, -10,10), xlim = c(-70, 150), ylim = c(-80, 0), proj = proj_dec)
bba_grat <- sf_to_df(st_as_sf(bba_grat), fill = TRUE, unlist = NULL)

p_bba02 <- ggplot() + geom_path(data = bba_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = country_dec, aes(x = x , y = y, group = polygon_id), fill = "grey80",colour = "grey30", size = 0.2) +
  geom_polygon(data = ind_df[[3]], aes(x=x, y=y, group = id_polygon), fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
  geom_polygon(data = pop_df[[3]], aes(x = x, y=y, group = polygon_id), fill = NA, col = "black", linetype = "dashed", size = 0.65)+
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),legend.position = "none", axis.ticks=element_blank())  +
  coord_fixed(ratio = 1, xlim = c(-60, 141), ylim = c(-58, -11))+
  geom_point(data = col_dec, aes(x = x, y = y),  pch = 24, size = 2.8, fill = "yellow") 
p_bba02

out.name <- paste0(plot_direx, "4_bba02.png")
ggsave(out.name, dpi = 600, width = 5500, height = 5500, units = "px")
dev.off()

##### 4D. BBA 2003  #####

p_bba03 <- ggplot() + geom_path(data = bba_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = country_dec, aes(x = x , y = y, group = polygon_id), fill = "grey80",colour = "grey30", size = 0.2) +
  geom_polygon(data = ind_df[[4]], aes(x=x, y=y, group = id_polygon), fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
  geom_polygon(data = pop_df[[4]], aes(x = x, y=y, group = polygon_id), fill = NA, col = "black", linetype = "dashed", size = 0.65)+
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text = element_blank(), axis.title = element_blank(), legend.position = "none", axis.ticks=element_blank())  +
  coord_fixed(ratio = 1, xlim = c(-60, 141), ylim = c(-58, -11))+
  geom_point(data = col_dec, aes(x = x, y = y),  pch = 24, size = 2.8, fill = "yellow")
p_bba03

out.name <- paste0(plot_direx, "4_bba03.png")
ggsave(out.name, dpi = 600, width = 5500, height = 5500, units = "px")
dev.off()

##### 4E. BP 2011 #####

bp11_grat <- graticule(lons = seq(-180, 60, 10), lats = seq(-80, -20,10), xlim = c(-180, 60), ylim = c(-80, 0), proj = proj_dec)
bp11_grat <- sf_to_df(st_as_sf(bp11_grat), fill = TRUE, unlist = NULL)

p_bp11 <- ggplot() + geom_path(data = bp11_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = country_dec, aes(x = x , y = y, group = polygon_id), fill = "grey80",colour = "grey30", size = 0.2) +
  geom_polygon(data = ind_df[[5]], aes(x=x, y=y, group = id_polygon), fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
  geom_polygon(data = pop_df[[5]], aes(x = x, y=y, group = polygon_id), fill = NA, col = "black", linetype = "dashed", size = 0.65)+
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                              axis.text = element_blank(), axis.title = element_blank(),legend.position = "none", axis.ticks=element_blank())  +
  coord_fixed(ratio = 1, xlim = c(-170, 10), ylim = c(-70, -38))+
  geom_point(data = col_dec, aes(x = x, y = y),  pch = 24, size = 3, fill = "yellow") 
p_bp11

out.name <- paste0(plot_direx, "4_bp11.png")
ggsave(out.name, dpi = 600, width = 5500, height = 5500, units = "px")
dev.off()

##### 4F. BS 2002 #####

bs02_grat <- graticule(lons = seq(-80, 60, 10), lats = seq(-60, -20,10), xlim = c(-80, 60), ylim = c(-80, 0), proj = proj_dec)
bs02_grat <- sf_to_df(st_as_sf(bs02_grat), fill = TRUE, unlist = NULL)

p_bs02 <- ggplot() + geom_path(data = bs02_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = country_dec, aes(x = x , y = y, group = polygon_id), fill = "grey80", colour = "grey30", size = 0.2) +
  geom_polygon(data = ind_df[[6]], aes(x=x, y=y, group = id_polygon), fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
  geom_polygon(data = pop_df[[6]], aes(x = x, y=y, group = polygon_id), fill = NA, col = "black", linetype = "dashed", size = 1)+
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(), legend.position = "none", axis.ticks=element_blank())  +
  coord_fixed(ratio = 1, xlim = c(-66, -19), ylim = c(-67, -24))+
  geom_point(data = col_dec, aes(x = x, y = y),  pch = 24, size = 4.8, fill = "yellow") 
p_bs02

out.name <- paste0(plot_direx, "4_bs02.png")
ggsave(out.name, dpi = 600, width = 3000, height = 3000, units = "px")
dev.off()

##### 4G. BS 2012 #####

bs12_grat <- graticule(lons = seq(-80, 60, 10), lats = seq(-60, -20,10), xlim = c(-80, 60), ylim = c(-80, 0), proj = proj_dec)
bs12_grat <- sf_to_df(st_as_sf(bs12_grat), fill = TRUE, unlist = NULL)

p_bs12 <- ggplot() + geom_path(data = bs12_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = country_dec, aes(x = x , y = y, group = polygon_id), fill = "grey80",colour = "grey30", size = 0.2) +
  geom_polygon(data = ind_df[[7]], aes(x=x, y=y, group = id_polygon), fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
  geom_polygon(data = pop_df[[7]], aes(x = x, y=y, group = polygon_id), fill = NA, col = "black", linetype = "dashed", size = 0.55)+
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(), legend.position = "none", axis.ticks=element_blank())  +
  coord_fixed(ratio = 1, xlim = c(-66, 27), ylim = c(-67, -24))+
  geom_point(data = col_dec, aes(x = x, y = y),  pch = 24, size = 2.8, fill = "yellow") 
p_bs12

out.name <- paste0(plot_direx, "4_bs12.png")
ggsave(out.name, dpi = 600, width = 3500, height = 3500, units = "px")
dev.off()

##### 4H. GHA 1999 #####

gha99_grat <- graticule(lons = seq(-160, 140, 10), lats = seq(-80, -20,10), xlim = c(-160, 140), ylim = c(-80, 0), proj = proj_dec)
gha99_grat <- sf_to_df(st_as_sf(gha99_grat), fill = TRUE, unlist = NULL)

p_gha99 <- ggplot() + geom_path(data = gha99_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = country_dec, aes(x = x , y = y, group = polygon_id), fill = "grey80",colour = "grey30", size = 0.05) +
  geom_polygon(data = ind_df[[8]], aes(x=x, y=y, group = id_polygon), fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
  geom_polygon(data = pop_df[[8]], aes(x = x, y=y, group = polygon_id), fill = NA, col = "black", linetype = "dashed", size = 0.65)+
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),legend.position = "none", axis.ticks=element_blank())  +
  coord_fixed(ratio = 1, xlim = c(-146, 120), ylim = c(-70, -28))+
  geom_point(data = col_dec, aes(x = x, y = y),  pch = 24, size = 3.3, fill = "yellow") 
p_gha99

out.name <- paste0(plot_direx, "4_gha99.png")
ggsave(out.name, dpi = 600, width = 6000, height = 6000, units = "px")
dev.off()

##### 4I. LMA 2004 #####

lma_grat <- graticule(lons = seq(-90, 120, 10), lats = seq(-80, -20,10), xlim = c(-90, 120), ylim = c(-80, 0), proj = proj_dec)
lma_grat <- sf_to_df(st_as_sf(lma_grat), fill = TRUE, unlist = NULL)

p_lma04 <- ggplot() + geom_path(data = lma_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = country_dec, aes(x = x , y = y, group = polygon_id), fill = "grey80", colour = "grey30", size = 0.2) +
  geom_polygon(data = ind_df[[9]], aes(x=x, y=y, group = id_polygon), fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
  geom_polygon(data = pop_df[[9]], aes(x = x, y=y, group = polygon_id), fill = NA, col = "black", linetype = "dashed", size = 0.65)+
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),legend.position = "none", axis.ticks=element_blank())  +
  coord_fixed(ratio = 1, xlim = c(-70, 110), ylim = c(-70, -38))+
  geom_point(data = col_dec, aes(x = x, y = y),  pch = 24, size = 3, fill = "yellow") 
p_lma04

out.name <- paste0(plot_direx, "4_lma04.png")
ggsave(out.name, dpi = 600, width = 5500, height = 5500, units = "px")
dev.off()

##### 4J. LMA 2005 #####

p_lma05 <- ggplot() + geom_path(data = lma_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = country_dec, aes(x = x , y = y, group = polygon_id), fill = "grey80",colour = "grey30", size = 0.2) +
  geom_polygon(data = ind_df[[10]], aes(x=x, y=y, group = id_polygon), fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
  geom_polygon(data = pop_df[[10]], aes(x = x, y=y, group = polygon_id), fill = NA, col = "black", linetype = "dashed", size = 0.65)+
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),legend.position = "none", axis.ticks=element_blank())  +
  coord_fixed(ratio = 1, xlim = c(-70, 110), ylim = c(-70, -38))+
  geom_point(data = col_dec, aes(x = x, y = y),  pch = 24, size = 3, fill = "yellow") 
p_lma05

out.name <- paste0(plot_direx, "4_lma05.png")
ggsave(out.name, dpi = 600, width = 5500, height = 5500, units = "px")
dev.off()

##### 4K. NGP 2000 #####

p_ngp <- ggplot() + geom_path(data = bs02_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = country_dec, aes(x = x , y = y, group = polygon_id), fill = "grey80",colour = "grey30", size = 0.2) +
  geom_polygon(data = ind_df[[11]], aes(x = x, y = y, group = id_polygon), fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
  geom_polygon(data = pop_df[[11]], aes(x = x, y  =y, group = polygon_id), fill = NA, col = "black", linetype = "dashed", size = 1)+
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),legend.position = "none", axis.ticks=element_blank())  +
  coord_fixed(ratio = 1, xlim = c(-66, -19), ylim = c(-67, -24))+
  geom_point(data = col_dec, aes(x = x, y = y),  pch = 24, size = 4.8, fill = "yellow") 
p_ngp

out.name <- paste0(plot_direx, "4_ngp00.png")
ggsave(out.name, dpi = 600, width = 3000, height = 3000, units = "px")
dev.off()

##### 4L. SGP 2000 #####

p_sgp <- ggplot() + geom_path(data = bs12_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = country_dec, aes(x = x , y = y, group = polygon_id), fill = "grey80", colour = "grey30", size = 0.2) +
  geom_polygon(data = ind_df[[12]], aes(x = x, y = y, group = id_polygon),  fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
 geom_polygon(data = pop_df[[12]], aes(x = x, y = y, group = polygon_id), fill = NA, col = "black", linetype = "dashed", size = 0.55)+
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(), legend.position = "none", axis.ticks = element_blank())  +
  coord_fixed(ratio = 1, xlim = c(-66, 27), ylim = c(-67, -24))+
  geom_point(data = col_dec, aes(x = x, y = y),  pch = 24, size = 2.8, fill = "yellow") 
p_sgp

out.name <- paste0(plot_direx, "4_sgp00.png")
ggsave(out.name, dpi = 600, width = 3500, height = 3500, units = "px")
dev.off()

##### 4M. WA 2004 #####

# setting up graticule, plus 2nd graticule as there is a gap between 160 and 180
wa_grat <- graticule(lons = seq(-180, 180, 20), lats = c(-60, -40,-20), xlim = c(-180, 160), ylim = c(-90, 0), proj = proj_utm)
wa_grat <- sf_to_df(st_as_sf(wa_grat), fill = TRUE, unlist = NULL)
wa_grat2 <- graticule(lons = seq(160, 180, 20), lats = c(-60, -40,-20), xlim = c(160, 180), ylim = c(-90, 0), proj = proj_utm)
wa_grat2 <- sf_to_df(st_as_sf(wa_grat2), fill = TRUE, unlist = NULL)

p_wa04 <- ggplot() + geom_path(data = wa_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_path(data = wa_grat2, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = ind_df[[13]], aes(x=x, y=y, group = id_polygon), fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
  geom_polygon(data = pop_df[[13]], aes(x = x, y=y, group = polygon_id), fill = NA, col = "black", linetype = "dashed", size = 0.65)+
  geom_polygon(data = country_proj, aes(x = x , y = y, group = polygon_id), fill = "grey80",colour = "grey30", size = 0.2) +
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),legend.position = "none", axis.ticks=element_blank()) +
  coord_fixed(ratio = 1, xlim = c(-5500000, 5800000), ylim = c(-5300000, 5500000))+
  geom_point(data = col_proj, aes(x = x, y = y),  pch = 24, size = 3, fill = "yellow")
p_wa04

out.name <- paste0(plot_direx, "4_wa04.png")
ggsave(out.name, dpi = 600, width = 3000, height = 3000, units = "px")
dev.off()

##### 4N. WA 2008 #####

p_wa08 <- ggplot() + geom_path(data = wa_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_path(data = wa_grat2, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = ind_df[[14]], aes(x=x, y=y, group = id_polygon), fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
 geom_polygon(data = pop_df[[14]], aes(x = x, y=y, group = polygon_id), fill = NA, col = "black", linetype = "dashed", size = 0.65)+
  geom_polygon(data = country_proj, aes(x = x , y = y, group = polygon_id), fill = "grey80", colour = "grey30", size = 0.2) +
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(), legend.position = "none", axis.ticks=element_blank()) +
  coord_fixed(ratio = 1, xlim = c(-5500000, 5800000), ylim = c(-5300000, 5500000))+
  geom_point(data = col_proj, aes(x = x, y = y),  pch = 24, size = 3, fill = "yellow")
p_wa08

out.name <- paste0(plot_direx, "4_wa08.png")
ggsave(out.name, dpi = 600, width = 3000, height = 3000, units = "px")
dev.off()

##### 4O. WA 2009 #####

p_wa09 <- ggplot() + geom_path(data = wa_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_path(data = wa_grat2, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = ind_df[[15]], aes(x=x, y=y, group = id_polygon), fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
  geom_polygon(data = pop_df[[15]], aes(x = x, y=y, group = polygon_id), fill = NA, col = "black", linetype = "dashed", size = 0.65)+
  geom_polygon(data = country_proj, aes(x = x , y = y, group = polygon_id), fill = "grey80",colour = "grey30", size = 0.2) +
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),legend.position = "none", axis.ticks=element_blank()) +
  coord_fixed(ratio = 1, xlim = c(-5500000, 5800000), ylim = c(-5300000, 5500000))+
  geom_point(data = col_proj, aes(x = x, y = y),  pch = 24, size = 3, fill = "yellow")
p_wa09

out.name <- paste0(plot_direx, "4_wa09.png")
ggsave(out.name, dpi = 600, width = 3000, height = 3000, units = "px")
dev.off()

#####  4P. WCP 2003 #####

p_wcp03 <- ggplot() + geom_path(data = ap_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = country_dec, aes(x = x , y = y, group = polygon_id), fill = "grey80",colour = "grey30", size = 0.2) +
  geom_polygon(data = ind_df[[16]], aes(x=x, y=y, group = id_polygon), fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
  geom_polygon(data = pop_df[[16]], aes(x = x, y=y), fill = NA, col = "black", linetype = "dashed", size = 0.7)+
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(), legend.position = "none", axis.ticks=element_blank())  +
  coord_fixed(ratio = 1, xlim = c(-101, -37), ylim = c(-57, -26))+
  geom_point(data = col_dec, aes(x = x, y = y),  pch = 24, size = 3, fill = "yellow") 
p_wcp03

out.name <- paste0(plot_direx, "4_wcp03.png")
ggsave(out.name, dpi = 600, width = 3000, height = 3000, units = "px")
dev.off()

##### 4Q. WCP 2013 #####

p_wcp13 <- ggplot() + geom_path(data = ap_grat, aes(x = x , y = y, group = sfg_id), colour = "grey50", size = 0.05) +
  geom_polygon(data = country_dec, aes(x = x , y = y, group = polygon_id), fill = "grey80",colour = "grey30", size = 0.2) +
  geom_polygon(data = ind_df[[17]], aes(x=x, y=y, group = id_polygon), fill = "royalblue", alpha = 0.1, size = 0.2, col = "grey30") +
  geom_polygon(data = pop_df[[17]], aes(x = x, y=y, group = polygon_id), fill = NA, col = "black", linetype = "dashed", size = 0.7)+
  theme_bw() + borders(fill = "gray") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),legend.position = "none", axis.ticks=element_blank())  +
  coord_fixed(ratio = 1, xlim = c(-101, -37), ylim = c(-57, -26))+
  geom_point(data = col_dec, aes(x = x, y = y),  pch = 24, size = 3, fill = "yellow") 
p_wcp13

out.name <- paste0(plot_direx, "4_wcp13.png")
ggsave(out.name, dpi = 600, width = 3000, height = 3000, units = "px")
dev.off()
