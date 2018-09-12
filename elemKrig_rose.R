## Initialization
rm(list = ls())
source("constant.R");source("anaTls_spatialView.R");source("uniTls_csv2latex.R")
pkgInitialization(c("dplyr","tidyr","sp","gstat","ggplot2","RColorBrewer","gridExtra"))
source("grid.R");source("grid_resamp.R")
dirInitialization(c("element","element/krig"))

## Functions 
distCalc <- function(lata, lona, latb, lonb, r = 6371){
  r*acos(cos(lata/180*pi)*cos(lona/180*pi)*cos(latb/180*pi)*cos(lonb/180*pi)+
           cos(lata/180*pi)*sin(lona/180*pi)*cos(latb/180*pi)*sin(lonb/180*pi)+
           sin(lata/180*pi)*sin(latb/180*pi))
}

azimuthCalc <- function(lata, lona, latb, lonb, latc =90, lonc = lonb){
  p <- (acos(cos(lata/180*pi)*cos(lona/180*pi)*cos(latb/180*pi)*cos(lonb/180*pi)+
               cos(lata/180*pi)*sin(lona/180*pi)*cos(latb/180*pi)*sin(lonb/180*pi)+
               sin(lata/180*pi)*sin(latb/180*pi))+
          (90-lata)/180*pi +
          (90-latb)/180*pi )/2.0
  cos.azi <- sqrt(sin(p)*sin(p-(90-latb)/180*pi)/sin((90-lata)/180*pi)/
                    sin(acos(cos(lata/180*pi)*cos(lona/180*pi)*cos(latb/180*pi)*cos(lonb/180*pi)+
                               cos(lata/180*pi)*sin(lona/180*pi)*cos(latb/180*pi)*sin(lonb/180*pi)+
                               sin(lata/180*pi)*sin(latb/180*pi))))
   
  if(lonb >= lona) {
    if(latb >= lata) {0.5*pi - 2*acos(cos.azi)} else {2.5*pi - 2*acos(cos.azi)}
  } else {
    if(latb >= lata) {0.5*pi - 2*acos(cos.azi)} else {2.5*pi - 2*acos(cos.azi)}
  }
}

azimuthCalc.rough <- azimuthCalc <- function(lata, lona, latb, lonb){
  tan.azi <- (latb-lata)/(lonb-lona)
  
  if(lonb >= lona) {
    if(latb >= lata) {atan(tan.azi)} else {2*pi + atan(tan.azi)}
  } else {
    pi + atan(tan.azi)
  }
}


## data calculation
grid.perm.tot <- read.csv("data/result_element_perm.csv") 

grid.content.tot <- NULL

seldat <- grid.perm.tot %>%
  dplyr::filter(trait == "Pb") %>%
  dplyr::mutate(delta = var.samp-var.site)
coordinates(seldat) <- ~lon+lat
grid.content.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean",
                                          cutoff = 1,
                                          modsel = vgm(0.02,"Sph",0.5),
                                          quietmode = T)) %>%
  dplyr::select(lon, lat, mean = var1.pred)
grid.content.tot <- rbind(grid.content.tot,
                          as.data.frame(cbind(grid.content.mean,
                                              trait = "Pb")))

seldat <- grid.perm.tot %>% filter(trait == "Cr") 
coordinates(seldat) <- ~lon+lat
grid.content.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                          cutoff = 1.5, 
                                          modsel = vgm(0.04,"Sph",0.5), 
                                          quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
grid.content.tot <- rbind(grid.content.tot, 
                          as.data.frame(cbind(grid.content.mean, 
                                              trait = "Cr")))

seldat <- grid.perm.tot %>% filter(trait == "Ni") 
coordinates(seldat) <- ~lon+lat
grid.content.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                          cutoff = 1.5, 
                                          modsel = vgm(0.05, "Mat", 1, kappa = 2), 
                                          quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
grid.content.tot <- rbind(grid.content.tot, 
                          as.data.frame(cbind(grid.content.mean,
                                              trait = "Ni")))

seldat <- grid.perm.tot %>% filter(trait == "Cu")
coordinates(seldat) <- ~lon+lat
grid.content.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                          cutoff = 1, 
                                          modsel = vgm(0.15, "Mat", 0.6, kappa = 3), 
                                          quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
grid.content.tot <- rbind(grid.content.tot, 
                          as.data.frame(cbind(grid.content.mean,
                                              trait = "Cu")))

seldat <- grid.perm.tot %>% filter(trait == "Zn") 
coordinates(seldat) <- ~lon+lat
grid.content.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                          cutoff = 0.7, 
                                          modsel = vgm(0.15, "Mat", 0.6, kappa = 3), 
                                          quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
grid.content.tot <- rbind(grid.content.tot, 
                          as.data.frame(cbind(grid.content.mean,
                                              trait = "Zn")))

seldat <- grid.perm.tot %>% filter(trait == "Cd") 
coordinates(seldat) <- ~lon+lat
grid.content.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean",
                                          cutoff = 0.8,
                                          modsel = vgm(0.15, "Mat", 0.6, kappa = 3.5),
                                          quietmode = T)) %>%
  dplyr::select(lon, lat, mean = var1.pred)
grid.content.tot <- rbind(grid.content.tot, 
                          as.data.frame(cbind(grid.content.mean,
                                              trait = "Cd")))

grid.content.tot$value <- exp(grid.content.tot$mean) 
grid.content.tot <- grid.content.tot %>% dplyr::select(-mean)
remove(grid.content.mean, grid.perm.tot, seldat)

grid.content.tot$dist <- distCalc(coo.oldmouth$lat, coo.oldmouth$lon, 
                                  grid.content.tot$lat, grid.content.tot$lon)

for(i in 1:length(grid.content.tot$lon)) {
  if (!i%%1000) print(paste("calculating:",i,"/",length(grid.content.tot$lon)))
  grid.content.tot$azi[i] <- azimuthCalc.rough(coo.oldmouth$lat, coo.oldmouth$lon, 
                                      grid.content.tot$lat[i], grid.content.tot$lon[i])
}

grid.content.tot$azideg <- ceiling(grid.content.tot$azi/3/pi*180) *3

grid.content.tot$azideggroup <- ceiling((grid.content.tot$azideg)/30) 
#grid.content.tot$azideggroup[grid.content.tot$azideg > 337.5] <- 9

p <- 
  ggplot() +
  geom_point(aes(x = azideg, y = value, col = dist), alpha = 0.05, 
               data = grid.content.tot) +
  scale_color_gradient(low = "yellow", high = "blue") +
  #coord_polar(start = -pi/2, direction = -1) +
  facet_wrap(~trait, scale = "free") +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.position = "none")
ggsave(plot = p, filename = "element/gather_pointrose_element.png", width = 8, height = 6, dpi = 300)
rm(p)

p <- 
  ggplot() +
  geom_point(aes(x = dist, y = value), alpha = 0.01, data = grid.content.tot) +
  #geom_smooth(aes(x = dist, y = value), method = "loess", data = grid.content.tot) +
  #coord_polar(start = -pi/2, direction = -1) +
  facet_grid(trait~azideggroup, scale = "free_y") +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.position = "none")
ggsave(plot = p, filename = "element/gather_point_element_dist.png", width = 12, height = 6, dpi = 300)
rm(p)
