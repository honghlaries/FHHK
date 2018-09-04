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

grid.content.tot$dist <- distCalc(grid.content.tot$lat, grid.content.tot$lon,
                                  coo.oldmouth$lat, coo.oldmouth$lon)
