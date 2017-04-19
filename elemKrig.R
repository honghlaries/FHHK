## Initialization
rm(list = ls())
source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R");source("anaTls_spatialView.R");
pkgInitialization(c("dplyr","tidyr","sp","gstat"))

## Functions 
datareadln <- function() { ## data readln
  pkgLoad("dplyr");pkgLoad("tidyr")
  read.csv("data/result_element.csv") %>%
    dplyr::inner_join(read.csv("data/result_grainSize.csv"), by = c("sampleID" = "sampleID")) %>%
    dplyr::right_join(read.csv("data/meta_splList.csv"), by = c("sampleID" = "sampleID")) %>%
    dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
    dplyr::mutate(orgC = C.ac, AVS = S - S.ac, 
                  isComplete = complete.cases(Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd,C,N,S,orgC,AVS,clay,silt,sand)) %>%
    dplyr::select(siteID:depth,isComplete,Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd,C,N,S,orgC,AVS,clay,silt,sand)
}

## Example
dat <- datareadln() %>% 
  gather(trait, value, Al:sand) %>%
  select(siteID, trait, value) %>%
  group_by(siteID, trait) %>%
  summarise(value = mean(value)) %>%
  spread(trait, value) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  dplyr::select(siteID:depth,Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd,C,N,S,orgC,AVS,clay,silt,sand) %>%
  dplyr::mutate(AvsRatio = AVS / orgC)
dat <- as.data.frame(dat)
coordinates(dat) <- ~lon+lat

longiRange <-  seq(from = 119.9, to = 121.8, length.out = 125)
latiRange <-  seq(from = 33.7, to = 34.9, length.out = 125)
dat.grid <- data.frame(lat = c(1), lon = c(1))
for (i in 1:125) {
  for (j in 1:125) {
    if((longiRange[j] - 119.9) * (latiRange[i] - 33.7) - (longiRange[j] - 120.6) * (latiRange[i] - 34.5) > 0) {
      dat.grid <- rbind(dat.grid, c(latiRange[i],longiRange[j]))
    }
  }
}
dat.grid <- dat.grid[-1,]
coordinates(dat.grid) <- ~lon+lat

dat.grid.silt <- doKrig(dat, dat.grid, log(silt) ~ 1, tag = "silt", modsel = vgm(800,"Sph",1.0), outTpe = "data")
dat.grid.silt <- as.data.frame(dat.grid.silt) %>% select(lon, lat, silt = var1.pred)
coordinates(dat.grid.silt) <- ~lon+lat

grid.value.tot <- NULL
grid.value <- as.data.frame(doKrig(dat, dat.grid, log(Al) ~ 1, tag = "Al", 
                                   outTpe = "data",modsel = vgm(0.15,"Sph",0.5))) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Al")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, log(Fe) ~ 1, tag = "Fe", 
                                   outTpe = "data",modsel = vgm(0.01,"Mat",0.5,0.01,kappa = 1))) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Fe")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, log(Mn) ~ 1, tag = "Mn", 
                                   outTpe = "data",modsel = vgm(0.035,"Sph",0.3))) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Mn")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, log(Pb) ~ 1, tag = "Pb", 
                                   outTpe = "data",modsel = vgm(0.35,"Sph",0.5))) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Pb")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, log(Cr) ~ 1, tag = "Cr", 
                                   outTpe = "data",modsel = vgm(0.06,"Mat",0.7,kappa = 1))) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cr")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, log(Ni) ~ 1, tag = "Ni", 
                                   outTpe = "data",modsel = vgm(0.06,"Sph",1))) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Ni")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, log(Cu) ~ 1, tag = "Cu", 
                                   outTpe = "data",modsel = vgm(0.25,"Sph",0.75))) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cu")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, log(Zn) ~ 1, tag = "Zn", 
                                   outTpe = "data",modsel = vgm(0.1,"Sph",0.5))) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Zn")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, log(As) ~ 1, tag = "As", 
                                   outTpe = "data",modsel = vgm(0.3,"Sph",0.75))) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "As")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, log(Cd) ~ 1, tag = "Cd", 
                                   outTpe = "data",modsel = vgm(0.2,"Sph",0.7))) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cd")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, log(C) ~ 1, tag = "C", 
                                   outTpe = "data",modsel = vgm(0.15,"Sph",0.7))) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "C")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, log(orgC) ~ 1, tag = "orgC", 
                                   outTpe = "data",modsel = vgm(0.06,"Mat",0.7,kappa = 1))) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "orgC")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, log(N) ~ 1, tag = "N", 
                                   outTpe = "data",modsel = vgm(0.25,"Sph",1))) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "N")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, log(S) ~ 1, tag = "S", 
                                   outTpe = "data",modsel = vgm(0.4,"Sph",1))) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "S")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, log(AVS) ~ 1, tag = "AVS", 
                                   outTpe = "data",modsel = vgm(0.3,"Sph",0.7))) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "AVS")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, log(silt) ~ 1, tag = "silt", 
                                   outTpe = "data",modsel = vgm(200,"Sph",1.0))) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "silt")))

spView.grid(dat = grid.value.tot,file = "tmp.png",col.gardient,
            lonC = 120.8, latC = 34.2, lonRange = c(119.9,121.8),latRange = c(33.7,34.9), zoom = 9,
            maptype = "hybrid", mapsource = "google")


doKrig(dat, dat.grid, tag = "Al", cutoff = 1.5, krigFormula = log(Al) ~ 1,
       outTpe = "draft")
