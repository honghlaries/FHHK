## clean ----
rm(list = ls())
source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R");source("anaTls_spatialView.R");
pkgInitialization(c("dplyr","tidyr","sp","gstat"))

## Functions ----
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

longiRange <-  seq(from = 119.9, to = 121.8, length.out = 125)
latiRange <-  seq(from = 33.7, to = 34.9, length.out = 125)
dat.grid <- data.frame(lat = c(1), lon = c(1))
for (i in 1:125) {
  for (j in 1:125) {
    #if((longiRange[j] - 119.9) * (latiRange[i] - 33.7) - (longiRange[j] - 120.6) * (latiRange[i] - 34.5) > 0) {
    dat.grid <- rbind(dat.grid, c(latiRange[i],longiRange[j]))
    #}
  }
}
dat.grid <- dat.grid[-1,]
coordinates(dat.grid) <- ~lon+lat

# for Igeo
background <- datareadln() %>% 
  gather(trait, bk, Al:orgC) %>%
  select(siteID, trait, bk) %>%
  group_by(trait) %>%
  summarise(bk = mean(bk)) %>% 
  mutate(bk2 = bk / 17092.2)

dat <- datareadln() %>% 
  gather(trait, value, Al:orgC) %>% 
  inner_join(background, by = c("trait" = "trait")) %>%
  mutate(value = log2(value / bk /1.5)) %>%
  select(siteID, trait, value) %>%
  group_by(siteID, trait) %>%
  summarise(value = mean(value)) %>%
  spread(trait, value) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  dplyr::select(siteID:depth,Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd,C,N,S,orgC) 
dat <- as.data.frame(dat)
coordinates(dat) <- ~lon+lat

grid.value.tot <- NULL
grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Fe", cutoff = 2,
                                   modsel = vgm(0.1,"Mat",1,0.02,kappa = 1), dir = "riskAssment/krig/igeo/igeo")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Fe")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Mn", cutoff = 2,
                                   modsel = vgm(0.02,"Mat",1,0.04,kappa = 1), dir = "riskAssment/krig/igeo")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Mn")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Pb", cutoff = 2,
                                   modsel = vgm(0.06,"Sph",1), dir = "riskAssment/krig/igeo")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Pb")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Cr", cutoff = 1.7,
                                   modsel = vgm(0.1,"Mat",0.7,kappa = 1), dir = "riskAssment/krig/igeo")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cr")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Ni", cutoff = 1.5,
                                   modsel = vgm(0.15,"Sph",1), dir = "riskAssment/krig/igeo")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Ni")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Cu", cutoff = 1.5,
                                   modsel = vgm(0.35,"Sph",0.75), dir = "riskAssment/krig/igeo")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cu")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Zn", cutoff = 1.5,
                                   modsel = vgm(0.4,"Sph",0.5), dir = "riskAssment/krig/igeo")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Zn")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "As", cutoff = 1.5,
                                   modsel = vgm(0.5,"Sph",0.5), dir = "riskAssment/krig/igeo")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "As")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Cd", cutoff = 1.5,
                                   modsel = vgm(0.4,"Sph",0.5), dir = "riskAssment/krig/igeo")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cd")))

spView.grid(dat = grid.value.tot %>% 
              filter(trait == "Fe"|trait == "Mn"|trait == "Cu"|
                       trait == "Zn"|trait == "Pb"|trait == "Cr"|
                       trait == "Ni"|trait == "As"|trait == "Cd"),
            leg.name = "Geo-accumulation Index",grad.value = c(-2,-1,0,1), 
            grad.tag = c(-2,-1,0,1), dir = "riskAssment/krig",file = "krig_igeo.png",
            lonRange = c(119.2,121.8),latRange = c(33.7,35))

spView.grid(dat = grid.value.tot %>% 
              filter(trait == "Cu"|trait == "Zn"|trait == "Pb"|
                       trait == "Ni"|trait == "As"|trait == "Cd"),
            leg.name = "Geo-accumulation Index",grad.value = c(-2,-1,0,1), 
            grad.tag = c(-2,-1,0,1), dir = "riskAssment/krig",file = "krig_igeo_sel.png",
            lonRange = c(119.2,121.8),latRange = c(33.7,35))

# for enrichment factor
dat <- datareadln() %>% 
  gather(trait, value, Al:Cd) %>% 
  inner_join(datareadln() %>%
               select(siteID, Al), by = c("siteID" = "siteID")) %>%
  mutate(value = value / Al) %>%
  inner_join(background, by = c("trait" = "trait")) %>%
  mutate(value = value / bk2) %>%
  select(siteID, trait, value) %>%
  group_by(siteID, trait) %>%
  summarise(value = mean(value)) %>%
  spread(trait, value) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  dplyr::select(siteID:depth,Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd) %>%
  dplyr::mutate(latitudes = conveySp(latitudes),
                longitude = conveySp(longitude))
dat <- as.data.frame(dat)
coordinates(dat) <- ~longitude+latitudes

doKrig(dat, dat.grid, tag = "Al", suffix = "_Ef", modsel = vgm(0.2,"Sph",0.7))
doKrig(dat, dat.grid, tag = "Fe", suffix = "_Ef", modsel = vgm(0.04,"Lin",0,0.003))
doKrig(dat, dat.grid, tag = "Mn", suffix = "_Ef", modsel = vgm(0.03,"Sph",0.7))
doKrig(dat, dat.grid, tag = "Pb", suffix = "_Ef", modsel = vgm(0.025,"Lin",0,0.015))
doKrig(dat, dat.grid, tag = "Cr", suffix = "_Ef", modsel = vgm(0.03,"Lin",0))
doKrig(dat, dat.grid, tag = "Ni", suffix = "_Ef", modsel = vgm(0.05,"Lin",0))
doKrig(dat, dat.grid, tag = "Cu", suffix = "_Ef", modsel = vgm(0.2,"Lin",0))
doKrig(dat, dat.grid, tag = "Zn", suffix = "_Ef", modsel = vgm(0.4,"Lin",0))
doKrig(dat, dat.grid, tag = "As", suffix = "_Ef", modsel = vgm(0.4,"Lin",0))
doKrig(dat, dat.grid, tag = "Cd", suffix = "_Ef", modsel = vgm(0.3,"Lin",0))

