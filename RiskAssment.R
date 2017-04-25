## clean 
rm(list = ls())
source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R");source("anaTls_spatialView.R");
pkgInitialization(c("dplyr","tidyr","sp","gstat"))
source("grid.R")
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
            lonRange = c(119.2,121.8),latRange = c(33.7,35))  -> plot.igeo

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
  dplyr::select(siteID:depth,Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd) 
dat <- as.data.frame(dat)
coordinates(dat) <- ~lon+lat

grid.value.tot <- NULL
grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Pb", cutoff = 1.5,
                                   modsel = vgm(0.1,"Sph",0.75), dir = "riskAssment/krig/ef")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Pb")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Ni", cutoff = 1.5,
                                   modsel = vgm(0.06,"Sph",1), dir = "riskAssment/krig/ef")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Ni")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Cu", cutoff = 1.5,
                                   modsel = vgm(0.25,"Sph",0.75), dir = "riskAssment/krig/ef")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cu")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Zn", cutoff = 1.5,
                                   modsel = vgm(0.2,"Sph",0.5), dir = "riskAssment/krig/ef")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Zn")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "As", cutoff = 1.5,
                                   modsel = vgm(0.6,"Sph",0.5), dir = "riskAssment/krig/ef")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "As")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Cd", cutoff = 1.5,
                                   modsel = vgm(0.15,"Sph",0.5), dir = "riskAssment/krig/ef")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cd")))

spView.grid(dat = grid.value.tot %>% 
              filter(trait == "Cu"|trait == "Zn"|trait == "Pb"|
                       trait == "Ni"|trait == "As"|trait == "Cd"),
            leg.name = "Enrichment Factor",grad.value = c(0,0.67,1,1.5,2,3,4), 
            grad.tag = c(0,0.67,1,1.5,2,3,4), dir = "riskAssment/krig/",file = "krig_enrichmentfactor.png",
            lonRange = c(119.2,121.8),latRange = c(33.7,35))  -> plot.ef

ggsave(filename = "riskAssment/krig/gather_krig_risk.png",  
       plot = grid.arrange(plot.igeo,plot.ef, ncol = 1, widths = c(15), heights = c(10,10)))
