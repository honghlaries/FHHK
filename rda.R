# Initialization
rm(list = ls())
source("constant.R");source("anaTls_spatialView.R");source("anaTls_multivariate.R");
pkgInitialization(c("dplyr","tidyr","sp","gstat","gridExtra","ColorPalette"))
source("grid.R")

# Functions 

# Example 
env <- datareadln() %>% 
  mutate(proton = 10^(-pH)) %>%
  select(depth,distance,salinity,proton,Al,Fe,Mn,orgC,AVS,clay,silt,sand) 
trait <- datareadln() %>% select(Pb:Cd)
samptag <- datareadln() %>% select(siteID)
rda <- rdaLoadingCal(env, trait, samptag, 3, dir = "rda", log = F)
rdaLoadingPlot(rda,3) -> plot.load


dat <- rda$sampload %>%
  group_by(siteID) %>%
  summarise(RDA1 = mean(RDA1), RDA2 = mean(RDA2), RDA3 = mean(RDA3)) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  dplyr::select(lat,lon,RDA1,RDA2,RDA3)
dat <- as.data.frame(dat)
coordinates(dat) <- ~lon+lat

doKrig(dat1, dat.grid, tag = "group1", cutoff = 1.7,
       modsel = vgm(0.25,"Sph",1,0), quietmode = T)

grid.value.tot <- NULL
grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "RDA1", cutoff = 1.5, 
                                   modsel = vgm(0.015,"Sph",0.7), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "RDA1")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "RDA2", cutoff = 2.3, 
                                   modsel = vgm(0.05,"Mat",1, kappa = 0.25), quietmode = T)) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "RDA2")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "RDA3", cutoff = 1.5, 
                                   modsel = vgm(0.05,"Sph",1), quietmode = T)) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "RDA3")))

spView.grid(dat = grid.value.tot, leg.name = "Loading",
            grad.value = c(-0.25,-0.2,-0.1,0,0.1,0.2,0.25), 
            grad.tag = c(-0.25,-0.2,-0.1,0,0.1,0.2,0.25), 
            #grad.col = c("#D2E9FF","#97CBFF","#66B3FF","#2894FF","#0072E3","#005AB5","#003D79"),
            lonRange = lonRange,
            latRange = latRange, pncol = 1) -> plot.sp

grid.arrange(plot.load, plot.sp, ncol = 2, widths = c(5,5), heights = 10) -> plot.gather

# saving plot
ggsave(plot = plot.load, filename = "rda/rdaloading.png", dpi = 600)
ggsave(plot = plot.sp, filename = "rda/rdaspatialView.png", dpi = 600)
ggsave(plot = plot.gather, filename = "rda/gather_rdaPlot.png", dpi = 600)