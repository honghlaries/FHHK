# Initialization
rm(list = ls())
source("constant.R");source("anaTls_spatialView.R");source("anaTls_multivariate.R");
pkgInitialization(c("dplyr","tidyr","sp","gstat","gridExtra"))
source("grid.R")

# Functions 

# Example 
env <- datareadln() %>% select(depth,Al,Fe,Mn,orgC,AVS,clay,silt,sand)
trait <- datareadln() %>% select(Pb:Cd)
samptag <- datareadln() %>% select(siteID)
rda <- rdaLoadingCal(env, trait, samptag, dir = "rda", log = T)
rdaLoadingPlot(rda) -> plot.load


dat <- rda$sampload %>%
  group_by(siteID) %>%
  summarise(RDA1 = mean(RDA1), RDA2 = mean(RDA2)) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  dplyr::select(lat,lon,RDA1,RDA2)
dat <- as.data.frame(dat)
coordinates(dat) <- ~lon+lat

grid.value.tot <- NULL
doKrig(dat, dat.grid, tag = "RDA1", cutoff = 1.5, modsel = vgm(0.015,"Sph",0.7), dir = "rda/krig") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "RDA1")))

doKrig(dat, dat.grid, tag = "RDA2", cutoff = 2.3, modsel = vgm(0.05,"Mat",1, kappa = 0.25), dir = "rda/krig") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "RDA2")))

spView.grid(dat = grid.value.tot, leg.name = "Loading",
            grad.value = c(-0.25,-0.2,-0.1,0,0.1,0.2,0.3), 
            grad.tag = c(-0.25,-0.2,-0.1,0,0.1,0.2,0.3), 
            binwidth = 0.05, lonRange = c(119.2,121.8),
            latRange = c(33.7,35), pncol = 1) -> plot.sp

grid.arrange(plot.load, plot.sp, ncol = 2, widths = c(5,5), heights = 5) -> plot.gather

# saving plot
ggsave(plot = plot.load, filename = "rda/rdaloading.png", dpi = 600)
ggsave(plot = plot.sp, filename = "rda/rdaspatialView.png", dpi = 600)
ggsave(plot = plot.gather, filename = "rda/gather_rdaPlot.png", dpi = 600)