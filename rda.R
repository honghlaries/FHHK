# Initialization
rm(list = ls())
source("constant.R");source("anaTls_spatialView.R");source("anaTls_multivariate.R");
pkgInitialization(c("dplyr","tidyr","sp","gstat","gridExtra","ColorPalette"))
source("grid.R")
dirInitialization(c("rda","rda/krig"))

# Functions 

# Example 
env <- datareadln() %>% 
  select(depth,distance,salinity,pH,Fe,Mn,orgC,AVS,clay,silt,sand) 
trait <- datareadln() %>% select(Pb:Cd)
samptag <- datareadln() %>% select(siteID)
rda <- rdaLoadingCal(env, trait, samptag, 3, dir = "rda", log = F)

rdaLoadingPlot(rda,3)+
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) -> plot.load


dat <- rda$sampload %>%
  group_by(siteID) %>%
  summarise(RDA1 = mean(RDA1), RDA2 = mean(RDA2), RDA3 = mean(RDA3)) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  dplyr::select(lat,lon,RDA1,RDA2,RDA3)
dat <- as.data.frame(dat)
coordinates(dat) <- ~lon+lat

grid.value.tot <- NULL
grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "RDA1", cutoff = 1.5, 
                                   modsel = vgm(0.015,"Sph",0.7), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "RDA1")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "RDA2", cutoff = 2, 
                                   modsel = vgm(0.02,"Mat",1.5, kappa = 0.25), quietmode = T)) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "RDA2")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "RDA3", cutoff = 1.2, 
                                   modsel = vgm(0.05,"Sph",0.5), quietmode = T)) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "RDA3")))

spView.grid.interval(dat = grid.value.tot, leg.name = "Loading",
            grad.value = c(-0.2,-0.1,-0.05,0.05,0.1,0.25), 
            grad.col = blues9[2:8],
            lonRange = lonRange,
            latRange = latRange, pncol = 1) + 
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", data = fortify(readShapePoly("data/bou2_4p.shp"))) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) -> plot.sp

grid.arrange(plot.load, plot.sp, ncol = 2, widths = c(5,6), heights = 10) -> plot.gather

# saving plot
ggsave(plot = plot.load, filename = "rda/rdaloading.png", dpi = 600)
ggsave(plot = plot.sp, filename = "rda/rdaspatialView.png", dpi = 600)
ggsave(plot = plot.gather, filename = "rda/gather_rdaPlot.png", dpi = 600)