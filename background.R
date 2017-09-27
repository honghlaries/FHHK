## Initialization
rm(list = ls())
source("constant.R");source("anaTls_spatialView.R");
pkgInitialization(c("dplyr","tidyr","sp","gstat","ggplot2"))
source("grid.R")

## Functions 
spView.trait <- function(traittag,unit,...) {
  spView(dat = grid.value.tot %>% 
           filter(trait == traittag),
         leg.name = paste(traittag,unit,sep = " "),
         lonRange = lonRange, latRange = latRange) +
    geom_polygon(aes(x = long, y = lat, group = group), 
                 colour = "black", fill = "grey80", data = fortify(readShapePoly("data/bou2_4p.shp"))) +
    geom_point(aes(x = lon, y = lat), color = "black", data = sites) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
}

dat <- datareadln() %>% 
  dplyr::select(siteID,lon,lat,salinity,pH,orgC,AVS,clay,silt,sand) %>%
  tidyr::gather(trait,value,salinity:sand) %>%
  dplyr::group_by(siteID,trait,lon,lat) %>%
  dplyr::summarise(value = mean(value)) %>%
  tidyr::spread(trait,value)

dat <- as.data.frame(dat)
coordinates(dat) <- ~lon+lat

# sand silt clay
grid.value.tot <- NULL
grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "sand", cutoff = 1.5, 
                                   modsel = vgm(800,"Sph",0.7), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "sand")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "silt", cutoff = 1.5, 
                                   modsel = vgm(800,"Sph",0.7), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "silt")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "clay", cutoff = 1.5, 
                                   modsel = vgm(5,"Sph",0.7), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "clay")))

spView.trait("silt","") + 
  scale_fill_gradientn("silt",
                       breaks = c(10,30,50,70,85), 
                       labels = paste(c(10,30,50,70,85),"%",sep = ""),
                       colors = rainbow(7)) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 6)) -> p.silt

spView.trait("sand","") + 
  scale_fill_gradientn("sand",guide = "colourbar",
                       breaks = c(10,30,50,70,85), 
                       labels = paste(c(10,30,50,70,85),"%",sep = ""),
                       colors = rainbow(7)) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 6)) -> p.sand

spView.trait("clay","") + 
  scale_fill_gradientn("clay",guide = "colourbar",
                       breaks = c(0.5,1.5,3,4.5,6,7.5), 
                       labels = paste(c(0.5,1.5,3,4.5,6,7.5),"%",sep = ""),
                       colors = rainbow(7)) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 6)) -> p.clay

ggsave(filename = "map/sand.silt.clay.png", dpi = 600, width = 4,height = 5.5,
       plot = grid.arrange(p.clay,p.silt,p.sand,ncol=1,widths = 15))

# salinity pH orgC AVS
grid.value.tot <- NULL
grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "salinity", cutoff = 1.5, 
                                   modsel = vgm(10,"Lin",0,10), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "salinity")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "pH", cutoff = 1.5, 
                                   modsel = vgm(0.15,"Lin",0), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "pH")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "orgC", cutoff = 1.2, 
                                   modsel = vgm(1000000,"Sph",0.7), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "orgC")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "AVS", cutoff = 1.2, 
                                   modsel = vgm(30000,"Sph",0.7,20000), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "AVS")))

spView.trait("salinity","") + 
  scale_fill_gradientn("salinity(ppt)",
                       breaks = 9:15, 
                       colors = rainbow(7)) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 6)) -> p.salinity

spView.trait("pH","") + 
  scale_fill_gradientn("pH",guide = "colourbar",
                       colors = rainbow(7)) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 6)) -> p.pH

spView.trait("orgC","") + 
  scale_fill_gradientn("orgC(g/kg)",guide = "colourbar",
                       breaks = 300+1:5*1000, 
                       labels = 0.3+1:5,
                       colors = rainbow(7)) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 6)) -> p.orgC

spView.trait("AVS","") + 
  scale_fill_gradientn("AVS(mg/kg)",guide = "colourbar",
                       breaks = 170+1:7*50, 
                       labels = 170+1:7*50,
                       colors = rainbow(7)) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 6)) -> p.AVS

ggsave(filename = "map/sal.ph.orgc.avs.png", dpi = 600, width = 7.5,height = 5,
       plot = grid.arrange(p.salinity,p.pH,p.orgC,p.AVS,ncol=2,widths = c(15,15)))