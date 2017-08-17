# clean 
rm(list = ls())
source("constant.R");source("anaTls_spatialView.R");
pkgInitialization(c("dplyr","tidyr","sp","gstat","ggplot2","directlabels","maptools"))
source("grid.R")


# data
dat <- datareadln() %>%
  tidyr::gather(trait,value,depth) %>%
  dplyr::group_by(siteID,trait,lon,lat) %>%
  dplyr::summarise(value = mean(value)) %>%
  tidyr::spread(trait,value) 

dat <- as.data.frame(dat)
coordinates(dat) <- ~lon+lat

grid.value.tot <- NULL
grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = depth~1, 
                                   tag = "depth", cutoff = 1.2, 
                                   modsel = vgm(80,"Lin",0,0.01), 
                                   dir = "map/")) %>% 
  select(lon, lat, value = var1.pred) 
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, 
                                                            trait = "depth")))

# map
break.lon <- findInt.coord(lonRange)
label.lon <- conveySpCoord(break.lon,conveyFrom = "num.degree.nosuffix",
                           conveyTo = "char.minute.suffix",type = "lon")
break.lat <- findInt.coord(latRange) 
label.lat <- conveySpCoord(break.lat,conveyFrom = "num.degree.nosuffix",
                           conveyTo = "char.minute.suffix",type = "lat")

plot.map <- ggplot() + 
  geom_raster(aes(x = lon, y = lat, fill = value),
              show.legend = T,  data = grid.value.tot) +
  geom_contour(aes(x = lon, y = lat, z = value,  col = ..level..), 
               breaks = c(17+5*0:5,18,19), show.legend = F, size = 0.8, data = grid.value.tot) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", 
               data = fortify(readShapePoly("data/bou2_4p.shp"))) +
  #geom_path(aes(x = long, y = lat, group = group), 
  #          colour = "blue", 
  #          data = fortify(rivers)) +
  geom_point(aes(x = lon, y = lat), color = "black", data = sites) +
  geom_label(aes(x = lon, y = lat, label = siteID), vjust = -0.3, hjust = 0.7, 
            color = "black", data = sites) +
  scale_x_continuous(name = "", breaks = break.lon, labels = label.lon, 
                     expand = c(0,0)) +
  scale_y_continuous(name = "", breaks = break.lat, labels = label.lat, 
                     expand = c(0,0)) +
  scale_fill_gradient("Depth",low = blues9[5], high = blues9[8]) +
  scale_color_gradient(low = "black", high = "black") +
  coord_quickmap(xlim = lonRange, ylim = latRange) +
  theme_bw() + 
  theme(aspect.ratio = (latiRange[2]-latiRange[1])/(longiRange[2]-longiRange[1]),
        panel.grid = element_blank(),
        legend.box.margin = margin(t = 2, r = 2, b = 2, l = 2),
        legend.box.spacing = unit(1,units = "pt"),
        axis.text.y = element_text(angle = 30),
        axis.ticks = element_blank())

p <- direct.label(plot.map,list("bottom.pieces",vjust = -0.2))

ggsave("map/samplemap.png",plot.map,dpi = 600)
ggsave("map/samplemap_label.png",p,dpi = 600)



