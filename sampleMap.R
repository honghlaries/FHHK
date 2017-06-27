# clean 
rm(list = ls())
source("constant.R");source("anaTls_spatialView.R");
pkgInitialization(c("dplyr","tidyr","sp","gstat","ggplot2","directlabels"))
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
grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = depth~1, tag = "depth", cutoff = 1.2, 
                                   modsel = vgm(80,"Lin",0,0.01), dir = "map/")) %>% 
  select(lon, lat, value = var1.pred) %>%
  dplyr::mutate(var1.pred = format(var1.pred,digits = 3))
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "depth")))

# map
break.lon <- findInt.coord(lonRange)
label.lon <- conveySpCoord(break.lon,conveyFrom = "num.degree.nosuffix",conveyTo = "char.minute.suffix",type = "lon")
break.lat <- findInt.coord(latRange) 
label.lat <- conveySpCoord(break.lat,conveyFrom = "num.degree.nosuffix",conveyTo = "char.minute.suffix",type = "lat")

ggplot() + 
  geom_raster(aes(x = lon, y = lat, fill = value),
              show.legend = T,  data = grid.value.tot) +
  geom_contour(aes(x = lon, y = lat, z = value, col = ..level..), 
              show.legend = F, size = 0.8,  bins = 7, data = grid.value.tot) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", data = fortify(readShapePoly("data/bou2_4p.shp"))) +
  geom_point(aes(x = lon, y = lat), color = "black", data = sites) +
  scale_x_continuous(name = "", breaks = break.lon, labels = label.lon, expand = c(0,0)) +
  scale_y_continuous(name = "", breaks = break.lat, labels = label.lat, expand = c(0,0)) +
  scale_fill_gradient(low = blues9[5], high = blues9[9]) +
  scale_color_gradient(low = "grey90", high = "grey30") +
  coord_quickmap(xlim = lonRange, ylim = latRange) +
  theme_bw() + 
  theme(aspect.ratio = (latiRange[2]-latiRange[1])/(longiRange[2]-longiRange[1]),
        panel.grid = element_blank(),
        legend.box.margin = margin(t = 2, r = 2, b = 2, l = 2),
        legend.box.spacing = unit(1,units = "pt"),
        plot.margin = margin(0,0,0,0))
  
ggsave("map/samplemap.png",dpi = 1200)

#plot <- direct.label(plot, method="top.pieces")
plot


