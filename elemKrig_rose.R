## Initialization
rm(list = ls())
source("constant.R");source("anaTls_spatialView.R");source("anaTls_spatialCalc.R")
pkgInitialization(c("dplyr","tidyr","gstat","ggplot2","RColorBrewer","gridExtra"))
dirInitialization(c("element","element/krig"))

## Functions 

## data calculation

### sample dataset

grid.samp.tot <- datareadln() %>%
  dplyr::select(lon, lat, Pb:Cd) %>%
  tidyr::gather(trait, value, Pb:Cd)
grid.samp.tot <- grid.samp.tot %>% dplyr::select(lon, lat, value, trait)

grid.samp.tot$dist <- distCalc(coo.oldmouth$lat, coo.oldmouth$lon, 
                              grid.samp.tot$lat, grid.samp.tot$lon)

for(i in 1:length(grid.samp.tot$lon)) {
  if (!i%%1000) print(paste("calculating:",i,"/",length(grid.samp.tot$lon)))
  grid.samp.tot$azi[i] <- azimuthCalc.rough(coo.oldmouth$lat, coo.oldmouth$lon, 
                                            grid.samp.tot$lat[i], grid.samp.tot$lon[i])
}

grid.samp.tot$azideg <- floor(grid.samp.tot$azi/pi*180) 

grid.samp.tot$azideggroup <- ceiling((grid.samp.tot$azideg)/30) 
#grid.samp.tot$azideggroup[grid.samp.tot$azideg > 337.5] <- 1
grid.samp.tot$azideggroup <- factor(grid.samp.tot$azideggroup, 
                                    levels = 1:12,
                                    labels = c("NEE","NE","NEN",
                                               "NWN","NW","NWW",
                                               "SWW","SW","SWE",
                                               "SES","SE","SEE"))

### perm dataset
grid.perm.tot <- read.csv("data/result_element_perm.csv") %>%
  dplyr::filter(trait != "As")

grid.perm.tot$value <- exp(grid.perm.tot$mean) 
grid.perm.tot <- grid.perm.tot %>% dplyr::select(lon, lat, value, trait)

grid.perm.tot$dist <- distCalc(coo.oldmouth$lat, coo.oldmouth$lon, 
                                  grid.perm.tot$lat, grid.perm.tot$lon)

for(i in 1:length(grid.perm.tot$lon)) {
  if (!i%%1000) print(paste("calculating:",i,"/",length(grid.perm.tot$lon)))
  grid.perm.tot$azi[i] <- azimuthCalc.rough(coo.oldmouth$lat, coo.oldmouth$lon, 
                                      grid.perm.tot$lat[i], grid.perm.tot$lon[i])
}

grid.perm.tot$azideg <- floor(grid.perm.tot$azi/pi*180) 

grid.perm.tot$azideggroup <- ceiling((grid.perm.tot$azideg)/30) 
#grid.perm.tot$azideggroup[grid.perm.tot$azideg > 337.5] <- 1
grid.perm.tot$azideggroup <- factor(grid.perm.tot$azideggroup, 
                                    levels = 1:12,
                                    labels = c("NEE","NE","NEN",
                                               "NWN","NW","NWW",
                                               "SWW","SW","SWE",
                                               "SES","SE","SEE"))

p <- 
  ggplot() +
  geom_point(aes(x = azideg, y = value, col = dist), alpha = 0.5, 
               data = grid.perm.tot) +
  scale_x_continuous(limits = c(0,360), breaks = 45* 0:8, labels = c("E","NE","N","NW","W","SW","S","SE","E")) + 
  scale_color_gradient(low = "yellow", high = "blue") +
  #coord_polar(start = -pi/2, direction = -1) +
  facet_wrap(~trait, scale = "free") +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.position = "right")
ggsave(plot = p, filename = "element/gather_pointrose_element.png", width = 8, height = 6, dpi = 300)
rm(p)

p <- 
  ggplot() +
  geom_point(aes(x = dist, y = value),
             alpha = 0.25, data = grid.perm.tot) +
  geom_smooth(aes(x = dist, y = value), 
              method = "loess", 
              #method.args = list(family = "possion", link = "logit"),
              data = grid.perm.tot) +
  #coord_polar(start = -pi/2, direction = -1) +
  facet_grid(trait~azideggroup, scale = "free_y") +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.position = "none")
ggsave(plot = p, filename = "element/gather_point_element_dist.png", width = 12, height = 6, dpi = 300)
rm(p)
