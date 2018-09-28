## Initialization
rm(list = ls())
source("constant.R");source("anaTls_spatialView.R")
pkgInitialization(c("dplyr","tidyr","gstat","ggplot2","RColorBrewer","gridExtra"))
dirInitialization(c("element","element/krig"))

## Functions 
distCalc <- function(lata, lona, latb, lonb, r = 6371){
  r*acos(cos(lata/180*pi)*cos(lona/180*pi)*cos(latb/180*pi)*cos(lonb/180*pi)+
           cos(lata/180*pi)*sin(lona/180*pi)*cos(latb/180*pi)*sin(lonb/180*pi)+
           sin(lata/180*pi)*sin(latb/180*pi))
}

azimuthCalc <- function(lata, lona, latb, lonb, latc =90, lonc = lonb){
  p <- (acos(cos(lata/180*pi)*cos(lona/180*pi)*cos(latb/180*pi)*cos(lonb/180*pi)+
               cos(lata/180*pi)*sin(lona/180*pi)*cos(latb/180*pi)*sin(lonb/180*pi)+
               sin(lata/180*pi)*sin(latb/180*pi))+
          (90-lata)/180*pi +
          (90-latb)/180*pi )/2.0
  cos.azi <- sqrt(sin(p)*sin(p-(90-latb)/180*pi)/sin((90-lata)/180*pi)/
                    sin(acos(cos(lata/180*pi)*cos(lona/180*pi)*cos(latb/180*pi)*cos(lonb/180*pi)+
                               cos(lata/180*pi)*sin(lona/180*pi)*cos(latb/180*pi)*sin(lonb/180*pi)+
                               sin(lata/180*pi)*sin(latb/180*pi))))
   
  if(lonb >= lona) {
    if(latb >= lata) {0.5*pi - 2*acos(cos.azi)} else {2.5*pi - 2*acos(cos.azi)}
  } else {
    if(latb >= lata) {0.5*pi - 2*acos(cos.azi)} else {2.5*pi - 2*acos(cos.azi)}
  }
}

azimuthCalc.rough <- azimuthCalc <- function(lata, lona, latb, lonb){
  tan.azi <- (latb-lata)/(lonb-lona)
  
  if(lonb >= lona) {
    if(latb >= lata) {atan(tan.azi)} else {2*pi + atan(tan.azi)}
  } else {
    pi + atan(tan.azi)
  }
}


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
