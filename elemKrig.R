## Initialization
rm(list = ls())
source("constant.R");source("anaTls_spatialView.R");
pkgInitialization(c("dplyr","tidyr","sp","gstat","ggplot2"))
source("grid.R")

## Functions 
spView.elem <- function(elem,...) {
  spView(dat = grid.value.tot %>% 
           filter(trait == elem),
         leg.name = paste(elem,"(mg/kg)",sep = " "),
         lonRange = lonRange, latRange = latRange) +
    geom_polygon(aes(x = long, y = lat, group = group), 
                 colour = "black", fill = "grey80", data = fortify(readShapePoly("data/bou2_4p.shp"))) +
    theme(#legend.key.height = unit(5, "mm"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(0,0,0,0))
}

## Example
dat <- datareadln() %>%
  dplyr::mutate(AvsRatio = AVS / orgC) %>%
  tidyr::gather(trait,value,lon,lat,depth,Al:AvsRatio) %>%
  dplyr::group_by(siteID,trait) %>%
  dplyr::summarise(value = mean(value)) %>%
  tidyr::spread(trait,value)

dat <- as.data.frame(dat)
coordinates(dat) <- ~lon+lat

grid.value.tot <- NULL
grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Al)~1, tag = "Al", cutoff = 1.2, 
                                   modsel = vgm(0.15,"Sph",0.5), dir = "element/krig")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Al")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Fe) ~ 1, tag = "Fe", cutoff = 1.5,
                                   modsel = vgm(0.01,"Mat",0.5,0.01,kappa = 1), dir = "element/krig")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Fe")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Mn) ~ 1, tag = "Mn", cutoff = 1.5,
                                   modsel = vgm(0.035,"Sph",0.3), dir = "element/krig")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Mn")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Pb) ~ 1, tag = "Pb", cutoff = 1.5,
                                   modsel = vgm(0.35,"Sph",0.5), dir = "element/krig")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Pb")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Cr) ~ 1, tag = "Cr", cutoff = 1.5,
                                   modsel = vgm(0.06,"Mat",0.7,kappa = 1), dir = "element/krig")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cr")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Ni) ~ 1, tag = "Ni", cutoff = 1.5,
                                   modsel = vgm(0.06,"Sph",1), dir = "element/krig")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Ni")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Cu) ~ 1, tag = "Cu", cutoff = 1.5,
                                   modsel = vgm(0.25,"Sph",0.75), dir = "element/krig")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cu")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Zn) ~ 1, tag = "Zn", cutoff = 1.5,
                                   modsel = vgm(0.1,"Sph",0.5), dir = "element/krig")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Zn")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(As) ~ 1, tag = "As", cutoff = 1.5,
                                   modsel = vgm(0.3,"Sph",0.75), dir = "element/krig")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "As")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Cd) ~ 1, tag = "Cd", cutoff = 1.5,
                                   modsel = vgm(0.2,"Sph",0.7), dir = "element/krig")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cd")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(C) ~ 1, tag = "C", cutoff = 1.5,
                                   modsel = vgm(0.15,"Sph",0.7), dir = "element/krig")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "C")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(orgC) ~ 1, tag = "orgC", cutoff = 1.5,
                                   modsel = vgm(0.06,"Mat",0.7,kappa = 1), dir = "element/krig")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "orgC")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(N) ~ 1, tag = "N", cutoff = 1.5,
                                   modsel = vgm(0.25,"Sph",1), dir = "element/krig")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "N")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(S) ~ 1, tag = "S", cutoff = 1.5,
                                   modsel = vgm(0.4,"Sph",1), dir = "element/krig")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "S")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(S) ~ 1, tag = "AVS", cutoff = 1.5,
                                   modsel = vgm(0.4,"Sph",1), dir = "element/krig")) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "AVS")))

grid.value.tot$value <- exp(grid.value.tot$value)


plot.cu <- spView.elem("Cu") 
plot.zn <- spView.elem("Zn")
plot.pb <- spView.elem("Pb")
plot.cr <- spView.elem("Cr")
plot.ni <- spView.elem("Ni")
plot.as <- spView.elem("As")
plot.cd <- spView.elem("Cd")

p <- grid.arrange(plot.cu,plot.zn,plot.pb,plot.cr,plot.ni,plot.as,plot.cd,
             nrow = 4, ncol = 2, widths = c(11,11), heights = rep(7,4))

ggsave(filename = "element/krig/gather_krig_element.png", plot = p)
