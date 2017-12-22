## Initialization
rm(list = ls())
source("constant.R");source("anaTls_spatialView.R");source("uniTls_csv2latex.R")
pkgInitialization(c("dplyr","tidyr","sp","gstat","ggplot2"))
source("grid.R");source("grid_resamp.R")

## Functions 
spView.elem <- function(elem,...) {
  spView(dat = grid.value.tot %>% 
           filter(trait == elem),
         leg.name = elem,
         lonRange = lonRange, latRange = latRange) +
    geom_polygon(aes(x = long, y = lat, group = group), 
                 colour = "black", fill = "grey80", 
                 data = fortify(readShapePoly("data/bou2_4p.shp"))) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
}

spView.delta <- function(elem,...) {
  spView(dat = grid.value.tot %>% 
           filter(trait == elem),
         leg.name = elem,
         lonRange = lonRange, latRange = latRange) +
    geom_contour(aes(x = lon, y = lat,  z = value),col= "black",
                 show.legend = F, size = 0.8, breaks = 0, linetype = 2,
                 data = grid.value.tot %>% filter(trait == elem)) +
    geom_polygon(aes(x = long, y = lat, group = group), 
                 colour = "black", fill = "grey80", 
                 data = fortify(readShapePoly("data/bou2_4p.shp"))) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
}


## Example
set.seed(20171216.082111)

dat <- datareadln() 

grid.perm.tot = NULL
nsamp = 99; nsite = 99

grid.perm <- as.data.frame(doKrig.resamp(dat, dat.grid.resample, 
                                          krigFormula = log(Pb) ~ 1, tag = "Pb", 
                                          cutoff = 1.5,
                                          modsel = vgm(0.35,"Sph",0.5), 
                                          nsamp = nsamp, nsite = nsite, 
                                          group = "siteID") )
grid.perm.tot <- rbind(grid.perm.tot, 
                        as.data.frame(cbind(grid.perm, trait = "Pb")))

grid.perm <- as.data.frame(doKrig.resamp(dat, dat.grid.resample, 
                                          krigFormula = log(Cr) ~ 1, tag = "Cr", 
                                          cutoff = 1.5,
                                          modsel = vgm(0.06,"Mat",0.7,kappa = 1), 
                                          nsamp = nsamp, nsite = nsite, 
                                          group = "siteID") )
grid.perm.tot <- rbind(grid.perm.tot, 
                        as.data.frame(cbind(grid.perm, trait = "Cr")))

grid.perm <- as.data.frame(doKrig.resamp(dat, dat.grid.resample, 
                                          krigFormula = log(Ni) ~ 1, tag = "Ni", 
                                          cutoff = 1.5,
                                          modsel = vgm(0.06,"Sph",1), 
                                          nsamp = nsamp, nsite = nsite, 
                                          group = "siteID") )
grid.perm.tot <- rbind(grid.perm.tot, 
                        as.data.frame(cbind(grid.perm, trait = "Ni")))

grid.perm <- as.data.frame(doKrig.resamp(dat, dat.grid.resample, 
                                          krigFormula = log(Cu) ~ 1, tag = "Cu", 
                                          cutoff = 1.5,
                                          modsel = vgm(0.25,"Sph",0.75), 
                                          nsamp = nsamp, nsite = nsite, 
                                          group = "siteID") )
grid.perm.tot <- rbind(grid.perm.tot, 
                        as.data.frame(cbind(grid.perm, trait = "Cu")))

grid.perm <- as.data.frame(doKrig.resamp(dat, dat.grid.resample, 
                                          krigFormula = log(Zn) ~ 1, tag = "Zn", 
                                          cutoff = 1.5,
                                          modsel = vgm(0.1,"Sph",0.5), 
                                          nsamp = nsamp, nsite = nsite, 
                                          group = "siteID") )
grid.perm.tot <- rbind(grid.perm.tot, 
                        as.data.frame(cbind(grid.perm, trait = "Zn")))

grid.perm <- as.data.frame(doKrig.resamp(dat, dat.grid.resample, 
                                          krigFormula = log(As) ~ 1, tag = "As", 
                                          cutoff = 1.5,
                                          modsel = vgm(0.3,"Sph",0.75), 
                                          nsamp = nsamp, nsite = nsite, 
                                          group = "siteID") )
grid.perm.tot <- rbind(grid.perm.tot, 
                        as.data.frame(cbind(grid.perm, trait = "As")))

grid.perm <- as.data.frame(doKrig.resamp(dat, dat.grid.resample, 
                                          krigFormula = log(Cd) ~ 1, tag = "Cd", 
                                          cutoff = 1.5,
                                          modsel = vgm(0.2,"Sph",0.7), 
                                          nsamp = nsamp, nsite = nsite, 
                                          group = "siteID") )
grid.perm.tot <- rbind(grid.perm.tot, 
                        as.data.frame(cbind(grid.perm, trait = "Cd")))


grid.value.tot <- NULL

seldat <- grid.perm.tot %>% 
  dplyr::filter(trait == "Pb") %>% 
  dplyr::mutate(delta = var.samp-var.site)
subdat <- as.data.frame(seldat)
coordinates(seldat) <- ~lon+lat
grid.value.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                        cutoff = 1, 
                                        modsel = vgm(0.02,"Sph",0.5), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
grid.value.mod <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.mod", 
                                       cutoff = 1.5, 
                                       modsel = vgm(0.02,"Sph",0.5), 
                                       quietmode = T)) %>% 
  dplyr::select(lon, lat, var.mod = var1.pred)
grid.value.samp <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.samp", 
                                        cutoff = 0.7, 
                                        modsel = vgm(0.04,"Sph",0.4),
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, var.samp = var1.pred)
grid.value.site <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.site", 
                                        cutoff = 1, 
                                        modsel = vgm(0.04,"Sph",0.4), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, var.site = var1.pred)
grid.value.delta <- as.data.frame(doKrig(seldat, dat.grid, tag = "delta", 
                                         cutoff = 0.8, 
                                         modsel = vgm(0.003,"Sph",0.4),
                                         quietmode = T)) %>% 
  dplyr::select(lon, lat, delta = var1.pred)
grid.value.tot <- rbind(grid.value.tot, 
                        as.data.frame(cbind(grid.value.mean %>%
                                              inner_join(grid.value.mod) %>%
                                              inner_join(grid.value.samp) %>%
                                              inner_join(grid.value.site) %>%
                                              inner_join(grid.value.delta), 
                                            trait = "Pb")))

seldat <- grid.perm.tot %>% filter(trait == "Cr") %>% mutate(delta = var.samp-var.site)
subdat <- as.data.frame(seldat)
coordinates(seldat) <- ~lon+lat
grid.value.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                        cutoff = 1.5, 
                                        modsel = vgm(0.04,"Sph",0.5), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
grid.value.mod <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.mod", 
                                       cutoff = 1.5, 
                                       modsel = vgm(0.02,"Sph",0.5), 
                                       quietmode = T)) %>% 
  dplyr::select(lon, lat, var.mod = var1.pred)
grid.value.samp <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.samp", 
                                        cutoff = 1.1, 
                                        modsel = vgm(0.0015,"Sph",1.0),
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, var.samp = var1.pred)
grid.value.site <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.site", 
                                        cutoff = 1.5, 
                                        modsel = vgm(0.0015,"Sph",1), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, var.site = var1.pred)
grid.value.delta <- as.data.frame(doKrig(seldat, dat.grid, tag = "delta", 
                                         cutoff = 1.2, 
                                         modsel = vgm(0.0012,"Sph",0.8),
                                         quietmode = T)) %>% 
  dplyr::select(lon, lat, delta = var1.pred)
grid.value.tot <- rbind(grid.value.tot, 
                        as.data.frame(cbind(grid.value.mean %>%
                                              inner_join(grid.value.mod) %>%
                                              inner_join(grid.value.samp) %>%
                                              inner_join(grid.value.site) %>%
                                              inner_join(grid.value.delta), 
                                            trait = "Cr")))

seldat <- grid.perm.tot %>% filter(trait == "Ni") %>% mutate(delta = var.samp-var.site)
subdat <- as.data.frame(seldat)
coordinates(seldat) <- ~lon+lat
grid.value.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                        cutoff = 1.5, 
                                        modsel = vgm(0.05, "Mat", 1, kappa = 2), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
grid.value.mod <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.mod", 
                                       cutoff = 1.5, 
                                       modsel = vgm(0.0003,"Lin",0), 
                                       quietmode = T)) %>% 
  dplyr::select(lon, lat, var.mod = var1.pred)
grid.value.samp <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.samp", 
                                        cutoff = 0.8, 
                                        modsel = vgm(0.0006,"Sph",0.4),
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, var.samp = var1.pred)
grid.value.site <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.site", 
                                        cutoff = 1.5, 
                                        modsel = vgm(0.0005,"Sph",1.2), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, var.site = var1.pred)
grid.value.delta <- as.data.frame(doKrig(seldat, dat.grid, tag = "delta", 
                                         cutoff = 1.7, 
                                         modsel = vgm(0.001,"Sph",1.2),
                                         quietmode = T)) %>% 
  dplyr::select(lon, lat, delta = var1.pred)
grid.value.tot <- rbind(grid.value.tot, 
                        as.data.frame(cbind(grid.value.mean %>%
                                              inner_join(grid.value.mod) %>%
                                              inner_join(grid.value.samp) %>%
                                              inner_join(grid.value.site) %>%
                                              inner_join(grid.value.delta), 
                                            trait = "Ni")))

seldat <- grid.perm.tot %>% filter(trait == "Cu") %>% mutate(delta = var.samp-var.site)
subdat <- as.data.frame(seldat)
coordinates(seldat) <- ~lon+lat
grid.value.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                        cutoff = 1, 
                                        modsel = vgm(0.15, "Mat", 0.6, kappa = 3), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
grid.value.mod <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.mod", 
                                       cutoff = 1.5, 
                                       modsel = vgm(0.006,"Lin",0), 
                                       quietmode = T)) %>% 
  dplyr::select(lon, lat, var.mod = var1.pred)
grid.value.samp <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.samp", 
                                        cutoff = 0.8, 
                                        modsel = vgm(0.006,"Mat",0.6,kappa = 3),
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, var.samp = var1.pred)
grid.value.site <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.site", 
                                        cutoff = 1, 
                                        modsel = vgm(0.0015,"Sph",0.8), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, var.site = var1.pred)
grid.value.delta <- as.data.frame(doKrig(seldat, dat.grid, tag = "delta", 
                                         cutoff = 1, 
                                         modsel = vgm(0.003,"Sph",0.4),
                                         quietmode = T)) %>% 
  dplyr::select(lon, lat, delta = var1.pred)
grid.value.tot <- rbind(grid.value.tot, 
                        as.data.frame(cbind(grid.value.mean %>%
                                              inner_join(grid.value.mod) %>%
                                              inner_join(grid.value.samp) %>%
                                              inner_join(grid.value.site) %>%
                                              inner_join(grid.value.delta), 
                                            trait = "Cu")))

seldat <- grid.perm.tot %>% filter(trait == "Zn") %>% mutate(delta = var.samp-var.site)
subdat <- as.data.frame(seldat)
coordinates(seldat) <- ~lon+lat
grid.value.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                        cutoff = 0.7, 
                                        modsel = vgm(0.15, "Mat", 0.6, kappa = 3), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
grid.value.mod <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.mod", 
                                       cutoff = 1.5, 
                                       modsel = vgm(0.0008,"Sph",1), 
                                       quietmode = T)) %>% 
  dplyr::select(lon, lat, var.mod = var1.pred)
grid.value.samp <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.samp", 
                                        cutoff = 0.6, 
                                        modsel = vgm(0.02,"Mat",0.4,kappa = 3),
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, var.samp = var1.pred)
grid.value.site <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.site", 
                                        cutoff = 0.7, 
                                        modsel = vgm(0.004,"Mat",0.5, kappa = 3), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, var.site = var1.pred)
grid.value.delta <- as.data.frame(doKrig(seldat, dat.grid, tag = "delta", 
                                         cutoff = 0.7, 
                                         modsel = vgm(0.004,"Mat",0.5, kappa = 3),
                                         quietmode = T)) %>% 
  dplyr::select(lon, lat, delta = var1.pred)
grid.value.tot <- rbind(grid.value.tot, 
                        as.data.frame(cbind(grid.value.mean %>%
                                              inner_join(grid.value.mod) %>%
                                              inner_join(grid.value.samp) %>%
                                              inner_join(grid.value.site) %>%
                                              inner_join(grid.value.delta), 
                                            trait = "Zn")))

seldat <- grid.perm.tot %>% filter(trait == "As") %>% mutate(delta = var.samp-var.site)
subdat <- as.data.frame(seldat)
coordinates(seldat) <- ~lon+lat
grid.value.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                        cutoff = 0.8, 
                                        modsel = vgm(0.25, "Mat", 0.6, kappa = 3.5), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
grid.value.mod <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.mod", 
                                       cutoff = 1.5, 
                                       modsel = vgm(0.0008,"Lin",0), 
                                       quietmode = T)) %>% 
  dplyr::select(lon, lat, var.mod = var1.pred)
grid.value.samp <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.samp", 
                                        cutoff = 0.6, 
                                        modsel = vgm(0.02,"Mat",0.4,kappa = 3),
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, var.samp = var1.pred)
grid.value.site <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.site", 
                                        cutoff = 0.7, 
                                        modsel = vgm(0.0012,"Mat",0.5, kappa = 3), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, var.site = var1.pred)
grid.value.delta <- as.data.frame(doKrig(seldat, dat.grid, tag = "delta", 
                                         cutoff = 0.7, 
                                         modsel = vgm(0.004,"Mat",0.5, kappa = 3),
                                         quietmode = T)) %>% 
  dplyr::select(lon, lat, delta = var1.pred)
grid.value.tot <- rbind(grid.value.tot, 
                        as.data.frame(cbind(grid.value.mean %>%
                                              inner_join(grid.value.mod) %>%
                                              inner_join(grid.value.samp) %>%
                                              inner_join(grid.value.site) %>%
                                              inner_join(grid.value.delta), 
                                            trait = "As")))

seldat <- grid.perm.tot %>% filter(trait == "Cd") %>% mutate(delta = var.samp-var.site)
subdat <- as.data.frame(seldat)
coordinates(seldat) <- ~lon+lat
grid.value.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                        cutoff = 0.8, 
                                        modsel = vgm(0.15, "Mat", 0.6, kappa = 3.5), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
grid.value.mod <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.mod", 
                                       cutoff = 1.5, 
                                       modsel = vgm(0.0008,"Lin",0), 
                                       quietmode = T)) %>% 
  dplyr::select(lon, lat, var.mod = var1.pred)
grid.value.samp <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.samp", 
                                        cutoff = 0.6, 
                                        modsel = vgm(0.004,"Mat",0.4,kappa = 3),
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, var.samp = var1.pred)
grid.value.site <- as.data.frame(doKrig(seldat, dat.grid, tag = "var.site", 
                                        cutoff = 0.7, 
                                        modsel = vgm(0.0015,"Mat",0.5, kappa = 3), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, var.site = var1.pred)
grid.value.delta <- as.data.frame(doKrig(seldat, dat.grid, tag = "delta", 
                                         cutoff = 0.7, 
                                         modsel = vgm(0.0015,"Mat",0.5, kappa = 3),
                                         quietmode = T)) %>% 
  dplyr::select(lon, lat, delta = var1.pred)
grid.value.tot <- rbind(grid.value.tot, 
                        as.data.frame(cbind(grid.value.mean %>%
                                              inner_join(grid.value.mod) %>%
                                              inner_join(grid.value.samp) %>%
                                              inner_join(grid.value.site) %>%
                                              inner_join(grid.value.delta), 
                                            trait = "Cd")))

grid.value.tot$value <- exp(grid.value.tot$mean)

dat.site <- datareadln() %>%
  tidyr::gather(trait,value,depth,Pb:Cd) %>%
  dplyr::group_by(siteID,trait,lon,lat) %>%
  dplyr::summarise(value = mean(value)) %>%
  tidyr::spread(trait,value)

guides.elem <- guides(fill = guide_colourbar(barwidth = 1, barheight = 6))

plot.cu <- spView.elem("Cu") + guides.elem +
  geom_point(aes(x = lon, y = lat), color = "white", size = 2,
             data = as.data.frame(dat.site[dat.site$Cu >= quantile(read.csv("data/result_element.csv")$Cu, 
                                                         probs = 0.75),]))

plot.zn <- spView.elem("Zn") + guides.elem +
  geom_point(aes(x = lon, y = lat), color = "white", size = 2,
             data = as.data.frame(dat.site[dat.site$Zn >= quantile(read.csv("data/result_element.csv")$Zn, 
                                                         probs = 0.75),]))

plot.pb <- spView.elem("Pb") + guides.elem +
  geom_point(aes(x = lon, y = lat), color = "white", size = 2,
             data = as.data.frame(dat.site[dat.site$Pb >= quantile(read.csv("data/result_element.csv")$Pb, 
                                                         probs = 0.75),]))

plot.cr <- spView.elem("Cr") + guides.elem +
  geom_point(aes(x = lon, y = lat), color = "white", size = 2,
             data = as.data.frame(dat.site[dat.site$Cr >= quantile(read.csv("data/result_element.csv")$Cr, 
                                                         probs = 0.75),]))

plot.ni <- spView.elem("Ni") + guides.elem +
  geom_point(aes(x = lon, y = lat), color = "white", size = 2,
             data = as.data.frame(dat.site[dat.site$Ni >= quantile(read.csv("data/result_element.csv")$Ni, 
                                                         probs = 0.75),]))

plot.as <- spView.elem("As") + guides.elem +
  geom_point(aes(x = lon, y = lat), color = "white", size = 2,
             data = as.data.frame(dat.site[dat.site$As >= quantile(read.csv("data/result_element.csv")$As, 
                                                         probs = 0.75),]))

plot.cd <- spView.elem("Cd") + guides.elem +
  geom_point(aes(x = lon, y = lat), color = "white", size = 2,
             data = as.data.frame(dat.site[dat.site$Cd >= quantile(read.csv("data/result_element.csv")$Cd, 
                                                         probs = 0.75),]))

p <- grid.arrange(plot.pb,plot.ni,plot.cu,plot.zn,plot.cr,plot.as,plot.cd,
                  ncol = 2, widths = c(11,11))

ggsave(filename = "element/krig/gather_krig_perm_element.png", plot = p, 
       dpi = 600, width = 8, height = 8)

grid.value.tot <- grid.value.tot %>% 
  mutate(value = 100*(exp(mean+var.mod)/exp(mean)-1))
plot.pb.mod <- spView.elem("Pb") + guides.elem
plot.cr.mod <- spView.elem("Cr") + guides.elem
plot.ni.mod <- spView.elem("Ni") + guides.elem
plot.cu.mod <- spView.elem("Cu") + guides.elem
plot.zn.mod <- spView.elem("Zn") + guides.elem
plot.as.mod <- spView.elem("As") + guides.elem
plot.cd.mod <- spView.elem("Cd") + guides.elem

p <- grid.arrange(plot.pb.mod,plot.ni.mod,plot.cu.mod,plot.zn.mod,
                  plot.cr.mod,plot.as.mod,plot.cd.mod,
                  ncol = 2, widths = c(11,11))

ggsave(filename = "element/krig/gather_krig_perm_element_mod.png", plot = p, 
       dpi = 600, width = 8, height = 8)

grid.value.tot <- grid.value.tot %>% 
  mutate(value = 100*(exp(mean+var.samp)/exp(mean)-1))
plot.pb.samp <- spView.elem("Pb") + guides.elem
plot.cr.samp <- spView.elem("Cr") + guides.elem
plot.ni.samp <- spView.elem("Ni") + guides.elem
plot.cu.samp <- spView.elem("Cu") + guides.elem
plot.zn.samp <- spView.elem("Zn") + guides.elem
plot.as.samp <- spView.elem("As") + guides.elem
plot.cd.samp <- spView.elem("Cd") + guides.elem

p <- grid.arrange(plot.pb.samp,plot.ni.samp,plot.cu.samp,plot.zn.samp,
                  plot.cr.samp,plot.as.samp,plot.cd.samp,
                  ncol = 2, widths = c(11,11))

ggsave(filename = "element/krig/gather_krig_perm_element_samp.png", plot = p, 
       dpi = 600, width = 8, height = 8)

grid.value.tot <- grid.value.tot %>% 
  mutate(value = 100*(exp(mean+var.site)/exp(mean)-1))
plot.pb.site <- spView.elem("Pb") + guides.elem
plot.cr.site <- spView.elem("Cr") + guides.elem
plot.ni.site <- spView.elem("Ni") + guides.elem
plot.cu.site <- spView.elem("Cu") + guides.elem
plot.zn.site <- spView.elem("Zn") + guides.elem
plot.as.site <- spView.elem("As") + guides.elem
plot.cd.site <- spView.elem("Cd") + guides.elem

p <- grid.arrange(plot.pb.site,plot.ni.site,plot.cu.site,plot.zn.site,
                  plot.cr.site,plot.as.site,plot.cd.site,
                  ncol = 2, widths = c(11,11))
ggsave(filename = "element/krig/gather_krig_perm_element_site.png", plot = p, 
       dpi = 600, width = 8, height = 8)

grid.value.tot <- grid.value.tot %>% 
  mutate(value = 100*((exp(mean+var.samp)/exp(mean))/
                        (exp(mean+var.site)/exp(mean))-1))
spView.grid.interval(dat = grid.value.tot, leg.name = "Vsamp vs Vsite(%)",
                     grad.value = c(-20,-5,-1,1,5,20), 
                     grad.col = c("#3A5FCD","#436EEE","#4876FF",
                                  "grey90",
                                  "#FF83FA","#EE7AE9","#CD69C9"),
                     lonRange = lonRange,
                   latRange = latRange, pncol = 3) + 
  geom_contour(aes(x = lon, y = lat,  z = value),col= "black",
               show.legend = F, size = 0.8, breaks = 0, linetype = 2,
               data = grid.value.tot) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", data = fortify(readShapePoly("data/bou2_4p.shp"))) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(1,0),
        legend.direction = "horizontal",
        legend.justification = c(1.1,-1.15)) -> p

ggsave(filename = "element/krig/gather_krig_perm_element_delta.png", plot = p, 
       dpi = 600, width = 8, height = 6)