# clean 
rm(list = ls())
source("constant.R");source("anaTls_spatialView.R");
pkgInitialization(c("dplyr","tidyr","sp","gstat","ggplot2","directlabels"))
source("grid.R")
dirInitialization(c("riskAssment","riskAssment/krig","riskAssment/krig/ef","riskAssment/krig/igeo"))

# Functions 
spView.igeo <- function(elem,...) {
  spView(dat = grid.value.tot %>% 
           filter(trait == elem),
         leg.name = paste("Igeo(",elem,")",sep = " "),
         lonRange = lonRange, latRange = latRange) +
    geom_polygon(aes(x = long, y = lat, group = group), 
                 colour = "black", fill = "grey80", data = fortify(readShapePoly("data/bou2_4p.shp"))) +
    #scale_color_gradient(low = "black", high = "black") +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(0,0,0,0))
}

spView.ef <- function(elem,...) {
  spView(dat = grid.value.tot %>% 
           filter(trait == elem),
         leg.name = paste("EF(",elem,")",sep = " "),
         lonRange = lonRange, latRange = latRange) +
    geom_polygon(aes(x = long, y = lat, group = group), 
                 colour = "black", fill = "grey80", data = fortify(readShapePoly("data/bou2_4p.shp"))) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(0,0,0,0))
}
# Example

background <- read.csv("data/meta_baseline.csv") %>%
  mutate(bk2 = bk / 27786)

## for Igeo
dat <- datareadln() %>% 
  gather(trait, value, Fe:Cd) %>% 
  inner_join(background, by = c("trait" = "trait")) %>% 
  mutate(value = log2(value / bk /1.5)) %>%
  select(siteID, trait, value) 

ggplot(data = dat %>% filter(trait %in% c("Pb","Cr","Ni","Cu","Zn","As","Cd")) %>%
         mutate(trait = factor(trait, levels = c("Pb","Cr","Ni","Cu","Zn","As","Cd")))) + 
  geom_hline(yintercept = c(0:2), col = "red", linetype = 2) +
  geom_boxplot(aes(x = trait, y = value), fill = "grey80") +
  scale_x_discrete("Element") +
  scale_y_continuous("Geo-accumulation Index",
                     breaks = c(0:2)) +
  coord_flip() +
  theme_bw() + 
  theme(aspect.ratio = 1,
        panel.grid = element_blank())  -> plot.igeo.box

dat <- dat %>%
  group_by(siteID, trait) %>%
  summarise(value = mean(value)) %>%
  spread(trait, value) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  dplyr::select(siteID:depth,Pb,Cr,Ni,Cu,Zn,As,Cd) 

dat <- as.data.frame(dat)
coordinates(dat) <- ~lon+lat

grid.value.tot <- NULL

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Pb", cutoff = 2, 
                                   modsel = vgm(0.06,"Sph",1), quietmode = T)) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Pb")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Cr", cutoff = 1.7, 
                                   modsel = vgm(0.1,"Mat",0.7,kappa = 1), quietmode = T)) %>%
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cr")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Ni", cutoff = 1.5, 
                                   modsel = vgm(0.15,"Sph",1), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Ni")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Cu", cutoff = 1.5, 
                                   modsel = vgm(0.35,"Sph",0.75), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cu")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Zn", cutoff = 1.5, 
                                   modsel = vgm(0.4,"Sph",0.5), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Zn")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "As", cutoff = 1.5, 
                                   modsel = vgm(0.5,"Sph",0.5), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "As")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Cd", cutoff = 1.5, 
                                   modsel = vgm(0.4,"Sph",0.5), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cd")))

grid.value.tot <- grid.value.tot %>%
  mutate(trait = factor(trait, levels = c("Pb","Cr","Ni","Cu","Zn","As","Cd")))

plot.cu <- spView.igeo("Cu");ggsave(filename = "riskAssment/krig/igeo/Cu_Igeo.png")
plot.zn <- spView.igeo("Zn");ggsave(filename = "riskAssment/krig/igeo/Zn_Igeo.png")
plot.pb <- spView.igeo("Pb");ggsave(filename = "riskAssment/krig/igeo/Pb_Igeo.png")
plot.cr <- spView.igeo("Cr");ggsave(filename = "riskAssment/krig/igeo/Cr_Igeo.png")
plot.ni <- spView.igeo("Ni");ggsave(filename = "riskAssment/krig/igeo/Ni_Igeo.png")
plot.as <- spView.igeo("As");ggsave(filename = "riskAssment/krig/igeo/As_Igeo.png")
plot.cd <- spView.igeo("Cd");ggsave(filename = "riskAssment/krig/igeo/Cd_Igeo.png")

spView.grid.interval(dat = grid.value.tot, leg.name = "Igeo",
                     grad.value = c(-2,-1,0,1), 
                     grad.col = blues9[3:7],
                     lonRange = lonRange,
                     latRange = latRange, pncol = 3) + 
  geom_contour(aes(x = lon, y = lat,  z = value),col= "black",
               show.legend = F, size = 0.8, breaks = 0, linetype = 2,
               data = grid.value.tot) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", data = fortify(readShapePoly("data/bou2_4p.shp"))) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) -> plot.igeo.sp.all

spView.grid.interval(dat = grid.value.tot%>% 
                       filter(trait %in% c("Zn", "As", "Cd")), 
                     leg.name = "Igeo",
                     grad.value = c(-2,-1,0,1), 
                     grad.col = blues9[3:7],
                     lonRange = lonRange,
                     latRange = latRange, pncol = 2) + 
  geom_contour(aes(x = lon, y = lat,  z = value), col= "black",
               show.legend = F, size = 0.8, breaks = 0, linetype = 2,
               data = grid.value.tot %>% 
                 filter(trait %in% c("Zn", "As", "Cd"))) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", data = fortify(readShapePoly("data/bou2_4p.shp"))) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) -> plot.igeo.sp.sel


grid.arrange(plot.igeo.box, plot.igeo.sp.sel, ncol = 2, widths = c(5,7), heights = 5) -> plot.igeo.gather

## for enrichment factor
dat <- datareadln() %>% 
  gather(trait, value, As,Cd,Cr,Cu,Ni,Pb,Zn,Fe) %>% 
  inner_join(datareadln() %>%
               select(siteID, Fe), by = c("siteID" = "siteID")) %>%
  mutate(value = value / Fe) %>%
  inner_join(background, by = c("trait" = "trait")) %>%
  mutate(value = value / bk2) %>%
  select(siteID, trait, value) 

ggplot(data = dat %>% filter(trait %in% c("Pb","Cr","Ni","Cu","Zn","As","Cd")) %>% 
         mutate(trait = factor(trait, levels = c("Pb","Cr","Ni","Cu","Zn","As","Cd")))) + 
  geom_hline(yintercept = c(0.5,1,1.5,2,3,4), col = "red", linetype = 2) +
  geom_boxplot(aes(x = trait, y = value),fill = "grey80") +
  scale_x_discrete("Element") +
  scale_y_continuous("Enrichment Factor",labels = c("0.5","1.0","1.5","2.0","3.0","4.0"),
                     breaks = c(0.667,1,1.5,2,3,4)) +
  coord_flip() +
  theme_bw() + 
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) -> plot.ef.box

dat <- dat %>%
  group_by(siteID, trait) %>%
  summarise(value = mean(value)) %>%
  spread(trait, value) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  dplyr::select(siteID:depth,Pb,Cr,Ni,Cu,Zn,As,Cd) 

dat <- as.data.frame(dat)
coordinates(dat) <- ~lon+lat

grid.value.tot <- NULL

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Pb", cutoff = 1.5, 
                                   modsel = vgm(0.1,"Sph",0.75), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Pb")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Cr", cutoff = 1.5, 
                                   modsel = vgm(0.05,"Lin",2,0.01), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cr")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Ni", cutoff = 1.5,
                                   modsel = vgm(0.06,"Sph",1), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Ni")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Cu", cutoff = 1.5, 
                                   modsel = vgm(0.25,"Sph",0.75), quietmode = T)) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cu")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Zn", cutoff = 1.5, 
                                   modsel = vgm(0.2,"Sph",0.5), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Zn")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "As", cutoff = 1.5, 
                                   modsel = vgm(0.6,"Sph",0.5), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "As")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Cd", cutoff = 1.5, 
                                   modsel = vgm(0.15,"Sph",0.5), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cd")))

grid.value.tot <- grid.value.tot %>%
  mutate(trait = factor(trait, levels = c("Pb","Cr","Ni","Cu","Zn","As","Cd")))


plot.cu <- spView.ef("Cu");ggsave(filename = "riskAssment/krig/EF/Cu_EF.png")
plot.zn <- spView.ef("Zn");ggsave(filename = "riskAssment/krig/EF/Zn_EF.png")
plot.pb <- spView.ef("Pb");ggsave(filename = "riskAssment/krig/EF/Pb_EF.png")
plot.cr <- spView.ef("Cr");ggsave(filename = "riskAssment/krig/EF/Cr_EF.png")
plot.ni <- spView.ef("Ni");ggsave(filename = "riskAssment/krig/EF/Ni_EF.png")
plot.as <- spView.ef("As");ggsave(filename = "riskAssment/krig/EF/As_EF.png")
plot.cd <- spView.ef("Cd");ggsave(filename = "riskAssment/krig/EF/Cd_EF.png")

spView.grid.interval(dat = grid.value.tot, leg.name = "EF",
                     grad.value = c(0.5,1,1.5,2,3), 
                     grad.col = blues9[3:8],
                     lonRange = lonRange,
                     latRange = latRange, pncol = 3) + 
  geom_contour(aes(x = lon, y = lat,  z = value),col= "black",
               show.legend = F, size = 0.8, breaks = 1.5, linetype = 2,
               data = grid.value.tot) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", data = fortify(readShapePoly("data/bou2_4p.shp"))) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) -> plot.ef.sp.all

spView.grid.interval(dat = grid.value.tot%>% 
                       filter(trait %in% c("Zn", "As", "Cd")), 
                     leg.name = "EF",
                     grad.value = c(0.5,1,1.5,2,3),  
                     grad.col = blues9[3:8],
                     lonRange = lonRange,
                     latRange = latRange, pncol = 2) + 
  geom_contour(aes(x = lon, y = lat,  z = value),col= "black",
               show.legend = F, size = 0.8, breaks = 1.5, linetype = 2,
               data = grid.value.tot %>% 
                 filter(trait %in% c("Zn", "As", "Cd"))) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", data = fortify(readShapePoly("data/bou2_4p.shp"))) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) -> plot.ef.sp.sel

grid.arrange(plot.ef.box, plot.ef.sp.sel, ncol = 2, widths = c(5,7), heights = 5) -> plot.ef.gather

# saving plot 
ggsave(plot = plot.igeo.box, filename = "riskAssment/box_igeo.png", dpi = 600)
ggsave(plot = plot.igeo.sp.all, filename = "riskAssment/map_igeo_all.png", dpi = 600)
ggsave(plot = plot.igeo.sp.sel, filename = "riskAssment/map_igeo_sel.png", dpi = 600)

ggsave(plot = plot.ef.box, filename = "riskAssment/box_Ef.png", dpi = 600)
ggsave(plot = plot.ef.sp.all, filename = "riskAssment/map_Ef_all.png", dpi = 600)
ggsave(plot = plot.ef.sp.sel, filename = "riskAssment/map_Ef_sel.png", dpi = 600)

grid.arrange(plot.igeo.box, plot.igeo.sp.sel,plot.ef.box, plot.ef.sp.sel, 
             ncol = 2, widths = c(5,7), heights = c(5,5)) -> plot.risk.gather
ggsave(plot = plot.risk.gather, filename = "riskAssment/gather_risk.png", dpi = 600)
