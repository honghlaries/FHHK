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

outlierTest <- function(dat) {
  dat.sel <- dat[dat$outlier == F,]
  mod <- lm(value~orgC, data = dat.sel)
  cooksd <- cooks.distance(mod)
  dat.sel$outlier <- (cooksd > 4* mean(cooksd))
  dat <- rbind(dat.sel, dat[dat$outlier == T,])
  if (sum(dat.sel$outlier)) {outlierTest(dat)} else dat
}

# Example

## for Igeo
dat <- datareadln() %>% 
  gather(trait, value, Fe:Cd) %>% 
  inner_join(background, by = c("trait" = "trait")) %>% 
  mutate(value = log2(value / bk /1.5)) %>%
  select(siteID, trait, value) 

plot.igeo.box <- 
  ggplot(data = dat %>% filter(trait %in% c("Pb","Cr","Ni","Cu","Zn","Cd")) %>%
           mutate(trait = factor(trait, levels = c("Pb","Cr","Ni","Cu","Zn","Cd")))) + 
  geom_hline(yintercept = c(0:2), col = "red", linetype = 2) +
  geom_boxplot(aes(x = trait, y = value), fill = "grey80") +
  labs(title = "Geo-accumulation Index")+
  scale_x_discrete("") +
  scale_y_continuous("",breaks = c(0:2)) +
  coord_flip() +
  theme_bw() + 
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) 


dat <- dat %>%
  group_by(siteID, trait) %>%
  summarise(value = mean(value)) %>%
  spread(trait, value) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  dplyr::select(siteID:depth,Pb,Cr,Ni,Cu,Zn,Cd) 

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

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Cd", cutoff = 1.5, 
                                   modsel = vgm(0.4,"Sph",0.5), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cd")))

grid.value.tot <- grid.value.tot %>%
  mutate(trait = factor(trait, levels = c("Pb","Cr","Ni","Cu","Zn","Cd")))

plot.cu <- spView.igeo("Cu");ggsave(filename = "riskAssment/krig/Cu_Igeo.png")
plot.zn <- spView.igeo("Zn");ggsave(filename = "riskAssment/krig/Zn_Igeo.png")
plot.pb <- spView.igeo("Pb");ggsave(filename = "riskAssment/krig/Pb_Igeo.png")
plot.cr <- spView.igeo("Cr");ggsave(filename = "riskAssment/krig/Cr_Igeo.png")
plot.ni <- spView.igeo("Ni");ggsave(filename = "riskAssment/krig/Ni_Igeo.png")
plot.cd <- spView.igeo("Cd");ggsave(filename = "riskAssment/krig/Cd_Igeo.png")

plot.igeo.sp <- 
  spView.grid.interval(dat = grid.value.tot, leg.name = "Igeo",
                       grad.value = c(-2,-1,0,1), 
                       grad.col = blues9[3:7],
                       lonRange = lonRange,
                       latRange = latRange, pncol = 3) + 
  geom_contour(aes(x = lon, y = lat,  z = value),col= "black",
               show.legend = F, size = 0.8, breaks = 0, linetype = 1,
               data = grid.value.tot) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", data = fortify(readShapePoly("data/bou2_4p.shp"))) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) 

plot.igeo.gather <- 
  grid.arrange(plot.igeo.box, plot.igeo.sp, 
               ncol = 2, widths = c(5,10), heights = 5) 

## for enrichment factor
### original
dat <- datareadln() %>% 
  data.frame(n=1:108) %>%
  gather(trait, value,Cd,Cr,Cu,Ni,Pb,Zn) %>% 
  select(trait,value,orgC,siteID,n) %>%
  inner_join(datareadln() %>%
               data.frame(n=1:108) %>%
               select(n, Fe),  by = c("n" = "n")) %>%
  mutate(value = value / Fe) %>%
  inner_join(background, by = c("trait" = "trait")) %>%
  mutate(value = value / bk2) %>%
  select(siteID, trait, value, orgC) 

plot.ef.box <-
  ggplot(data = dat %>% filter(trait %in% c("Pb","Cr","Ni","Cu","Zn","Cd")) %>% 
         mutate(trait = factor(trait, levels = c("Pb","Cr","Ni","Cu","Zn","Cd")))) + 
  geom_hline(yintercept = c(0.5,1,1.5,2,3,4), col = "red", linetype = 2) +
  geom_boxplot(aes(x = trait, y = value),fill = "grey80") +
  labs(title = "Enrichment Factor") +
  scale_x_discrete("") +
  scale_y_continuous("",labels = c("0.5","1.0","1.5","2.0","3.0","4.0"),
                     breaks = c(0.667,1,1.5,2,3,4)) +
  coord_flip() +
  theme_bw() + 
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) 

dat <- dat %>%
  group_by(siteID, trait) %>%
  summarise(value = mean(value)) %>%
  spread(trait, value) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  dplyr::select(siteID:depth,Pb,Cr,Ni,Cu,Zn,Cd) 

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

grid.value <- as.data.frame(doKrig(dat, dat.grid, tag = "Cd", cutoff = 1.5, 
                                   modsel = vgm(0.15,"Sph",0.5), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cd")))

grid.value.tot <- grid.value.tot %>%
  mutate(trait = factor(trait, levels = c("Pb","Cr","Ni","Cu","Zn","Cd")))


plot.cu <- spView.ef("Cu");ggsave(filename = "riskAssment/krig/Cu_EF.png")
plot.zn <- spView.ef("Zn");ggsave(filename = "riskAssment/krig/Zn_EF.png")
plot.pb <- spView.ef("Pb");ggsave(filename = "riskAssment/krig/Pb_EF.png")
plot.cr <- spView.ef("Cr");ggsave(filename = "riskAssment/krig/Cr_EF.png")
plot.ni <- spView.ef("Ni");ggsave(filename = "riskAssment/krig/Ni_EF.png")
plot.cd <- spView.ef("Cd");ggsave(filename = "riskAssment/krig/Cd_EF.png")

plot.ef.sp <- 
  spView.grid.interval(dat = grid.value.tot, leg.name = "EF",
                       grad.value = c(0.66,1,1.5,2,3), 
                       grad.col = blues9[3:8],
                       lonRange = lonRange,
                       latRange = latRange, pncol = 3) + 
  geom_contour(aes(x = lon, y = lat,  z = value),col= "black",
               show.legend = F, size = 0.8, breaks = 1.5, linetype = 1,
               data = grid.value.tot) +
  geom_contour(aes(x = lon, y = lat,  z = value),col= "black",
               show.legend = F, size = 0.8, breaks = 0.66, linetype = 2,
               data = grid.value.tot) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", data = fortify(readShapePoly("data/bou2_4p.shp"))) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())

plot.ef.gather <- 
  grid.arrange(plot.ef.box, plot.ef.sp, 
               ncol = 2, widths = c(5,10), heights = 5)

### organic-carbon-adjusted Enrichment factor
dat <- datareadln() %>% 
  data.frame(n=1:108) %>%
  gather(trait, value, Cd,Cr,Cu,Ni,Pb,Zn) %>% 
  select(trait,value,orgC,siteID,n) %>%
  inner_join(datareadln() %>%
               data.frame(n=1:108) %>%
               select(n, Fe),  by = c("n" = "n")) %>%
  mutate(value = value / Fe) %>%
  inner_join(background, by = c("trait" = "trait")) %>%
  mutate(value = value / bk2) %>%
  select(siteID, trait, value, orgC) %>%
  mutate(outlier = F)

dat.tmp <- NULL
for(i in unique(dat$trait)) {dat.tmp <- rbind(dat.tmp,outlierTest(dat[dat$trait == i,]))}

plot.ef.orgc.lm <- 
  ggplot() + 
  geom_hline(yintercept = 1,col = "red", linetype = 2) + 
  geom_vline(xintercept = 3.250,col = "red", linetype = 2) + 
  geom_point(aes(x = orgC/1000, y = value, shape = outlier), 
             col = "grey50", show.legend = F, data = dat.tmp)+
  geom_smooth(aes(x = orgC/1000, y = value), 
              col = "black",
              method = "lm", se = T, 
              data = dat.tmp%>% filter(outlier == F))+ 
  scale_x_continuous("orgC (g/kg)") +
  scale_y_continuous("Enrichment Factor",
                     labels = c(".67","1.0","1.5","2.0","3.0","4.0"),
                     breaks = c(0.667,1,1.5,2,3,4)) +
  scale_shape_manual(values = c(2,3))+
  coord_cartesian(xlim = c(0.5,5.5)) +
  facet_wrap(~trait, nrow = 2) +
  theme_bw() +
  theme(aspect.ratio = 0.672, 
        panel.grid = element_blank())

dat.tmp <- dat.tmp %>% 
  dplyr::filter(outlier == F) 

regFactor <- NULL
for(i in unique(dat.tmp$trait)) {
  mod <- lm(value ~ orgC, data = dat.tmp[dat.tmp$trait == i,c("value","orgC")])
  regFactor <- rbind(regFactor, 
                     data.frame(trait = i, 
                                intercept = mod$coefficients[1], 
                                coefficients = mod$coefficients[2]))
}

dat <- dat %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>% 
  dplyr::inner_join(regFactor, by = c("trait" = "trait")) %>%
  dplyr::mutate(value_adjusted = value - intercept - coefficients * orgC,
                value_orgC = value - value_adjusted) %>%
  dplyr::select(siteID, trait, value, value_adjusted, value_orgC) 

plot.ef.orgc.resid <-
  ggplot() +
  geom_abline(slope =1, intercept = 0, linetype = 2, col = "black") +
  geom_point(aes(x =  value, y = value_adjusted), 
             col = "grey20", shape = 3, data = dat) + 
  scale_x_continuous("Enrichment Factor",limits = c(-1,4)) + 
  scale_y_continuous("Residual",limits = c(-1,4)) + 
  facet_wrap(~trait, nrow = 2) +
  theme_bw() +
  theme(aspect.ratio = 0.672,
        panel.grid = element_blank())

plot.ef.orgc.eff <-
  ggplot() +
  geom_abline(slope =1, intercept = 0, linetype = 2, col = "black") +
  geom_point(aes(x =  value, y = value_orgC), 
             col = "grey20", shape = 4, data = dat) + 
  scale_x_continuous("Enrichment Factor",limits = c(-1,4)) + 
  scale_y_continuous("Effect",limits = c(-1,4)) + 
  facet_wrap(~trait, nrow = 2) +
  theme_bw() +
  theme(aspect.ratio = 0.672,
        panel.grid = element_blank()) 

plot.ef.orgc.mod <-
  ggplot() +
  geom_point(aes(x =  value, y = value_adjusted/value_orgC), 
             col = "grey20", shape = 1, data = dat) + 
  scale_x_continuous("Enrichment Factor",limits = c(-1,4)) + 
  scale_y_continuous("Residual/Effect",limits = c(-1,3.3),
                     labels = paste(-1:3, ".0", sep = ""),
                     breaks = -1:3) + 
  facet_wrap(~trait, nrow = 2) +
  theme_bw() +
  theme(aspect.ratio = 0.672,
        panel.grid = element_blank()) 

plot.ef.orgc <-
  grid.arrange(plot.ef.orgc.lm, plot.ef.orgc.mod, 
               ncol = 1, heights = c(10,10)) 

plot.ef.orgc.effresid <-
  grid.arrange(plot.ef.orgc.eff, plot.ef.orgc.resid, 
               ncol = 1, heights = c(10,10))

# saving plot 
ggsave(plot = plot.igeo.box, filename = "riskAssment/box_igeo.png", dpi = 600)
ggsave(plot = plot.igeo.sp, filename = "riskAssment/map_igeo_all.png", dpi = 600)

ggsave(plot = plot.ef.box, filename = "riskAssment/box_Ef.png", dpi = 600)
ggsave(plot = plot.ef.sp, filename = "riskAssment/map_Ef_all.png", dpi = 600)

ggsave(plot = plot.ef.orgc, 
       filename = "riskAssment/scatter_Ef_orgC.png", 
       height = 8, width = 9, dpi = 600)
ggsave(plot = plot.ef.orgc.effresid, 
       filename = "riskAssment/scatter_Ef_EffResid.png", 
       height = 8, width = 9, dpi = 600)

plot.risk.gather <- 
  grid.arrange(plot.igeo.box, plot.igeo.sp,plot.ef.box, plot.ef.sp, 
               ncol = 2, widths = c(5,10), heights = c(5,5)) 
ggsave(plot = plot.risk.gather, 
       filename = "riskAssment/gather_risk.png",
       height = 6.72, width = 10, dpi = 600)
