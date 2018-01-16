## Initialization
rm(list = ls())
source("constant.R");source("anaTls_spatialView.R");source("uniTls_csv2latex.R")
pkgInitialization(c("dplyr","tidyr","sp","gstat","ggplot2","RColorBrewer","gridExtra"))
source("grid.R");source("grid_resamp.R")
dirInitialization(c("element","element/krig"))

## Functions 
spView.elem <- function(elem,...) {
  spView(dat = grid.value.tot %>% 
           filter(trait == elem) %>% 
           mutate(value = exp(mean)),
         leg.name = elem,
         lonRange = lonRange, latRange = latRange) +
    geom_polygon(aes(x = long, y = lat, group = group), 
                 colour = "black", fill = "grey80", 
                 data = fortify(readShapePoly("data/bou2_4p.shp"))) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
}

spView.delta <- function(elem,...) {
  spView.interval(dat = dplyr::filter(grid.value.tot, trait %in% elem) %>%
                    dplyr::mutate(value = var.samp/var.site), 
                       leg.name = "", grad.value = c(0.1, 0.3, 0.5, 0.7, 0.9,
                                                     1.1, 1.3, 1.5, 2, 2.5), 
                       grad.col = RColorBrewer::brewer.pal(11,"RdBu")[11:1],
                       lonRange = lonRange,
                       latRange = latRange) + 
    geom_polygon(aes(x = long, y = lat, group = group), 
                 colour = "black", fill = "grey80", 
                 data = fortify(readShapePoly("data/bou2_4p.shp"))) +
    theme(legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

view.delta <- function(elem,range,...) {
  dat <- dplyr::filter(grid.perm.tot, trait %in% elem) %>%
    dplyr::mutate(group = factor((var.samp > 2.5*var.site) + 
                                   (var.samp > 2*var.site) + 
                                   (var.samp > 1.5*var.site) + 
                                   (var.samp > 1.3*var.site) + 
                                   (var.samp > 1.1*var.site) + 
                                   (var.samp > 0.9*var.site) + 
                                   (var.samp > 0.7*var.site) + 
                                   (var.samp > 0.5*var.site) + 
                                   (var.samp > 0.3*var.site) + 
                                   (var.samp > 0.1*var.site)))
  ggplot()+
    geom_abline(slope = 1, intercept = 0, linetype = 2)+
    geom_point(aes(x = abs(var.site/mean), y = abs(var.samp/mean), col =  group), 
               alpha  = 1, data = dat)+
    geom_text(aes(x = 0.05*range, y = range, label = elem), size = 10)+ 
    scale_color_manual(breaks = c(0.1, 0.3, 0.5, 0.7, 0.9,
                                  1.1, 1.3, 1.5, 2, 2.5),
                       values = RColorBrewer::brewer.pal(11,"RdBu")[11:1])+
    scale_x_continuous("Variation from missing sites",breaks = 0.05*0:8, labels = paste(5*0:8,"%",sep =""))+
    scale_y_continuous("Variation from deficient subsample",breaks = 0.05*0:8, labels = paste(5*0:8,"%",sep =""))+
    coord_equal(xlim = c(0,range), ylim = c(0,range))+
    theme_bw()+
    theme(legend.position = "none",
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 16))
}


## submaple v.s. site
grid.perm.tot <- read.csv("data/result_element_perm.csv") 
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
grid.value.tot <- rbind(grid.value.tot,
                        as.data.frame(cbind(grid.value.mean %>%
                                              inner_join(grid.value.samp) %>%
                                              inner_join(grid.value.site),
                                            trait = "Pb")))

seldat <- grid.perm.tot %>% filter(trait == "Cr") %>% mutate(delta = var.samp-var.site)
subdat <- as.data.frame(seldat)
coordinates(seldat) <- ~lon+lat
grid.value.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                        cutoff = 1.5, 
                                        modsel = vgm(0.04,"Sph",0.5), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
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
grid.value.tot <- rbind(grid.value.tot, 
                        as.data.frame(cbind(grid.value.mean %>%
                                              inner_join(grid.value.samp) %>%
                                              inner_join(grid.value.site), 
                                            trait = "Cr")))

seldat <- grid.perm.tot %>% filter(trait == "Ni") %>% mutate(delta = var.samp-var.site)
subdat <- as.data.frame(seldat)
coordinates(seldat) <- ~lon+lat
grid.value.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                        cutoff = 1.5, 
                                        modsel = vgm(0.05, "Mat", 1, kappa = 2), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
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
grid.value.tot <- rbind(grid.value.tot, 
                        as.data.frame(cbind(grid.value.mean %>%
                                              inner_join(grid.value.samp) %>%
                                              inner_join(grid.value.site),
                                            trait = "Ni")))

seldat <- grid.perm.tot %>% filter(trait == "Cu") %>% mutate(delta = var.samp-var.site)
subdat <- as.data.frame(seldat)
coordinates(seldat) <- ~lon+lat
grid.value.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                        cutoff = 1, 
                                        modsel = vgm(0.15, "Mat", 0.6, kappa = 3), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
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
grid.value.tot <- rbind(grid.value.tot, 
                        as.data.frame(cbind(grid.value.mean %>%
                                              inner_join(grid.value.samp) %>%
                                              inner_join(grid.value.site),
                                            trait = "Cu")))

seldat <- grid.perm.tot %>% filter(trait == "Zn") %>% mutate(delta = var.samp-var.site)
subdat <- as.data.frame(seldat)
coordinates(seldat) <- ~lon+lat
grid.value.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                        cutoff = 0.7, 
                                        modsel = vgm(0.15, "Mat", 0.6, kappa = 3), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
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
grid.value.tot <- rbind(grid.value.tot, 
                        as.data.frame(cbind(grid.value.mean %>%
                                              inner_join(grid.value.samp) %>%
                                              inner_join(grid.value.site), 
                                            trait = "Zn")))

seldat <- grid.perm.tot %>% filter(trait == "As") %>% mutate(delta = var.samp-var.site)
subdat <- as.data.frame(seldat)
coordinates(seldat) <- ~lon+lat
grid.value.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                        cutoff = 0.8, 
                                        modsel = vgm(0.25, "Mat", 0.6, kappa = 3.5), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
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
grid.value.tot <- rbind(grid.value.tot, 
                        as.data.frame(cbind(grid.value.mean %>%
                                              inner_join(grid.value.samp) %>%
                                              inner_join(grid.value.site),
                                            trait = "As")))

seldat <- grid.perm.tot %>% filter(trait == "Cd") %>% mutate(delta = var.samp-var.site)
subdat <- as.data.frame(seldat)
coordinates(seldat) <- ~lon+lat
grid.value.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean",
                                        cutoff = 0.8,
                                        modsel = vgm(0.15, "Mat", 0.6, kappa = 3.5),
                                        quietmode = T)) %>%
  dplyr::select(lon, lat, mean = var1.pred)
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
grid.value.tot <- rbind(grid.value.tot,
                        as.data.frame(cbind(grid.value.mean %>%
                                               inner_join(grid.value.samp) %>%
                                               inner_join(grid.value.site),
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


nym <- function(elem,range) {
  p <- view.delta(elem,range)
  ggsave(plot = p, height = 6, width = 6,
         filename = paste("element/krig/",elem,"_delta.png", sep = ""))
  p <- spView.delta(elem)
  ggsave(plot = p, height = 5, width = 7.5,
         filename = paste("element/krig/",elem,"_deltaSp.png", sep = ""))
}
  
nym("Pb",0.15)
nym("Cr",0.1)
nym("Ni",0.1)
nym("Cu",0.2)
nym("Zn",0.1)
nym("As",0.4)
nym("Cd",0.3)
# view.delta("",0.5)

## 
rm(list = ls())
grid.perm.tot <- read.csv("data/result_element_perm_siteImp.csv") 

ggplot() +
  # geom_path(aes(x = factor(sitelv), y = var, group = 1), size = 1.5,
  #           data = grid.perm.tot%>%
  #             dplyr::select(var.site1:var.site5,trait)%>%
  #             tidyr::gather(sitelv,var,var.site1:var.site5)%>%
  #             dplyr::group_by(sitelv,trait)%>%
  #             dplyr::summarise(var = mean(var)))+
  geom_hline(aes(yintercept = var,col = sitelv), linetype = 2, size =1,
             data = grid.perm.tot%>%
               dplyr::select(var.site1:var.site5,trait)%>%
               tidyr::gather(sitelv,var,var.site1:var.site5)%>%
               dplyr::group_by(sitelv,trait)%>%
               dplyr::summarise(var = mean(var)))+
  geom_hline(aes(yintercept = var), linetype = 1, size =1.5, col = "red",
             data = grid.perm.tot%>%
               dplyr::select(var.samp,trait)%>%
               dplyr::group_by(trait)%>%
               dplyr::summarise(var = mean(var.samp)))+
  geom_boxplot(aes(x = sitelv, y = var, col = sitelv),
               data = grid.perm.tot%>%
                 dplyr::select(var.samp:var.site5,trait)%>%
                 tidyr::gather(sitelv,var,var.samp:var.site5))+
  scale_y_continuous(limits = c(0,0.5), breaks = 0.1*0:5, 
                     labels = paste(10*0:5,"%",sep = ""))+
  scale_color_manual(values = c("red","grey50","grey40","grey30","grey20","black"))+
  facet_wrap(~trait, ncol = 4)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
ggsave("element/krig/Imp.png")
  