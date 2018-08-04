## Initialization
rm(list = ls())
source("constant.R");source("anaTls_spatialView.R");source("uniTls_csv2latex.R")
pkgInitialization(c("dplyr","tidyr","sp","gstat","ggplot2","RColorBrewer","gridExtra"))
source("grid.R");source("grid_resamp.R")
dirInitialization(c("element","element/krig"))

## Functions 
spView.elem <- function(elem, grad.value, dat.sample) {
  
  dat.top <- as.data.frame(dat.sample[dat.sample[,elem] >=  
                                        quantile(dat.sample[,elem],
                                                 probs = 0.80),])
  colnames(dat.top)[colnames(dat.top) == elem] <- "tag"
  
  dat.bot <- as.data.frame(dat.sample[dat.sample[,elem] <
                                        quantile(dat.sample[,elem], 
                                                 probs = 0.80),])
  colnames(dat.bot)[colnames(dat.bot) == elem] <- "tag"
  
  plot.map <- 
  spView.interval(dat = grid.value.tot %>% 
                    filter(trait == elem) %>% 
                    mutate(value = exp(mean)),
                  leg.name = elem, 
                  grad.col = blues9[3:8], grad.value = grad.value,
                  lonRange = lonRange, latRange = latRange) +
    geom_polygon(aes(x = long, y = lat, group = group), 
                 colour = "black", fill = "grey80", 
                 data = fortify(readShapePoly("data/bou2_4p.shp"))) +
    geom_path(aes(x = lon, y = lat), col = "red", size = 0.8, linetype = 2,
              data = coo.1855) +
    geom_point(aes(x = lon, y = lat), color = "black", size = 1, shape = 3,
               data = dat.bot) +
    geom_point(aes(x = lon, y = lat), color = "red", size = 2, shape = 4,
               data = dat.top) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank()) 
  ggsave(filename = paste("element/map_",elem,".png",sep = ""),
         plot = plot.map, dpi = 600, width = 4, height = 2.1)
  
  plot.den <-
    ggplot() + 
    geom_density(aes(x = tag), fill = "red", alpha = 0.5, 
                 data = dat.top) +
    geom_density(aes(x = tag), fill = "black", alpha = 0.5, 
                 data = dat.bot) +
    geom_point(aes(x = mean,y = I(-0.02)), col = "red", alpha = 0.5, size = 5,
               data = dat.top %>% summarise(mean = mean(tag))) +
    geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd, x = mean, y = I(-0.02)), 
                   col = "red", alpha = 0.5, height = 0, size = 1.5,
                   data = dat.top %>% summarise(mean = mean(tag), sd = sd(tag))) +
    geom_point(aes(x = mean,y = I(-0.02)), col = "black", alpha = 0.5, size = 5,
               data = dat.bot %>% summarise(mean = mean(tag))) +
    geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd, x = mean, y = I(-0.02)), 
                   col = "black", alpha = 0.5, height = 0, size = 1.5,
                   data = dat.bot %>% summarise(mean = mean(tag), sd = sd(tag))) +
    #geom_hline(yintercept = -0.02) +
    coord_flip(ylim = c(-0.04,0.2)) + 
    theme_bw() + 
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "grey80"))
  ggsave(filename = paste("element/den_",elem,".png",sep = ""), 
         plot = plot.den, dpi = 600, width = 2, height = 3)
  
  plot.map
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

dat.sample <- datareadln()
dat.site <- datareadln() %>%
  tidyr::gather(trait,value,depth,Pb:Cd) %>%
  dplyr::group_by(siteID,trait,lon,lat) %>%
  dplyr::summarise(value = mean(value)) %>%
  tidyr::spread(trait,value)

guides.elem <- guides(fill = guide_colourbar(barwidth = 1, barheight = 6))

plot.cu <- spView.elem("Cu", c(5,10,15,20,30),dat.sample) 

plot.zn <- spView.elem("Zn", c(20,30,40,60,100),dat.sample) 

plot.pb <- spView.elem("Pb", c(25,30,35,40,50),dat.sample) 

plot.cr <- spView.elem("Cr", c(20,25,30,40,50),dat.sample) 

plot.ni <- spView.elem("Ni", c(10,15,20,25,30),dat.sample) 
  
plot.cd <- spView.elem("Cd", c(0.05,0.10,0.15,0.20,0.30),dat.sample) 

p <- grid.arrange(plot.pb,plot.ni,plot.cu,plot.zn,plot.cr,plot.cd,
                  ncol = 2, widths = c(11,11))

ggsave(filename = "element/gather_krig_perm_element.png", plot = p, 
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
  