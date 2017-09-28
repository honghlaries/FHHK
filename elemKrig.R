## Initialization
rm(list = ls())
source("constant.R");source("anaTls_spatialView.R");source("uniTls_csv2latex.R")
pkgInitialization(c("dplyr","tidyr","sp","gstat","ggplot2"))
source("grid.R")


## Functions 
spView.elem <- function(elem,...) {
  spView(dat = grid.value.tot %>% 
           filter(trait == elem),
         leg.name = elem,
         lonRange = lonRange, latRange = latRange) +
    geom_polygon(aes(x = long, y = lat, group = group), 
                 colour = "black", fill = "grey80", data = fortify(readShapePoly("data/bou2_4p.shp"))) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
}

compare_lv <- function(dat,tag) {
  ratio <- mean(unlist(dat[dat[,tag] >=
                           quantile(read.csv("data/result_element.csv")[,tag],
                                               probs = 0.75),][,tag])) /
    mean(unlist(dat[dat[,tag] < 
                      quantile(read.csv("data/result_element.csv")[,tag],
                               probs = 0.75),][,tag])) 
  data.frame(tag,ratio)
}

summary.tab <- function(dat,tag,digit=3) {
  dat <- dat[,tag] 
  min <- format(min(unlist(dat)),digit = digit)
  max <- format(max(unlist(dat)),digit = digit)
  mean <- format(mean(unlist(dat),na.rm = T),digit = digit)
  sd <- format(sd(unlist(dat),na.rm = T),digit = digit)
  paste("$",mean,"\\","pm",sd,"(",min,"-",max,")","$",sep="")
}

## Example
dat <- datareadln() %>%
  tidyr::gather(trait,value,depth,Pb:Cd) %>%
  dplyr::group_by(siteID,trait,lon,lat) %>%
  dplyr::summarise(value = mean(value)) %>%
  tidyr::spread(trait,value)

dat <- as.data.frame(dat)
coordinates(dat) <- ~lon+lat

grid.value.tot <- NULL
#grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Al)~1, tag = "Al", cutoff = 1.2, 
#                                   modsel = vgm(0.15,"Sph",0.5), quietmode = T)) %>% 
#  select(lon, lat, value = var1.pred)
#grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Al")))

#grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Fe) ~ 1, tag = "Fe", cutoff = 1.5,
#                                   modsel = vgm(0.01,"Mat",0.5,0.01,kappa = 1), quietmode = T)) %>% 
#  select(lon, lat, value = var1.pred)
#grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Fe")))
#
#grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Mn) ~ 1, tag = "Mn", cutoff = 1.5,
#                                   modsel = vgm(0.035,"Sph",0.3), quietmode = T)) %>% 
#  select(lon, lat, value = var1.pred)
#grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Mn")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Pb) ~ 1, tag = "Pb", cutoff = 1.5,
                                   modsel = vgm(0.35,"Sph",0.5), quietmode = T)) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Pb")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Cr) ~ 1, tag = "Cr", cutoff = 1.5,
                                   modsel = vgm(0.06,"Mat",0.7,kappa = 1), quietmode = T)) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cr")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Ni) ~ 1, tag = "Ni", cutoff = 1.5,
                                   modsel = vgm(0.06,"Sph",1), quietmode = T)) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Ni")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Cu) ~ 1, tag = "Cu", cutoff = 1.5,
                                   modsel = vgm(0.25,"Sph",0.75), quietmode = T)) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cu")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Zn) ~ 1, tag = "Zn", cutoff = 1.5,
                                   modsel = vgm(0.1,"Sph",0.5), quietmode = T)) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Zn")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(As) ~ 1, tag = "As", cutoff = 1.5,
                                   modsel = vgm(0.3,"Sph",0.75), quietmode = T)) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "As")))

grid.value <- as.data.frame(doKrig(dat, dat.grid, krigFormula = log(Cd) ~ 1, tag = "Cd", cutoff = 1.5,
                                   modsel = vgm(0.2,"Sph",0.7), quietmode = T)) %>% 
  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cd")))


grid.value.tot$value <- exp(grid.value.tot$value)


guides.elem <- guides(fill = guide_colourbar(barwidth = 1, barheight = 6))

plot.cu <- spView.elem("Cu") + guides.elem +
  geom_point(aes(x = lon, y = lat), color = "white", size = 2,
             data = as.data.frame(dat[dat$Cu >= quantile(read.csv("data/result_element.csv")$Cu, 
                                                         probs = 0.75),]))

plot.zn <- spView.elem("Zn") + guides.elem +
  geom_point(aes(x = lon, y = lat), color = "white", size = 2,
             data = as.data.frame(dat[dat$Zn >= quantile(read.csv("data/result_element.csv")$Zn, 
                                                         probs = 0.75),]))

plot.pb <- spView.elem("Pb") + guides.elem +
  geom_point(aes(x = lon, y = lat), color = "white", size = 2,
             data = as.data.frame(dat[dat$Pb >= quantile(read.csv("data/result_element.csv")$Pb, 
                                                         probs = 0.75),]))
  
plot.cr <- spView.elem("Cr") + guides.elem +
  geom_point(aes(x = lon, y = lat), color = "white", size = 2,
             data = as.data.frame(dat[dat$Cr >= quantile(read.csv("data/result_element.csv")$Cr, 
                                                         probs = 0.75),]))

plot.ni <- spView.elem("Ni") + guides.elem +
  geom_point(aes(x = lon, y = lat), color = "white", size = 2,
             data = as.data.frame(dat[dat$Ni >= quantile(read.csv("data/result_element.csv")$Ni, 
                                                         probs = 0.75),]))

plot.as <- spView.elem("As") + guides.elem +
  geom_point(aes(x = lon, y = lat), color = "white", size = 2,
             data = as.data.frame(dat[dat$As >= quantile(read.csv("data/result_element.csv")$As, 
                                                         probs = 0.75),]))

plot.cd <- spView.elem("Cd") + guides.elem +
  geom_point(aes(x = lon, y = lat), color = "white", size = 2,
             data = as.data.frame(dat[dat$Cd >= quantile(read.csv("data/result_element.csv")$Cd, 
                                                         probs = 0.75),]))

p <- grid.arrange(plot.pb,plot.ni,plot.cu,plot.zn,plot.cr,plot.as,plot.cd,
             ncol = 2, widths = c(11,11))

ggsave(filename = "element/krig/gather_krig_element.png", plot = p, 
       dpi = 600, width = 8, height = 8)

## estimate uneven
dat <- datareadln() %>%
  tidyr::gather(trait,value,depth,Pb:Cd) %>%
  dplyr::group_by(siteID,trait,lon,lat) %>%
  dplyr::summarise(value = mean(value)) %>%
  tidyr::spread(trait,value)

compare_lv(dat,"Pb") %>%
  rbind(compare_lv(dat,"Cr")) %>%
  rbind(compare_lv(dat,"Ni")) %>%
  rbind(compare_lv(dat,"Cu")) %>%
  rbind(compare_lv(dat,"Zn")) %>%
  rbind(compare_lv(dat,"As")) %>%
  rbind(compare_lv(dat,"Cd")) 

## latex table for basic info
108 %>% 
  cbind(summary.tab(dat,"Pb")) %>%
  cbind(summary.tab(dat,"Ni")) %>%
  cbind(summary.tab(dat,"Cu")) %>%
  cbind(summary.tab(dat,"Zn")) %>%
  cbind(summary.tab(dat,"Cr")) %>%
  cbind(summary.tab(dat,"As")) %>%
  cbind(summary.tab(dat,"Cd")) %>%
  cbind("This study") %>%
  as.data.frame()-> tab.sum

colnames(tab.sum) <- c("Sample size","Pb","Ni","Cu","Zn","Cr","As","Cd","Source") 

for(i in 1:9) {
  tab.sum[,i]<- as.character(tab.sum[,i])
}

conveyLaTex(tab.sum,"element/summary.txt")
  


  

