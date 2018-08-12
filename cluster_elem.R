# Initialization
rm(list = ls())
source("constant.R");source("anaTls_multivariate.R");source("anaTls_spatialView.R");
pkgInitialization(c("dplyr","tidyr","sp","gdata","gridExtra","ggplot2","maptools"))
source("grid.R")

# Constants 
bulkDensity <- 0.7
area.tile <- 110.95 * (latiRange[2] - latiRange[1]) * 111.314 * cos(34.3/180*pi) * (longiRange[2] - longiRange[1])
weightunit.tile <- bulkDensity*10*bulkDensity*100000*100000
# Functions 


# Processing

## contents
grid.perm.tot <- read.csv("data/result_element_perm.csv") 

grid.content.tot <- NULL

seldat <- grid.perm.tot %>%
  dplyr::filter(trait == "Pb") %>%
  dplyr::mutate(delta = var.samp-var.site)
coordinates(seldat) <- ~lon+lat
grid.content.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean",
                                        cutoff = 1,
                                        modsel = vgm(0.02,"Sph",0.5),
                                        quietmode = T)) %>%
  dplyr::select(lon, lat, mean = var1.pred)
grid.content.tot <- rbind(grid.content.tot,
                        as.data.frame(cbind(grid.content.mean,
                                            trait = "Pb")))

seldat <- grid.perm.tot %>% filter(trait == "Cr") 
coordinates(seldat) <- ~lon+lat
grid.content.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                        cutoff = 1.5, 
                                        modsel = vgm(0.04,"Sph",0.5), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
grid.content.tot <- rbind(grid.content.tot, 
                        as.data.frame(cbind(grid.content.mean, 
                                            trait = "Cr")))

seldat <- grid.perm.tot %>% filter(trait == "Ni") 
coordinates(seldat) <- ~lon+lat
grid.content.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                        cutoff = 1.5, 
                                        modsel = vgm(0.05, "Mat", 1, kappa = 2), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
grid.content.tot <- rbind(grid.content.tot, 
                        as.data.frame(cbind(grid.content.mean,
                                            trait = "Ni")))

seldat <- grid.perm.tot %>% filter(trait == "Cu")
coordinates(seldat) <- ~lon+lat
grid.content.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                        cutoff = 1, 
                                        modsel = vgm(0.15, "Mat", 0.6, kappa = 3), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
grid.content.tot <- rbind(grid.content.tot, 
                        as.data.frame(cbind(grid.content.mean,
                                            trait = "Cu")))

seldat <- grid.perm.tot %>% filter(trait == "Zn") 
coordinates(seldat) <- ~lon+lat
grid.content.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean", 
                                        cutoff = 0.7, 
                                        modsel = vgm(0.15, "Mat", 0.6, kappa = 3), 
                                        quietmode = T)) %>% 
  dplyr::select(lon, lat, mean = var1.pred)
grid.content.tot <- rbind(grid.content.tot, 
                         as.data.frame(cbind(grid.content.mean,
                                            trait = "Zn")))

seldat <- grid.perm.tot %>% filter(trait == "Cd") 
coordinates(seldat) <- ~lon+lat
grid.content.mean <- as.data.frame(doKrig(seldat, dat.grid, tag = "mean",
                                        cutoff = 0.8,
                                        modsel = vgm(0.15, "Mat", 0.6, kappa = 3.5),
                                        quietmode = T)) %>%
  dplyr::select(lon, lat, mean = var1.pred)
grid.content.tot <- rbind(grid.content.tot, 
                          as.data.frame(cbind(grid.content.mean,
                                            trait = "Cd")))

grid.content.tot$value <- exp(grid.content.tot$mean) 
grid.content.tot <- grid.content.tot %>% dplyr::select(-mean)
remove(grid.content.mean, grid.perm.tot, seldat)

## cluster
dat <- datareadln() %>%
  dplyr::select(Fe:Cd,orgC:sand,depth,siteID)

cluster.samp <- hcluster(dat %>% select(Pb:Cd),
                         rname = dat$siteID)

dat1 <- data.frame(dat,class = cutree(cluster.samp,4)) %>%
  group_by(siteID) %>%
  summarise(group1 = sum(class == 1)/n(),group2 = sum(class == 2)/n(),
            group3 = sum(class == 3)/n(),group4 = sum(class == 4)/n()) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  select(lon,lat,siteID,group1:group4)

dat1 <- as.data.frame(dat1)
coordinates(dat1) <- ~lon+lat

grid.cluster.tot <- NULL
grid.cluster <- as.data.frame(doKrig(dat1, dat.grid, tag = "group1", cutoff = 1.7,
                                   modsel = vgm(0.15,"Sph",1), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred) 
grid.cluster.tot <- rbind(grid.cluster.tot, as.data.frame(cbind(grid.cluster, class = 1)))

grid.cluster <- as.data.frame(doKrig(dat1, dat.grid, tag = "group2", cutoff = 2,
                                   modsel = vgm(0.25,"Sph",1,0), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred) 
grid.cluster.tot <- rbind(grid.cluster.tot, as.data.frame(cbind(grid.cluster, class = 2)))

grid.cluster <- as.data.frame(doKrig(dat1, dat.grid, tag = "group3", cutoff = 1,
                                   modsel = vgm(0.15,"Lin",0,0.15), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred) 
grid.cluster.tot <- rbind(grid.cluster.tot, as.data.frame(cbind(grid.cluster, class = 3)))

grid.cluster <- as.data.frame(doKrig(dat1, dat.grid, tag = "group4", cutoff = 1.6,
                                   modsel = vgm(0.10,"Sph",0.5), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred) 
grid.cluster.tot <- rbind(grid.cluster.tot, as.data.frame(cbind(grid.cluster, class = 4)))

nym <- function(a,b) {
  if(abs(a - b) <= 0.00001) a else 0
}

grid.cluster.tot1 <- grid.cluster.tot %>%
  dplyr::mutate(class = paste("c",class,sep = "")) %>%
  tidyr::spread(class,value) 
for (i in 1: length(grid.cluster.tot1$c1)) {
  tag <- with(grid.cluster.tot1[i,], max(c1,c2,c3,c4))
  grid.cluster.tot1[i,"c1"] <- nym(grid.cluster.tot1[i,"c1"], tag)
  grid.cluster.tot1[i,"c2"] <- nym(grid.cluster.tot1[i,"c2"], tag)
  grid.cluster.tot1[i,"c3"] <- nym(grid.cluster.tot1[i,"c3"], tag)
  grid.cluster.tot1[i,"c4"] <- nym(grid.cluster.tot1[i,"c4"], tag)
}
grid.cluster.tot1 <- grid.cluster.tot1 %>% 
  tidyr::gather(class,value,c1:c4) %>%
  dplyr::filter(value != 0) %>%
  dplyr::select(-value)
remove(grid.cluster, grid.cluster.tot , nym, dat, dat1, cluster.samp)

## gather
grid.tot <- inner_join(grid.cluster.tot1, grid.content.tot)
remove(grid.cluster.tot1, grid.content.tot)

## calculation
nymlow <- function(dat) {
  mod <- t.test(dat) 
  mod$conf.int[1]
}

nymhigh <- function(dat) {
  mod <- t.test(dat) 
  mod$conf.int[2]
}

bk <- grid.tot %>% dplyr::filter(class == "c1") %>%
  dplyr::group_by(trait) %>%
  dplyr::summarise(bk = mean(value))

cont <- grid.tot %>% dplyr::inner_join(bk) %>% 
  dplyr::mutate(con = value - bk)

bkint <- cont %>%
  dplyr::group_by(trait) %>%
  dplyr::summarise(intlow = nymlow(con),
                   inthigh = nymhigh(con))


stock <- cont %>% 
  dplyr::inner_join(bkint) %>%
  dplyr::mutate(stock = con * weightunit.tile / 1E6,
                stocklow = intlow * weightunit.tile / 1E6,
                stockhigh = inthigh * weightunit.tile / 1E6) %>%
  dplyr::group_by(trait,class) %>%
  dplyr::summarise(stock = sum(stock), 
                   #stocklow = sum(stocklow),
                   #stockhigh = sum(stockhigh),
                   conmean = mean(con)) %>%
  tidyr::gather(type,value,stock:conmean) %>%
  tidyr::spread(trait,value)


