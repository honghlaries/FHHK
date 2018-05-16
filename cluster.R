# Initialization
rm(list = ls())
source("constant.R");source("anaTls_multivariate.R");source("anaTls_spatialView.R");
pkgInitialization(c("dplyr","tidyr","sp","gdata","gridExtra","ggplot2","maptools"))
dirInitialization("hca")
bkmap <- readShapePoly("data/bou2_4p.shp")
source("grid.R")

# Functions 

# Examples
dat <- datareadln() %>%
  group_by(siteID) %>%
  summarise(depth = mean(depth),
            Fe = mean(Fe),
            Mn = mean(Mn),
            Pb = mean(Pb),
            Cr = mean(Cr),
            Ni = mean(Ni),
            Cu = mean(Cu),
            Zn = mean(Zn),
            Cd = mean(Cd),
            orgC = mean(orgC),
            AVS = mean(AVS),
            clay = mean(clay),
            silt = mean(silt),
            sand = mean(sand)) %>%
  dplyr::select(Fe:Cd,orgC:sand,depth,siteID)

cluster.site <- hcluster(dat %>% select(Pb:Cd),
                         rname = dat$siteID)

plot.ca.tree <- plot(cluster.site) 


dat <- datareadln() %>%
  dplyr::select(Fe:Cd,orgC:sand,depth,siteID)

cluster.samp <- hcluster(dat %>% select(Pb:Cd),
                         rname = dat$siteID)

plot.ca.tree <- plot(cluster.samp) 

dat1 <- data.frame(dat,class = cutree(cluster.samp,4)) %>%
  group_by(siteID) %>%
  summarise(group1 = sum(class == 1)/n(),group2 = sum(class == 2)/n(),
            group3 = sum(class == 3)/n(),group4 = sum(class == 4)/n()) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  select(lon,lat,siteID,group1:group4)
  
dat1 <- as.data.frame(dat1)
coordinates(dat1) <- ~lon+lat

grid.value.tot <- NULL
grid.value <- as.data.frame(doKrig(dat1, dat.grid, tag = "group1", cutoff = 1.7,
                                   modsel = vgm(0.15,"Sph",1), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred) 
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, class = 1)))

grid.value <- as.data.frame(doKrig(dat1, dat.grid, tag = "group2", cutoff = 2,
                                   modsel = vgm(0.25,"Sph",1,0), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred) 
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, class = 2)))

grid.value <- as.data.frame(doKrig(dat1, dat.grid, tag = "group3", cutoff = 1,
                                   modsel = vgm(0.15,"Lin",0,0.15), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred) 
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, class = 3)))

grid.value <- as.data.frame(doKrig(dat1, dat.grid, tag = "group4", cutoff = 1.6,
                                   modsel = vgm(0.10,"Sph",0.5), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred) 
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, class = 4)))

nym <- function(a,b) {
  if(abs(a - b) <= 0.00001) a else 0
}

grid.value.tot1 <- grid.value.tot %>%
  dplyr::mutate(class = paste("c",class,sep = "")) %>%
  tidyr::spread(class,value) 
for (i in 1: length(grid.value.tot1$c1)) {
  tag <- with(grid.value.tot1[i,], max(c1,c2,c3,c4))
  grid.value.tot1[i,"c1"] <- nym(grid.value.tot1[i,"c1"], tag)
  grid.value.tot1[i,"c2"] <- nym(grid.value.tot1[i,"c2"], tag)
  grid.value.tot1[i,"c3"] <- nym(grid.value.tot1[i,"c3"], tag)
  grid.value.tot1[i,"c4"] <- nym(grid.value.tot1[i,"c4"], tag)
}
grid.value.tot1 <- grid.value.tot1 %>%
  tidyr::gather(class,value,c1:c4) %>%
  dplyr::filter(value != 0) 
  
plot.ca.indicator1 <- 
  ggplot() + 
  geom_raster(aes(x = lon, y = lat, fill = class),
              interpolate = T, show.legend = F, data = grid.value.tot1) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", data = fortify(bkmap)) +
  geom_text(aes(x = c(120.48,121.0,120.38,120.58,121.38),
                y = c(34.36,34.6,34.57,34.20,34.00),
                label = c("Class 4","Class 1","Class 3","Class 3","Class 2"),
                angle = c(0,-25,25,-5,65),
                size = c(20,20,20,20,20))) + 
  scale_fill_manual(values = blues9[1:4*2]) +
  xlab("") + ylab("") +
  coord_quickmap(xlim = lonRange, ylim = latRange) +
  theme_bw() + 
  theme(aspect.ratio = (latiRange[2]-latiRange[1])/(longiRange[2]-longiRange[1]),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

dat1 <-data.frame(dat,class = cutree(cluster.samp,4)) %>%  
  gather(trait, value, Fe:depth) %>%
  group_by(trait,class) %>%
  summarise(mean = mean(value,na.rm = T),
            se = sd(value,na.rm= T)/sqrt(n()-sum(is.na(value)))) %>%
  group_by()%>% 
  dplyr::filter(trait %in% c("silt","depth","Fe",
                             "Mn","orgC","AVS","Pb","Cr","Ni","Cu",
                             "Zn","Cd")) %>%
  mutate(trait = factor(trait,levels = c("sand","silt","clay","depth","Fe",
                                         "Mn","orgC","AVS","Pb","Cr","Ni","Cu",
                                         "Zn","Cd")),
         class = factor(class))%>%
  arrange(trait)

cache <- rep(NA,length(row.names(dat1)))
tmp <- c("sand","silt","clay","depth","Fe","Mn","orgC","AVS")
for(i in 1:length(tmp)) {
  cache[dat1$trait == tmp[i]] <- "Environment Factor" 
}
tmp <- c("Pb","Cr","Ni","Cu","Zn","Cd")
for(i in 1:length(tmp)) {
  cache[dat1$trait == tmp[i]] <- "Trageted Heavy Metal" 
}
dat1 <- data.frame(dat1, taggroup = cache)

plot.ca.bar <- 
  ggplot(data = dat1) + 
  geom_bar(aes(x = trait, y = mean, group = class, fill = class), stat = "identity",
           position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(x = trait, ymin = mean - se, ymax = mean + se, group = class),
                width = 0.5 ,size = 0.7, col = "black",
                position = position_dodge(width = 0.9)) + 
  facet_wrap(~trait, scale = "free") + 
  scale_fill_manual(values = blues9[1:4*2]) +
  xlab("") + ylab("") +
  theme_bw() + 
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right") 

plot.gather <- 
  grid.arrange(plot.ca.indicator1,plot.ca.bar,
               ncol = 1, heights = c(9.75,12))

ggsave(filename = "hca/gather_hcaPlot.png", 
       plot = plot.gather, 
       dpi = 600, width = 5.75, height = 9)
