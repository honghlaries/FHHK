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
  
plot.ca.indicator <- 
ggplot() + 
  geom_raster(aes(x = lon, y = lat, fill = class),
              interpolate = T, show.legend = F, data = grid.value.tot1) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", data = fortify(bkmap)) +
  geom_path(aes(x = lon, y = lat), col = "red", size = 0.8, linetype = 2,
            data = coo.1855) +
  geom_text(aes(x = c(120.6,121.0,120.38,120.58,121.38),
                y = c(34.36,34.6,34.57,34.20,34.00),
                label = c("Class 4","Class 1","Class 3","Class 3","Class 2"),
                angle = c(0,-25,25,-5,65),
                size = c(20,20,20,20,20))) + 
  scale_fill_manual(values = blues9[1:4*2]) +
  xlab("") + ylab("") +
  coord_quickmap(xlim = lonRange + c(0.08,-0.08), ylim = latRange+ c(0.05,-0.05)) +
  theme_bw() + 
  theme(aspect.ratio = (latRange[2]-latRange[1]+0.16)/(lonRange[2]-lonRange[1]+0.10),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

ggsave(filename = "hca/map_hcaPlot.png", plot = plot.ca.indicator, 
       dpi = 600, height = 5, width = 7.3)

area.tile <- 110.95 * (latiRange[2] - latiRange[1]) * 111.314 * cos(34.3/180*pi) * (longiRange[2] - longiRange[1])
area <- data.frame(class = paste("c",1:4,sep = ""),area = NA)
for (i in 1:4) area[i,"area"] <- sum(grid.value.tot1$class == area[i,"class"]) * area.tile 
area$area <- round(area$area/100,1)*100
area <- rbind(area,data.frame(class = c("LON",
                                        "NY",
                                        "SIN",
                                        "HK",
                                        "SHH"),
                              area = c(1737.9,
                                       1213.37,
                                       721.5,
                                       2755,
                                       4000)))
area$class <- as.character(area$class)
area$class[1:4] <- paste("C",1:4, sep = "")
area$class <- factor(area$class, levels = c("SHH",
                                            "HK",
                                            "SIN",
                                            "NY",
                                            "LON",
                                            paste("C",4:1, sep = "")))
area$area <- area$area / 1000

plot.ca.area <- 
ggplot(data = area) +
  geom_bar(aes(x = class, y = area, fill = class), 
           stat = "identity", position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c(terrain.colors(12)[11:7],blues9[4:1*2])) +
  ylab(expression(Area~(10^3~km^2))) +
  coord_flip() +
  theme_bw() +
  theme(plot.background = element_rect(fill = "grey80"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

ggsave(filename = "hca/bar_hcaArea.png", plot = plot.ca.area, 
       dpi = 600, height = 3, width = 2)


dat1 <-data.frame(dat,class = cutree(cluster.samp,4)) %>%  
  gather(trait, value, Fe:depth) %>%
  group_by(trait,class) %>%
  summarise(mean = mean(value,na.rm = T),
            se = sd(value,na.rm= T)/sqrt(n()-sum(is.na(value)))) %>%
  group_by()%>% 
  dplyr::filter(trait %in% c("Pb","Cr","Ni","Cu","Zn","Cd")) %>%
  mutate(class = factor(class))%>%
  arrange(trait)

plot.ca.elem <- 
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
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") 
ggsave(filename = "hca/bar_hcaCont.png", plot = plot.ca.elem, 
       dpi = 600, height = 4, width = 3)

dat1 <-data.frame(dat,class = cutree(cluster.samp,4)) %>%  
  gather(trait, value, Fe:depth) %>%
  group_by(trait,class) %>%
  summarise(mean = mean(value,na.rm = T),
            se = sd(value,na.rm= T)/sqrt(n()-sum(is.na(value)))) %>%
  group_by()%>% 
  dplyr::filter(trait %in% c("silt","depth","Fe","Mn","orgC","AVS")) %>%
  mutate(class = factor(class))%>%
  arrange(trait)

plot.ca.env <- 
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
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") 
ggsave(filename = "hca/bar_hcaEnv.png", plot = plot.ca.env, 
       dpi = 600, height = 4, width = 3)

plot.elem_env <- 
grid.arrange(plot.ca.elem,plot.ca.env,
             ncol = 2, heights = 4)

ggsave(filename = "hca/gather_hcaElemEnv.png", 
       plot = plot.elem_env, 
       dpi = 600, width = 6, height = 4)
