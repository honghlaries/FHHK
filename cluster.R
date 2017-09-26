# Initialization
rm(list = ls())
source("constant.R");source("anaTls_multivariate.R");source("anaTls_spatialView.R");
pkgInitialization(c("dplyr","tidyr","sp","gdata","gridExtra","ggplot2","maptools"))
dirInitialization("hca")
bkmap <- readShapePoly("data/bou2_4p.shp")
lonRange = c(119.2,121.8);latRange = c(33.7,35)
source("grid.R")

# Functions 

# Examples
dat <- datareadln() %>%
  group_by(siteID) %>%
  summarise(depth = mean(depth),
            Al = mean(Al),
            Fe = mean(Fe),
            Mn = mean(Mn),
            Pb = mean(Pb),
            Cr = mean(Cr),
            Ni = mean(Ni),
            Cu = mean(Cu),
            Zn = mean(Zn),
            As = mean(As),
            Cd = mean(Cd),
            orgC = mean(orgC),
            AVS = mean(AVS),
            clay = mean(clay),
            silt = mean(silt),
            sand = mean(sand)) %>%
  dplyr::select(Al:Cd,orgC:sand,depth,siteID)

cluster.site <- hcluster(dat %>% select(Pb:Cd),
                         rname = dat$siteID)

plot.ca.tree <- plot(cluster.site) 

dat1 <- data.frame(dat,class = cutree(cluster.site,3)) %>%
  gather(trait, value, Al:sand,depth) %>%
  spread(trait, value) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  dplyr::select(lon, lat, class) %>%
  mutate(group1 = (class == 1),group2 = (class == 2),group3 = (class == 3),
         group4 = (class == 4),class = factor(class))

dat1 <- as.data.frame(dat1)
coordinates(dat1) <- ~lon+lat


plot.ca.sp <- ggplot() + 
  geom_point(aes(x = lon, y = lat, col = class), size = 2, data = as.data.frame(dat1)) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", data = fortify(bkmap)) +
  coord_quickmap(xlim = lonRange, ylim = latRange) +
  theme_bw() + 
  theme(aspect.ratio = 1/2,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")

grid.value.tot <- NULL
grid.value <- as.data.frame(doKrig(dat1, dat.grid, tag = "group1", cutoff = 1.7,
                                   modsel = vgm(0.25,"Sph",1,0), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred) 
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, class = 1)))

grid.value <- as.data.frame(doKrig(dat1, dat.grid, tag = "group2", cutoff = 2,
                                   modsel = vgm(0.25,"Sph",1,0), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred) 
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, class = 2)))

grid.value <- as.data.frame(doKrig(dat1, dat.grid, tag = "group3", cutoff = 1.5,
                                   modsel = vgm(0.10,"Sph",0.5), quietmode = T)) %>%  
  select(lon, lat, value = var1.pred) 
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, class = 3)))


plot.ca.group <- ggplot() + 
  geom_raster(aes(x = lon, y = lat, fill = value),
              interpolate = T, show.legend = F, data = grid.value.tot) +
  geom_contour(aes(x = lon, y = lat, z = value), breaks = 0.5, col = "black", 
               size = 1, show.legend = F, data = grid.value.tot) +
  geom_point(aes(x = lon, y = lat), col = "black", shape = 2,
             size = 2, data = as.data.frame(dat1)) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", data = fortify(bkmap)) +
  scale_fill_gradient(low = "white", high = "grey20") +
  facet_wrap(~class,ncol = 1) + 
  coord_quickmap(xlim = lonRange, ylim = latRange) +
  theme_bw() + 
  theme(aspect.ratio = 1/2,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")

dat <- dat %>%
  group_by(trait,class) %>%
  summarise(mean = mean(value,na.rm=T), 
            se = sd(value,na.rm=T)/sqrt(sum(1-is.na(value))))

dat <-dat %>%  
  group_by() %>%
  mutate(trait = factor(trait,levels = c("sand","silt","clay","depth","Al","Fe","Mn","orgC","AVS","Pb","Cr","Ni","Cu","Zn","As","Cd")),
         class = factor(class))%>%
  arrange(trait)

cache <- rep(NA,length(row.names(dat)))
tmp <- c("sand","silt","clay","depth","Al","Fe","Mn","orgC","AVS")
for(i in 1:length(tmp)) {
  cache[dat$trait == tmp[i]] <- "Environment Factor" 
}
tmp <- c("Pb","Cr","Ni","Cu","Zn","As","Cd")
for(i in 1:length(tmp)) {
  cache[dat$trait == tmp[i]] <- "Traget Heavy Metal" 
}
dat <- data.frame(dat, taggroup = cache)

ggplot(data = dat) + 
  geom_bar(aes(x = trait, y = mean, group = class, fill = class), stat = "identity",
           position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(x = trait, ymin = mean - se, ymax = mean + se, group = class), width = 0.5 ,size = 0.7, col = "black",
                position = position_dodge(width = 0.9)) + 
  facet_wrap(~trait, scale = "free") + 
  theme_bw() + 
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right") -> plot.ca.bar

## sample layer
dat <- datareadln() %>%
  dplyr::select(Al:Cd,orgC:sand,depth,siteID)

cluster.samp <- hcluster(dat %>% select(Pb:Cd),
                         rname = dat$siteID)

plot.ca.tree <- plot(cluster.samp) 

dat1 <- data.frame(dat,class = cutree(cluster.samp,5)) %>%
  group_by(siteID) %>%
  summarise(group1 = sum(class == 1)/n(),group2 = sum(class == 2)/n(),
            group3 = sum(class == 3)/n(),group4 = sum(class == 4)/n(),
            group5 = sum(class == 5)/n()) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  select(lon,lat,group1:group5)
  
dat1 <- as.data.frame(dat1)
coordinates(dat1) <- ~lon+lat

plot.ca.sp <- ggplot() + 
  geom_point(aes(x = lon, y = lat, col = class), 
             size = 2, data = as.data.frame(dat1)) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", data = fortify(bkmap)) +
  coord_quickmap(xlim = lonRange, ylim = latRange) +
  theme_bw() + 
  theme(aspect.ratio = 1/2,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")

doKrig(dat1, dat.grid, tag = "group1", cutoff = 2, 
       modsel = vgm(0.1,"Mat",1,0.02,kappa = 1), 
       dir = "hca/samp/krig/group1") -> tmp
grid.value <- as.data.frame(tmp) %>%  
  select(lon, lat, value = var1.pred) %>%
  mutate(value = (value>=0.5)) %>%
  as.data.frame()

plot.ca.group1 <- ggplot() + 
  geom_raster(aes(x = lon, y = lat, fill = value),
              interpolate = T, show.legend = T, data = grid.value) +
  geom_point(aes(x = lon, y = lat), 
             size = 2, data = as.data.frame(dat1)) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", data = fortify(bkmap)) +
  scale_color_grey(start = 0.8, end = 0.2) +
  scale_alpha_continuous(range = c(0,0.8)) +
  coord_quickmap(xlim = lonRange, ylim = latRange) +
  theme_bw() + 
  theme(aspect.ratio = 1/2,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")

## saving plot
ggsave(plot = plot.ca.sp, filename = "hca/casp.png", dpi = 600)
ggsave(plot = plot.ca.bar, filename = "hca/cabar.png", dpi = 600)

grid.arrange(plot.ca.sp, plot.ca.bar, ncol = 2, widths = c(10,6), heights = 3) -> plot.gather
ggsave(plot = plot.gather, filename = "hca/gather_hcaPlot.png", dpi = 600)

