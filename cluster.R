## Initialization
rm(list = ls())
source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R");
source("anaTls_spatialView.R");source("anaTls_multivariate.R");
pkgInitialization(c("dplyr","tidyr","sp","gdata","gridExtra","ggplot2","maptools"))
#source("grid.R")

## Functions 
datareadln <- function() { 
  pkgLoad("dplyr");pkgLoad("tidyr")
  read.csv("data/result_element.csv") %>%
    dplyr::inner_join(read.csv("data/result_grainSize.csv"), by = c("sampleID" = "sampleID")) %>%
    dplyr::right_join(read.csv("data/meta_splList.csv"), by = c("sampleID" = "sampleID")) %>%
    dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
    dplyr::mutate(orgC = C.ac, AVS = S - S.ac, 
                  isComplete = complete.cases(Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd,C,N,S,orgC,AVS,clay,silt,sand)) %>%
    dplyr::select(siteID:depth,isComplete,Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd,C,N,S,orgC,AVS,clay,silt,sand)
}

## Examples
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

dat <- data.frame(dat,class = cutree(cluster.site,5)) %>%
  gather(trait, value, Al:sand,depth) 

dat1 <- dat %>%
  spread(trait, value) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  dplyr::select(lon, lat, class) %>%
  mutate(class = factor(class))

dat1 <- as.data.frame(dat1)

bkmap <- readShapePoly("data/bou2_4p.shp")
lonRange = c(119.2,121.8);latRange = c(33.7,35)

ggplot() + 
  geom_point(aes(x = lon, y = lat, col = class), size = 2, data = dat1) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80", data = fortify(bkmap)) +
  coord_quickmap(xlim = lonRange, ylim = latRange) +
  theme_bw() + 
  theme(aspect.ratio = 1/2,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> plot.ca.sp


dat <- dat %>%
  group_by(trait,class) %>%
  summarise(mean = mean(value,na.rm=T), se = sd(value,na.rm=T)/sqrt(sum(1-is.na(value))))

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

## saving plot
ggsave(plot = plot.ca.sp, filename = "hca/casp.png", dpi = 600)
ggsave(plot = plot.ca.bar, filename = "hca/cabar.png", dpi = 600)

grid.arrange(plot.ca.sp, plot.ca.bar, ncol = 2, widths = c(10,6), heights = 3) -> plot.gather
ggsave(plot = plot.gather, filename = "hca/gather_rdaPlot.png", dpi = 600)



  


  
 



