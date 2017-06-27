# clean 
rm(list = ls())
source("constant.R");source("anaTls_spatialView.R");
pkgInitialization(c("dplyr","tidyr","sp","gstat","ggplot2"))
source("grid.R")

# Functions 
spView.elem <- function(elem,...) {
  spView(dat = grid.value.tot %>% 
           filter(trait == elem),
         leg.name = paste(elem,"(mg/kg)",sep = " "),
         lonRange = lonRange, latRange = latRange) +
    geom_polygon(aes(x = long, y = lat, group = group), 
                 colour = "black", fill = "grey80", data = fortify(readShapePoly("data/bou2_4p.shp"))) +
    theme(#legend.key.height = unit(5, "mm"),
      #axis.text = element_blank(),
      #axis.ticks = element_blank(),
      plot.margin = margin(0,0,0,0))
}
# Example

## for Igeo
background <- datareadln() %>% 
  gather(trait, bk, Al:orgC) %>%
  select(siteID, trait, bk) %>%
  group_by(trait) %>%
  summarise(bk = mean(bk)) %>% 
  mutate(bk2 = bk / 28066.7)

dat <- datareadln() %>% 
  gather(trait, value, Al:orgC) %>% 
  inner_join(background, by = c("trait" = "trait")) %>%
  mutate(value = log2(value / bk /1.5)) %>%
  select(siteID, trait, value) 

ggplot(data = dat %>% 
         mutate(trait = factor(trait, levels = c("Al","Fe","Mn","Pb","Cr","Ni","Cu","Zn","As","Cd")))) + 
  geom_hline(yintercept = c(0:2), col = "red", linetype = 2) +
  geom_boxplot(aes(x = trait, y = value), fill = "grey80") +
  scale_x_discrete("Element") +
  scale_y_continuous("Geo-accumulation Index",
                     breaks = c(0:2)) +
  coord_flip() +
  theme_bw() + 
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) -> plot.igeo.box

dat <- dat %>%
  group_by(siteID, trait) %>%
  summarise(value = mean(value)) %>%
  spread(trait, value) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  dplyr::select(siteID:depth,Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd) 

dat <- as.data.frame(dat)
coordinates(dat) <- ~lon+lat

grid.value.tot <- NULL
doKrig(dat, dat.grid, tag = "Fe", cutoff = 2, modsel = vgm(0.1,"Mat",1,0.02,kappa = 1), dir = "riskAssment/krig/igeo") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Fe")))

doKrig(dat, dat.grid, tag = "Mn", cutoff = 2, modsel = vgm(0.02,"Mat",1,0.04,kappa = 1), dir = "riskAssment/krig/igeo") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Mn")))

doKrig(dat, dat.grid, tag = "Pb", cutoff = 2, modsel = vgm(0.06,"Sph",1), dir = "riskAssment/krig/igeo") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Pb")))

doKrig(dat, dat.grid, tag = "Cr", cutoff = 1.7, modsel = vgm(0.1,"Mat",0.7,kappa = 1), dir = "riskAssment/krig/igeo") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cr")))

doKrig(dat, dat.grid, tag = "Ni", cutoff = 1.5, modsel = vgm(0.15,"Sph",1), dir = "riskAssment/krig/igeo") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Ni")))

doKrig(dat, dat.grid, tag = "Cu", cutoff = 1.5, modsel = vgm(0.35,"Sph",0.75), dir = "riskAssment/krig/igeo") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cu")))

doKrig(dat, dat.grid, tag = "Zn", cutoff = 1.5, modsel = vgm(0.4,"Sph",0.5), dir = "riskAssment/krig/igeo") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Zn")))

doKrig(dat, dat.grid, tag = "As", cutoff = 1.5, modsel = vgm(0.5,"Sph",0.5), dir = "riskAssment/krig/igeo") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "As")))

doKrig(dat, dat.grid, tag = "Cd", cutoff = 1.5, modsel = vgm(0.4,"Sph",0.5), dir = "riskAssment/krig/igeo") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cd")))

grid.value.tot <- grid.value.tot %>%
  mutate(element = factor(element, levels = c("Al","Fe","Mn","Pb","Cr","Ni","Cu","Zn","As","Cd")))

spView.grid(dat = grid.value.tot %>% 
              filter(trait == "Al"|trait == "Fe"|trait == "Mn"|trait == "Cu"|
                       trait == "Zn"|trait == "Pb"|trait == "Cr"|
                       trait == "Ni"|trait == "As"|trait == "Cd"),
            leg.name = "Igeo",grad.value = c(-2,-1,0,1), 
            grad.tag = c(-2,-1,0,1), lonRange = lonRange,
            latRange = latRange, pncol = 3) -> plot.igeo.sp.all

spView.grid(dat = grid.value.tot %>% 
              filter(trait == "Cu"|trait == "Zn"|trait == "As"|trait == "Cd"),
            leg.name = "Igeo",grad.value = c(-2,-1,0,1), 
            grad.tag = c(-2,-1,0,1),lonRange = lonRange,
            latRange = latRange, pncol = 2) -> plot.igeo.sp.sel

grid.arrange(plot.igeo.box, plot.igeo.sp.sel, ncol = 2, widths = c(5,10), heights = 5) -> plot.igeo.gather

## for enrichment factor
dat <- datareadln() %>% 
  gather(trait, value, Al:Cd) %>% 
  inner_join(datareadln() %>%
               select(siteID, Fe), by = c("siteID" = "siteID")) %>%
  mutate(value = value / Fe) %>%
  inner_join(background, by = c("trait" = "trait")) %>%
  mutate(value = value / bk2) %>%
  select(siteID, trait, value) 

ggplot(data = dat %>% 
         filter(trait != "Fe") %>%
         mutate(trait = factor(trait, levels = c("Al","Mn","Pb","Cr","Ni","Cu","Zn","As","Cd")))) + 
  geom_hline(yintercept = c(0.667,1,1.5,2,3,4), col = "red", linetype = 2) +
  geom_boxplot(aes(x = trait, y = value),fill = "grey80") +
  scale_x_discrete("Element") +
  scale_y_continuous("Enrichment Factor",labels = c("0.67","1.0","1.5","2.0","3.0","4.0"),
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
  dplyr::select(siteID:depth,Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd) 

dat <- as.data.frame(dat)
coordinates(dat) <- ~lon+lat

grid.value.tot <- NULL
doKrig(dat, dat.grid, tag = "Al", cutoff = 1.5, modsel = vgm(0.05,"Sph",0.75), dir = "riskAssment/krig/ef") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Al")))

doKrig(dat, dat.grid, tag = "Mn", cutoff = 1.7, modsel = vgm(0.05,"Sph",1.1), dir = "riskAssment/krig/ef") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Mn")))

doKrig(dat, dat.grid, tag = "Pb", cutoff = 1.5, modsel = vgm(0.1,"Sph",0.75), dir = "riskAssment/krig/ef") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Pb")))

doKrig(dat, dat.grid, tag = "Cr", cutoff = 1.5, modsel = vgm(0.05,"Lin",2,0.01), dir = "riskAssment/krig/ef") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cr")))

doKrig(dat, dat.grid, tag = "Ni", cutoff = 1.5, modsel = vgm(0.06,"Sph",1), dir = "riskAssment/krig/ef") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Ni")))

doKrig(dat, dat.grid, tag = "Cu", cutoff = 1.5, modsel = vgm(0.25,"Sph",0.75), dir = "riskAssment/krig/ef") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cu")))

doKrig(dat, dat.grid, tag = "Zn", cutoff = 1.5, modsel = vgm(0.2,"Sph",0.5), dir = "riskAssment/krig/ef") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Zn")))

doKrig(dat, dat.grid, tag = "As", cutoff = 1.5, modsel = vgm(0.6,"Sph",0.5), dir = "riskAssment/krig/ef") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "As")))

doKrig(dat, dat.grid, tag = "Cd", cutoff = 1.5, modsel = vgm(0.15,"Sph",0.5), dir = "riskAssment/krig/ef") -> tmp;tmp
grid.value <- as.data.frame(tmp) %>%  select(lon, lat, value = var1.pred)
grid.value.tot <- rbind(grid.value.tot, as.data.frame(cbind(grid.value, trait = "Cd")))

grid.value.tot <- grid.value.tot %>%
  mutate(element = factor(element, levels = c("Al","Fe","Mn","Pb","Cr","Ni","Cu","Zn","As","Cd")))

spView.grid(dat = grid.value.tot %>% 
              filter(trait == "Al"|trait == "Mn"|trait == "Cr"|
                       trait == "Cu"|trait == "Zn"|trait == "Pb"|
                       trait == "Ni"|trait == "As"|trait == "Cd"),
            leg.name = "EF",grad.value = c(0,0.67,1,1.5,2,2.5), 
            grad.tag = c(0,0.67,1,1.5,2,2.5), lonRange = lonRange,
            latRange = latRange, pncol = 3)  -> plot.ef.sp.all

spView.grid(dat = grid.value.tot %>% 
              filter(trait == "Cu"|trait == "Zn"|trait == "As"|trait == "Cd"),
            leg.name = "EF",grad.value = c(0,0.67,1,1.5,2,2.5), 
            grad.tag = c(0,0.67,1,1.5,2,2.5), lonRange = lonRange,
            latRange = latRange, pncol = 2)  -> plot.ef.sp.sel

grid.arrange(plot.ef.box, plot.ef.sp.sel, ncol = 2, widths = c(5,10), heights = 5) -> plot.ef.gather

# saving plot 
ggsave(plot = plot.igeo.box, filename = "riskAssment/box_igeo.png", dpi = 600)
ggsave(plot = plot.igeo.sp.all, filename = "riskAssment/map_igeo_all.png", dpi = 600)
ggsave(plot = plot.igeo.sp.sel, filename = "riskAssment/map_igeo_sel.png", dpi = 600)

ggsave(plot = plot.ef.box, filename = "riskAssment/box_Ef.png", dpi = 600)
ggsave(plot = plot.ef.sp.all, filename = "riskAssment/map_Ef_all.png", dpi = 600)
ggsave(plot = plot.ef.sp.sel, filename = "riskAssment/map_Ef_sel.png", dpi = 600)

grid.arrange(plot.igeo.box, plot.igeo.sp.sel,plot.ef.box, plot.ef.sp.sel, 
             ncol = 2, widths = c(5,10), heights = c(5,5)) -> plot.risk.gather
ggsave(plot = plot.risk.gather, filename = "riskAssment/gather_risk.png", dpi = 600)
