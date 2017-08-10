## Initialization
rm(list = ls())
source("constant.R");
pkgInitialization(c("dplyr","tidyr","ggplot2","gridExtra"))

## Functions 
relationPlot <- function(dat, fact, resp, col) {
  ggplot() +
    geom_point(aes_string(x = fact, y = resp), col = "black", data = dat) +
    geom_smooth(aes_string(x = fact, y = resp), method = "loess", col = col, fill = col, data = dat) +
    theme_bw()
}

relationPlot.gather <- function(dat,tag) {
  plot.dist <- relationPlot(dat, "distance", tag, "green")
  plot.clay <- relationPlot(dat, "clay", tag, "black")
  plot.dep <- relationPlot(dat, "depth", tag, "blue")
  plot.al <- relationPlot(dat, "Al", tag, "purple")
  plot.h <- relationPlot(dat, "proton", tag, "orange")
  plot.sal <- relationPlot(dat, "salinity", tag, "brown")
  p <- grid.arrange(plot.dist,plot.dep,plot.clay,plot.al,plot.h,plot.sal, nrow=2, ncol=3)
  p
}

stepFitting <- function(dat,tag) {
  null <- lm(as.formula(paste(tag,"~1",sep = "")), data = dat) 
  full <- lm(as.formula(paste(tag,"~distance*depth*clay*Al*proton*salinity",sep = "")), data = dat) 
  step(null, scope = formula(full), test = "F")
}

## Example
dat <- datareadln() %>%
  dplyr::select(depth,distance,salinity:sand) %>% 
  dplyr::mutate(proton = 10^(-pH))

taglist <- c("Cr","As","Ni","Cu","Pb","Zn","Cd")

for(i in 1:length(taglist)) {
  relationPlot.gather(dat,taglist[i]) -> p
  ggsave(plot = p,
         filename = paste(dirPreset("relation/driver"),"/",taglist[i],".png",sep = ""),dpi = 600)
  p
}

plot.dist.Cr <- relationPlot(dat, "distance", "Cr", "green")
plot.clay.Cr <- relationPlot(dat, "clay", "Cr", "black")
plot.dep.Cr <- relationPlot(dat, "depth", "Cr", "blue")
plot.al.Cr <- relationPlot(dat, "Al", "Cr", "purple")
plot.h.Cr <- relationPlot(dat, "proton", "Cr", "orange")
plot.sal.Cr <- relationPlot(dat, "salinity", "Cr", "brown")

plot.dist.As <- relationPlot(dat, "distance", "As", "green")
plot.clay.As <- relationPlot(dat, "clay", "As", "black")
plot.dep.As <- relationPlot(dat, "depth", "As", "blue")
plot.al.As <- relationPlot(dat, "Al", "As", "purple")
plot.h.As <- relationPlot(dat, "proton", "As", "orange")
plot.sal.As <- relationPlot(dat, "salinity", "As", "brown")

plot.dist.Ni <- relationPlot(dat, "distance", "Ni", "green")
plot.clay.Ni <- relationPlot(dat, "clay", "Ni", "black")
plot.dep.Ni <- relationPlot(dat, "depth", "Ni", "blue")
plot.al.Ni <- relationPlot(dat, "Al", "Ni", "purple")
plot.h.Ni <- relationPlot(dat, "proton", "Ni", "orange")
plot.sal.Ni <- relationPlot(dat, "salinity", "Ni", "brown")

plot.dist.Cu <- relationPlot(dat, "distance", "Cu", "green")
plot.clay.Cu <- relationPlot(dat, "clay", "Cu", "black")
plot.dep.Cu <- relationPlot(dat, "depth", "Cu", "blue")
plot.al.Cu <- relationPlot(dat, "Al", "Cu", "purple")
plot.h.Cu <- relationPlot(dat, "proton", "Cu", "orange")
plot.sal.Cu <- relationPlot(dat, "salinity", "Cu", "brown")

plot.dist.Pb <- relationPlot(dat, "distance", "Pb", "green")
plot.clay.Pb <- relationPlot(dat, "clay", "Pb", "black")
plot.dep.Pb <- relationPlot(dat, "depth", "Pb", "blue")
plot.al.Pb <- relationPlot(dat, "Al", "Pb", "purple")
plot.h.Pb <- relationPlot(dat, "proton", "Pb", "orange")
plot.sal.Pb <- relationPlot(dat, "salinity", "Pb", "brown")

plot.dist.Zn <- relationPlot(dat, "distance", "Zn", "green")
plot.clay.Zn <- relationPlot(dat, "clay", "Zn", "black")
plot.dep.Zn <- relationPlot(dat, "depth", "Zn", "blue")
plot.al.Zn <- relationPlot(dat, "Al", "Zn", "purple")
plot.h.Zn <- relationPlot(dat, "proton", "Zn", "orange")
plot.sal.Zn <- relationPlot(dat, "salinity", "Zn", "brown")

plot.dist.Cd <- relationPlot(dat, "distance", "Cd", "green")
plot.clay.Cd <- relationPlot(dat, "clay", "Cd", "black")
plot.dep.Cd <- relationPlot(dat, "depth", "Cd", "blue")
plot.al.Cd <- relationPlot(dat, "Al", "Cd", "purple")
plot.h.Cd <- relationPlot(dat, "proton", "Cd", "orange")
plot.sal.Cd <- relationPlot(dat, "salinity", "Cd", "brown")

grid.arrange(plot.dist.Pb,plot.dep.Pb,plot.clay.Pb,plot.al.Pb,plot.h.Pb,plot.sal.Pb, 
             plot.dist.Cr,plot.dep.Cr,plot.clay.Cr,plot.al.Cr,plot.h.Cr,plot.sal.Cr,
             plot.dist.Ni,plot.dep.Ni,plot.clay.Ni,plot.al.Ni,plot.h.Ni,plot.sal.Ni,
             plot.dist.Cu,plot.dep.Cu,plot.clay.Cu,plot.al.Cu,plot.h.Cu,plot.sal.Cu,
             plot.dist.Zn,plot.dep.Zn,plot.clay.Zn,plot.al.Zn,plot.h.Zn,plot.sal.Zn,
             plot.dist.As,plot.dep.As,plot.clay.As,plot.al.As,plot.h.As,plot.sal.As,
             plot.dist.Cd,plot.dep.Cd,plot.clay.Cd,plot.al.Cd,plot.h.Cd,plot.sal.Cd,
             nrow=7, ncol=6) -> p.gather
ggsave(plot = p.gather,
       filename = paste(dirPreset("relation/driver"),"/gather_relation.png",sep = ""),
       dpi = 300, width = 12, height = 14)


