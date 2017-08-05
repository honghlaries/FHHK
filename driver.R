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
