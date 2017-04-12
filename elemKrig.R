## clean ----
rm(list = ls())
source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R")
pkgInitialization(c("dplyr","tidyr","sp"))

## Functions ----
datareadln <- function() { ## data readln
  pkgLoad("dplyr");pkgLoad("tidyr")
  read.csv("data/result_element.csv") %>%
    dplyr::inner_join(read.csv("data/result_grainSize.csv"), by = c("sampleID" = "sampleID")) %>%
    dplyr::right_join(read.csv("data/meta_splList.csv"), by = c("sampleID" = "sampleID")) %>%
    dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
    dplyr::mutate(orgC = C.ac, AVS = S - S.ac, 
                  isComplete = complete.cases(Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd,C,N,S,orgC,AVS,clay,silt,sand)) %>%
    dplyr::select(siteID:depth,isComplete,Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd,C,N,S,orgC,AVS,clay,silt,sand)
}

conveySp <- function(x) {
  tmp <- strsplit(as.character(x),"°")[[1]]
  as.numeric(tmp[1]) + as.numeric(substr(tmp[2],1,nchar(tmp[2])-1))/60
  }

doKrig <- function(dat, dat.grid, tag, suffix = "", modsel) {
  pkgLoad("sp");pkgLoad("gstat");pkgLoad("gridExtra")
  p1 <- spplot(dat, tag, do.log = F, main = tag, xlab = "Longi", ylab = "Lati") 
  krigFormal <- as.formula(paste("log(",tag,")~1"))
  mod <- variogram(krigFormal,dat, alpha = c(-50 + 90 * 0:1, 90 * 0:1))
  fit <- fit.variogram(mod, model = modsel)
  p2 <- plot(mod,fit, main = tag)
  krig <- krige(krigFormal, dat, dat.grid, model = modsel)
  p3 <- spplot(krig["var1.pred"], main = tag, xlab = "Longi", ylab = "Lati")
  png(paste(dirPreset("element/krig"),"/",tag,"_modelfix.png",sep = ""))
  grid.arrange(p1,p3,p2, ncol = 2, widths = c(15,15), heights = c(5,5))
  dev.off()
}

dat <- datareadln() %>% 
  gather(trait, value, Al:sand) %>%
  select(siteID, trait, value) %>%
  group_by(siteID, trait) %>%
  summarise(value = mean(value)) %>%
  spread(trait, value) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  dplyr::select(siteID:depth,Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd,C,N,S,orgC,AVS,clay,silt,sand) %>%
  dplyr::mutate(latitudes = conveySp(latitudes),
                longitude = conveySp(longitude),
                AvsRatio = AVS / orgC)
dat <- as.data.frame(dat)
coordinates(dat) <- ~longitude+latitudes

longiRange <-  seq(from = 119.9, to = 121.8, length.out = 125)
latiRange <-  seq(from = 33.7, to = 34.9, length.out = 125)
dat.grid <- data.frame(latitudes = c(1), longitude = c(1))
for (i in 1:125) {
  for (j in 1:125) {
    if((longiRange[j] - 119.9) * (latiRange[i] - 33.7) - (longiRange[j] - 120.6) * (latiRange[i] - 34.5) > 0) {
      dat.grid <- rbind(dat.grid, c(latiRange[i],longiRange[j]))
    }
  }
}
dat.grid <- dat.grid[-1,]
coordinates(dat.grid) <- ~longitude+latitudes

doKrig(dat, dat.grid, tag = "Al", modsel = vgm(0.2,"Sph",0.7))
doKrig(dat, dat.grid, tag = "Fe", modsel = vgm(0.04,"Lin",0,0.003))
doKrig(dat, dat.grid, tag = "Mn", modsel = vgm(0.03,"Sph",0.7))
doKrig(dat, dat.grid, tag = "Pb", modsel = vgm(0.025,"Lin",0,0.015))
doKrig(dat, dat.grid, tag = "Cr", modsel = vgm(0.03,"Lin",0))
doKrig(dat, dat.grid, tag = "Ni", modsel = vgm(0.05,"Lin",0))
doKrig(dat, dat.grid, tag = "Cu", modsel = vgm(0.2,"Lin",0))
doKrig(dat, dat.grid, tag = "Zn", modsel = vgm(0.4,"Lin",0))
doKrig(dat, dat.grid, tag = "As", modsel = vgm(0.4,"Lin",0))
doKrig(dat, dat.grid, tag = "Cd", modsel = vgm(0.3,"Lin",0))
doKrig(dat, dat.grid, tag = "C", modsel = vgm(0.15,"Lin",0))
doKrig(dat, dat.grid, tag = "orgC", modsel = vgm(0.2,"Sph",0.7))
doKrig(dat, dat.grid, tag = "N", modsel = vgm(0.2,"Sph",0.7))
doKrig(dat, dat.grid, tag = "S", modsel = vgm(0.2,"Sph",0.7))
doKrig(dat, dat.grid, tag = "AVS", modsel = vgm(0.2,"Sph",0.7))


doKrig(dat, dat.grid, tag = "AvsRatio", modsel = vgm(0.04,"Wav",0.7,0.02))
