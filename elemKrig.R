## clean ----
rm(list = ls())
source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R")
pkgInitialization(c("dplyr","tidyr","sp", "gstat"))

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

doKrig <- function(dat, dat.grid, krigFormula, tag, suffix = "", modsel, isData = F) {
  pkgLoad("sp");pkgLoad("gstat");pkgLoad("gridExtra")
  p1 <- spplot(dat, tag, do.log = F, main = tag, xlab = "Longi", ylab = "Lati") 
  mod <- variogram(krigFormula,dat, cutoff = 1.5)
  fit <- fit.variogram(mod, model = modsel)
  p2 <- plot(mod,fit, main = tag)
  krig <- krige(krigFormula, dat, dat.grid, model = modsel)
  p3 <- spplot(krig["var1.pred"], main = tag, xlab = "Longi", ylab = "Lati")
  png(paste(dirPreset("element/krig"),"/",tag,suffix,".png",sep = ""))
  grid.arrange(p1,p3,p2, ncol = 2, widths = c(15,15), heights = c(5,5))
  dev.off()
  if(isData) krig else p3
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

dat.grid.silt <- doKrig(dat, dat.grid, log(silt) ~ 1, tag = "silt", suffix = "_ord_", modsel = vgm(800,"Sph",1.0), isData = T)
dat.grid.silt <- as.data.frame(dat.grid.silt) %>% select(longitude, latitudes, silt = var1.pred)
coordinates(dat.grid.silt) <- ~longitude+latitudes

doKrig(dat, dat.grid.silt, log(Al) ~ silt, tag = "Al", suffix = "_rm_silt_",modsel = vgm(0.15,"Sph",0.5))
doKrig(dat, dat.grid.silt, log(Fe) ~ silt, tag = "Fe", suffix = "_rm_silt_",modsel = vgm(0.01,"Mat",0.5,0.01,kappa = 1))
doKrig(dat, dat.grid.silt, log(Mn) ~ silt, tag = "Mn", suffix = "_rm_silt_", modsel = vgm(0.035,"Sph",0.3))
doKrig(dat, dat.grid.silt, log(Pb) ~ silt, tag = "Pb", suffix = "_rm_silt_",modsel = vgm(0.35,"Sph",0.5))
doKrig(dat, dat.grid.silt, log(Cr) ~ silt, tag = "Cr", suffix = "_rm_silt_",modsel = vgm(0.06,"Mat",0.7,kappa = 1))
doKrig(dat, dat.grid.silt, log(Ni) ~ silt, tag = "Ni", suffix = "_rm_silt_",modsel = vgm(0.06,"Sph",1))
doKrig(dat, dat.grid.silt, log(Cu) ~ silt, tag = "Cu", suffix = "_rm_silt_",modsel = vgm(0.25,"Sph",0.75))
doKrig(dat, dat.grid.silt, log(Zm) ~ silt, tag = "Zn", suffix = "_rm_silt_",modsel = vgm(0.1,"Sph",0.5))
doKrig(dat, dat.grid.silt, log(As) ~ silt, tag = "As", suffix = "_rm_silt_",modsel = vgm(0.3,"Sph",0.75))
doKrig(dat, dat.grid.silt, log(Cd) ~ silt, tag = "Cd", suffix = "_rm_silt_",modsel = vgm(0.2,"Sph",0.7))
doKrig(dat, dat.grid.silt, log(C) ~ silt, tag = "C", suffix = "_rm_silt_",modsel = vgm(0.15,"Sph",0.7))
doKrig(dat, dat.grid.silt, log(orgC) ~ silt, tag = "orgC", suffix = "_rm_silt_",modsel = vgm(0.06,"Mat",0.7,kappa = 1))
doKrig(dat, dat.grid.silt, log(N) ~ silt, tag = "N", suffix = "_rm_silt_",modsel = vgm(0.25,"Sph",1))
doKrig(dat, dat.grid.silt, log(S) ~ silt, tag = "S", suffix = "_rm_silt_",modsel = vgm(0.4,"Sph",1))
doKrig(dat, dat.grid.silt, log(AVS) ~ silt, tag = "AVS", suffix = "_rm_silt_",modsel = vgm(0.3,"Sph",0.7))
doKrig(dat, dat.grid.silt, log(silt) ~ silt, tag = "silt", suffix = "_rm_silt_",modsel = vgm(200,"Sph",1.0))
#mod <- variogram(log(Cd) ~ depth, dat,  cutoff = 1.5)
#fit <- fit.variogram(mod, model = vgm(0.2,"Sph",0.7))
#plot(mod,fit)



g <- gstat(NULL, "logAl", log(Al) ~ 1, dat)
g <- gstat(g, "logAVS", log(AVS) ~ 1, dat)

mod <- variogram(g , cutoff = 1.5)
fit <- fit.lmc(mod, g, vgm(0.2, "Sph", 0.7))
plot(mod,fit)
maps <- predict(fit, dat.grid)
spplot.vcov(maps)