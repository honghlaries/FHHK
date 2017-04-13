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
  p1 <- spplot(dat, tag, do.log = F, main = paste(tag,suffix,sep = ""), xlab = "Longi", ylab = "Lati") 
  krigFormal <- as.formula(paste("(",tag,")~1"))
  mod <- variogram(krigFormal,dat, alpha = c(-50 + 45 * 0:4))
  fit <- fit.variogram(mod, model = modsel)
  p2 <- plot(mod,fit, main = tag)
  krig <- krige(krigFormal, dat, dat.grid, model = modsel)
  p3 <- spplot(krig["var1.pred"], main = paste(tag,suffix,sep = ""), xlab = "Longi", ylab = "Lati")
  png(paste(dirPreset("riskAssment/krig"),"/",tag,suffix,"_modelfix.png",sep = ""))
  grid.arrange(p1,p3,p2, ncol = 2, widths = c(15,15), heights = c(5,5))
  dev.off()
}

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

background <- datareadln() %>% 
  gather(trait, bk, Al:orgC) %>%
  select(siteID, trait, bk) %>%
  group_by(trait) %>%
  summarise(bk = mean(bk)) %>% 
  mutate(bk2 = bk / 17092.2)

dat <- datareadln() %>% 
  gather(trait, value, Al:orgC) %>% 
  inner_join(background, by = c("trait" = "trait")) %>%
  mutate(value = log2(value / bk /1.5)) %>%
  select(siteID, trait, value) %>%
  group_by(siteID, trait) %>%
  summarise(value = mean(value)) %>%
  spread(trait, value) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  dplyr::select(siteID:depth,Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd,C,N,S,orgC) %>%
  dplyr::mutate(latitudes = conveySp(latitudes),
                longitude = conveySp(longitude))
dat <- as.data.frame(dat)
coordinates(dat) <- ~longitude+latitudes

doKrig(dat, dat.grid, tag = "Al", suffix = "_Igeo", modsel = vgm(0.2,"Sph",0.7))
doKrig(dat, dat.grid, tag = "Fe", suffix = "_Igeo", modsel = vgm(0.04,"Lin",0,0.003))
doKrig(dat, dat.grid, tag = "Mn", suffix = "_Igeo", modsel = vgm(0.03,"Sph",0.7))
doKrig(dat, dat.grid, tag = "Pb", suffix = "_Igeo", modsel = vgm(0.025,"Lin",0,0.015))
doKrig(dat, dat.grid, tag = "Cr", suffix = "_Igeo", modsel = vgm(0.03,"Lin",0))
doKrig(dat, dat.grid, tag = "Ni", suffix = "_Igeo", modsel = vgm(0.05,"Lin",0))
doKrig(dat, dat.grid, tag = "Cu", suffix = "_Igeo", modsel = vgm(0.2,"Lin",0))
doKrig(dat, dat.grid, tag = "Zn", suffix = "_Igeo", modsel = vgm(0.4,"Lin",0))
doKrig(dat, dat.grid, tag = "As", suffix = "_Igeo", modsel = vgm(0.4,"Lin",0))
doKrig(dat, dat.grid, tag = "Cd", suffix = "_Igeo", modsel = vgm(0.3,"Lin",0))
doKrig(dat, dat.grid, tag = "C", suffix = "_Igeo", modsel = vgm(0.15,"Lin",0))
doKrig(dat, dat.grid, tag = "orgC", suffix = "_Igeo", modsel = vgm(0.2,"Sph",0.7))
doKrig(dat, dat.grid, tag = "N", suffix = "_Igeo", modsel = vgm(0.2,"Sph",0.7))
doKrig(dat, dat.grid, tag = "S", suffix = "_Igeo", modsel = vgm(0.2,"Sph",0.7))



dat <- datareadln() %>% 
  gather(trait, value, Al:Cd) %>% 
  inner_join(datareadln() %>%
               select(siteID, Al), by = c("siteID" = "siteID")) %>%
  mutate(value = value / Al) %>%
  inner_join(background, by = c("trait" = "trait")) %>%
  mutate(value = value / bk2) %>%
  select(siteID, trait, value) %>%
  group_by(siteID, trait) %>%
  summarise(value = mean(value)) %>%
  spread(trait, value) %>%
  dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
  dplyr::select(siteID:depth,Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd) %>%
  dplyr::mutate(latitudes = conveySp(latitudes),
                longitude = conveySp(longitude))
dat <- as.data.frame(dat)
coordinates(dat) <- ~longitude+latitudes

doKrig(dat, dat.grid, tag = "Al", suffix = "_Ef", modsel = vgm(0.2,"Sph",0.7))
doKrig(dat, dat.grid, tag = "Fe", suffix = "_Ef", modsel = vgm(0.04,"Lin",0,0.003))
doKrig(dat, dat.grid, tag = "Mn", suffix = "_Ef", modsel = vgm(0.03,"Sph",0.7))
doKrig(dat, dat.grid, tag = "Pb", suffix = "_Ef", modsel = vgm(0.025,"Lin",0,0.015))
doKrig(dat, dat.grid, tag = "Cr", suffix = "_Ef", modsel = vgm(0.03,"Lin",0))
doKrig(dat, dat.grid, tag = "Ni", suffix = "_Ef", modsel = vgm(0.05,"Lin",0))
doKrig(dat, dat.grid, tag = "Cu", suffix = "_Ef", modsel = vgm(0.2,"Lin",0))
doKrig(dat, dat.grid, tag = "Zn", suffix = "_Ef", modsel = vgm(0.4,"Lin",0))
doKrig(dat, dat.grid, tag = "As", suffix = "_Ef", modsel = vgm(0.4,"Lin",0))
doKrig(dat, dat.grid, tag = "Cd", suffix = "_Ef", modsel = vgm(0.3,"Lin",0))

