# Initialization
source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R");
pkgInitialization(c("dplyr","tidyr","sp","maptools"))

# constant
#lonRange = c(119.2,121.8);latRange = c(33.7,35)
lonRange = c(119.9,121.8);latRange = c(33.7,34.9)
alphalevel = 0.05

sites <- read.csv("data/meta_sites.csv")

background <- read.csv("data/meta_baseline.csv") %>%
  mutate(bk2 = bk / 27786)

#rivers <- readShapeLines("data/item_rivers.shp")
#proj4string(rivers) <- CRS("+init=epsg:4030")
#rivers <- spTransform(rivers, 
#                      CRS("+init=EPSG:4326"))
#ggplot() + 
#  geom_polygon(aes(x = long, y = lat, group = group), 
#               colour = "black", fill = "grey80", 
#               data = fortify(readShapePoly("data/bou2_4p.shp"))) +
#  geom_path(aes(x = long, y = lat, group = group), 
#            colour = "blue", data = fortify(rivers)) 

# data readln
datareadln <- function() { 
  pkgLoad("dplyr");pkgLoad("tidyr")
  read.csv("data/result_element.csv") %>%
    dplyr::inner_join(read.csv("data/result_basic.csv"), by = c("sampleID" = "sampleID")) %>%
    dplyr::inner_join(read.csv("data/result_grainSize.csv"), by = c("sampleID" = "sampleID")) %>%
    dplyr::right_join(read.csv("data/meta_splList.csv"), by = c("sampleID" = "sampleID")) %>%
    dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
    dplyr::mutate(orgC = C.ac, AVS = S - S.ac, 
                  isComplete = complete.cases(Fe,Mn,Pb,Cr,Ni,Cu,Zn,Cd,C,N,S,orgC,AVS,clay,silt,sand)) %>%
    dplyr::select(siteID:depth,distance,
                  isComplete,salinity,pH,
                  Fe,Mn,Pb,Cr,Ni,Cu,Zn,Cd,C,N,S,orgC,AVS,
                  clay,silt,sand)
}

cname <- c("lon","lat")
coo.river1 <- as.data.frame(getKMLcoordinates("data/river chan 1.kml", ignoreAltitude = T)[[1]])
colnames(coo.river1) <- cname
coo.river2_1 <- as.data.frame(getKMLcoordinates("data/river chan 2-1.kml", ignoreAltitude = T)[[1]])
colnames(coo.river2_1) <- cname
coo.river2_2 <- as.data.frame(getKMLcoordinates("data/river chan 2-2.kml", ignoreAltitude = T)[[1]])
colnames(coo.river2_2) <- cname
coo.river2_3 <- as.data.frame(getKMLcoordinates("data/river chan 2-3.kml", ignoreAltitude = T)[[1]])
colnames(coo.river2_3) <- cname
coo.river2_4 <- as.data.frame(getKMLcoordinates("data/river chan 2-4.kml", ignoreAltitude = T)[[1]])
colnames(coo.river2_4) <- cname
coo.river2_5 <- as.data.frame(getKMLcoordinates("data/river chan 2-5.kml", ignoreAltitude = T)[[1]])
colnames(coo.river2_5) <- cname
coo.river2_6 <- as.data.frame(getKMLcoordinates("data/river chan 2-6.kml", ignoreAltitude = T)[[1]])
colnames(coo.river2_6) <- cname
coo.river2_7 <- as.data.frame(getKMLcoordinates("data/river chan 2-7.kml", ignoreAltitude = T)[[1]])
colnames(coo.river2_7) <- cname
coo.river2_8 <- as.data.frame(getKMLcoordinates("data/river chan 2-8.kml", ignoreAltitude = T)[[1]])
colnames(coo.river2_8) <- cname
coo.river2_9 <- as.data.frame(getKMLcoordinates("data/river chan 2-9.kml", ignoreAltitude = T)[[1]])
colnames(coo.river2_9) <- cname
coo.river3 <- as.data.frame(getKMLcoordinates("data/river chan 3.kml", ignoreAltitude = T)[[1]])
colnames(coo.river3) <- cname
coo.river4_1 <- as.data.frame(getKMLcoordinates("data/river chan 4-1.kml", ignoreAltitude = T)[[1]])
colnames(coo.river4_1) <- cname
coo.river4_2 <- as.data.frame(getKMLcoordinates("data/river chan 4-2.kml", ignoreAltitude = T)[[1]])
colnames(coo.river4_2) <- cname
coo.river4_3 <- as.data.frame(getKMLcoordinates("data/river chan 4-3.kml", ignoreAltitude = T)[[1]])
colnames(coo.river4_3) <- cname
coo.river5 <- as.data.frame(getKMLcoordinates("data/river chan 5.kml", ignoreAltitude = T)[[1]])
colnames(coo.river5) <- cname
coo.river6 <- as.data.frame(getKMLcoordinates("data/river chan 6.kml", ignoreAltitude = T)[[1]])
colnames(coo.river6) <- cname
coo.river7 <- as.data.frame(getKMLcoordinates("data/river chan 7.kml", ignoreAltitude = T)[[1]])
colnames(coo.river7) <- cname
coo.river8 <- as.data.frame(getKMLcoordinates("data/river chan 8.kml", ignoreAltitude = T)[[1]])
colnames(coo.river8) <- cname
coo.river9 <- as.data.frame(getKMLcoordinates("data/river chan 9.kml", ignoreAltitude = T)[[1]])
colnames(coo.river9) <- cname
coo.river10 <- as.data.frame(getKMLcoordinates("data/river chan 10.kml", ignoreAltitude = T)[[1]])
colnames(coo.river10) <- cname

coo.1855 <- read.csv("data/meta_shoreline_1855.csv")

coo.oldmouth <- data.frame(lon = 120.45, lat = 34.36)
