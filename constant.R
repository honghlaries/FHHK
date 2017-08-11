# Initialization
source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R");
pkgInitialization(c("dplyr","tidyr","sp","maptools"))

# constant
#lonRange = c(119.2,121.8);latRange = c(33.7,35)
lonRange = c(119.9,121.8);latRange = c(33.7,34.9)
alphalevel = 0.05

sites <- read.csv("data/meta_sites.csv")

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
                  isComplete = complete.cases(Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd,C,N,S,orgC,AVS,clay,silt,sand)) %>%
    dplyr::select(siteID:depth,distance,
                  isComplete,salinity,pH,
                  Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd,C,N,S,orgC,AVS,
                  clay,silt,sand)
}
