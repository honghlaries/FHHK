## Initialization
rm(list = ls())
source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R");
source("anaTls_spatialView.R");source("anaTls_multivariate.R");
pkgInitialization(c("dplyr","tidyr","sp","gstat","gridExtra"))
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
hcluster(datareadln() %>%
           select(Al:Cd,orgC:sand,depth),
         rname = datareadln() %>%
           select(Al:Cd,orgC:sand,depth) %>% colnames(),
         flop = T) -> a

hcluster(datareadln() %>%
           select(Al:Cd),
         rname = datareadln()$siteID) -> a


