## Initialization
rm(list = ls())
source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R");source("anaTls_spatialPerm.R");
pkgInitialization(c("dplyr","tidyr","sp","gstat"))
source("grid_resamp.R")


## Readln fit list
fitlist <- read.csv("data/meta_fitlist.csv")
fitlist$tag <- as.character(fitlist$tag)
fitlist$krigFormula <- as.character(fitlist$krigFormula)
fitlist$cutoff <- as.numeric(fitlist$cutoff)
fitlist$psill <- as.numeric(fitlist$psill)
fitlist$mod <- as.character(fitlist$mod)
fitlist$range <- as.numeric(fitlist$range)
fitlist$kappa <- as.numeric(fitlist$kappa)
fitlist$nsamp <- as.numeric(fitlist$nsamp)
fitlist$nsite <- as.numeric(fitlist$nsite)
fitlist$group <- as.character(fitlist$group)
fitlist$processor <- as.numeric(fitlist$processor)

## Functions 
doPerm.array <- function(fitlist, processor, seed = 31415926) {
  
  dat <- datareadln()
  set.seed(seed)
 
  fitlist <- fitlist[fitlist$processor == processor,]
  
  for (i in 1:length(fitlist$tag)) { 
    
    print(paste("Apply Perm On:",fitlist$tag[i]))
    
    start_time <- Sys.time()
    
    dat.perm <- 
      with(data = fitlist,
           doKrig.resamp(dat, dat.grid.resample, krigFormula = as.formula(krigFormula[i]), 
                         tag = tag[i], 
                         cutoff = cutoff[i],
                         modsel = vgm(psill = psill[i], 
                                      model = mod[i], 
                                      range = range[i], 
                                      kappa = kappa[i]), 
                         nsamp = nsamp[i], nsite = nsite[i], 
                         group = group[i]))
    
    print(paste("Complete Perm:",fitlist$tag[i]))
    
    write.csv(dat.perm, paste("data/perm_",tag[i],".csv", sep = ""), row.names = F)
    rm(dat.perm)
    
    print(paste("Time in Perm:",timeLog(start_time)))
  }
}

datareadln <- function() { 
  pkgLoad("dplyr");pkgLoad("tidyr")
  read.csv("data/result_element.csv") %>%
    dplyr::inner_join(read.csv("data/result_basic.csv"), by = c("sampleID" = "sampleID")) %>%
    dplyr::inner_join(read.csv("data/result_grainSize.csv"), by = c("sampleID" = "sampleID")) %>%
    dplyr::right_join(read.csv("data/meta_splList.csv"), by = c("sampleID" = "sampleID")) %>%
    dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
    dplyr::mutate(orgC = C.ac, AVS = S - S.ac) %>%
    dplyr::select(siteID,lon,lat,C,N,S,orgC,AVS,
                  salinity,pH,
                  Pb,Cr,Ni,Cu,Zn,Cd,
                  clay,silt,sand)
}



