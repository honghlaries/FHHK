## Initialization
rm(list = ls())
source("constant.R"); source("anaTls_spatialView.R"); source("grid_resamp.R")
pkgInitialization(c("dplyr","tidyr","sp","gstat","ggplot2","RColorBrewer","gridExtra"))

## Functions 
doPerm.array <- function(fitlist, processor, seed = 20171216.082111) {
  dat <- datareadln()
  set.seed(seed)
  grid.perm.tot = NULL
  
  fitlist <- fitlist[fitlist$processor == processor,]
  
  for (i in 1:length(fitlist$tag)) { 
    grid.perm <- as.data.frame(with(data = fitlist,
                                    doKrig.resamp(dat, dat.grid.resample, 
                                             krigFormula = as.formula(krigFormula[i]), 
                                             tag = tag[i], 
                                             cutoff = cutoff[i],
                                             modsel = vgm(psill = psill[i], 
                                                          model = mod[i], 
                                                          range = range[i], 
                                                          kappa = kappa[i]), 
                                             nsamp = nsamp[i], nsite = nsite[i], 
                                             group = group[i]) ) )
    grid.perm.tot <- rbind(grid.perm.tot, 
                           as.data.frame(cbind(grid.perm, trait = fitlist$tag)))
  }
  
  grid.perm.tot
}

## Example
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









