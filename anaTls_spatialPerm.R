## Initialization
rm(list = ls())
source("constant.R"); source("anaTls_spatialView.R"); source("grid_resamp.R")
pkgInitialization(c("dplyr","tidyr","sp","gstat"))

## Functions 
doKrig.resamp <- function(dat, dat.grid.resample, tag, cutoff, 
                          krigFormula = as.formula(paste(tag,"~1",sep="")), 
                          modsel, nsamp, nsite, group) {
  pkgLoad("sp");pkgLoad("gstat");pkgLoad("gridExtra");pkgLoad("dplyr")
  
  dat <- dat %>% select_("lon","lat",tag,group)
  krig.tot <- NULL
  for(i in 1:nsamp) {
    if (i %in% c(1,ceiling(1:9*0.1*nsamp),nsamp)) print(paste("Calculating on",tag,"(",i,"/",nsamp,")"))
    subdat <- dat %>% dplyr::group_by_(group) %>% dplyr::sample_n(1)
    subdat <- as.data.frame(subdat)
    coordinates(subdat) <- ~lon+lat
    mod <- variogram(krigFormula, subdat, cutoff = cutoff)
    #fit <- fit.variogram(mod, model = modsel)
    krig <- krige(krigFormula, subdat, dat.grid.resample, 
                  model = modsel, debug.level = 0)
    krig <- as.data.frame(krig) %>% dplyr::rename(value = var1.pred, var.mod =  var1.var)
    
    subdat <- as.data.frame(subdat)
    kirg.siteresamp.tot <- NULL
    for(j in 1:nsite) {
      if ((i %in% c(1,ceiling(1:9*0.1*nsamp),nsamp)) && 
          (j %in% c(1,ceiling(1:3*0.25*nsite),nsite))) print(paste("Site resampling:","(",j,"/",nsite,")"))
      subsubdat <- subdat %>% dplyr::filter(siteID %in% sample(levels(siteID), ceiling(0.75*length(levels(siteID)))))
      subsubdat <- as.data.frame(subsubdat)
      coordinates(subsubdat) <- ~lon+lat
      mod <- variogram(krigFormula, subsubdat, cutoff = cutoff)
      #fit <- fit.variogram(mod, model = modsel)
      kirg.siteresamp <- krige(krigFormula, subsubdat, dat.grid.resample, 
                               model = modsel, debug.level = 0)
      kirg.siteresamp <- as.data.frame(kirg.siteresamp) %>% dplyr::select(lon, lat, value = var1.pred) 
      kirg.siteresamp.tot <- rbind(kirg.siteresamp.tot, kirg.siteresamp)
    }
    krig <- kirg.siteresamp.tot %>% dplyr::group_by(lon,lat) %>% 
      dplyr::summarise(var.site = sd(value)) %>% dplyr::inner_join(krig)
    
    krig.tot <- rbind(krig.tot,krig)
  }
  
  krig.tot %>% dplyr::group_by(lon,lat) %>% 
    dplyr::summarise(mean = mean(value), var.mod = mean(var.mod), 
                     var.samp = sd(value), var.site = mean(var.site))
}

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
                           as.data.frame(cbind(grid.perm, trait = fitlist$tag[i])))
  }
  
  grid.perm.tot
}

doKrig.siteImpor <- function(dat, dat.grid.resample, tag, cutoff, 
                          krigFormula = as.formula(paste(tag,"~1",sep="")), 
                          modsel, nsamp, nsite, group) {
  pkgLoad("sp");pkgLoad("gstat");pkgLoad("gridExtra");pkgLoad("dplyr")
  
  dat <- dat %>% select_("lon","lat",tag,group)
  krig.tot <- NULL
  for(i in 1:nsamp) {
    print(paste("Calculating on",tag,"(",i,"/",nsamp,")"))
    subdat <- dat %>% dplyr::group_by_(group) %>% dplyr::sample_n(1)
    subdat <- as.data.frame(subdat)
    coordinates(subdat) <- ~lon+lat
    mod <- variogram(krigFormula, subdat, cutoff = cutoff)
    #fit <- fit.variogram(mod, model = modsel)
    krig <- krige(krigFormula, subdat, dat.grid.resample, 
                  model = modsel, debug.level = 0)
    krig <- as.data.frame(krig) %>% dplyr::rename(value = var1.pred, var.mod =  var1.var)
    
    subdat <- as.data.frame(subdat)
    kirg.siteresamp.tot <- NULL
    reaminSiteSize <- floor(length(levels(subdat$siteID)) - 
                              c(0.01,0.05,0.10,0.25,0.50) * length(levels(subdat$siteID)))
    print(paste("Five Discard level(",1:5,"):", 1 - reaminSiteSize/length(levels(subdat$siteID)) ))
    for(j in 1:nsite) {
      print(paste("Site resampling:","(",j,"/",nsite,")"))
      subsubdat <- subdat %>% dplyr::filter(siteID %in% sample(levels(siteID), reaminSiteSize[1]))
      subsubdat <- as.data.frame(subsubdat)
      coordinates(subsubdat) <- ~lon+lat
      mod <- variogram(krigFormula, subsubdat, cutoff = cutoff)
      #fit <- fit.variogram(mod, model = modsel)
      kirg.siteresamp <- krige(krigFormula, subsubdat, dat.grid.resample, 
                               model = modsel, debug.level = 0)
      kirg.siteresamp1 <- as.data.frame(kirg.siteresamp) %>% dplyr::select(lon, lat, value1 = var1.pred) 
      subsubdat <- subdat %>% dplyr::filter(siteID %in% sample(levels(siteID), reaminSiteSize[2]))
      subsubdat <- as.data.frame(subsubdat)
      coordinates(subsubdat) <- ~lon+lat
      mod <- variogram(krigFormula, subsubdat, cutoff = cutoff)
      #fit <- fit.variogram(mod, model = modsel)
      kirg.siteresamp <- krige(krigFormula, subsubdat, dat.grid.resample, 
                               model = modsel, debug.level = 0)
      kirg.siteresamp2 <- as.data.frame(kirg.siteresamp) %>% dplyr::select(lon, lat, value2 = var1.pred) 
      subsubdat <- subdat %>% dplyr::filter(siteID %in% sample(levels(siteID), reaminSiteSize[3]))
      subsubdat <- as.data.frame(subsubdat)
      coordinates(subsubdat) <- ~lon+lat
      mod <- variogram(krigFormula, subsubdat, cutoff = cutoff)
      #fit <- fit.variogram(mod, model = modsel)
      kirg.siteresamp <- krige(krigFormula, subsubdat, dat.grid.resample, 
                               model = modsel, debug.level = 0)
      kirg.siteresamp3 <- as.data.frame(kirg.siteresamp) %>% dplyr::select(lon, lat, value3 = var1.pred) 
      subsubdat <- subdat %>% dplyr::filter(siteID %in% sample(levels(siteID), reaminSiteSize[4]))
      subsubdat <- as.data.frame(subsubdat)
      coordinates(subsubdat) <- ~lon+lat
      mod <- variogram(krigFormula, subsubdat, cutoff = cutoff)
      #fit <- fit.variogram(mod, model = modsel)
      kirg.siteresamp <- krige(krigFormula, subsubdat, dat.grid.resample, 
                               model = modsel, debug.level = 0)
      kirg.siteresamp4 <- as.data.frame(kirg.siteresamp) %>% dplyr::select(lon, lat, value4 = var1.pred) 
      subsubdat <- subdat %>% dplyr::filter(siteID %in% sample(levels(siteID), reaminSiteSize[5]))
      subsubdat <- as.data.frame(subsubdat)
      coordinates(subsubdat) <- ~lon+lat
      mod <- variogram(krigFormula, subsubdat, cutoff = cutoff)
      #fit <- fit.variogram(mod, model = modsel)
      kirg.siteresamp <- krige(krigFormula, subsubdat, dat.grid.resample, 
                               model = modsel, debug.level = 0)
      kirg.siteresamp5 <- as.data.frame(kirg.siteresamp) %>% dplyr::select(lon, lat, value5 = var1.pred) 
      kirg.siteresamp <- kirg.siteresamp1 %>%
        dplyr::inner_join(kirg.siteresamp2, by = c("lon", "lat")) %>%
        dplyr::inner_join(kirg.siteresamp3, by = c("lon", "lat")) %>%
        dplyr::inner_join(kirg.siteresamp4, by = c("lon", "lat")) %>%
        dplyr::inner_join(kirg.siteresamp5, by = c("lon", "lat"))
      kirg.siteresamp.tot <- rbind(kirg.siteresamp.tot, kirg.siteresamp)
    }
    krig <- kirg.siteresamp.tot %>% dplyr::group_by(lon,lat) %>% 
      dplyr::summarise(var.site1 = sd(value1),
                       var.site2 = sd(value2),
                       var.site3 = sd(value3),
                       var.site4 = sd(value4),
                       var.site5 = sd(value5)) %>% 
      dplyr::inner_join(krig, by = c("lon", "lat"))
    
    krig.tot <- rbind(krig.tot,krig)
  }
  
  krig.tot %>% dplyr::group_by(lon,lat) %>% 
    dplyr::summarise(mean = mean(value), var.mod = mean(var.mod), 
                     var.samp = sd(value), var.site1 = mean(var.site1),
                     var.site2 = mean(var.site2), var.site3 = mean(var.site3),
                     var.site4 = mean(var.site4), var.site5 = mean(var.site5))
}

siteImpor.array <- function(fitlist, processor, seed = 20171216.082111) {
  dat <- datareadln()
  set.seed(seed)
  grid.perm.tot = NULL
  
  fitlist <- fitlist[fitlist$processor == processor,]
  
  for (i in 1:length(fitlist$tag)) { 
    grid.perm <- as.data.frame(with(data = fitlist,
                                    doKrig.siteImpor(dat, dat.grid.resample, 
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
                           as.data.frame(cbind(grid.perm, trait = fitlist$tag[i])))
  }
  
  grid.perm.tot
}


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

