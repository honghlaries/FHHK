source("uniTls_pkgInstall.R")

doKrig.resamp <- function(dat, dat.grid.resample, tag, cutoff, 
                          krigFormula = as.formula(paste(tag,"~1",sep="")), 
                          modsel, nsamp, nsite, group) {
  pkgLoad("dplyr");pkgLoad("sp");pkgLoad("gstat");pkgLoad("foreach");pkgLoad("doParallel")
  
  cl<-makeCluster(4)
  registerDoParallel(cl)
  
  dat <- dat %>% select_("lon","lat",tag,group)
  
  reaminSiteSize <- floor(length(levels(dat$siteID)) - 
                            c(0.01,0.05,0.10,0.25,0.50) * length(levels(dat$siteID)))
  
  kirg.sampresamp.tot <- 
    foreach(i = 1:nsamp, .combine="rbind", .packages = c("dplyr","sp","gstat","foreach")) %dopar% {
      
      if (i %in% c(1,ceiling(1:9*0.1*nsamp),nsamp)) print(paste("Calculating on",tag,"(",i,"/",nsamp,")"))
      
      subdat <- dat %>% dplyr::group_by_(group) %>% dplyr::sample_n(1)
      subdat <- as.data.frame(subdat)
      coordinates(subdat) <- ~lon+lat
      
      mod <- variogram(krigFormula, subdat, cutoff = cutoff)
      #fit <- fit.variogram(mod, model = modsel)
      krig <- krige(krigFormula, subdat, dat.grid.resample, 
                    model = modsel, debug.level = 0)
      krig <- as.data.frame(krig) %>% 
        dplyr::select(lon, lat, value = var1.pred, var.mod =  var1.var) %>%
        dplyr::mutate(trait = tag, nsamp = i)
      
      subdat <- as.data.frame(subdat)
      
      kirg.siteresamp <-
        foreach(j = 1:nsite, .combine="rbind", .packages = c("dplyr","sp","gstat")) %do% {
          if (j %in% c(1,ceiling(1:9*0.1*nsite),nsite)) print(paste("Site resampling:","(",j,"/",nsite,")"))
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
          
          kirg.siteresamp1 %>%
            dplyr::inner_join(kirg.siteresamp2, by = c("lon", "lat")) %>%
            dplyr::inner_join(kirg.siteresamp3, by = c("lon", "lat")) %>%
            dplyr::inner_join(kirg.siteresamp4, by = c("lon", "lat")) %>%
            dplyr::inner_join(kirg.siteresamp5, by = c("lon", "lat")) %>%
            dplyr::mutate(nsite = j)
          
        }
      kirg.siteresamp %>% dplyr::inner_join(krig, by = c("lon", "lat")) 
    }

  stopImplicitCluster()
  stopCluster(cl)
  
  kirg.sampresamp.tot
}




timeLog <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}

meanExtract <- function(path) {
  
  pkgLoad("dplyr");pkgLoad("foreach");pkgLoad("doParallel")

  files <- Sys.glob(paste(path,"perm_*.csv",sep = '')) 
  
  cl<-makeCluster(4)
  registerDoParallel(cl)
  
  mean <- 
    foreach(i = files, .combine="rbind", .packages = c("dplyr")) %dopar% {
      read.csv(i) %>%
        dplyr::mutate(lon = as.character(lon),
                      lat = as.character(lat)) %>%
        dplyr::group_by(lon,lat,trait) %>%
        dplyr::summarise(value = mean(value, na.rm = T))
    }
  
  stopImplicitCluster()
  stopCluster(cl)
  
  mean
  
}





