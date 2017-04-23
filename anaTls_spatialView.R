source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R")

spView.grid <- function(dat,dir,file,col.gardient,
                        lonC,latC,lonRange,latRange,zoom,
                        maptype,mapsource) {
  pkgLoad("dplyr");pkgLoad("tidyr");pkgLoad("ggplot2");pkgLoad("ggmap");
  plot <- ggmap(get_map(location = c(lat = latC,lon = lonC), zoom = zoom, 
                        maptype = maptype, source = mapsource)) + 
    geom_tile(aes(x = lon, y = lat, fill = value), data = dat) +
    scale_fill_gradient(low = "blue", high = "orange") +
    coord_quickmap(xlim = lonRange, ylim = latRange) + 
    facet_wrap(~trait)
  ggsave(filename = file, plot = plot)
  plot
}

spView <- function(dat,dir,file,col.gardient,
                   lonC,latC,lonRange,latRange,zoom,
                   maptype,mapsource) {
  pkgLoad("dplyr");pkgLoad("tidyr");pkgLoad("ggplot2");pkgLoad("ggmap");pkgLoad("maptools")
  bkmap <- readShapePoly("data/bou2_4p.shp")
  
  plot <- ggplot() + 
    geom_tile(aes(x = lon, y = lat, fill = value), data = dat) +
    geom_polygon(aes(x = long, y = lat, group = id), 
                 colour = "black", fill = "orange", data = fortify(bkmap)) +
    scale_fill_gradient(low = "blue", high = "red") +
    coord_quickmap(xlim = lonRange, ylim = latRange) 
  ggsave(filename = file, plot = plot)
  plot
}

conveySp <- function(x) {
  tmp <- strsplit(as.character(x),"°")[[1]]
  as.numeric(tmp[1]) + as.numeric(substr(tmp[2],1,nchar(tmp[2])-1))/60
}

doKrig <- function(dat, dat.grid, tag, cutoff, 
                   krigFormula = as.formula(paste(tag,"~1",sep="")), 
                   modsel, outTpe = "draft") {
  pkgLoad("sp");pkgLoad("gstat");pkgLoad("gridExtra");pkgLoad("geoR")
  
  p1 <- spplot(dat, tag, do.log = F, main = tag, xlab = "Longi", ylab = "Lati") 
  mod <- variogram(krigFormula,  dat, cutoff = cutoff)
  fit <- fit.variogram(mod, model = modsel)
  p2 <- plot(mod,fit, main = tag)
  krig <- krige(krigFormula, dat, dat.grid, model = modsel)
  col <- brewer.pal(9,"OrRd")
  p3 <- spplot(krig["var1.pred"], main = tag, xlab = "Longi", ylab = "Lati", 
               col.regions = col)
  
  if(outTpe == "data") return(krig)
  if(outTpe == "draft") return(grid.arrange(p1,p3,p2, ncol = 2, widths = c(15,15), heights = c(5,5)))
  if(outTpe == "hemivar") return(p2)
  if(outTpe == "map") return(p3)
}

doKrig.grid <- function(tag) {
  grid.value.tot <- NULL
  for(i in 1:length(tag)) {
    grid.value <- as.data.frame(doKrig(dat, dat.grid, log(Al) ~ 1, tag = tag[i], 
                                       outTpe = "data",modsel = vgm(0.15,"Sph",0.5))) %>% 
      select(lon, lat, value = var1.pred)
    grid.value.tot <- cbind(grid.value.tot, as.data.frame(grid.value, trait = tag[i]))
  }
  
}
