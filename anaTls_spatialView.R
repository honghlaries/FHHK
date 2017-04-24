source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R")

spView.grid <- function(dat,dir,file,leg.name, grad.value, grad.tag,
                        lonRange,latRange) {
  pkgLoad("dplyr");pkgLoad("tidyr");pkgLoad("ggplot2");pkgLoad("maptools");
  pkgLoad("rgdal");pkgLoad("directlabels")
  #CRS("+proj=tmerc +lat_0=0 +lon_0=108 +k=1 +x_0=36500000 +y_0=0 +ellps=krass +towgs84=15")
  #CRS("+proj=longlat +ellps=WGS84")
  #bkmap <- readOGR("data/zone_province.shp")
  #proj4string(bkmap) <- CRS("+proj=tmerc +lat_0=0 +lon_0=108 +k=1 +x_0=36500000 +y_0=0 +ellps=krass +towgs84=15")
  #bkmap <- spTransform(bkmap,  CRS("+proj=longlat +ellps=WGS84"))
  #rivers <- readShapeLines("data/item_rivers.shp")
  #railway <- readShapeLines("data/item_railway.shp")
  #roads <- readShapeLines("data/item_roads.shp")
  bkmap <- readShapePoly("data/bou2_4p.shp")
  plot <- ggplot() + 
    geom_tile(aes(x = lon, y = lat, fill = value), data = dat) +
    geom_contour(aes(x = lon, y = lat, z = value), col = "black", binwidth = 0.25, data = dat) +
    geom_polygon(aes(x = long, y = lat, group = group), 
                 colour = "black", fill = "orange", data = fortify(bkmap)) +
    scale_fill_gradientn(leg.name, breaks = grad.value, labels = grad.tag,colours = rainbow(15)) +
    
    coord_quickmap(xlim = lonRange, ylim = latRange) +
    facet_wrap(~trait,ncol = 3) + 
    theme_bw() + 
    theme(aspect.ratio = 1/2,
          panel.grid = element_blank(),
          strip.background = element_blank()
          )
  
  #plot <- direct.label(plot, method="top.pieces")
  ggsave(filename = paste(dirPreset(dir),"/", file, sep = ""), plot = plot)
  plot
}

spView <- function(dat,dir,file,col.gardient,
                   lonRange,latRange) {
  pkgLoad("dplyr");pkgLoad("tidyr");pkgLoad("ggplot2");pkgLoad("maptools");
  pkgLoad("rgdal");pkgLoad("directlabels")
  #CRS("+proj=tmerc +lat_0=0 +lon_0=108 +k=1 +x_0=36500000 +y_0=0 +ellps=krass +towgs84=15")
  #CRS("+proj=longlat +ellps=WGS84")
  #bkmap <- readOGR("data/zone_province.shp")
  #proj4string(bkmap) <- CRS("+proj=tmerc +lat_0=0 +lon_0=108 +k=1 +x_0=36500000 +y_0=0 +ellps=krass +towgs84=15")
  #bkmap <- spTransform(bkmap,  CRS("+proj=longlat +ellps=WGS84"))
  #rivers <- readShapeLines("data/item_rivers.shp")
  #railway <- readShapeLines("data/item_railway.shp")
  #roads <- readShapeLines("data/item_roads.shp")
  bkmap <- readShapePoly("data/bou2_4p.shp")
  plot <- ggplot() + 
    geom_tile(aes(x = lon, y = lat, fill = value), data = dat) +
    geom_contour(aes(x = lon, y = lat, z = value, colour = ..level..), binwidth = 0.2, data = dat) +
    geom_polygon(aes(x = long, y = lat, group = group), 
                 colour = "black", fill = "orange", data = fortify(bkmap)) +
    scale_fill_continuous(low = "blue", high = "red") +
    scale_colour_gradient(low = "black", high = "black") +
    coord_quickmap(xlim = lonRange, ylim = latRange) +
    theme_bw()
  
  plot <- direct.label(plot, method="top.pieces")
  ggsave(filename = paste(dirPreset(dir),"/", file, sep = ""), plot = plot)
  plot
}

conveySp <- function(x) {
  tmp <- strsplit(as.character(x),"°")
  out <- NULL
  for(i in 1:length(tmp)) {
    out <- c(out,as.numeric(tmp[[i]][1]) + as.numeric(substr(tmp[[i]][2],1,nchar(tmp[[i]][2])-1))/60)
  }
  out
}

doKrig <- function(dat, dat.grid, tag, cutoff, 
                   krigFormula = as.formula(paste(tag,"~1",sep="")), 
                   modsel, dir) {
  pkgLoad("sp");pkgLoad("gstat");pkgLoad("gridExtra");pkgLoad("RColorBrewer")
  
  p1 <- spplot(dat, tag, do.log = F, main = tag, xlab = "Longi", ylab = "Lati") 
  mod <- variogram(krigFormula,  dat, cutoff = cutoff)
  fit <- fit.variogram(mod, model = modsel)
  p2 <- plot(mod,fit, main = tag)
  krig <- krige(krigFormula, dat, dat.grid, model = modsel)
  col <- brewer.pal(9,"OrRd")
  p3 <- spplot(krig["var1.pred"], main = tag, xlab = "Longi", ylab = "Lati", 
               col.regions = col)
  ggsave(filename = paste(dirPreset(dir),"/", "modelfit_",tag,".png", sep = ""), 
         plot = grid.arrange(p1,p3,p2, ncol = 2, widths = c(15,15), heights = c(5,5)))
  krig
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
