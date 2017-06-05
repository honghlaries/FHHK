source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R")

conveySp <- function(x) {
  tmp <- strsplit(as.character(x),"°")
  out <- NULL
  for(i in 1:length(tmp)) {
    out <- c(out,as.numeric(tmp[[i]][1]) + as.numeric(substr(tmp[[i]][2],1,nchar(tmp[[i]][2])-1))/60)
  }
  out
}

spView.grid <- function(dat,leg.name, grad.value, grad.tag, grad.col = rainbow(15), 
                        binwidth, lonRange, latRange, pncol) {
  pkgLoad("dplyr");pkgLoad("tidyr");pkgLoad("ggplot2");pkgLoad("maptools");
  pkgLoad("rgdal");pkgLoad("directlabels")
  bkmap <- readShapePoly("data/bou2_4p.shp")
  ggplot() + 
    geom_raster(aes(x = lon, y = lat, fill = value), data = dat) +
    geom_contour(aes(x = lon, y = lat, z = value), col = "black", binwidth = binwidth, data = dat) +
    geom_polygon(aes(x = long, y = lat, group = group), 
                 colour = "black", fill = "grey80", data = fortify(bkmap)) +
    scale_fill_gradientn(leg.name, breaks = grad.value, labels = grad.tag, colours = grad.col) +
    coord_quickmap(xlim = lonRange, ylim = latRange) +
    facet_wrap(~trait,ncol = pncol) + 
    theme_bw() + 
    theme(aspect.ratio = (latiRange[2]-latiRange[1])/(longiRange[2]-longiRange[1]),
          panel.grid = element_blank(),
          strip.background = element_blank())
}

spView <- function(dat,
                   leg.name, grad.value, grad.tag, grad.col = rainbow(15), 
                   binwidth, lonRange, latRange) {
  pkgLoad("dplyr");pkgLoad("tidyr");pkgLoad("ggplot2");pkgLoad("maptools");
  pkgLoad("rgdal");pkgLoad("directlabels")
  bkmap <- readShapePoly("data/bou2_4p.shp")
  latiLab <- 
    
  plot <- ggplot() + 
    geom_raster(aes(x = lon, y = lat, fill = value), data = dat) +
    geom_contour(aes(x = lon, y = lat,  z = value), 
                 size = 0.8, col = "black", bins = 9, data = dat) +
    geom_polygon(aes(x = long, y = lat, group = group), 
                 colour = "black", fill = "grey80", data = fortify(bkmap)) +
    #scale_fill_gradientn(leg.name, guide = "colourbar",
    #                     breaks = grad.value, labels = grad.tag, colours = grad.col) +
    #scale_fill_gradient(leg.name, guide = "colourbar",
    #                     breaks = grad.value[c(1,7)], labels = grad.tag[c(1,7)], low = "white", high = "blue") +
    scale_colour_gradient(low = "black", high = "black") +
    scale_x_continuous(name = "") +
    scale_y_continuous(name = "") +
    coord_quickmap(xlim = lonRange, ylim = latRange) +
    theme_bw() + 
    theme(aspect.ratio = (latiRange[2]-latiRange[1])/(longiRange[2]-longiRange[1]),
          panel.grid = element_blank())
  
  #plot <- direct.label(plot, method="top.pieces")
  plot
}

doKrig <- function(dat, dat.grid, tag, cutoff, 
                   krigFormula = as.formula(paste(tag,"~1",sep="")), 
                   modsel, dir = NULL) {
  pkgLoad("sp");pkgLoad("gstat");pkgLoad("gridExtra");pkgLoad("RColorBrewer")
  
  mod <- variogram(krigFormula,  dat, cutoff = cutoff)
  fit <- fit.variogram(mod, model = modsel)
  p2 <- plot(mod,fit, main = tag)
  if (is.null(dir)) return(p2)
  krig <- krige(krigFormula, dat, dat.grid, model = modsel)
  col <- brewer.pal(9,"OrRd")
  p1 <- spplot(dat, tag, do.log = F, main = tag, xlab = "Longi", ylab = "Lati") 
  p3 <- spplot(krig["var1.pred"], main = tag, xlab = "Longi", ylab = "Lati", 
               col.regions = col)
  ggsave(filename = paste(dirPreset(dir),"/", "modelfit_",tag,".png", sep = ""), 
         plot = grid.arrange(p1,p3,p2, ncol = 2, widths = c(15,15), heights = c(5,5)))
  krig
}

