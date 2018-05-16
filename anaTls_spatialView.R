source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R")

conveySpCoord <- function(x, conveyFrom, conveyTo, type) {
  if(conveyFrom == "char.degree.nosuffix") {
    tmp <- strsplit(as.character(x),"°")
    out <- NULL
    for(i in 1:length(tmp)) {
      out <- c(out,as.numeric(tmp[[i]][1]) + as.numeric(substr(tmp[[i]][2],1,nchar(tmp[[i]][2])-1))/60)
    }
    out
  }
  
  if(conveyFrom == "num.degree.nosuffix") {
    out <- x
  }
  
  if(conveyTo == "char.minute.suffix") {
    if(type == "lat") {suffix <- if(out >= 0) "N" else "S"}
    if(type == "lon") {suffix <- if(out >= 0) "E" else "W"}
    out <- paste(floor(abs(out)),"°",format(60*(abs(out) %% 1),3),"'",suffix,sep = "")
  }
  out
}

findInt.coord <- function(range) {
  if(length(range) != 2) stop("error range input")
  low <- 0.9*min(range) + 0.1*max(range); high <- 0.9*max(range) + 0.1*min(range)
  
  r <- high - low
  res <- rep((r / (5:3)) %% 1, 5) - rep(c(0,1/2,1/3,1/4,1/6),each = 3)
  res[res >= 0.5] = res[res >= 0.5] - 1
  res <- matrix(abs(res),3,5)
  fac <-(5:3)[as.vector((res == min(res)) %*% matrix(rep(1,5),5,1))][1]
  seq(low,high,r/fac)
}

spView.grid <- function(dat, lonRange, latRange, leg.name, 
                        bins = 7, grad.col = rainbow(15),
                        grad.value = quantile(dat$value, probs = (1:bins)/(bins+1)), 
                        grad.tag = grad.value, pncol) {
  pkgLoad("dplyr");pkgLoad("tidyr");pkgLoad("ggplot2");pkgLoad("maptools");
  pkgLoad("rgdal")
  
  ggplot() + 
    geom_raster(aes(x = lon, y = lat, fill = value), data = dat) +
    scale_fill_gradientn(leg.name, breaks = grad.value, labels = grad.tag, colours = grad.col) +
    scale_x_continuous(name = "", expand = c(0,0)) +
    scale_y_continuous(name = "", expand = c(0,0)) +
    coord_quickmap(xlim = lonRange, ylim = latRange) +
    facet_wrap(~trait,ncol = pncol) + 
    theme_bw() + 
    theme(aspect.ratio = (latiRange[2]-latiRange[1])/(longiRange[2]-longiRange[1]),
          panel.grid = element_blank(),
          strip.background = element_blank())
}

spView.grid.interval <- function(dat, lonRange, latRange, leg.name, 
                                 bins = 7, grad.col = rainbow(15),
                                 grad.value = quantile(dat$value, probs = (1:bins)/(bins+1)),
                                 grad.tag = grad.value, pncol) {
  pkgLoad("dplyr");pkgLoad("tidyr");pkgLoad("ggplot2");pkgLoad("maptools");
  pkgLoad("rgdal")
  
  dat$class <- "c0"
  labels <- paste("<", grad.value[1])
  
  for(i in 1:(length(grad.value)-1)) {
    dat[(dat$value - grad.value[i])> 0.0001,"class"] <- paste("c", i, sep = "")
    labels <- c(labels, paste(grad.value[i], "~", grad.value[i+1]))
  }
  
  dat[(dat$value - grad.value[length(grad.value)])> 0.0001,"class"] <- paste("c", length(grad.value), sep = "")
  labels <- c(labels, paste(">", grad.value[length(grad.value)]))
  
  ggplot() + 
    geom_raster(aes(x = lon, y = lat, fill = class), data = dat) +
    scale_fill_manual(leg.name, labels = labels, 
                      limits = paste("c",0:length(grad.value),sep = ""), values = grad.col) +
    scale_x_continuous(name = "", expand = c(0,0)) +
    scale_y_continuous(name = "", expand = c(0,0)) +
    coord_quickmap(xlim = lonRange, ylim = latRange) +
    facet_wrap(~trait,ncol = pncol) + 
    theme_bw() + 
    theme(aspect.ratio = (latiRange[2]-latiRange[1])/(longiRange[2]-longiRange[1]),
          panel.grid = element_blank(),
          strip.background = element_blank())
}

spView <- function(dat, lonRange, latRange, leg.name, 
                   bins = 7, grad.col = rainbow(bins+1),
                   grad.value = NULL, grad.tag = NULL) {
  pkgLoad("dplyr");pkgLoad("tidyr");pkgLoad("ggplot2");pkgLoad("maptools");
  pkgLoad("rgdal")
  
  if(is.null(grad.value)) grad.value <- (0:bins)*(max(dat$value) - min(dat$value))/bins + min(dat$value)
  if(is.null(grad.tag)) grad.tag <- format(grad.value,digit = 2)
  
  break.lon <- findInt.coord(lonRange)
  label.lon <- conveySpCoord(break.lon,conveyFrom = "num.degree.nosuffix",
                             conveyTo = "char.minute.suffix",type = "lon")
  break.lat <- findInt.coord(latRange)
  label.lat <- conveySpCoord(break.lat,conveyFrom = "num.degree.nosuffix",
                             conveyTo = "char.minute.suffix",type = "lat")
  
  ggplot() + 
    geom_raster(aes(x = lon, y = lat, fill = value), 
                interpolate = T, show.legend = T, data = dat) +
    #geom_contour(aes(x = lon, y = lat,  z = value),
    #             col= "black", show.legend = T, size = 0.8,  bins = bins, data = dat) +
    scale_fill_gradientn(leg.name, guide = "colourbar",
                         breaks = grad.value, labels = grad.tag, 
                         colours = grad.col) +
    #scale_color_gradient2("leg.name", guide = "colourbar", low = "black", high = "white") +
    scale_x_continuous(name = "", breaks = break.lon, labels = label.lon, expand = c(0,0)) +
    scale_y_continuous(name = "", breaks = break.lat, labels = label.lat, expand = c(0,0)) +
    coord_quickmap(xlim = lonRange, ylim = latRange) +
    theme_bw() + 
    theme(aspect.ratio = (latiRange[2]-latiRange[1])/(longiRange[2]-longiRange[1]),
          panel.grid = element_blank(),
          legend.box.margin = margin(t = 2, r = 2, b = 2, l = 2),
          legend.box.spacing = unit(1,units = "pt"))
}

spView.interval <- function(dat, lonRange, latRange, leg.name, 
                            bins = 7, grad.col = rainbow(bins+1),
                            grad.value = NULL, grad.tag = NULL) {
  pkgLoad("dplyr"); pkgLoad("tidyr"); pkgLoad("ggplot2"); pkgLoad("maptools");
  pkgLoad("rgdal")
  
  dat$class <- "c0"
  labels <- paste("<", grad.value[1])
  
  for(i in 1:(length(grad.value)-1)) {
    dat[(dat$value - grad.value[i])> 0.0001,"class"] <- paste("c", i, sep = "")
    labels <- c(labels, paste(grad.value[i], "~", grad.value[i+1]))
  }
  
  dat[(dat$value - grad.value[length(grad.value)])> 0.0001,"class"] <- paste("c", length(grad.value), sep = "")
  labels <- c(labels, paste(">", grad.value[length(grad.value)]))
  
  if(is.null(grad.value)) grad.value <- (0:bins)*(max(dat$value) - min(dat$value))/bins + min(dat$value)
  if(is.null(grad.tag)) grad.tag <- format(grad.value,digit = 2)
  
  break.lon <- findInt.coord(lonRange)
  label.lon <- conveySpCoord(break.lon,conveyFrom = "num.degree.nosuffix",
                             conveyTo = "char.minute.suffix",type = "lon")
  break.lat <- findInt.coord(latRange)
  label.lat <- conveySpCoord(break.lat,conveyFrom = "num.degree.nosuffix",
                             conveyTo = "char.minute.suffix",type = "lat")
  
  ggplot() + 
    geom_raster(aes(x = lon, y = lat, fill = class), 
                interpolate = T, show.legend = T, data = dat) +
    scale_fill_manual(leg.name, labels = labels, 
                      limits = paste("c",0:length(grad.value),sep = ""), values = grad.col) +
    scale_x_continuous(name = "", breaks = break.lon, labels = label.lon, expand = c(0,0)) +
    scale_y_continuous(name = "", breaks = break.lat, labels = label.lat, expand = c(0,0)) +
    coord_quickmap(xlim = lonRange, ylim = latRange) +
    theme_bw() + 
    theme(aspect.ratio = (latiRange[2]-latiRange[1])/(longiRange[2]-longiRange[1]),
          panel.grid = element_blank(),
          legend.box.margin = margin(t = 2, r = 2, b = 2, l = 2),
          legend.box.spacing = unit(1,units = "pt"))
}


doKrig <- function(dat, dat.grid, tag, cutoff, 
                   krigFormula = as.formula(paste(tag,"~1",sep="")), 
                   modsel, quietmode = FALSE, addlog = NULL) {
  pkgLoad("sp");pkgLoad("gstat");pkgLoad("gridExtra");pkgLoad("RColorBrewer")
  
  fitting <- FALSE
  while (!fitting) {
    mod <- variogram(krigFormula,  dat, cutoff = cutoff)
    fit <- fit.variogram(mod, model = modsel)
    krig <- krige(krigFormula, dat, dat.grid, model = modsel, debug.level = 0)
    if(quietmode) {return(krig)}
    col <- brewer.pal(9,"OrRd")
    p1 <- spplot(dat, tag, do.log = F, main = tag, xlab = "Longi", ylab = "Lati") 
    p2 <- plot(mod,fit, main = tag)
    p3 <- spplot(krig["var1.pred"], main = tag, xlab = "Longi", ylab = "Lati", 
                 col.regions = col)
    grid.arrange(p1,p3,p2, ncol = 2, widths = c(15,15), heights = c(5,5))
    print("Is fitting correct? 1:YES; 0:NO.")
    fitting <- scan(n=1)
    if(!fitting) {
      print("Change condition? 1:YES; 0:NO."); if(!scan(n=1)) stop("fitting stopped") 
      print("Change formula? 1:YES; 0:NO.")
      if(scan(n=1)) {print("Input New Formula:"); krigFormula = as.formula(scan(n=1,what = character()))}
      print("Change cutoff? 1:YES; 0:NO.")
      if(scan(n=1)) {print("Input New Cutoff:"); krigFormula = as.formula(scan(n=1))}
      #print("Change model? 1:YES; 0:NO.")
      #if(scan(n=1)) {print("Input New Model:"); modsel = as.formula(scan(n=1))}
    }  
  }
  
  if(is.null(addlog)) {print("Add log file? 1:YES; 0:NO."); addlog <- scan(n=1)}
  if(addlog) ggsave(filename = scan(what = character()), 
                    plot = grid.arrange(p1,p3,p2, ncol = 2, 
                                        widths = c(15,15), heights = c(5,5)))
  krig
}



