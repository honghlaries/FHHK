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

distCalc <- function(lata, lona, latb, lonb, r = 6371){
  r*acos(cos(lata/180*pi)*cos(lona/180*pi)*cos(latb/180*pi)*cos(lonb/180*pi)+
           cos(lata/180*pi)*sin(lona/180*pi)*cos(latb/180*pi)*sin(lonb/180*pi)+
           sin(lata/180*pi)*sin(latb/180*pi))
}

azimuthCalc <- function(lata, lona, latb, lonb, latc =90, lonc = lonb){
  p <- (acos(cos(lata/180*pi)*cos(lona/180*pi)*cos(latb/180*pi)*cos(lonb/180*pi)+
               cos(lata/180*pi)*sin(lona/180*pi)*cos(latb/180*pi)*sin(lonb/180*pi)+
               sin(lata/180*pi)*sin(latb/180*pi))+
          (90-lata)/180*pi +
          (90-latb)/180*pi )/2.0
  cos.azi <- sqrt(sin(p)*sin(p-(90-latb)/180*pi)/sin((90-lata)/180*pi)/
                    sin(acos(cos(lata/180*pi)*cos(lona/180*pi)*cos(latb/180*pi)*cos(lonb/180*pi)+
                               cos(lata/180*pi)*sin(lona/180*pi)*cos(latb/180*pi)*sin(lonb/180*pi)+
                               sin(lata/180*pi)*sin(latb/180*pi))))
  
  if(lonb >= lona) {
    if(latb >= lata) {0.5*pi - 2*acos(cos.azi)} else {2.5*pi - 2*acos(cos.azi)}
  } else {
    if(latb >= lata) {0.5*pi - 2*acos(cos.azi)} else {2.5*pi - 2*acos(cos.azi)}
  }
}

azimuthCalc.rough <- azimuthCalc <- function(lata, lona, latb, lonb){
  tan.azi <- (latb-lata)/(lonb-lona)
  
  if(lonb >= lona) {
    if(latb >= lata) {atan(tan.azi)} else {2*pi + atan(tan.azi)}
  } else {
    pi + atan(tan.azi)
  }
}