## Initialization
rm(list = ls())
source("constant.R");
pkgInitialization(c("dplyr","tidyr","ggplot2","gridExtra","MASS"))

## Functions 
relationPlot <- function(dat, fact, resp, col, family) {
 
  ggplot(aes_string(x = fact, y = resp),data = dat) +
    geom_point(col = "grey80") +
    geom_smooth(method="loess", se = F, linetype = 2, col = "black") +
    stat_smooth(method="glm", method.args=list(family = family), 
                col = col, fill = col, alpha = 0.3) +
    theme_bw() + 
    theme(panel.grid = element_blank())
}

relationPlot.gather <- function(dat, tag, modfamilies = rep("poisson",6)) {
 
  plot.clay <- relationPlot(dat, "clay", tag, "black",family = modfamilies[1])
  plot.dep <- relationPlot(dat, "depth", tag, "blue",family = modfamilies[2])
  plot.h <- relationPlot(dat, "pH", tag, "orange",family = modfamilies[3])
  plot.fe <- relationPlot(dat, "Fe", tag, "purple",family = modfamilies[4])
  plot.sal <- relationPlot(dat, "salinity", tag, "brown",family = modfamilies[5])
  plot.avs <- relationPlot(dat, "AVS", tag, "green",family = modfamilies[6])
  
  p <- grid.arrange(plot.clay, plot.dep, plot.h, plot.fe, plot.sal, plot.avs, 
                    nrow=2, ncol=3)
  p
}

stepFitting <- function(dat, tag, method = "lm", family) {
  
  if (method == "lm") {
    null <- lm(as.formula(paste(tag,"~1",sep = "")), data = dat) 
    # full mod was based on RDA result
    full <- lm(as.formula(paste(tag,"~ clay * depth * pH * Fe * salinity * AVS", 
                                sep = "")), data = dat) 
    return(step(null, scope = formula(full), test = "F"))
  }
 
  if (method == "glm") {
    #glm <- glm(as.formula(paste(tag,"~clay + depth + pH + Fe + salinity + AVS",
    #                            sep = "")), family = family, data = dat)
    return(MASS::stepAIC(glm, trace = T))
  }
  stop("Wrong Fitting Method.")
  
}

modExamine <- function(mod) {
  
  mod.est <- summary(mod)$coefficients
  mod.est <- cbind (parameter = rownames(mod.est),mod.est)
  colnames(mod.est)[5] <- "p.Estimate."
  mod.est <- as.data.frame(mod.est)
  
  mod.aov <- anova(mod)
  mod.aov <- cbind (parameter = rownames(mod.aov),mod.aov)
  colnames(mod.aov)[6] <- "p.ANOVA."
  mod.aov <- as.data.frame(mod.aov)
  
  full_join(mod.est,mod.aov,by = c("parameter" = "parameter"))
}

## Example
dat <- datareadln() %>%
  dplyr::select(depth,distance,salinity:sand) %>%
  dplyr::mutate(Fe = Fe/10000)

taglist <- c("Cr","As","Ni","Cu","Pb","Zn","Cd")

for(i in 1:length(taglist)) {
  relationPlot.gather(dat,taglist[i]) -> p
  ggsave(plot = p,
         filename = paste(dirPreset("relation/driver"),"/",taglist[i],".png",sep = ""),dpi = 600)
  p
}

plot.avs.Cr <- relationPlot(dat, "AVS", "Cr", "green", "poisson")
plot.clay.Cr <- relationPlot(dat, "clay", "Cr", "black", "poisson")
plot.dep.Cr <- relationPlot(dat, "depth", "Cr", "blue", "poisson")
plot.fe.Cr <- relationPlot(dat, "Fe", "Cr", "purple", "poisson")
plot.h.Cr <- relationPlot(dat, "pH", "Cr", "orange", "poisson")
plot.sal.Cr <- relationPlot(dat, "salinity", "Cr", "brown", "poisson")

plot.avs.As <- relationPlot(dat, "AVS", "As", "green", "poisson")
plot.clay.As <- relationPlot(dat, "clay", "As", "black", "poisson")
plot.dep.As <- relationPlot(dat, "depth", "As", "blue", "poisson")
plot.fe.As <- relationPlot(dat, "Fe", "As", "purple", "poisson")
plot.h.As <- relationPlot(dat, "pH", "As", "orange", "poisson")
plot.sal.As <- relationPlot(dat, "salinity", "As", "brown", "poisson")

plot.avs.Ni <- relationPlot(dat, "AVS", "Ni", "green", "poisson")
plot.clay.Ni <- relationPlot(dat, "clay", "Ni", "black", "poisson")
plot.dep.Ni <- relationPlot(dat, "depth", "Ni", "blue", "poisson")
plot.fe.Ni <- relationPlot(dat, "Fe", "Ni", "purple", "poisson")
plot.h.Ni <- relationPlot(dat, "pH", "Ni", "orange", "poisson")
plot.sal.Ni <- relationPlot(dat, "salinity", "Ni", "brown", "poisson")

plot.avs.Cu <- relationPlot(dat, "AVS", "Cu", "green", "poisson")
plot.clay.Cu <- relationPlot(dat, "clay", "Cu", "black", "poisson")
plot.dep.Cu <- relationPlot(dat, "depth", "Cu", "blue", "poisson")
plot.fe.Cu <- relationPlot(dat, "Fe", "Cu", "purple", "poisson")
plot.h.Cu <- relationPlot(dat, "pH", "Cu", "orange", "poisson")
plot.sal.Cu <- relationPlot(dat, "salinity", "Cu", "brown", "poisson")

plot.avs.Pb <- relationPlot(dat, "AVS", "Pb", "green", "poisson")
plot.clay.Pb <- relationPlot(dat, "clay", "Pb", "black", "poisson")
plot.dep.Pb <- relationPlot(dat, "depth", "Pb", "blue", "poisson")
plot.fe.Pb <- relationPlot(dat, "Fe", "Pb", "purple", "poisson")
plot.h.Pb <- relationPlot(dat, "pH", "Pb", "orange", "poisson")
plot.sal.Pb <- relationPlot(dat, "salinity", "Pb", "brown", "poisson")

plot.avs.Zn <- relationPlot(dat, "AVS", "Zn", "green", "poisson")
plot.clay.Zn <- relationPlot(dat, "clay", "Zn", "black", "poisson")
plot.dep.Zn <- relationPlot(dat, "depth", "Zn", "blue", "poisson")
plot.fe.Zn <- relationPlot(dat, "Fe", "Zn", "purple", "poisson")
plot.h.Zn <- relationPlot(dat, "pH", "Zn", "orange", "poisson")
plot.sal.Zn <- relationPlot(dat, "salinity", "Zn", "brown", "poisson")

plot.avs.Cd <- relationPlot(dat, "AVS", "Cd", "green", "poisson")
plot.clay.Cd <- relationPlot(dat, "clay", "Cd", "black", "poisson")
plot.dep.Cd <- relationPlot(dat, "depth", "Cd", "blue", "poisson")
plot.fe.Cd <- relationPlot(dat, "Fe", "Cd", "purple", "poisson")
plot.h.Cd <- relationPlot(dat, "pH", "Cd", "orange", "poisson")
plot.sal.Cd <- relationPlot(dat, "salinity", "Cd", "brown", "poisson")

p.gather <- grid.arrange(plot.clay.Pb, plot.dep.Pb, plot.h.Pb, plot.fe.Pb, plot.sal.Pb, plot.avs.Pb, 
                         plot.clay.Cr, plot.dep.Cr, plot.h.Cr, plot.fe.Cr, plot.sal.Cr, plot.avs.Cr, 
                         plot.clay.Ni, plot.dep.Ni, plot.h.Ni, plot.fe.Ni, plot.sal.Ni, plot.avs.Ni, 
                         plot.clay.Cu, plot.dep.Cu, plot.h.Cu, plot.fe.Cu, plot.sal.Cu, plot.avs.Cu, 
                         plot.clay.Zn, plot.dep.Zn, plot.h.Zn, plot.fe.Zn, plot.sal.Zn, plot.avs.Zn, 
                         plot.clay.As, plot.dep.As, plot.h.As, plot.fe.As, plot.sal.As, plot.avs.As, 
                         plot.clay.Cd, plot.dep.Cd, plot.h.Cd, plot.fe.Cd, plot.sal.Cd, plot.avs.Cd,
                         nrow=7, ncol=6)
ggsave(plot = p.gather,
       filename = paste(dirPreset("relation/driver"),"/gather_relation.png",sep = ""),
       dpi = 300, width = 12, height = 14)

taglist <- c("Cr","As","Ni","Cu","Pb","Zn","Cd")

for(i in 1:length(taglist)) {
  mod <- stepFitting(dat, taglist[i])
  mod.sum <- modExamine(mod)
  write.csv(mod.sum, 
            paste(dirPreset("relation/driver"),"/multiRegExplo_",taglist[i],".csv",sep = ""),
            row.names = F)
  mod.sum.refined <- mod.sum 
  for(k in 1:10) {
    tag <- mod.sum.refined %>%
      dplyr::filter(parameter != "(Intercept)",
             parameter != "Residuals",
             as.numeric(as.character(p.Estimate.)) < alphalevel,
             as.numeric(as.character(p.ANOVA.)) < alphalevel) %>%
      dplyr::select(parameter)  %>% unlist() %>% as.vector()
    if (length(tag) ==0) next
    formulas <- paste(taglist[i],"~",sep = "")
    for (j in 1:length(tag)) {
      formulas <- paste(formulas,tag[j],sep = "")
      if (j < length(tag)) formulas <- paste(formulas,"+",sep = "")
    }
    formulas <- as.formula(formulas)
    mod.refined <- lm(formulas, data = dat)
    mod.sum.refined <- modExamine(mod.refined)
    
  }
  if(length(tag) !=0) {
    write.csv(mod.sum.refined, 
              paste(dirPreset("relation/driver"),"/multiRegRefined_",taglist[i],".csv",sep = ""),
              row.names = F)
    }
}
