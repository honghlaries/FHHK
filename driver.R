## Initialization
rm(list = ls())
source("constant.R");
pkgInitialization(c("dplyr","tidyr","ggplot2","gridExtra","MASS"))
dirInitialization(c("driver"))

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

taglist <- c("Cr","Ni","Cu","Pb","Zn","Cd")

for(i in 1:length(taglist)) {
  relationPlot.gather(dat,taglist[i]) -> p
  ggsave(plot = p,
         filename = paste(dirPreset("driver"),"/",taglist[i],".png",sep = ""),dpi = 600)
  p
}

# factor reg
plot.orgC.Cr <- relationPlot(dat, "orgC", "Cr", "green", "poisson")
plot.clay.Cr <- relationPlot(dat, "clay", "Cr", "black", "poisson")
plot.dep.Cr <- relationPlot(dat, "depth", "Cr", "blue", "poisson")
plot.fe.Cr <- relationPlot(dat, "Fe", "Cr", "purple", "poisson")
plot.h.Cr <- relationPlot(dat, "pH", "Cr", "orange", "poisson")
plot.dist.Cr <- relationPlot(dat, "distance", "Cr", "brown", "poisson")

plot.orgC.Ni <- relationPlot(dat, "orgC", "Ni", "green", "poisson")
plot.clay.Ni <- relationPlot(dat, "clay", "Ni", "black", "poisson")
plot.dep.Ni <- relationPlot(dat, "depth", "Ni", "blue", "poisson")
plot.fe.Ni <- relationPlot(dat, "Fe", "Ni", "purple", "poisson")
plot.h.Ni <- relationPlot(dat, "pH", "Ni", "orange", "poisson")
plot.dist.Ni <- relationPlot(dat, "distance", "Ni", "brown", "poisson")

plot.orgC.Cu <- relationPlot(dat, "orgC", "Cu", "green", "poisson")
plot.clay.Cu <- relationPlot(dat, "clay", "Cu", "black", "poisson")
plot.dep.Cu <- relationPlot(dat, "depth", "Cu", "blue", "poisson")
plot.fe.Cu <- relationPlot(dat, "Fe", "Cu", "purple", "poisson")
plot.h.Cu <- relationPlot(dat, "pH", "Cu", "orange", "poisson")
plot.dist.Cu <- relationPlot(dat, "distance", "Cu", "brown", "poisson")

plot.orgC.Pb <- relationPlot(dat, "orgC", "Pb", "green", "poisson")
plot.clay.Pb <- relationPlot(dat, "clay", "Pb", "black", "poisson")
plot.dep.Pb <- relationPlot(dat, "depth", "Pb", "blue", "poisson")
plot.fe.Pb <- relationPlot(dat, "Fe", "Pb", "purple", "poisson")
plot.h.Pb <- relationPlot(dat, "pH", "Pb", "orange", "poisson")
plot.dist.Pb <- relationPlot(dat, "distance", "Pb", "brown", "poisson")

plot.orgC.Zn <- relationPlot(dat, "orgC", "Zn", "green", "poisson")
plot.clay.Zn <- relationPlot(dat, "clay", "Zn", "black", "poisson")
plot.dep.Zn <- relationPlot(dat, "depth", "Zn", "blue", "poisson")
plot.fe.Zn <- relationPlot(dat, "Fe", "Zn", "purple", "poisson")
plot.h.Zn <- relationPlot(dat, "pH", "Zn", "orange", "poisson")
plot.dist.Zn <- relationPlot(dat, "distance", "Zn", "brown", "poisson")

plot.orgC.Cd <- relationPlot(dat, "orgC", "Cd", "green", "poisson")
plot.clay.Cd <- relationPlot(dat, "clay", "Cd", "black", "poisson")
plot.dep.Cd <- relationPlot(dat, "depth", "Cd", "blue", "poisson")
plot.fe.Cd <- relationPlot(dat, "Fe", "Cd", "purple", "poisson")
plot.h.Cd <- relationPlot(dat, "pH", "Cd", "orange", "poisson")
plot.dist.Cd <- relationPlot(dat, "distance", "Cd", "brown", "poisson")

p.gather <- grid.arrange(plot.clay.Pb, plot.dep.Pb, plot.h.Pb, plot.fe.Pb, plot.dist.Pb, plot.orgC.Pb, 
                         plot.clay.Cr, plot.dep.Cr, plot.h.Cr, plot.fe.Cr, plot.dist.Cr, plot.orgC.Cr, 
                         plot.clay.Ni, plot.dep.Ni, plot.h.Ni, plot.fe.Ni, plot.dist.Ni, plot.orgC.Ni, 
                         plot.clay.Cu, plot.dep.Cu, plot.h.Cu, plot.fe.Cu, plot.dist.Cu, plot.orgC.Cu, 
                         plot.clay.Zn, plot.dep.Zn, plot.h.Zn, plot.fe.Zn, plot.dist.Zn, plot.orgC.Zn, 
                         plot.clay.Cd, plot.dep.Cd, plot.h.Cd, plot.fe.Cd, plot.dist.Cd, plot.orgC.Cd,
                         nrow=6, ncol=6)
ggsave(plot = p.gather,
       filename = paste(dirPreset("relation/driver"),"/gather_relation.png",sep = ""),
       dpi = 300, width = 12, height = 12)

# factor reg
taglist <- c("Cr","Ni","Cu","Pb","Zn","Cd")

for(i in 1:length(taglist)) {
  mod <- stepFitting(dat, taglist[i])
  mod.sum <- modExamine(mod)
  write.csv(mod.sum, 
            paste(dirPreset("driver"),"/multiRegExplo_",taglist[i],".csv",sep = ""),
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

# sem
pkgLoad("sem");pkgLoad("semPlot");pkgLoad("lavaan");

dat <- datareadln()%>% 
  dplyr::select(depth:sand) %>% 
  dplyr::mutate(acc_depth = 1/depth,
                acc_dist = 1/distance,
                ad_pH = (pH-8.2)^2)

model <- ' 
            # latent variable definitions
              accessibility =~ acc_depth + acc_dist
              adsorbability =~ orgC + AVS + Fe + clay + ad_pH

            # regressions
              Pb ~ accessibility +  adsorbability
              Zn ~ accessibility +  adsorbability
              Cd ~ accessibility +  adsorbability
              Cu ~ accessibility +  adsorbability
              Ni ~ accessibility +  adsorbability
              Cr ~ accessibility +  adsorbability         

           # residual correlations

        '

fit <- lavaan::sem(model, data = scale(dat))

png("driver/SEM.png",width = 1400, height = 800)
semPaths(fit,"std",  style = "lisrel",
         edge.label.cex = 0.5, exoVar = FALSE, exoCov = FALSE)
dev.off()