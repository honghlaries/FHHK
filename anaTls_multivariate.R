source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R");

circleFun <- function(center = c(0,0), r = 1, npoints = 100){
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  data.frame(x = xx, y = yy)
}

pcaLoadingCal <- function(dat, grouped = T, log = T) { ## cal and log the pca loading
  pkgLoad("dplyr");pkgLoad("tidyr");pkgLoad("psych")
  dat1 <- dat %>%
    filter(isComplete) %>%
    select(Al:sand,depth)
  png(paste(dirPreset("pca"),"/screenPlot_","all",".png",sep = ""))
  fa.parallel(dat1,fa="pc", n.iter=10, show.legend=T,
              main = paste("Screen plot with parallel analysis"))
  dev.off()
  pca <- principal(dat1, nfactors=3, scores = F, rotate="varimax")
  tag <- dimnames(pca$loadings)[[1]]
  loadings <- matrix(as.numeric(pca$loadings), ncol = 3, 
                     dimnames = list(tag, c("PC1","PC2","PC3"))) 
  rst <- data.frame(loadings,tag, group = "all")
  if (grouped) {
    grouptag <- levels(dat$group)
    for(i in 1:length(grouptag)) {
      dat1 <- dat %>% 
        filter(group == grouptag[i]) %>%
        filter(isComplete) %>%
        select(Al:sand,depth)
      png(paste(dirPreset("pca"),,"/screenPlot_",grouptag[i],".png",sep = ""))
      fa.parallel(dat1,fa="pc", n.iter=10, show.legend=T,
                  main = paste("Screen plot with parallel analysis\n, group =",grouptag[i]))
      dev.off()
      pca <- principal(dat1, nfactors=2, scores = T, rotate="varimax")
      tag <- dimnames(pca$loadings)[[1]]
      loadings <- matrix(as.numeric(pca$loadings), ncol = 2, 
                         dimnames = list(tag, c("PC1","PC2"))) 
      rst <- rbind(rst, data.frame(loadings,tag, group = grouptag[i]))
    }
  }
  if(log) {
    write.csv(rst, paste(dirPreset("pca"),"/pcaLoading.csv",sep = ""))
  }
  rst
}

pcaLoadingPlot <- function(dat, grouped = T, themeset, suffix = "", emphtag = NA){ ## output the pca loading plot
  pkgLoad("dplyr");pkgLoad("tidyr");pkgLoad("ggplot2")
  if (grouped) {
    dat <- dat%>% filter(group != "all")
  } else {
    dat <- dat%>% filter(group == "all")
  }
  zeros <- data.frame(PC1 = 0, PC2 = 0, tag = dat$tag, group = dat$group)
  dat2 <- rbind(zeros,dat)
  p <- ggplot() +
    geom_path(aes(x = PC1, y = PC2, group = tag, col = tag, alpha = tag, size = tag),
              arrow = arrow(angle = 15, length = unit(0.10, "inches"),
                            ends = "last", type = "open"),
              data = dat2) +
    geom_path(aes(x = x, y = y),col = "black", size = 0.7, linetype = 2, data = circleFun())+
    geom_text(aes(x = 1.1 * PC1, y = 1.1 * PC2, label = tag, col = tag, alpha = tag),
              check_overlap = F,data = dat) +
    facet_wrap(~ group) +
    xlim(-1.1, 1.1) +ylim(-1.1, 1.1) + 
    themeset
  taglv <- levels(dat$tag); ntaglv <- length(taglv)
  alphaSet <- rep(if(grouped) 0.3 else 1, ntaglv); sizeSet <- rep(0.5, ntaglv); colSet <- rep("black", ntaglv)
  if(!is.na(emphtag)) {
    for(i in 1:length(emphtag)) {
      alphaSet[taglv == emphtag[i]] <- 1
      sizeSet[taglv == emphtag[i]] <- 0.75
      colSet[taglv == emphtag[i]] <- "red"
    }
  }
  p + scale_alpha_manual(breaks = taglv, values = alphaSet) +
    scale_size_manual(breaks = taglv, values = sizeSet) +
    scale_color_manual(breaks = taglv, values = colSet) 
}

pairPlot <- function(dat, xtag, ytag, themeset, suffix = "") {
  p <- ggplot(aes_string(x = xtag, y = ytag, col = "group"), data = dat) + 
    geom_point() + 
    geom_smooth(aes(fill = group), method = "lm", col = "black") + 
    facet_wrap(~ group, scales = "free") + 
    scale_color_manual("Location", breaks = c("CL","EA","NV","WE"), 
                       values = c("#B45F04","#31B404","grey50","#013ADF"))+
    scale_fill_manual("Location", breaks = c("CL","EA","NV","WE"), 
                      values = c("#B45F04","#31B404","grey50","#013ADF"))+
    themeset
  ggsave(plot = p,
         filename = paste("pca/plot/pair_", ytag, "_", xtag, suffix, ".png", sep = ""),
         width = 6, height = 6, dpi = 600)
}

rdaLoadingCal <- function(env, trait, samptag, dir, log = T) { 
  pkgLoad("dplyr");pkgLoad("tidyr");pkgLoad("vegan")
  
  null <- rda(trait~1, env, scale = T) 
  full <- rda(trait~., env, scale = T) 
  mod <- step(null, scope = formula(full), test = "perm")
  
  aov <- anova.cca(mod,permutations = how(nperm=9999))
  print(aov)
  # plot(mod)
  traitload <- mod$CCA$v[,1:2] 
  traitload <- data.frame(traitload,traittag = row.names(traitload))
  
  envload <- mod$CCA$biplot[,1:2] 
  envload <- data.frame(envload,envtag = row.names(envload))
  
  sampload <- mod$CCA$wa[,1:2]
  sampload <- data.frame(sampload,samptag)
  
  if(log) {
    write.csv(format(aov), paste(dirPreset(dir), "/rdaANOVA.csv", sep = ""))
    write.csv(traitload, paste(dirPreset(dir), "/rdaTraitLoad.csv", sep = ""))
    write.csv(envload, paste(dirPreset(dir), "/rdaEnvLoad.csv", sep = ""))
    write.csv(sampload, paste(dirPreset(dir), "/rdaSampLoad.csv", sep = ""))
  }
  list(traitload = traitload,
       envload = envload,
       sampload = sampload) -> dat
}

rdaLoadingPlot <- function(dat) {
  pkgLoad("ggplot2")
  
  traitload <- dat$traitload
  traitload <- rbind(data.frame(RDA1 = 0, RDA2 = 0, traittag = dat$traitload$traittag),traitload)
  ntrait <- length(row.names(traitload))
  
  envload <- dat$envload
  envload <- rbind(data.frame(RDA1 = 0, RDA2 = 0, envtag = dat$envload$envtag),envload)
  nenv <- length(row.names(envload))
  
  sampload <- dat$sampload
  
  ggplot() +
    geom_path(aes(x = x, y = y),col = "black", size = 0.7, linetype = 2, 
              data = circleFun()) +
    geom_point(aes(x = RDA1,y = RDA2), size = 1.5, color = "grey50", shape = 1,
               data = sampload) +
    geom_path(aes(x = RDA1,y = RDA2, group = envtag), size = 0.7,
              data = envload, col = "black") +
    geom_path(aes(x = RDA1,y = RDA2, group = traittag), size = 0.7,
              data = traitload, col = "blue") +
    geom_label(aes(x = RDA1,y = RDA2, label = envtag), size = 3.5, 
               data = envload[nenv/2+1:nenv,], col = "black") + 
    geom_label(aes(x = RDA1,y = RDA2, label = traittag), size = 3.5, 
               data = traitload[ntrait/2+1:ntrait,], col = "blue") +
    scale_x_continuous(name = "RDA1", limits = c(-1.1,1.1)) +
    scale_y_continuous(name = "RDA2", limits = c(-1.1,1.1)) +
    theme_bw() + 
    theme(aspect.ratio = 1,
          legend.position = "none",
          panel.grid = element_blank())
    
}