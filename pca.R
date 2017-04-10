## clean ----
rm(list = ls())
source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R")

## Functions ----
datareadln <- function() { ## data readln
  pkgLoad("dplyr");pkgLoad("tidyr")
  read.csv("data/result_element.csv") %>%
    dplyr::inner_join(read.csv("data/result_grainSize.csv"), by = c("sampleID" = "sampleID")) %>%
    dplyr::right_join(read.csv("data/meta_splList.csv"), by = c("sampleID" = "sampleID")) %>%
    dplyr::inner_join(read.csv("data/meta_sites.csv"), by = c("siteID" = "siteID")) %>%
    dplyr::mutate(orgC = C.ac, AVS = S - S.ac, 
                  isComplete = complete.cases(Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd,C,N,S,orgC,AVS,clay,silt,sand)) %>%
    dplyr::select(siteID:depth,isComplete,Al,Fe,Mn,Pb,Cr,Ni,Cu,Zn,As,Cd,C,N,S,orgC,AVS,clay,silt,sand)
}

pcaLoadingCal <- function(dat, grouped = T, log = T) { ## cal and log the pca loading
  pkgLoad("dplyr");pkgLoad("tidyr");pkgLoad("psych")
  dat1 <- dat %>%
    filter(isComplete) %>%
    select(Al:sand,depth)
  png(paste(dirPreset("pca"),"/screenPlot_","all",".png",sep = ""))
  fa.parallel(dat1,fa="pc", n.iter=10, show.legend=T,
              main = paste("Screen plot with parallel analysis\n, group =","all"))
  dev.off()
  pca <- principal(dat1, nfactors=2, scores = F, rotate="varimax")
  tag <- dimnames(pca$loadings)[[1]]
  loadings <- matrix(as.numeric(pca$loadings), ncol = 2, 
                     dimnames = list(tag, c("PC1","PC2"))) 
  rst.tot <- data.frame(loadings,tag, group = "all")
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
      rst.tot <- rbind(rst.tot, data.frame(loadings,tag, group = grouptag[i]))
    }
  }
  if(log) {
    write.csv(paste(dirPreset("pca"),"/pcaLoading.csv",sep = ""))
  }
  rst.tot
}

pcaLoadingPlot <- function(dat, grouped = T, themeset, suffix = "", emphtag = NA){ ## output the pca loading plot
  pkgLoad("dplyr");pkgLoad("tidyr");pkgLoad("ggplot2")
  circleFun <- function(center = c(0,0), r = 1, npoints = 100){
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
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
  p <- p + scale_alpha_manual(breaks = taglv, values = alphaSet) +
    scale_size_manual(breaks = taglv, values = sizeSet) +
    scale_color_manual(breaks = taglv, values = colSet) 
  ggsave(plot = p,
         filename = paste(dirPreset("pca"),"/pcaLoading", suffix, ".png", sep = ""),
         width = 6, height = 6, dpi = 600)
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

## Example ----
pcaLoading <- pcaLoadingCal(datareadln(), grouped = F)
library(ggplot2)
themeset <- theme_bw() + theme(aspect.ratio = 1, legend.position = "none")
pcaLoadingPlot(pcaLoading, grouped = F, themeset = themeset, suffix = "_a", emphtag = c("depth","silt","AVS"))

#pairPlot(dat = datareadln(), xtag = "orgC", ytag = "Cd", themeset = themeset)
#pairPlot(dat = datareadln(), xtag = "Al", ytag = "Cd", themeset = themeset)
#pairPlot(dat = datareadln(), xtag = "orgC", ytag = "S", themeset = themeset)

