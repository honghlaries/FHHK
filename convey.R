## Initialization
rm(list = ls())
source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R");source("uniTls_csv2latex.R");

pkgInitialization(c("dplyr","tidyr"))

## Functions 
sigNotation <- function(x) {
  if (is.na(x)) NA else if (x < 0.001) "***" else if (x < 0.01) "**" else if (x < 0.05) "*" else x
}

## Example
taglist <- c("Cr","As","Ni","Cu","Pb","Zn","Cd")

for(i in 1:length(taglist)) {
  dat <- read.csv(paste("relation/driver/polyRegression_",taglist[i],".csv",sep = ""))
  dat <- dat[,-6:-8]
  dat <- dplyr::transmute_(dat, "parameter" = "parameter",
                           "Estimate" = "signif(Estimate,3)",
                           "Std..Error" = "signif(Std..Error,3)",
                           "t.value" = "signif(t.value,3)",
                           "p.Estimate." = "signif(p.Estimate.,2)",
                           "F.value" = "signif(F.value,3)",
                           "p.ANOVA." = "signif(p.ANOVA.,2)")
  for(j in 1:length(dat$p.Estimate.)) dat$p.Estimate.[j] <-sigNotation(as.numeric(dat$p.Estimate.[j]) )
  for(j in 1:length(dat$p.Estimate.)) dat$p.ANOVA.[j] <- sigNotation(as.numeric(dat$p.ANOVA.[j]) )
  colnames(dat)[c(5,7)] <- c("p(Estimate)","p(ANOVA)")
  dat$parameter <- as.character(dat$parameter)
  conveyLaTex(dat,paste("relation/driver/polyRegression_",taglist[i],"_latexcode.txt",sep = ""))
}