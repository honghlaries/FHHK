## Initialization
rm(list = ls())
source("uniTls_pkgInstall.R");source("uniTls_presetPaths.R");source("uniTls_csv2latex.R");
pkgInitialization(c("dplyr","tidyr"))

## Functions 

## Example
taglist <- c("Cr","As","Ni","Cu","Pb","Zn","Cd")

for(i in 1:length(taglist)) {
  dat <- read.csv(paste("relation/driver/polyRegression_",taglist[i],".csv",sep = ""))
  conveyLaTex(dat,paste("relation/driver/polyRegression_",taglist[i],"_latexcode.txt",sep = ""))
}