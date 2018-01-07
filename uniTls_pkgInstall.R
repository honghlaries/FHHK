### =============================================
### Universal Tools: Install the required package
### HongHL 
###
### =============================================

pkgInitialization <- function(pkgs) {
  lapply(pkgs, FUN = function(pkg){
    if(!require(pkg, character.only =T)) {
      install.packages(pkg); paste("package installed:", pkg)
    } else {
      paste("package exists:", pkg)
    }
  }) -> out
  lapply(pkgs, FUN = function(pkg) {require(pkg, character.only =T)})
  cat(paste(unlist(out), "\n", sep = ""))
}

pkgLoad <- function(pkg,n.quit=3) {
  if(!require(pkg, character.only =T)) {
    for (i in 1:n.quit) {
      install.packages(pkg)
      if(require(pkg, character.only =T)) {
        cat(paste("package installed:", pkg));break
      }
    }
    if(!require(pkg, character.only =T)) stop("package installization failed")
  }
}