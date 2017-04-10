### =======================================
### Universal Tools: Preset the local paths
### HongHL 
###
### =======================================

dirInitialization <- function(dirs) {
  lapply(dirs, FUN = function(paths){
    if(!dir.exists(paths)) {
      dir.create(paths, recursive = T); paste("creating path:", paths)
    } else {
      paste("path exists:", paths)
    }
  }) -> out
  cat(paste(unlist(out), "\n", sep = ""))
}

dirPreset <- function(dirs, n.quit =3) {
  if(!dir.exists(dirs)) {
    for(i in 1:n.quit) {
      dir.create(dirs, recursive = T)
      if(dir.exists(dirs)) {
        cat(paste("path created:", dirs));break
      }
    }
  }
  if(!dir.exists(dirs)) stop("path creating failed")
  dirs
}

