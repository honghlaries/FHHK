conveyLaTex <- function(dat,filename) {
  dat <- as.data.frame(dat)
  cat("\\toprule \n", file = filename,append = F)
  cat(colnames(dat), file = filename, sep = " & ",append = T)
  cat(" \\\\ \n", file = filename,append = T)
  cat("\\midrule \n", file = filename,append = T)
  for(i in 1:length(rownames(dat))) {
    cat(unlist(dat[i,]), file = filename, sep = " & ",append = T)
    cat(" \\\\ \n", file = filename,append = T)
  }
  cat("\\bottomrule \n", file = filename,append = T)
}



