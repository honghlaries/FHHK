conveyLaTex <- function(dat,filename) {
  dat <- as.data.frame(dat)
  cat("\\toprule \n", file = "relation/test.txt",append = F)
  cat(colnames(dat), file = "relation/test.txt", sep = " & ",append = T)
  cat(" \\\\ \n", file = "relation/test.txt",append = T)
  cat("\\midrule \n", file = "relation/test.txt",append = T)
  for(i in 1:length(rownames(dat))) {
    cat(unlist(dat[i,]), file = "relation/test.txt", sep = " & ",append = T)
    cat(" \\\\ \n", file = "relation/test.txt",append = T)
  }
  cat("\\bottomrule \n", file = "relation/test.txt",append = T)
}



