landliplon <- 34.2762; landliplat <- 120.2609
land1lon <- 34.49081; land1lat <- 119.7353
land2lon <- 32.77474; land2lat <- 120.9608

longiRange <-  seq(from = 119.9, to = 121.8, length.out = 125)
latiRange <-  seq(from = 33.7, to = 34.9, length.out = 125)
dat.grid <- data.frame(lat = c(1), lon = c(1))
for (i in 1:125) {
  for (j in 1:125) {
    if(((longiRange[j] - landliplon) * (latiRange[i] - land1lat) - 
        (longiRange[j] - land1lon) * (latiRange[i] - landliplat)) > 0 ||
       ((longiRange[j] - landliplon) * (latiRange[i] - land2lat) - 
        (longiRange[j] - land2lon) * (latiRange[i] - landliplat)) > 0 ){
    dat.grid <- rbind(dat.grid, c(latiRange[i],longiRange[j]))
    }
  }
}
dat.grid <- dat.grid[-1,]
coordinates(dat.grid) <- ~lon+lat