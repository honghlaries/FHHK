landliplat <- 34.2762; landliplon <- 120.2609
land1lat <- 34.49081; land1lon <- 119.7353
land2lat <- 32.77474; land2lon <- 120.9608

longiRange <-  seq(from = 119.9, to = 121.8, length.out = 250)
latiRange <-  seq(from = 33.7, to = 34.9, length.out = 250)
dat.grid <- data.frame(lat = c(1), lon = c(1))
for (i in 1:250) {
  for (j in 1:250) {
    if((((longiRange[j] - landliplon) * (latiRange[i] - land1lat) - 
         (longiRange[j] - land1lon) * (latiRange[i] - landliplat)) < 0) ||
       (((longiRange[j] - landliplon) * (latiRange[i] - land2lat) -
         (longiRange[j] - land2lon) * (latiRange[i] - landliplat)) > 0) ){
    dat.grid <- rbind(dat.grid, c(latiRange[i],longiRange[j]))
    }
  }
}
dat.grid <- dat.grid[-1,]
coordinates(dat.grid) <- ~lon+lat