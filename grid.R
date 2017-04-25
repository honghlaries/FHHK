longiRange <-  seq(from = 119.9, to = 121.8, length.out = 125)
latiRange <-  seq(from = 33.7, to = 34.9, length.out = 125)
dat.grid <- data.frame(lat = c(1), lon = c(1))
for (i in 1:125) {
  for (j in 1:125) {
    #if((longiRange[j] - 119.9) * (latiRange[i] - 33.7) - (longiRange[j] - 120.6) * (latiRange[i] - 34.5) > 0) {
    dat.grid <- rbind(dat.grid, c(latiRange[i],longiRange[j]))
    #}
  }
}
dat.grid <- dat.grid[-1,]
coordinates(dat.grid) <- ~lon+lat