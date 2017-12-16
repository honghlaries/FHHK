longiRange <-  seq(from = 119.9, to = 121.8, length.out = 25)
latiRange <-  seq(from = 33.7, to = 34.9, length.out = 25)
dat.grid.resample <- data.frame(lat = c(1), lon = c(1))
for (i in 1:25) {
  for (j in 1:25) {
    #if((longiRange[j] - 119.9) * (latiRange[i] - 33.7) - (longiRange[j] - 120.6) * (latiRange[i] - 34.5) > 0) {
    dat.grid.resample <- rbind(dat.grid.resample, c(latiRange[i],longiRange[j]))
    #}
  }
}
dat.grid.resample <- dat.grid.resample[-1,]
coordinates(dat.grid.resample) <- ~lon+lat