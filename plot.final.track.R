# new type of plot.optim.results..
plot.final.track<-function(all.out, add=F, col="red", map.fill=T, pch=3) {  
  # add mean line:
  if (add) {
  lines(all.out$Final.Means$CENTRE.y~all.out$Final.Means$CENTRE.x, pch=3, col=col, lwd=2) } else {
  data("wrld_simpl", package="maptools")
  plot(all.out$Final.Means$CENTRE.y~all.out$Final.Means$CENTRE.x, pch=3, col=col, lwd=2, type="l")
  plot(wrld_simpl, add=T, col=ifelse(map.fill, "lightgray" ,0))
  lines(all.out$Final.Means$CENTRE.y~all.out$Final.Means$CENTRE.x, pch=3, col=col, lwd=2)
  }
  points(all.out$Final.Means$CENTRE.y~all.out$Final.Means$CENTRE.x, pch=pch, col=col)
}