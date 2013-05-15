
plot.optimisation.results<-function(object, all.arrays.object, add=F, col="red", map.fill=T, type="mean", pch=3) {
 # simple plotting function
  if (as.character(type)=="mean") {
    Track.coords<-t(apply(object,2, FUN=function(x) cbind(mean( all.arrays.object$Points.Land[x,1]), mean(all.arrays.object$Points.Land[x,2]))))}
  else {
    cat("plotting medians\n")
    Track.coords<-t(apply(object,2, FUN=function(x) cbind(median( all.arrays.object$Points.Land[x,1]), median(all.arrays.object$Points.Land[x,2]))))
  }
  if (add) {lines(Track.coords, type="l", col=col)} 
  else {
    plot(Track.coords, type="l", col=col)
   data("wrld_simpl", package="maptools")

    plot(wrld_simpl, add=T, col=ifelse(map.fill, "lightgray" ,0))
  }
  points(Track.coords, pch=pch, col=col)
  lines(Track.coords, type="l", col=col)
}

