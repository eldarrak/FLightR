
get.optimization.SD<-function(all.out.old, all.out, return.values=F) {
  require("sp")
  # this function will estimate SD between iterations.
  Coords<-cbind(all.out.old$Final.Means$CENTRE.x, all.out.old$Final.Means$CENTRE.y, all.out$Final.Means$CENTRE.x, all.out$Final.Means$CENTRE.y)
  Dists<-apply(Coords, 1, FUN=function(x) {sp:::spDistsN1(cbind(x[[1]], x[[2]]), cbind(x[[3]], x[[4]]), longlat=TRUE)})
  if(return.values) {return(Dists)
  } else {return(sd(Dists))}
}
