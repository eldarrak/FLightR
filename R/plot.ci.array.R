
plot.ci.array<-function(CI.array, xlim=c(-115, -94), ylim=c(16,38.5), all.out) {
  require(maptools)
  my.golden.colors <- colorRampPalette(
    c("white","#FF7100"),
    bias=2.0)	
  Col=my.golden.colors(1+length(CI.array))
  Col=Col[-1]
  #topo.colors(length(CI.array)+5)
  for (i in 1:length(CI.array)) {
  if(i==1) {
	add=F 
  } else { 
	add=T
  }
  if (class(CI.array[[i]])=="SpatialCollections") {
  
  suppressWarnings(plot(CI.array[[i]]@polyobj, col=Col[i], border=grey(0.5), xlim=xlim, ylim=ylim, add=add))
  suppressWarnings(plot(CI.array[[i]]@lineobj, col=Col[i], border=grey(0.5), xlim=xlim, ylim=ylim, add=add))
 } else  {
      suppressWarnings(plot(CI.array[[i]], col=Col[i], border=grey(0.5), xlim=xlim, ylim=ylim, add=add))
  }
  }
  data("wrld_simpl", package="maptools")
  plot(wrld_simpl, add=T, border=grey(0.8), lwd=2)
  
  # add mean line:
  points(all.out$Final.Means$CENTRE.y~all.out$Final.Means$CENTRE.x, pch=3, cex=0.5, col="red")
  lines(all.out$Final.Means$CENTRE.y~all.out$Final.Means$CENTRE.x, pch=3, col="red", lwd=2)
  
  # release position
  points(-98.7, 34.7, pch=1, cex=2, col="black", lwd=3)
  
  # add bars every month:
  Index<-which(as.POSIXlt(all.out$Matrix.Index.Table$Real.time)$mday==1 &  all.out$Main.Index$proposal.index=="Dusk")
  points(all.out$Final.Means$CENTRE.y[Index]~all.out$Final.Means$CENTRE.x[Index], pch=3, cex=3, col="black", lwd=2)
  
}

