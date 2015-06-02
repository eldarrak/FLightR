#plotting functions.R

# this functions are not tested yet for consistensy with new data structures..

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

# this function is needed to estimate CI
#################################################################
####
create.credible.intervals.array<-function(in.Data, intervals=c(0.8, 0.6, 0.4, 0.2) ) {
  require(rgeos)
  require(maptools)
  # 
  # intervals should be a vector like c(0.8, 0.6, 0.4)
  intervals=sort(intervals, decreasing=T)
  Polygons<-vector(mode="list", length=length(intervals))
  names(Polygons)<-paste("ci.", intervals, sep="")
  nParticles=sum(in.Data$Points.rle[[1]]$lengths) #dim(All.results.mat)[1]
  Seq<-intervals 
  Polygons<-lapply(Polygons, list)
  cat("estimating rle..")
  for (Iteration in 1:length(in.Data$Points.rle)) {
    #cat("Iteration", Iteration, "\n")
    points2read<-inverse.rle(in.Data$Points.rle[[Iteration]])
    Density<-cumsum(sort.int(table(points2read), decreasing=T))/nParticles
    for (int in 1:length(Seq)) {
      Pointstopoly<-as.integer(names(Density)[Density<=Seq[int]])
      if (length(Pointstopoly)==0) {
        Pointstopoly<-as.integer(names(Density)[1])}
      Pol<-Polygons(list(Polygon(in.Data$Points.Land[c(Pointstopoly, Pointstopoly[1]), ])), ID=Iteration)
      Polygons[[int]][[Iteration]]<-Pol
    }
  }
  cat("   Done!\n")
  Convexed<-list()
  ###
  # now we need to combine points by 2..
  All.cr<-list()
  for (cred in 1:length(intervals)) {
    Convexed<-list()
    cat(intervals[cred], "\n")
    Curr.cr<-Polygons[[cred]]
    for (i in 1:(length(Curr.cr)-1)) {
      Spat.Polygons<-lapply(Curr.cr, FUN=function(x) SpatialPolygons(list(x), proj4string=CRS("+proj=longlat +datum=WGS84")))
      if(nrow(Spat.Polygons[[i]]@polygons[[1]]@Polygons[[1]]@coords)<=3) 	{ Spat.Polygons[[i]]@polygons[[1]]@Polygons[[1]]@coords<-rbind(Spat.Polygons[[i]]@polygons[[1]]@Polygons[[1]]@coords, Spat.Polygons[[i]]@polygons[[1]]@Polygons[[1]]@coords, Spat.Polygons[[i]]@polygons[[1]]@Polygons[[1]]@coords)
      }
      if(nrow(Spat.Polygons[[i+1]]@polygons[[1]]@Polygons[[1]]@coords)<=3) 	{ Spat.Polygons[[i+1]]@polygons[[1]]@Polygons[[1]]@coords<-rbind(Spat.Polygons[[i+1]]@polygons[[1]]@Polygons[[1]]@coords, Spat.Polygons[[i+1]]@polygons[[1]]@Polygons[[1]]@coords, Spat.Polygons[[i+1]]@polygons[[1]]@Polygons[[1]]@coords,Spat.Polygons[[i+1]]@polygons[[1]]@Polygons[[1]]@coords)
      }
      if (i==1) {
        Convexed<-gConvexHull(gUnion(gConvexHull(Spat.Polygons[[i]]), gBoundary(Spat.Polygons[[i+1]])))
      }
      else {Convexed<-gUnion(Convexed, gConvexHull(gUnion(gConvexHull(Spat.Polygons[[i]]), gConvexHull(Spat.Polygons[[i+1]]))))}
    }
    All.cr[[cred]]<-Convexed
  }
  names(All.cr)<-names(Polygons)
  return(All.cr)
}

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


FLightR.r.save.kml<-function(all.out, file.name="result.kml") {
	require(plotKML)
	require(fossil)
	require(sp)
	require(spacetime)
	require(adehabitat)
	Data<-data.frame(Time=all.out$Matrix.Index.Table$Real.time, lon=all.out$Final.Means$CENTRE.x[-1], lat=all.out$Final.Means$CENTRE.y[-1])
	coordinates(Data) <- ~lon+lat
	proj4string(Data) <- CRS("+proj=longlat +datum=WGS84")
	xy <- as.list(data.frame(t(coordinates(Data))))
	Data$dist.km <- sapply(xy, function(x) { deg.dist(long1=x[1], lat1=x[2], long2=xy[[1]][1], lat2=xy[[1]][2]) } )
	Data.ltraj <- as.ltraj(coordinates(Data), Data$Time, id = "th")
	Data.st <- as(Data.ltraj, "STTDF")
	Data.st$speed <- Data$speed
	Data.st@sp@proj4string <- CRS("+proj=longlat +datum=WGS84")
	plotKML(Data.st, file.name=file.name)
}

