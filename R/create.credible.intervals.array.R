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
    #points2read<-inverse.rle(in.Data$Points.rle[[Iteration]])
    #Density<-cumsum(sort.int(table(points2read), decreasing=T))/nParticles
	Sorted<-sort.int(in.Data$Points.rle[[Iteration]]$length,  index.return=T, decreasing=T)
    Density<-cumsum(Sorted$x)/nParticles
	names(Density)<-in.Data$Points.rle[[Iteration]]$values[Sorted$ix]
	for (int in 1:length(Seq)) {
      Pointstopoly<-as.integer(names(Density)[Density<=Seq[int]])
      if (length(Pointstopoly)==0) {
        Pointstopoly<-as.integer(names(Density)[1])}
      Pol<-Polygons(list(Polygon(in.Data$Points.Land[c(Pointstopoly, Pointstopoly[1]), 1:2])), ID=Iteration)
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
