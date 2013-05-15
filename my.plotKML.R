
light.r.save.kml<-function(all.out, file.name="result.kml") {
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





