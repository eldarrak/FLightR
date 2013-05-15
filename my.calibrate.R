
my.calibrate<-function(Time.adjusted.stable, start.Lon, start.Lat, plot=T) {
	Start.Angles<-elevation(start.Lon, start.Lat, solar(Time.adjusted.stable))
	Angle<-c(Mean=mean(Start.Angles), Sd=sd(Start.Angles))
	if (plot) {
	#par(mfrow=c(2,1))
	plot(Start.Angles)
	hist(Start.Angles)
	#par(mfrow=c(1,1))
	}
	return(Angle)
}
