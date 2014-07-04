#===========
# i want to correct the function not to run any loess...
make.result.list<-function(Data, raw.X, raw.Y) {
	#Loess<-segmented.lnorm.loess(raw.X, raw.Y, Segments, span.correction=span.correction, dusk=dusk, window.size=window.size, cpus=cpus, maxiter=maxiter, use.first.interval=use.first.interval)
	#Int<-c(min(c(raw.Y, Loess$fit)), max(c(raw.Y, Loess$fit)))
	Int<-c(min(raw.Y), max(raw.Y))
	#plot(raw.Y.dusk, ylim=Int)
	#lines(Loess$fit, col="red", lwd=2)
	Index<-sapply(raw.X,  FUN=function(x) which.min(abs(as.numeric(Data$d$gmt-x))))
	Res<-Data$d[Index,]

	# now I want to adjust time and light
	Res$light<-approx(x=Data$d$gmt, y=Data$d$light, xout=raw.X)$y
	Res$gmt<-raw.X
	Res$Hour<-raw.Y
	
	Result<-list(Data=Res)
	#Result<-list(Data=Res, Loess.predict=Loess)
	#Adjust<-(raw.Y-Loess$fit)*3600
	#Result$Data$gmt.adj<-Result$Data$gmt-Adjust
	Result$Data$gmt.adj<-Result$Data$gmt
	Result$Data$gmt<-as.POSIXct(Result$Data$gmt, tz="UTC", origin="1970-01-01")
	Result$Data$gmt.adj<-as.POSIXct(Result$Data$gmt.adj, tz="UTC", origin="1970-01-01")
	return(Result)
}




# here I need a function that will correct results according to segemtned lnormloess..
run.segmented.lnorm.loess<-function(Data, raw.X, raw.Y,  Segments, span.correction=15, dusk=T, window.size=9, cpus=1, maxiter=100, use.first.interval=T, plot=T) {
	Loess<-segmented.lnorm.loess(raw.X, raw.Y, Segments, span.correction=span.correction, dusk=dusk, window.size=window.size, cpus=cpus, maxiter=maxiter, use.first.interval=use.first.interval)
	Int<-c(min(c(raw.Y, Loess$fit)), max(c(raw.Y, Loess$fit)))
	#plot(raw.Y.dusk, ylim=Int)
	#lines(Loess$fit, col="red", lwd=2)
	Index<-sapply(raw.X,  FUN=function(x) which.min(abs(as.numeric(Data$d$gmt-x))))
	Res<-Data$d[Index,]

	# now I want to adjust time and light
	Res$light<-approx(x=Data$d$gmt, y=Data$d$light, xout=raw.X)$y
	Res$gmt<-raw.X
	Res$Hour<-raw.Y
	
	if (plot) {
		plot(raw.Y~raw.X, col="orange", type="p", pch=3,) 
		lines(Loess$fit~Res$gmt, col="darkgreen", lwd=2)
	}
	Result<-list(Data=Res, Loess.predict=Loess)
	Adjust<-(raw.Y-Loess$fit)*3600
	Result$Data$gmt.adj<-Result$Data$gmt-Adjust
	Result$Data$gmt<-as.POSIXct(Result$Data$gmt, tz="UTC", origin="1970-01-01")
	Result$Data$gmt.adj<-as.POSIXct(Result$Data$gmt.adj, tz="UTC", origin="1970-01-01")
	return(Result)
}