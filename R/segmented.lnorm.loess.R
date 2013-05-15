#' this function performs lognormal segmented loess
#' @param raw.X vector of days

# this is the main log normal function
# some help file is need for it
segmented.lnorm.loess<-function(raw.X, raw.Y, discont.Segments, span.correction=15, dusk=T, window.size=5, cpus=1, maxiter=100, use.first.interval=F) {
	# this version will try to add weights depending on direction
	if (cpus>1) {
	require("parallel")
    mycl <- parallel:::makeCluster(cpus)
    parallel:::clusterSetRNGStream(mycl)
    ### we don' need to send all parameters to node. so keep it easy..
    parallel:::clusterExport(mycl,c("my.ralling.function", "ll.dlnorm3", "dlnorm3"))
    #parallel:::clusterEvalQ(mycl, library("truncnorm")) 
} else {mycl=NA}
	raw.X<-as.numeric(raw.X)
	raw.Y.old<-raw.Y
	if (dusk) raw.Y<--raw.Y+max(raw.Y)+1

	Loess<-list(shape=NULL, scale=NULL, fit=NULL, n=NULL, reversed=NULL)
	# part 1 estimating pure loess
		for (i in 1:nrow(discont.Segments)) {
		span.local<-min(1,span.correction/(discont.Segments[i,2]-discont.Segments[i,1]))
		if(use.first.interval & i>1) {
			Loess.local<-log.normal.loess(raw.X[discont.Segments[i,1]:discont.Segments[i,2]], raw.Y[discont.Segments[i,1]:discont.Segments[i,2]], span=span.local, window.size=window.size, cpus=cpus, cl=mycl, maxiter=maxiter, fixed.from.first=Fixed.from.first)
		} else {
			Loess.local<-log.normal.loess(raw.X[discont.Segments[i,1]:discont.Segments[i,2]], raw.Y[discont.Segments[i,1]:discont.Segments[i,2]], span=span.local, window.size=window.size, cpus=cpus, cl=mycl, maxiter=maxiter, fixed.from.first=NULL)
		}
		if (dusk) {
			Loess.local.fit<--Loess.local$thres+max(raw.Y.old)+1
		} else {Loess.local.fit<-Loess.local$thres}
		Loess<-list(shape=c(Loess$shape, Loess.local$shape), scale=c(Loess$scale, Loess.local$scale), fit=c(Loess$fit, Loess.local.fit), se.fit=c(Loess$se.fit, Loess.local$se.fit), n=c(Loess$n, rep(nrow(Loess.local), nrow(Loess.local))), reversed=c(Loess$reversed, rep(ifelse(dusk, T, F), nrow(Loess.local))))
		#print(str(Loess))
		if (use.first.interval & i==1) {Fixed.from.first=c(fixed.shape=unique(Loess.local$shape), fixed.scale=unique(Loess.local$scale))
		print(Fixed.from.first)
		}
	}
	if (cpus>1) stopCluster(mycl)
	Loess<-as.data.frame(Loess)
	Loess$Border<-0
	Loess$Border[as.vector(discont.Segments)]<-1
	return(Loess)
}	
