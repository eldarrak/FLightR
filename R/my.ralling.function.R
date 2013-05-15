# this internal function works inside segmented lognormal loess
my.ralling.function<-function(x) {
	Fixed=c(shape=x$shape$fixed, scale=x$scale$fixed, thres=x$thres$fixed)
	Initials<-c(shape=x$shape$initials, scale=x$scale$initials, thres=x$thres$initials)
	# autostart
	thres.init=min(x$Y.neighb)-0.01
	if (is.na(Initials[1])) Initials[1]<-log(sd(log(x$Y.neighb-thres.init)))
	if (is.na(Initials[2])) Initials[2]<-mean(log(x$Y.neighb-thres.init))
	if (is.na(Initials[3]) | Initials[3] > thres.init) Initials[3] <-thres.init
	#if (is.na(Initials[3])) Initials[3] <-thres.init
	#
	prm<-Initials[which(is.na(Fixed))]
	#cat("Fixed")
	#print(Fixed)
	#cat("prm")
	#print(prm)
	# I also need to make lower and upper dynamic!!!
	lower=c(-Inf,-Inf, min(x$Y.neighb)-1)[which(is.na(Fixed))]
	upper=c(Inf, Inf,min(x$Y.neighb))[which(is.na(Fixed))]
	Res=try(optim(prm, ll.dlnorm3, gr=NULL, x$Y.neighb, Fixed , FALSE, method="L-BFGS-B", control=list(trace=0, maxit=1e6), lower=lower, upper=upper))
	if (class(Res)!="try-error") {
	 return(Res$par)
	 } else {
	 return(rep(NA, length(prm)))
	 stop()
	 }
}