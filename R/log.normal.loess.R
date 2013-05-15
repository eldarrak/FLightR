# this function will run lnorm loess on the track segment	

log.normal.loess<-function(X, Y,  window.size=5, span=0.1, maxiter=100, cpus=1, cl=NA, fixed.from.first=NULL) {

	# ver.0.2 decided to use normalized by previous loess neighbours
	Input.List.Initial<-Input.List<-create.input.list(X, Y, window.size=window.size)
	#Lengths<-sapply(Input.List, FUN=function(x) length(x[["Y.neighb"]]))
	
	plot(Y, pch="+", ylim=c(min(Y)-1, max(Y)+1))

	
	if (is.null(fixed.from.first) ){
	if (cpus>1) {
	Res<-as.data.frame(t(parSapply(cl=cl, Input.List.Initial,  my.ralling.function)))
	} else {
	Res<-as.data.frame(t(sapply(Input.List.Initial, my.ralling.function)))
	}
	names(Res)<-c("shape", "scale", "thres")
	print(summary(Res))
    lines(Res$thres, col="red")
		# now I want to rerun the same thing but with new initials

	# Let's try to exclude parameters at the boundary!
	# they can be discovered as ones that have Res$thres as
	# min(x$Y.neighb)-1 for lower boundary
	# and min(x$Y.neighb) fro upper boundary
	At.the.boundary<-which((Y-1-Res$thres)>=0)
	cat(length(At.the.boundary), "parameters from", nrow(Res), "estimated at the boundary\n")
	Not.at.the.boundary<-unique(c(1, (1:nrow(Res))[!(1:nrow(Res)) %in% At.the.boundary], nrow(Res)))
	#Lengths=Lengths/max(Lengths)
	Mean.shape.old<-Mean.shape<-mean(Res$shape[Not.at.the.boundary], na.rm=T)
	Mean.scale.old<-Mean.scale<-mean(Res$scale[Not.at.the.boundary], na.rm=T)
	#Loess.thres.old<-Loess.thres<-predict(loess(Res$thres[Not.at.the.boundary]~X[Not.at.the.boundary], span=span, weights=Lengths[Not.at.the.boundary]), newdata=X, se=T)$fit
	Loess.thres.old<-Loess.thres<-predict(loess(Res$thres[Not.at.the.boundary]~X[Not.at.the.boundary]), span=span, newdata=X, se=T)$fit
	Res.old<-Res
	Diff.shape<-Inf
	Diff.scale<-Inf
	Diff.thres<-Inf
	Counter=1
	Stop=F
	cat("       estimating initial parameter values\n")

	# for bayesian means
	Mean.shape.Sum<-Mean.shape
	Mean.scale.Sum<-Mean.scale
	Mean.Loess.thres<-Loess.thres
	
	while ((abs(Diff.shape) > 0.05 | abs(Diff.shape)>0.05 |abs(Diff.thres)>0.01 ) & !Stop) {
	#if( Counter > maxiter ) stop("maximum amount of iterations reached!\n")
	cat("Iteration", Counter, "\n")
	# now we need to update paramters in the Input.List
	Input.List<-Input.List.Initial
	for (i in 1:length(Input.List)) {
		Input.List[[i]]$shape$initials<-Mean.shape.old
		Input.List[[i]]$scale$initials<-Mean.scale.old
		Input.List[[i]]$thres$initials<-Loess.thres.old[i]
		Index<-ceiling(i-(window.size/2)):floor(i+(window.size/2))
		# i want to extract thres from the neighbours.
		# and take it from the points..
		#Input.List[[i]]$Y.neighb<-Input.List[[i]]$Y.neighb-Loess.thres.old[which(1:length(Loess.thres.old) %in% Index)]	
		#Input.List[[i]]$thres$initials<-0
	}
	if (cpus>1) {
	Res<-as.data.frame(t(parSapply(cl=cl, Input.List,  my.ralling.function)))
	} else {
	Res<-as.data.frame(t(sapply(Input.List, my.ralling.function)))
	}
	names(Res)<-c("shape", "scale", "thres")
	print(summary(Res))

	#Res$thres<-Res$thres+Loess.thres.old
	
	if (is.null(fixed.from.first)) lines(Res.old$thres, col="grey")
	lines(Res$thres, col="red")
	lines(Loess.thres.old, col="yellow")
	# estimated shape:
	#Estimated.shape<-Res$shape[which(Res$shape>1e-0)]
	#cat("estimated ", length(Estimated.shape), "parameters from", length(Res$shape), "\n")
	At.the.boundary<-which((Y-1-Res$thres)>=0)
	cat(length(At.the.boundary), "parameters from", nrow(Res), "estimated at the boundary\n")
	Not.at.the.boundary<-unique(c(1, (1:nrow(Res))[!(1:nrow(Res)) %in% At.the.boundary], nrow(Res)))
	
	Mean.shape<-mean(Res$shape[Not.at.the.boundary], na.rm=T)
	#Mean.shape<-mean(Estimated.shape, na.rm=T)
	
	#Mean.scale<-mean(Res$scale[which(Res$shape>1e-0)], na.rm=T)
	Mean.scale<-mean(Res$scale[Not.at.the.boundary], na.rm=T)
	#Loess.thres<-predict(loess(Res$thres[Not.at.the.boundary]~X[Not.at.the.boundary], span=span, weights=Lengths[Not.at.the.boundary]), newdata=X, se=T)$fit
	Loess.thres<-predict(loess(Res$thres[Not.at.the.boundary]~X[Not.at.the.boundary]), span=span, newdata=X, se=T)$fit
	
	lines(Loess.thres, col="blue")
	Diff.shape<-Mean.shape-Mean.shape.old
	Diff.scale<-Mean.scale-Mean.scale.old
	Diff.thres<-sum(abs((Loess.thres-Loess.thres.old)/length(Loess.thres)))
	cat("Mean.shape: ", Mean.shape, "; Means scale ", Mean.scale, "\n" )

	cat("Diff.shape", Diff.shape, "Diff.scale" , Diff.scale, "Diff.thres", Diff.thres, "\n")

	Mean.shape.old<-Mean.shape
	Mean.scale.old<-Mean.scale
	Loess.thres.old<-Loess.thres
	Res.old<-Res
	
	Mean.shape.Sum<-Mean.shape.Sum+Mean.shape
	Mean.scale.Sum<-Mean.scale.Sum+Mean.scale
	Mean.Loess.thres<-Mean.Loess.thres+Loess.thres
	Counter=Counter+1
	if (Counter>maxiter) {
	Stop=T 
	Mean.shape.old<-Mean.shape.Sum/Counter
	Mean.scale.old<-Mean.scale.Sum/Counter
	Loess.thres.old<-Mean.Loess.thres/Counter
	}
	lines(Loess.thres.old, col="orange", lwd=2)
	}
} else {
	cat("running segments using external parameters from first iteration\n")
	Mean.shape.old<-fixed.from.first[1]
	Mean.scale.old<-fixed.from.first[2]
	}
	Diff.shape<-Inf
	Diff.scale<-Inf
	Diff.thres<-Inf
	Counter=1
	# now we want to fix parmaters and estimate loess
	# then we will fix loess and estimate parameters
	cat("          estimating final parameter values\n")
	while (abs(Diff.shape) > 0.05 | abs(Diff.shape)>0.05 |abs(Diff.thres)>0.01) {
	if( Counter> maxiter ) stop("maximum amount of iterations reached!\n")

	cat("Iteration", Counter, "\n")
	Input.List<-Input.List.Initial
	# here I'll fix shape and scale and estimate thres
	for (i in 1:length(Input.List)) {
		Input.List[[i]]$shape$fixed<-Mean.shape.old
		Input.List[[i]]$scale$fixed<-Mean.scale.old
		
		Index<-ceiling(i-(window.size/2)):floor(i+(window.size/2))
		# i want to extract thres from the neighbours.
		# and take it from the points..

		if (is.null(fixed.from.first) | Counter>1){
		Input.List[[i]]$Y.neighb<-Input.List[[i]]$Y.neighb-Loess.thres.old[which(1:length(Loess.thres.old) %in% Index)]	
		Input.List[[i]]$thres$initials<-0
		}
	}
	if (cpus>1) {
	Res<-as.data.frame(parSapply(cl=cl, Input.List,  my.ralling.function))
	} else {
	Res<-as.data.frame(sapply(Input.List, my.ralling.function))
	}
	names(Res)<-c("thres")
	if (is.null(fixed.from.first) | Counter>1) {
		Res$thres<-Res$thres +Loess.thres.old	
		lines(Loess.thres.old, col="yellow")
	}
	#Loess.pred<-predict(loess(Res$thres~X, span=span, weights=Lengths), se=T)
	Loess.pred<-predict(loess(Res$thres~X, span=span), se=T)
	Loess.thres<-Loess.pred$fit
	if (is.null(fixed.from.first)) lines(Res.old$thres, col="grey")
	lines(Res$thres, col="red")
	lines(Loess.thres, col="blue")
	if (is.null(fixed.from.first) | Counter>1) Diff.thres<-sum(abs((Loess.thres-Loess.thres.old)/length(Loess.thres)))
	
	Res.old<-Res

	#=======================================================
	# now I would like to run the same thing but with fixed loess...
	# 
	Input.List<-Input.List.Initial
	for (i in 1:length(Input.List)) {
		Input.List[[i]]$shape$initials<-Mean.shape.old
		Input.List[[i]]$scale$initials<-Mean.scale.old
		
		Index<-ceiling(i-(window.size/2)):floor(i+(window.size/2))
		# i want to extract thres from the neighbours.
		# and take it from the points..
		Input.List[[i]]$Y.neighb<-Input.List[[i]]$Y.neighb-Loess.thres[which(1:length(Loess.thres) %in% Index)]	
		Input.List[[i]]$thres$fixed<-0
	}
		if (cpus>1) {
	Res<-as.data.frame(t(parSapply(cl=cl, Input.List,  my.ralling.function)))
	} else {
	Res<-as.data.frame(t(sapply(Input.List, my.ralling.function)))
	}
	names(Res)<-c("shape","scale")

	#Mean.shape<-mean(Res$shape, na.rm=T)
	#Estimated.shape<-Res$shape[which(Res$shape>1e-0)]
	#cat("estimated ", length(Estimated.shape), "parameters from", length(Res$shape), "\n")
	if (is.null(fixed.from.first)) {
		Mean.shape<-weighted.mean(c(mean(Res$shape, na.rm=T), Mean.shape.old), w=c(0.6, 0.4))
	#Mean.shape<-mean(Estimated.shape, na.rm=T)
	#Mean.scale<-mean(Res$scale[which(Res$shape>1e-0)], na.rm=T)
	
		Mean.scale<-weighted.mean(c(mean(Res$scale, na.rm=T), Mean.scale.old), w=c(0.6, 0.4))
	} else { 
		Mean.shape=Mean.shape.old
		Mean.scale=Mean.scale.old
	}
	Diff.shape<-Mean.shape-Mean.shape.old
	Diff.scale<-Mean.scale-Mean.scale.old
	cat("Mean.shape: ", Mean.shape, "; Means scale ", Mean.scale, "\n" )
	cat( "Diff.shape", Diff.shape, "Diff.scale" , Diff.scale, "Diff.thres", Diff.thres, "\n")
	Mean.shape.old<-Mean.shape
	Mean.scale.old<-Mean.scale
	Loess.thres.old<-Loess.thres
	Counter=Counter+1
	}
	# and now I want to return results
	Final.res<-data.frame(shape=rep(Mean.shape, length(X)), scale=rep(Mean.scale, length(X)), thres=Loess.thres, se.fit=Loess.pred$se.fit)
	lines(Final.res$thres, col="red", lwd=2)
	return(Final.res)
}	
