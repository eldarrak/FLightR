EM.PF.safe.ram.wrapper<-function(all.out, iterations=10, save.Res=T, cpus=40, nParticles=1e6, known.last=T, precision.sd=25, sea.value=0.00, save.memory=T, k=NA, parallel=T, plot.each.iter=T, prefix="EMPF", max.kappa=100, min.SD=25, min.Prob=0.01, max.Prob=0.99, start.new.optimization.SD=T, save.points.distribution=T, max.attempts=3, fixed.parameters=NA, cluster.type="SOCK", sink2file=T, L=25,update.angle.drift=T, adaptive.resampling=F, plot2pdf=F, RStudio=F, a=45, b=500) {
  all.out$Call<-match.call()
  # fixed.parameters could be list(Kappa=1, M.sd=250)
  # in a new version it can also be M.mean..
  ###
  ### ok, first we need to initialize a new cluster with only 1 node!!!
  ##
  cat("creating first level cluster \n")
  require(parallel)
  WD<-getwd()
  if (start.new.optimization.SD) {
    all.out$Next<-"Going"
    all.out$EMIter<-1
	all.out$LL<-NULL
  }
  while (all.out$Next %in% c("Going", "Repeat") & all.out$EMIter<iterations) {
    cat("running iteration", all.out$EMIter, "\n")
	Start.time<-Sys.time()
    Curr.prefix<-paste(prefix, "iteration", all.out$EMIter, sep=".")
    maincl <- parallel:::makeCluster(1, type="SOCK")
    ### we don' need to send all parameters to node. so keep it easy..
    parallel:::clusterExport(maincl,c("generate.points.dirs", "pf.par.internal", "my.dvonmises", "pf.final.smoothing", "plot.optimisation.results", "pf.run.parallel.SO.resample", "return.matrix.from.char", "get.transition.rle", "update.proposal.PF", "get.coordinates.PF", "EM.PF.wrapper", "get.optimization.SD", "get.LL.PF", "mu.sigma.truncnorm", "get.angle.drift", "create.spatial.sun.matrices", "sun.matrix.internal", "node.run" ,"solar", "elevation"))
    parallel:::clusterExport(maincl, "WD", envir=environment())
    parallel:::clusterEvalQ(maincl, setwd(WD)) 
    parallel:::clusterEvalQ(maincl, library("circular")) 
    parallel:::clusterEvalQ(maincl, library("truncnorm")) 
    parallel:::clusterEvalQ(maincl, library("parallel")) 
    
    # now I want to run wrapper on cluster
    
    MainParameters<-list(all.out=all.out, iterations=1, save.Res=save.Res, cpus=cpus, nParticles=nParticles, known.last=known.last, precision.sd=precision.sd, sea.value=sea.value, save.memory=save.memory, k=k, parallel=parallel, plot.each.iter=FALSE, prefix=Curr.prefix, max.kappa=max.kappa, min.SD=min.SD, min.Prob=min.Prob, max.Prob=max.Prob, start.new.optimization.SD=start.new.optimization.SD, extend.prefix=FALSE, save.points.distribution=save.points.distribution, max.attempts=max.attempts, node.mode=TRUE, fixed.parameters=fixed.parameters, cluster.type=cluster.type, sink2file=sink2file, L=L, update.angle.drift=update.angle.drift, adaptive.resampling=adaptive.resampling, RStudio=RStudio, a=a, b=b)
    parallel:::clusterExport(maincl, "MainParameters" , envir=environment())
    all.out.old<-all.out
    all.out<-parallel:::clusterEvalQ(maincl, do.call(EM.PF.wrapper , MainParameters))
    all.out<-all.out[[1]]
    cat("cluster estimation finished\n")
		  
    rm(MainParameters)
    parallel:::stopCluster(cl = maincl)
    if (!is.null(all.out$SD)) print(all.out$SD)
    cat("Negative Log Likelihoods:\n")
    if (!is.null(all.out$LL)) print(all.out$LL)
    if ((all.out$Next)=="Going") cat("last iteration succesfully improved LL, new proposal accepted\n")
    if ((all.out$Next)=="Repeat") {
	cat("last iteration DID NOT  improve LL,rerunning estimation with old proposal, attempt", all.out$Curr.attempt, "\n")
	}
    if ((all.out$Next)=="Stop") {
	cat("with last", all.out$Curr.attempt, "attempts LL was not optimized any better...\n Assuming convergence\n")
	}
    
    start.new.optimization.SD<-F
    if (plot.each.iter & (all.out$Next == "Going")) {
	if (plot2pdf) pdf(file=paste(prefix, "iteration", all.out$EMIter-1, "pdf", sep="."))
      plot(CENTRE.y~CENTRE.x, type="p", pch=3, col="blue", main=paste("Iteration", all.out$EMIter-1, "(blue) and previous (red)"), data=all.out$Final.Means)
	  # new part
      if (all.out$EMIter-1>1 & !is.null(all.out.old$Final.Means)) {
	  Median.X<-sapply(all.out.old$Points.rle, FUN=function(x) median(all.out.old$Points.Land[inverse.rle(x),1]))
	  Median.Y<-sapply(all.out.old$Points.rle, FUN=function(x) median(all.out.old$Points.Land[inverse.rle(x),2]))
	  points(Median.Y~Median.X, type="p", data=all.out$Final.Means, pch=1, col="red")
	  lines(Median.Y~Median.X,  data=all.out$Final.Means, pch=1, col="red", lty=2)
	  
        points(CENTRE.y~CENTRE.x, type="p", data=all.out.old$Final.Means, pch=3, col="red")
        lines(CENTRE.y~CENTRE.x, data=all.out.old$Final.Means, col="red")
      }
	  Median.X<-sapply(all.out$Points.rle, FUN=function(x) median(all.out$Points.Land[inverse.rle(x),1]))
	  Median.Y<-sapply(all.out$Points.rle, FUN=function(x) median(all.out$Points.Land[inverse.rle(x),2]))
	  points(Median.Y~Median.X, type="p", data=all.out$Final.Means, pch=1, col="blue")
	  lines(Median.Y~Median.X,  data=all.out$Final.Means, pch=1, col="blue", lty=2)
 
      points(CENTRE.y~CENTRE.x, type="p", data=all.out$Final.Means, pch=3, col="blue")
      lines(CENTRE.y~CENTRE.x, data=all.out$Final.Means, col="blue")
      library("maptools")
   data("wrld_simpl", package="maptools")
      plot(wrld_simpl, add=T)
	  if (plot2pdf) dev.off()
    }
	if (all.out$Next != "Going") {
	# ok, I think if we are going to repeat then at this step we need to set all.out to old all.out
	all.out.old$Curr.attempt<-all.out$Curr.attempt
	all.out.old$Next<-all.out$Next
	all.out<-all.out.old
	}
  cat("Iteration took", difftime(Sys.time(), Start.time,     units = c("hours")), "hours", "\n")
    if (!save.Res) {
		save(all.out, file="all.out.last.iter.RData") 
	if (all.out$EMIter==2) {
				save(all.out, file=paste("all.out", prefix, "iteration.1.RData", sep="."))
	}
		#	}
		#}
  #if (!save.Res) {
          #save(all.out, file=paste("all.out", prefix, "iteration", all.out$EMIter, "RData", sep="."))
  }
  }
  return(all.out)
}
