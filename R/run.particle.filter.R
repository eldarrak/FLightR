run.particle.filter<-function(all.out, save.Res=T, cpus=NULL, nParticles=1e6, known.last=T, precision.sd=25, sea.value=0.00, save.memory=T, k=NA, parallel=T, plot.each.iter=T, prefix="pf", extend.prefix=T, max.kappa=100, min.SD=25, min.Prob=0.01, max.Prob=0.99, save.points.distribution=T, fixed.parameters=NA, cluster.type="SOCK", a=45, b=500, L=25, update.angle.drift=F, adaptive.resampling=0.5, save.transitions=T, check.outliers=F, sink2file=F) {

    all.out$SD<-vector(mode = "double")
    all.out$LL<-vector(mode = "double")
	
  if (parallel) {
    require("parallel")
	if (is.null(cpus)) cpus=detectCores()-1
    cat("creating cluster with", cpus, "threads\n")
    mycl <- parallel:::makeCluster(cpus, type=cluster.type)
    parallel:::clusterSetRNGStream(mycl)
    ### we don' need to send all parameters to node. so keep it easy..
    #parallel:::clusterExport(mycl,c("generate.points.dirs", "pf.par.internal", "my.dvonmises", "mu.sigma.truncnorm"))
    parallel:::clusterEvalQ(mycl, library("circular")) 
    parallel:::clusterEvalQ(mycl, library("truncnorm")) 
	parallel:::clusterEvalQ(mycl, library("FLightR")) 

  }	else mycl=NA

    Res<-pf.run.parallel.SO.resample(in.Data=all.out, cpus=cpus, nParticles=nParticles, known.last=known.last, precision.sd=precision.sd, sea.value=sea.value, k=k, parallel=parallel, plot=F, existing.cluster=mycl, cluster.type=cluster.type, a=a, b=b, L=L, sink2file=sink2file, adaptive.resampling=adaptive.resampling, RStudio=RStudio, check.outliers=check.outliers)
    # Part 2. Creating matrix of results.
    cat("creating results matrix \n")
    All.results.mat<-return.matrix.from.char(Res$All.results)
	all.out$outliers <- Res$outliers
	all.out$tmp<-Res$tmp
    # Part 2a. Estimating log likelihood
    LL<-get.LL.PF(all.out, All.results.mat)
    cat("+----------------------------------+\n")
    cat("|     estimated Log Likelihood is",  LL, "\n")
    cat("+----------------------------------+\n")
    #save(LL, file=paste(LL, "time", format(Sys.time(), "%H-%m"), ".RData"))
	# Part 2b comparing the likelihood with previous estimate
    # now we need to save it..
      all.out$LL<-LL

	  # Part 3. Updating proposal
      cat("estimating results object\n")
      all.out.old<-all.out
      all.out<-get.coordinates.PF(All.results.mat, all.out, save.points.distribution=save.points.distribution)
      all.out<-update.proposal.PF(All.results.mat, Res$Trans, all.out, fixed.parameters=fixed.parameters, a=a, b=b, parallel=parallel, existing.cluster=mycl, save.transitions=save.transitions)

    rm(Res)
    rm(All.results.mat)
    # plotting resuls
    if (plot.each.iter) {
      plot(CENTRE.y~CENTRE.x, type="p", data=all.out$Final.Means, pch=3, col="blue", main=paste("Iteration", all.out$EMIter, "(blue) and previous (red)"))
      if (all.out$EMIter>1 & !is.null(all.out$Final.Means)) {
        points(CENTRE.y~CENTRE.x, type="p", data=all.out.old$Final.Means, pch=3, col="red")
        lines(CENTRE.y~CENTRE.x, data=all.out.old$Final.Means, col="red")
      }
      points(CENTRE.y~CENTRE.x, type="p", data=all.out$Final.Means, pch=3, col="blue")
      lines(CENTRE.y~CENTRE.x, data=all.out$Final.Means, col="blue")
      data("wrld_simpl", package="maptools")
      plot(wrld_simpl, add=T)
    }
    gc()
  if (parallel) parallel:::stopCluster(cl =mycl)
  cat("DONE!\n")
  return(all.out)
}

