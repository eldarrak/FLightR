

EM.PF.wrapper<-function(all.out, iterations=10, save.Res=T, cpus=24, nParticles=1e6, known.last=T, precision.sd=25, sea.value=0.00, save.memory=T, k=NA, parallel=T, plot.each.iter=T, prefix="EMPF", extend.prefix=T, max.kappa=100, min.SD=25, min.Prob=0.01, max.Prob=0.99, start.new.optimization.SD=T, save.points.distribution=T, max.attempts=5, node.mode=F, fixed.parameters=NA, cluster.type="SOCK", a=45, b=500, sink2file=T, L=25, update.angle.drift=T, adaptive.resampling=F, RStudio=F) {
  if (sink2file & !(RStudio) ) sink(file=paste(prefix,"tmpsink.txt", sep="."))
  if (sink2file & (RStudio) ) sink()

  # this function is a wrapper around pf
  # it will take in.Data object and will go into several iterations updating in.Data by the results from the previous iteration..
  if((start.new.optimization.SD) | length(all.out$LL)==0) {
    all.out$SD<-vector(mode = "double")
    all.out$LL<-vector(mode = "double")
    all.out$Next="Going" # Next will show what is the next step
    all.out$Curr.attempt=1
    all.out$EMIter=1
  } else 	{
	all.out$EMIter=length(all.out$LL)+1 
	}
  # Part 0 starting cluster
  if (parallel) {
    cat("creating cluster\n")
    require("parallel")
    mycl <- parallel:::makeCluster(cpus, type=cluster.type)
    parallel:::clusterSetRNGStream(mycl)
    ### we don' need to send all parameters to node. so keep it easy..
    parallel:::clusterExport(mycl,c("generate.points.dirs", "pf.par.internal", "my.dvonmises", "mu.sigma.truncnorm"))
    parallel:::clusterEvalQ(mycl, library("circular")) 
    parallel:::clusterEvalQ(mycl, library("truncnorm")) 
  }	else mycl=NA
  
  while ((node.mode) | (all.out$Next %in% c("Going", "Repeat") & all.out$EMIter<=iterations )) {
    node.mode=F
    #exporting data
    cat("running iteration", all.out$EMIter,"\n")
    ###################
    # updating proposal for wrong parameter values
    #
    # part 4. updating stanard errors for the proposal
    ####
    
    #=================================================================
    # first part is update for point with migration probability ==0
    
    all.out$Matrix.Index.Table$M.mean[all.out$Matrix.Index.Table$Decision==0]<-mean(all.out$Matrix.Index.Table$M.mean, na.rm=T)
    
    all.out$Matrix.Index.Table$M.medians[all.out$Matrix.Index.Table$Decision==0]<-median(all.out$Matrix.Index.Table$M.medians, na.rm=T)
    
    all.out$Matrix.Index.Table$Kappa[all.out$Matrix.Index.Table$Decision==0]<-0
    
    all.out$Matrix.Index.Table$M.sd[all.out$Matrix.Index.Table$Decision==0]<-mean(all.out$Matrix.Index.Table$M.sd, na.rm=T)
    
    #==================================================================
    # now updating points with extremely low diversity
    
    all.out$Matrix.Index.Table$M.sd[all.out$Matrix.Index.Table$M.sd< min.SD]<-min.SD
    all.out$Matrix.Index.Table$M.sd[is.na(all.out$Matrix.Index.Table$M.sd)]<-min.SD
    all.out$Matrix.Index.Table$Kappa[all.out$Matrix.Index.Table$Kappa>max.kappa]<-max.kappa 
    all.out$Matrix.Index.Table$Kappa[is.na(all.out$Matrix.Index.Table$Kappa)]<-max.kappa 
    
    all.out$Matrix.Index.Table$Decision[all.out$Matrix.Index.Table$Decision<min.Prob]<-min.Prob
    
    all.out$Matrix.Index.Table$Decision[all.out$Matrix.Index.Table$Decision> max.Prob]<-max.Prob
    
    ####################################################
    # Part 1. first PF
    Res<-pf.run.parallel.SO.resample(in.Data=all.out, cpus=cpus, nParticles=nParticles, known.last=known.last, precision.sd=precision.sd, sea.value=sea.value, k=k, parallel=parallel, plot=F, existing.cluster=mycl, cluster.type=cluster.type, a=a, b=b, L=L, sink2file=sink2file, adaptive.resampling=adaptive.resampling, RStudio=RStudio)
    # Part 2. Creating matrix of results.
    cat("creating results matrix for iteration", all.out$EMIter, "\n")
    #All.results.mat<-return.matrix.from.char(Res, cluster=mycl)
    All.results.mat<-return.matrix.from.char(Res$All.results)
    # Part 2a. Estimating log likelihood
    LL<-get.LL.PF(all.out, All.results.mat)
	cat("+----------------------------------+\n")
	cat("|     new Log Likelihood is",  LL, "\n")
	cat("+----------------------------------+\n")
save(LL, file=paste(LL, "time", format(Sys.time(), "%Y-%m-%dT%H:%m"), ".RData", sep="."))
    # Part 2b comparing the likelihood with previous estimate
    # now we need to save it..
    if (length(all.out$LL)==0 ) {
      all.out$LL<-LL
      all.out$Next="Going"}
    else {
      # Part  compare LL with previous one..
      if ((rev(all.out$LL)[1]  <= LL)  & !(length(all.out$LL)==1 & (update.angle.drift))) {
        if (all.out$Curr.attempt < max.attempts ) {
          #try again
          all.out$Curr.attempt=all.out$Curr.attempt+1
          all.out$Next="Repeat"
          cat("LL (", LL ,") was higher then previous (",rev(all.out$LL)[1],") so going for rerun #", all.out$Curr.attempt, "\n" )
        } else {
          all.out$Next="Stop"
        }
      } else {
        all.out$Next="Going"
        all.out$Curr.attempt=1
        cat("delta LL was", rev(all.out$LL)[1] - LL, "\nmoving forward\n")
        all.out$LL=c(all.out$LL, LL)
      }
    }
    # Part 3. Updating proposal
    if (all.out$Next %in% c("Going")) {
      cat("updating proposal after iteration", all.out$EMIter, "\n")
      all.out.old<-all.out
	  all.out<-get.coordinates.PF(All.results.mat, all.out, save.rle=save.points.distribution)
      all.out<-update.proposal.PF(All.results.mat, Res$Trans, all.out, fixed.parameters=fixed.parameters, a=a, b=b, parallel=parallel, existing.cluster=mycl)
	    if (update.angle.drift) {  
	  cat("\n   estimating measured twilight angle drift\n")
save(all.out, file="all.out.tmp.beofre.angle.RData")
	  Angle.drift<-get.angle.drift(all.out)
	  all.out$Matrix.Index.Table$Angle<-Angle.drift$Mean
	  all.out$Matrix.Index.Table$Angle.sd<-Angle.drift$SD
	  cat("   estimating new sun matrices\n\n")
	  all.out$Phys.Mat<-create.spatial.sun.matrices(all.out, existing.cluster=mycl)
	  # for now cpus =1 as we can't send call from here to the nodes.. 
	  # or at least it is not obvious hiow to do it.
	  }
      # Part 4. Looking at iterations SD
      if (!is.null(all.out.old$Final.Means) & (!start.new.optimization.SD)) {
        SD<-get.optimization.SD(all.out.old, all.out)
        #if (extend.prefix) {
        #	names(SD)<-rev(unlist(strsplit(prefix, "\\.")))[1]
        #} else {
        names(SD)<-paste("Iter", all.out$EMIter, sep=".")
        #}
        all.out$SD<-c(all.out$SD, SD)
        start.new.optimization.SD<-T
        print(all.out$SD)
      }
      if (save.Res) {
        if(extend.prefix) {
          save( All.results.mat, file=paste("All.results.mat", prefix, "iteration", all.out$EMIter, "RData", sep="."))
        } else {
          save( All.results.mat, file=paste("All.results.mat", prefix,  "RData", sep="."))
        }
        if(extend.prefix) {
          save(all.out, file=paste("all.out", prefix, "iteration", all.out$EMIter, "RData", sep="."))
        } else {
          save(all.out, file=paste("all.out", prefix, "RData", sep="."))
        }
      }
    }
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
    if (all.out$Next=="Going") all.out$EMIter=length(all.out$LL)+1
    gc()
  }
  if (parallel) parallel:::stopCluster(cl =mycl)
  cat("all iters finished!\n")
  if (sink2file) sink()
  return(all.out)
}
