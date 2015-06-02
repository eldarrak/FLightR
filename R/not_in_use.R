# not_in_use.R
# this functions are currently not used but may be used at some point....
# I will not release them on CRAN and will archive locally.
stop("buildingore")

get.time.shift<-function(start, Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk,  log.light.borders=NA, diap=c(-600, 600), plot=T,  log.irrad.borders=c(-15, 50), verbose=F) {

	# THIS FUNCTION IS NOT CURRENTLY U
	# this version will try to make slope difference insignificant..
	# will try to make this criteria as the main.

	# for now - very simple function:
	# 3 loops - minutes - 10 seconds - seconds..
	# 	# we may want to do it with optim..
	x<-list(start[1], start[2])
	get.time.shift.internal<-function(start, Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk,  log.light.borders, diap, Step=1, plot,  log.irrad.borders=c(-8, 1.5)) {
		Final<-c()
		for (i in seq(diap[1], diap[2], Step)) {
#	cat("i", i, "Step", Step, "\n")	
			Calibration<-try(logger.template.calibration(Twilight.time.mat.Calib.dawn+i, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk+i, Twilight.log.light.mat.Calib.dusk, positions=start, log.light.borders=log.light.borders,  log.irrad.borders= log.irrad.borders))
			
			# and now I want to estimate LL that the bird was in the known location
			LL<-0
			All<- Inf
			
			# now combining all data together
			Twilight.time.mat.Calib<-cbind(Twilight.time.mat.Calib.dawn, Twilight.time.mat.Calib.dusk)
			Dusk<-c(rep(FALSE, dim(Twilight.time.mat.Calib.dawn)[2]), rep(TRUE, dim(Twilight.time.mat.Calib.dusk)[2]))
			Twilight.log.light.mat.Calib<-cbind(Twilight.log.light.mat.Calib.dawn, Twilight.log.light.mat.Calib.dusk)
			
			for (Twilight.ID in 1:(dim(Twilight.time.mat.Calib)[2])) {
			#print(Twilight.ID )
			Twilight.solar.vector<-solar(as.POSIXct(Twilight.time.mat.Calib[c(1:24, 26:49), Twilight.ID]+i, tz="gmt", origin="1970-01-01"))
			Twilight.log.light.vector<-Twilight.log.light.mat.Calib[c(1:24, 26:49), Twilight.ID]
			Res<-get.current.slope.prob(x, calibration=Calibration,  Twilight.solar.vector=Twilight.solar.vector, Twilight.log.light.vector=Twilight.log.light.vector, plot=F, verbose=verbose,  log.light.borders=log.light.borders,  log.irrad.borders= log.irrad.borders, dusk=Dusk[Twilight.ID])
			#cat("Twilight.ID", Twilight.ID, "Res", Res, "\n")
			if (is.na(Res) | Res==0 ) {
			Res<-get.current.slope.prob(x, calibration=Calibration,  Twilight.solar.vector=Twilight.solar.vector, Twilight.log.light.vector=Twilight.log.light.vector, plot=F, verbose=T,  log.light.borders=log.light.borders,  log.irrad.borders= log.irrad.borders, dusk=Dusk[Twilight.ID])
			
			break

			} 
			LL<-LL+log(Res)
			All<-c(All, Res)
		}
#print(length(All))
#print(All)
#print((dim(Twilight.time.mat.Calib)[2]))		
		if (length(All)<=(dim(Twilight.time.mat.Calib)[2])) {
		All=-Inf;
		LL=-Inf
		}
		Final=rbind(Final, cbind(i, LL, min(All),  Calibration$Significance.of.dusk.dawn.diff))
		if (TRUE) print(Final)
		}
		return(Final)
		}
		
	# minutes
	Final.min<-get.time.shift.internal(start, Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk,  log.light.borders, diap=diap, Step=60, plot=F, log.irrad.borders=log.irrad.borders)
	cat("\n")
	if (max(Final.min[,2])==-Inf) stop("something went wrong\n more likely that the logger was not stable during the calibration time\n exclude some of the twilights and try again\n")
	#Diap.10.sec<-c(Final.min[which.max(Final.min[,4]),1]-60, Final.min[which.max(Final.min[,4]),1]+60)
	Diap.10.sec<-c(Final.min[which.max(Final.min[,2]),1]-60, Final.min[which.max(Final.min[,2]),1]+60)
	Final.10secs.min<-get.time.shift.internal(start, Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk,  log.light.borders, diap=Diap.10.sec, Step=10,  log.irrad.borders= log.irrad.borders)
	
	Diap.1.sec<-c(Final.10secs.min[which.max(Final.10secs.min[,2]),1]-10, Final.10secs.min[which.max(Final.10secs.min[,2]),1]+10)
	Final.1sec.min<-get.time.shift.internal(start, Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk,  log.light.borders, diap=Diap.1.sec, Step=1,  log.irrad.borders= log.irrad.borders)
	
	return(Final.1sec.min[which.max(Final.1sec.min[,2]),1])
	}
	



get.angle.drift<-function(in.Data, plot=F) {
  #require(tripEstimation)
  # last.point.fun<-function(x) {
  #	in.Data$Points.Land[x%%1e5,]
  #}
  Mean.angle<-NULL
  SD.angle<-NULL
  Angles.rle<-in.Data$Points.rle
  for (i in 1:length(in.Data$Matrix.Index.Table$Real.time)) {
    Angles.rle[[i+1]]$values<-apply(matrix(in.Data$Points.Land[in.Data$Points.rle[[i+1]]$values,1:2], ncol=2), 1, FUN=function(x) elevation(x[[1]], x[[2]], solar(in.Data$Matrix.Index.Table$Real.time[i])))
    SD.angle<-c(SD.angle, sd(inverse.rle(Angles.rle[[i+1]])))
    Mean.angle<-c(Mean.angle, weighted.mean(Angles.rle[[i+1]]$values, Angles.rle[[i+1]]$lengths))
  }
  #Angles.rle<-Angles.rle[-1]
  if (plot) {
    plot(Mean.angle[in.Data$Matrix.Index.Table$Dusk]~as.numeric(in.Data$Matrix.Index.Table$Real.time[in.Data$Matrix.Index.Table$Dusk]), col="blue")
    lines(predict(loess(Mean.angle[in.Data$Matrix.Index.Table$Dusk]~as.numeric(in.Data$Matrix.Index.Table$Real.time[in.Data$Matrix.Index.Table$Dusk]), span=150/length(Mean.angle[in.Data$Matrix.Index.Table$Dusk])))~ as.numeric(in.Data$Matrix.Index.Table$Real.time[in.Data$Matrix.Index.Table$Dusk]), col="blue", type="l", lwd=2)
    
    points(Mean.angle[!in.Data$Matrix.Index.Table$Dusk]~as.numeric(in.Data$Matrix.Index.Table$Real.time[!in.Data$Matrix.Index.Table$Dusk]), pch="+", col="red")
    lines(predict(loess(Mean.angle[!in.Data$Matrix.Index.Table$Dusk]~as.numeric(in.Data$Matrix.Index.Table$Real.time[!in.Data$Matrix.Index.Table$Dusk]), span=150/length(Mean.angle[!in.Data$Matrix.Index.Table$Dusk])))~ as.numeric(in.Data$Matrix.Index.Table$Real.time[!in.Data$Matrix.Index.Table$Dusk]), col="red", type="l", lwd=2)
  }
  
  Mean.Dusk.angle.drift<-predict(loess(Mean.angle[in.Data$Matrix.Index.Table$Dusk]~as.numeric(in.Data$Matrix.Index.Table$Real.time[in.Data$Matrix.Index.Table$Dusk]), span=150/length(Mean.angle[in.Data$Matrix.Index.Table$Dusk])))
  SD.Dusk.angle.drift<-sd(Mean.angle[in.Data$Matrix.Index.Table$Dusk] - Mean.Dusk.angle.drift)
  Mean.Dawn.angle.drift<-predict(loess(Mean.angle[!in.Data$Matrix.Index.Table$Dusk]~as.numeric(in.Data$Matrix.Index.Table$Real.time[!in.Data$Matrix.Index.Table$Dusk]), span=150/length(Mean.angle[in.Data$Matrix.Index.Table$Dusk])))
  SD.Dawn.angle.drift<-sd(Mean.angle[!in.Data$Matrix.Index.Table$Dusk] - Mean.Dawn.angle.drift)
  # now we want to comibine them
  Mean.angle.drift<-rep(NA, length(Mean.angle))
  SD.angle.drift<-rep(NA, length(Mean.angle))
  Mean.angle.drift[in.Data$Matrix.Index.Table$Dusk]<-Mean.Dusk.angle.drift
  Mean.angle.drift[!in.Data$Matrix.Index.Table$Dusk]<-Mean.Dawn.angle.drift
  SD.angle.drift[in.Data$Matrix.Index.Table$Dusk]<-SD.Dusk.angle.drift
  SD.angle.drift[!in.Data$Matrix.Index.Table$Dusk]<-SD.Dawn.angle.drift
  return(data.frame(Mean=Mean.angle.drift, SD=SD.angle.drift))
}

# usage:
#Angle.drift<-get.angle.drift(in.Data)

# new functions will replace create.spatial.sun.matrix
# function 1 - get calibration angles

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

# now I want to get make write a new create.spatial.sun.matrix function that will work from in.Data
# the problem here is that it will still need loess data otherwise it will not work. how to incorporate loess?
# added into the create.proposal function..
# now Index.tab has got loess and potentially could get angle from the calibration..

# then we will 
# 1. run geologger.sampler.create.arrays
# 2. add calibration angle
# 3. and only after that we will create spatial sun matrices.

# now we need a function theat will estimate angles..
get.initial.calibration<-function(in.Data, days.stable=list(Dusk=c(1:20), Dawn=c(1:20))) {
  par(mfrow=c(2,2))
  Dusk.angle<-my.calibrate(in.Data$Matrix.Index.Table$Real.time[in.Data$Matrix.Index.Table$Dusk==T][days.stable$Dusk],  in.Data$Points.Land[in.Data$start.point,1],  in.Data$Points.Land[in.Data$start.point,2])
  Dawn.angle<-my.calibrate(in.Data$Matrix.Index.Table$Real.time[in.Data$Matrix.Index.Table$Dusk==F][days.stable$Dawn]	,  in.Data$Points.Land[in.Data$start.point,1],  in.Data$Points.Land[in.Data$start.point,2])
  in.Data$Matrix.Index.Table$Angle[in.Data$Matrix.Index.Table$Dusk==T]<-Dusk.angle[1]
  in.Data$Matrix.Index.Table$Angle.sd[in.Data$Matrix.Index.Table$Dusk==T]<-Dusk.angle[2]
  in.Data$Matrix.Index.Table$Angle[in.Data$Matrix.Index.Table$Dusk==F]<-Dawn.angle[1]
  in.Data$Matrix.Index.Table$Angle.sd[in.Data$Matrix.Index.Table$Dusk==F]<-Dawn.angle[2]
  par(mfrow=c(1,1))
  return(in.Data$Matrix.Index.Table)
}


sun.matrix.internal<-function(Matrix.Index.Table.Row, Points.Land) {
  # this function will preestimate SpSunM for one row of in.data
  # spatial.Points is an Rdata file that contains XXX?? object
  #require(tripEstimation)
  require(maptools)
  require(fields)
  
  correct.estimates<-function(Latitudes) {
    # this function estimates correction coef. for twilight matrices
    # it will look for a file "Correction.coefs.RData"
    Res<-suppressWarnings(try(load("Correction.coefs.RData"), silent=TRUE))
    if (class(Res)=="try-error") stop("Correction file: \"Correction.coefs.RData\" is needed in the working directory\n")
    if (any(abs(Latitudes)>84)) warning("there is no correction coefs for latitudes over 84 degree - NAs were returned\n")
    return(predict(A, Latitudes))
  }
  
  time.angle.integrate<-function(x) {
    sum(dnorm(elevation(x[[1]], x[[2]], solar(p.time.gmt)), mean=Angle[1], sd=Angle[2])*p.time.dens*p.dtime)
  }
  # integrating without the last point
  p.num = 100
  p.time.cumprob <- seq(0.001, 0.999, length.out=p.num)                     # in probabilities
  p.time.points  <- qnorm(p.time.cumprob,as.numeric(Matrix.Index.Table.Row$Real.time),Matrix.Index.Table.Row$Loess.se.fit*3600*sqrt(Matrix.Index.Table.Row$Loess.n))  # in time
  p.dtime        <- (p.time.points[-1] - p.time.points[-p.num])           # step in time - dTime 
  p.time.points1 <-p.time.points[-p.num]                                  # deleting last point
  p.time.dens    <- dnorm(p.time.points1,as.numeric(Matrix.Index.Table.Row$Real.time),Matrix.Index.Table.Row$Loess.se.fit*3600*sqrt(Matrix.Index.Table.Row$Loess.n))  # probability density - p(Time)
  p.time.gmt <- as.POSIXct(p.time.points1, tz="GMT", origin="1970-01-01")
  Angle<-as.matrix(Matrix.Index.Table.Row[c("Angle", "Angle.sd")])[1,]
  Res<-apply(Points.Land, 1, FUN=time.angle.integrate)
  Correction.coefs<-correct.estimates(Points.Land[,2])
  Res<-Res/Correction.coefs
  return(Res)
}

#Z<-sun.matrix.internal(in.Data$Matrix.Index.Table[1,], in.Data$Points.Land)
node.run<-function(x=x) {
  Points.Land<-get("Points.Land", globalenv())
  return(sun.matrix.internal(x, Points.Land))
}



#================================================
# and now we want to make a wrapper around this that will estimate Phys.Mat for all rows in the Matrix.Index.Table

create.spatial.sun.matrices<-function(in.Data, cpus=7, cluster.type="SOCK", existing.cluster=NULL) {	
  List<-suppressWarnings(split(in.Data$Matrix.Index.Table, 1:nrow(in.Data$Matrix.Index.Table)))
  Points.Land<-in.Data$Points.Land
  if (cpus >1) {
    if (is.null(existing.cluster)) {
      require("parallel")
      mycl <- parallel:::makeCluster(cpus, type=cluster.type)
      parallel:::clusterSetRNGStream(mycl)
      ### we don' need to send all parameters to node. so keep it easy..
    } else {mycl<-existing.cluster}
    WD<-getwd()
    parallel:::clusterExport(mycl, "WD", envir=environment())
    parallel:::clusterEvalQ(mycl, setwd(WD)) 
    parallel:::clusterExport(mycl,"Points.Land", envir=environment())
    #parallel:::clusterEvalQ(mycl, library("tripEstimation"))
    parallel:::clusterExport(mycl,c("sun.matrix.internal","node.run", "elevation", "solar"))
    Res<-simplify2array(parLapply(cl = mycl, List, node.run))
    if (is.null(existing.cluster)) stopCluster(mycl)
  } else {
    Res<-simplify2array(lapply(List, FUN=function(x) sun.matrix.internal(x, Points.Land)))
  }
  return(Res)
}



get.optimization.SD<-function(all.out.old, all.out, return.values=F) {
  require("sp")
  # this function will estimate SD between iterations.
  Coords<-cbind(all.out.old$Final.Means$CENTRE.x, all.out.old$Final.Means$CENTRE.y, all.out$Final.Means$CENTRE.x, all.out$Final.Means$CENTRE.y)
  Dists<-apply(Coords, 1, FUN=function(x) {sp:::spDistsN1(cbind(x[[1]], x[[2]]), cbind(x[[3]], x[[4]]), longlat=TRUE)})
  if(return.values) {return(Dists)
  } else {return(sd(Dists))}
}



EM.PF.wrapper<-function(all.out, iterations=10, save.Res=T, cpus=24, nParticles=1e6, known.last=T, precision.sd=25, sea.value=0.00, save.memory=T, k=NA, parallel=T, plot.each.iter=T, prefix="EMPF", extend.prefix=T, max.kappa=100, min.SD=25, min.Prob=0.01, max.Prob=0.99, start.new.optimization.SD=T, save.points.distribution=T, max.attempts=5, node.mode=F, fixed.parameters=NA, cluster.type="SOCK", a=45, b=500, sink2file=T, L=25, update.angle.drift=T, adaptive.resampling=0.5, RStudio=F, save.transitions=T, check.outliers=F) {
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
    Res<-pf.run.parallel.SO.resample(in.Data=all.out, cpus=cpus, nParticles=nParticles, known.last=known.last, precision.sd=precision.sd, sea.value=sea.value, k=k, parallel=parallel, plot=F, existing.cluster=mycl, cluster.type=cluster.type, a=a, b=b, L=L, sink2file=sink2file, adaptive.resampling=adaptive.resampling, RStudio=RStudio, check.outliers=check.outliers)
    # Part 2. Creating matrix of results.
    cat("creating results matrix for iteration", all.out$EMIter, "\n")
    #All.results.mat<-return.matrix.from.char(Res, cluster=mycl)
    All.results.mat<-return.matrix.from.char(Res$All.results)
	all.out$outliers <- Res$outliers
	all.out$tmp<-Res$tmp
    # Part 2a. Estimating log likelihood
    LL<-get.LL.PF(all.out, All.results.mat)
    cat("+----------------------------------+\n")
    cat("|     new Log Likelihood is",  LL, "\n")
    cat("+----------------------------------+\n")
    save(LL, file=paste(LL, "time", format(Sys.time(), "%H-%m"), ".RData"))
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
      all.out<-get.coordinates.PF(All.results.mat, all.out, save.points.distribution=save.points.distribution)
      all.out<-update.proposal.PF(All.results.mat, Res$Trans, all.out, fixed.parameters=fixed.parameters, a=a, b=b, parallel=parallel, existing.cluster=mycl, save.transitions=save.transitions, estimatetruncnorm=T)
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

###################################
# now we also need a function that will estimate
# SD of movement between iterations
# first we have to estimate differences between point positions.

####################################################################
# to prevent memory leak I decide to create a new wrapper around previous one.
# the idea is that it will restart main R node after every iteration.
# decided to change the name

EM.PF.safe.ram.wrapper<-function(all.out, iterations=10, save.Res=T, cpus=40, nParticles=1e6, known.last=T, precision.sd=25, sea.value=0.00, save.memory=T, k=NA, parallel=T, plot.each.iter=T, prefix="EMPF", max.kappa=100, min.SD=25, min.Prob=0.01, max.Prob=0.99, start.new.optimization.SD=T, save.points.distribution=T, save.transitions=T, max.attempts=3, fixed.parameters=NA, cluster.type="SOCK", sink2file=T, L=25,update.angle.drift=T, adaptive.resampling=0.5, plot2pdf=F, RStudio=F, a=45, b=500, check.outliers=F) {
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
    #parallel:::clusterExport(maincl,c("generate.points.dirs", "pf.par.internal", "my.dvonmises", "pf.final.smoothing", "plot.optimisation.results", "pf.run.parallel.SO.resample", "return.matrix.from.char", "get.transition.rle", "update.proposal.PF", "get.coordinates.PF", "EM.PF.wrapper", "get.optimization.SD", "get.LL.PF", "mu.sigma.truncnorm", "get.angle.drift", "create.spatial.sun.matrices", "sun.matrix.internal", "node.run" ,"solar", "elevation"))
    parallel:::clusterExport(maincl, "WD", envir=environment())
    parallel:::clusterEvalQ(maincl, setwd(WD)) 
    parallel:::clusterEvalQ(maincl, library("circular")) 
    parallel:::clusterEvalQ(maincl, library("truncnorm")) 
    parallel:::clusterEvalQ(maincl, library("parallel")) 
    parallel:::clusterEvalQ(maincl, library("FLightR")) 
    
    # now I want to run wrapper on cluster
    
    MainParameters<-list(all.out=all.out, iterations=1, save.Res=save.Res, cpus=cpus, nParticles=nParticles, known.last=known.last, precision.sd=precision.sd, sea.value=sea.value, save.memory=save.memory, k=k, parallel=parallel, plot.each.iter=FALSE, prefix=Curr.prefix, max.kappa=max.kappa, min.SD=min.SD, min.Prob=min.Prob, max.Prob=max.Prob, start.new.optimization.SD=start.new.optimization.SD, extend.prefix=FALSE, save.points.distribution=save.points.distribution, max.attempts=max.attempts, node.mode=TRUE, fixed.parameters=fixed.parameters, cluster.type=cluster.type, sink2file=sink2file, L=L, update.angle.drift=update.angle.drift, adaptive.resampling=adaptive.resampling, RStudio=RStudio, a=a, b=b, save.transitions=save.transitions, check.outliers=check.outliers)
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


## function to get optimization coefficient...
get.max.kappa<-function(max.sd, a=45, b=500, maxIter=10000, verbose=F) {
  require(truncnorm)
  # max(dtruncnorm())/min(dtruncnorm())==max(dvonmises/min(dvonmises))
  try.k<-100
  Run=T
  Iter=1
  while (Iter<maxIter & Run) {
    max.diff.dist<-truncnorm:::dtruncnorm(x=a, mean=a, sd=max.sd, a=a, b=b)/truncnorm:::dtruncnorm(x=b, mean=a, sd=max.sd, a=a, b=b)
    max.diff.dir<-suppressWarnings(dvonmises(x=0, mu=0, kappa=try.k)/dvonmises(x=pi, mu=0, kappa=try.k))
    if (verbose) {
      cat("Iteration:", Iter,  "\n")
      cat("max.diff.dist:", max.diff.dist, "\n")
      cat("max.diff.dir:", max.diff.dir, "\n")
      cat("ratio:", max.diff.dist/max.diff.dir, "\n") }
    Ratio=max.diff.dist/max.diff.dir
    if (abs(log(Ratio))>0.5) {
      if (Ratio>1) { try.k=try.k+0.5		}
      else {try.k=try.k-0.5}
      if (verbose) cat("try.k", try.k, "\n")
    } else {
      cat("OPTIMIZED! in", Iter, "iterations\n")
      Run=F
    }
    Iter=Iter+1
  }
  return(try.k)
}
# get.max.kappa(max.sd=25, a=1, b=500, maxIter=1000)

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


# this internal function works inside log normal loess...
create.input.list<-function( X, Y, window.size=5) {
	if (window.size%%2==0) {
		cat("window.size was corrected to",window.size+1, "\n" )
		window.size<-window.size+1
	}
	Input.List<-vector("list", length(X))
	names(Input.List)<-X	
	# now adding parts
	for (i in 1:length(Input.List)) {
	Input.List[[i]]$X<-X[i]
	Input.List[[i]]$Y<-Y[i]
	# now neighbors
	Index<-ceiling(i-(window.size/2)):floor(i+(window.size/2))
	Input.List[[i]]$Y.neighb<-Y[which(1:length(Y) %in% Index)]
	#
	Input.List[[i]]$scale<-list(initials=NA, fixed=NA)
	Input.List[[i]]$shape<-list(initials=NA, fixed=NA)
	Input.List[[i]]$thres<-list(initials=NA, fixed=NA)
	}
	return(Input.List)
}
	


# this function estimates likelihood of 3 paramter lognormal distribution..

ll.dlnorm3<-function(prm, x, fixed=c(NA, NA, NA), trace=F) {
#print(prm)
	if (!all(is.na((fixed)))) {
	if (!is.na(fixed[1])) {
		shape=fixed[1]
		} else {
		shape=prm[1]
		prm=prm[-1]}
	if (!is.na(fixed[2])) {
		scale=fixed[2]
		} else {
		scale=prm[1]
		prm=prm[-1]
	}
	if (!is.na(fixed[3])) {
		thres=fixed[3]
		} else {
		thres=prm[1]
	}
	} else {
		shape=prm[1]
		scale=prm[2]
		thres=prm[3]
	}
if (trace) {
	cat("prm\n")
	print(shape)
	print(scale)
	print(thres)
	cat("x\n")
	print(x)
	cat("lnorm\n")
	print(dlnorm3(x, shape=shape, scale=scale, thres=thres, log=T))
}
	Res<-min(-sum(dlnorm3(x, shape=shape, scale=scale, thres=thres, log=T)), 1e50)
	if (trace) print(Res)
	return(Res)
}

# density of 3 paramter lognormal distribution
# NOTE: it takes exp shape as an input to make it suitable for optimisation routines.
dlnorm3<-function (x, shape = 1, scale = 1, thres = 0, log = FALSE) {
    fx <- dlnorm(x - thres, scale, exp(shape))
    if (log) 
        return(log(fx))
    else return(fx)
}

