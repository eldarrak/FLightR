# run_particle.filter.R
# functions used during the main run

#' Run Particle Filter
#' 
#' Main function of FLightR, it takes fully prepared object created by \code{\link{make.prerun.object}} and produces a result object that can be used for plotiing etc.
#' @param all.out An object created by \code{\link{make.prerun.object}}.
#' @param threads An amount of threads to use while running in parallel. default is -1. if value 1 submitted package will run sequentially
#' @param cpus another way to specify  \code{threads}
#' @param nParticles  total amount of particles to be used with the run. 10 000 (1e4) is recommended for the preliminary run and 1 000 000 (1e6) for the final
#' @param known.last Set to FALSE if your bird was not at a known place during last twilight in the data
#' @param precision.sd if \code{known.last} then what is the precision of this information. Will be used to resample particles prportionally to their ditance from the known last point with probability \code{P = dnorm(0, precision.sd)}
#' @param behav.mask.low.value Probability value that will be used instead of 0 in the behavioural mask. If set to 1 behavioural mask will not be active anymore
#' @param k Kappa parameter from vonMises distribution. Default is NA, otherwise will generate particles in a direction of a previous transitions with kappa = k
#' @param plot Should function plot preliminary map in the end of the run?
#' @param cluster.type see help to package parallel for details
#' @param a minimum distance that is used in the movement model - left boundary for truncated normal distribtuon of ditances moved between twilights. Default is 45 for as default grid has a minimum ditance of 50 km.
#' @param b Maximum distance allowed to fly between two consequtive twilights
#' @param L how many consequitive particles to resample
#' @param adaptive.resampling Above what level of ESS resampling should be skipped
#' @param check.outliers switches ON the online outlier routine 
#' @param sink2file will write run details in a file instead of showing on the screen
#' @param add.jitter will add spatial jitter inside a grid cell for the median estiamtes
#' @return FLightR object, containing output and extracted results. It is a list with the following elements 
#' 
#'    \item{Indices}{List with prior information and indices}
#'    \item{Spatial}{Spatial data - Grid, Mask, spatial likelihood}
#'    \item{Calibration}{all calibration parameters}
#'    \item{Data}{original data}
#'    \item{Results}{The main results object. Main components of it are
#'       \describe{
#'       \item{Quantiles}{dataframe containing results on locations. Each line corresponds to a twilight}
#'       \item{Movement.results}{dataframe containing all the movement results, Note - time at line n means time of the end of transition between n and n-1}
#'       \item{outliers}{id of twilights excluded by online outlier detection tool}
#'       \item{LL}{-Log likelihood}
#'       \item{Points.rle}{run length encoding object with posterior distribution for every twilight. Note that numbers of points correspond to line numbers in \code{$Spatial$Grid}}
#'       \item{Transitions.rle}{run length encoding object with all the transitions}
#'        }
#'   }
#'
#' @examples
#' File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
#' # to run example fast we will cut the real data file by 2013 Aug 20
#' Proc.data<-get.tags.data(File, end.date=as.POSIXct('2013-07-02', tz='GMT'))
#' Calibration.periods<-data.frame(
#'        calibration.start=NA,
#'        calibration.stop=as.POSIXct("2013-08-20"),
#'        lon=5.43, lat=52.93) 
#'        #use c() also for the geographic coordinates, if you have more than one calibration location
#'        # (e. g.,  lon=c(5.43, 6.00), lat=c(52.93,52.94))
#' print(Calibration.periods)
#' 
#' # NB Below likelihood.correction is set to FALSE for fast run! 
#' # Leave it as default TRUE for real examples
#' Calibration<-make.calibration(Proc.data, Calibration.periods, likelihood.correction=FALSE)
#' 
#' Grid<-make.grid(left=0, bottom=50, right=10, top=56,
#'   distance.from.land.allowed.to.use=c(-Inf, Inf),
#'   distance.from.land.allowed.to.stay=c(-Inf, Inf))
#'
#' all.in<-make.prerun.object(Proc.data, Grid, start=c(5.43, 52.93),
#'                              Calibration=Calibration, threads=2)
#' # here we will run only 1e4 partilces for a very short track.
#' # One should use 1e6 particles for the full run.
#' Result<-run.particle.filter(all.in, threads=1,
#'            nParticles=1e3, known.last=TRUE,
#'            precision.sd=25, check.outliers=FALSE)
#'
#' @author Eldar Rakhimberdiev
#' @export
run.particle.filter<-function(all.out, cpus=NULL, threads=-1, nParticles=1e6, known.last=TRUE, precision.sd=25, behav.mask.low.value=0.00, k=NA, plot=TRUE, cluster.type="PSOCK", a=45, b=1500, L=90, adaptive.resampling=0.99, check.outliers=FALSE, sink2file=FALSE, add.jitter=FALSE) {
   if (!is.null(cpus)) {
      warning("use threads instead of cpus! cpus will be supressed in the newer versions\n")
      threads<-cpus
   }
   
   	all.out$Results<-list()
    all.out$Results$SD<-vector(mode = "double")
    all.out$Results$LL<-vector(mode = "double")

  
   if (threads!=1){
   parallel=TRUE
   Possible.threads=parallel::detectCores()
   if (threads<=0) Threads=max(Possible.threads+threads, 1)
   if (threads>0) Threads=min(Possible.threads,threads)
    cat("creating cluster with", Threads, "threads")
    mycl <- parallel::makeCluster(Threads, type=cluster.type)
    parallel::clusterSetRNGStream(mycl)
	parallel::clusterEvalQ(mycl, library("FLightR")) 
	cat('   Done\n')

    tryCatch(Res<- pf.run.parallel.SO.resample(in.Data=all.out, threads=Threads, nParticles=nParticles, known.last=known.last, precision.sd=precision.sd, behav.mask.low.value=behav.mask.low.value, k=k, parallel=parallel, plot=FALSE, existing.cluster=mycl, cluster.type=cluster.type, a=a, b=b, L=L, sink2file=sink2file, adaptive.resampling=adaptive.resampling, RStudio=FALSE, check.outliers=check.outliers), finally = parallel::stopCluster(mycl))
	} else {	
	mycl=NA
    parallel=FALSE
    Res<- pf.run.parallel.SO.resample(in.Data=all.out, threads=Threads, nParticles=nParticles, known.last=known.last, precision.sd=precision.sd, behav.mask.low.value=behav.mask.low.value, k=k, parallel=parallel, plot=FALSE, existing.cluster=mycl, cluster.type=cluster.type, a=a, b=b, L=L, sink2file=sink2file, adaptive.resampling=adaptive.resampling, RStudio=FALSE, check.outliers=check.outliers)
	}
    # Part 2. Creating matrix of results.
    #cat("creating results matrix \n")
    #All.results.mat<-return.matrix.from.char(Res$All.results)
	all.out$Results<-list()
	all.out$Results$outliers <- Res$Results$outliers
	all.out$Results$tmp.results<-Res$Results$tmp.results
    # Part 2a. Estimating log likelihood
    LL<- get.LL.PF(all.out, Res$Points)
    cat("+----------------------------------+\n")
    cat("|     estimated negative Log Likelihood is",  LL, "\n")
    cat("+----------------------------------+\n")
    #save(LL, file=paste(LL, "time", format(Sys.time(), "%H-%m"), ".RData"))
	# Part 2b comparing the likelihood with previous estimate
    # now we need to save it..
      

	  # Part 3. Updating proposal
      cat("estimating results object\n")
      all.out.old<-all.out
      all.out<- get.coordinates.PF(Res$Points, all.out, add.jitter=add.jitter)
      Movement.parameters<- estimate.movement.parameters(Res$Trans, all.out, fixed.parameters=NA, a=a, b=b, parallel=parallel, existing.cluster=mycl, nParticles=nParticles)
	  
	all.out$Results$Movement.results=Movement.parameters$Movement.results
	all.out$Results$Transitions.rle=Movement.parameters$Transitions.rle	
	
	  all.out$Results$LL<-LL
	  
	all.out$Results<-list(

        #Final.Means=cbind(all.out$Results$Final.Means[-1,],
		#time=all.out$Indices$Matrix.Index.Table$time),
		Quantiles=cbind(all.out$Results$Quantiles[-1,],
		time=all.out$Indices$Matrix.Index.Table$time),
		Movement.results=all.out$Results$Movement.results,
		outliers=all.out$Results$outliers,
		LL=all.out$Results$LL,
		SD=all.out$Results$SD,
		Points.rle=all.out$Results$Points.rle[-1],
		Transitions.rle=all.out$Results$Transitions.rle,
		tmp.results=all.out$Results$tmp.results)
	# [-1] was added in ver 0.3.5, and it makes schedules correct.
	# still there are problems with transitions - the first transition we have is from the starting point and it has no Dawn or Dusk! So one should remove that one if he/she wants to pair these transitions with the time
    rm(Res)
    #rm(All.results.mat)
    # plotting resuls
    if (plot) {
	   lazy.result.plot(all.out)
    }
    gc()
  
  all.out$Spatial$tmp<-NULL
  
  cat("DONE!\n")
  return(all.out)
}

generate.points.dirs<-function(x , in.Data, Current.Proposal, a=45, b=500) {
  # this function is needed to generate new points - it works as from input point Index and biological proposal
  ################
  # x has 3 columns
  # Index
  # number
  # nMoving
 # here is the important change for the mask started in vesion 3.0
 # the new idea is that we wnat to set probabilities of stop to 0 in case the bird is flying over the restricted habitat..
  
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  
  #if (!is.null(in.Data$Spatial$tmp$dDistance))in.Data$Spatial$tmp$Distance<-attach.big.matrix(in.Data$Spatial$tmp$dDistance)
  #if (!is.null(in.Data$Spatial$tmp$dAzimuths))in.Data$Spatial$tmp$Azimuths<-attach.big.matrix(in.Data$Spatial$tmp$dAzimuths)
  
  if (x[[3]]>0) {
	# here is the addition of clever mask
	#if (in.Data$Spatial$Grid[x[[1]],3]==0) x[[3]]=x[[2]]
	# end of addition
	#Dists.distr<-in.Data$Spatial$tmp$Distance[x[[1]],]	
	Dists.distr<- sp::spDists(in.Data$Spatial$Grid[x[[1]], c(1,2), drop=FALSE] ,  in.Data$Spatial$Grid[,c(1,2)], longlat=TRUE)
	
    Dists.probs<-truncnorm::dtruncnorm(as.numeric(Dists.distr), a=a, b=b, Current.Proposal$M.mean, Current.Proposal$M.sd)
    ###
    #fields::image.plot(fields::as.image(Dists.probs, x=in.Data$Spatial$Grid, nrow=50, ncol=50))
    #
    if (Current.Proposal$Kappa>0) {
      #Angles.dir<-in.Data$Spatial$tmp$Azimuths[x[[1]],]
      Angles.dir<-maptools::gzAzimuth(from=in.Data$Spatial$Grid[,c(1,2)], to=in.Data$Spatial$Grid[x[[1]], c(1,2), drop=FALSE], type="abdali")
	  
      Angles.probs<-as.numeric(suppressWarnings(circular::dvonmises(Angles.dir/180*pi, mu=Current.Proposal$Direction/180*pi, kappa=Current.Proposal$Kappa)))
      Angles.probs[is.na(Angles.probs)]<-0
    }
    else {
      Angles.probs<-0.1591549
    }
    # this is neede to catch an error 
    if (any(c(Angles.probs, Dists.probs)<0)) {
      cat("PF produced weird probs \n")
      tmp.out<-list(Angles.probs=Angles.probs, Dists.probs=Dists.probs, Current.Proposal=Current.Proposal, Data=x)
      save(tmp.out, file=paste("tmp.out", round(stats::runif(n=1,min=1, max=1000)), "RData", sep="."))
      Biol.proposal<-pmax.int(Dists.probs*Angles.probs, 0)
    } else {
      Biol.proposal<-Dists.probs*Angles.probs
    }
    pos.biol<-suppressWarnings(sample.int(length(Biol.proposal), size=x[[3]], replace=TRUE, prob=Biol.proposal))
    # ok, now we want to return more complicated stuff - final indices!
    return(resample(as.integer(c(pos.biol, rep(x[[1]],(x[[2]]-x[[3]]))))))
  } else {
    return(as.integer(rep(x[[1]],(x[[2]]))))
  }
}

pf.run.parallel.SO.resample<-function(in.Data, threads=2, nParticles=1e6, known.last=TRUE, precision.sd=25, behav.mask.low.value=0.01, k=NA, parallel=TRUE, plot=TRUE, existing.cluster=NA, cluster.type="PSOCK", a=45, b=500, sink2file=FALSE, L=25, adaptive.resampling=0.5, RStudio=FALSE, check.outliers=FALSE) {
  ### to make algorhythm work in a fast mode w/o directional proposal use k=NA
  if (sink2file & !RStudio)  sink(file=paste("pf.run.parallel.SO.resample", format(Sys.time(), "%H-%m"), "txt", sep="."))
  if (sink2file & RStudio) sink()
  #### if save.rle=TRUE function will save a out.rle object in out.rle.RData file in working directory BUT to make it work There have to be out.rle column in Main.Index that will contain true for the twilights where results needed.
  ####
  
  if (any(in.Data$Spatial$Behav.mask!=1)) { 
	smart.filter=TRUE
	cat("smart filter is ON\n")
	} else {
	smart.filter=FALSE
	cat("smart filter is OFF\n")
	}
  # so, the idea here will be that we don't need to create these complicated Indexes..
  in.Data.short<-list(Indices=in.Data$Indices,  Spatial=in.Data$Spatial)
  in.Data.short$Spatial$Behav.mask<-NULL
  in.Data.short$Spatial$Phys.Mat<-NULL
  
  Parameters<-list(in.Data=in.Data.short, a=a, b=b)
  #if (!is.null(Parameters$in.Data$Spatial$tmp$dDistance)) Parameters$in.Data$Spatial$tmp$Distance<-attach.big.matrix(Parameters$in.Data$Spatial$tmp$dDistance)
  #if (!is.null(Parameters$in.Data$Spatial$tmp$dAzimuths)) Parameters$in.Data$Spatial$tmp$Azimuths<-attach.big.matrix(Parameters$in.Data$Spatial$tmp$dAzimuths)
  if (parallel) {
    if (length(existing.cluster)==1) {
      mycl <- parallel::makeCluster(threads, type=cluster.type)
      parallel::clusterSetRNGStream(mycl)
      ### we don' need to send all parameters to node. so keep it easy..
      # cleaning dataset
    }
    else {
      mycl<-existing.cluster
    }
    parallel::clusterExport(mycl, "Parameters", envir=environment())
	#if (!is.null(in.Data$Spatial$tmp$dDistance))  {parallel::clusterEvalQ(mycl, {Parameters$in.Data$Spatial$tmp$Distance <- attach.big.matrix(Parameters$in.Data$Spatial$tmp$dDistance);1})
	#}
	#if (!is.null(in.Data$Spatial$tmp$dAzimuths))  parallel::clusterEvalQ(mycl, { Parameters$in.Data$Spatial$tmp$Azimuths <- attach.big.matrix(Parameters$in.Data$Spatial$tmp$dAzimuths);1})
  }
  
  
  #########################################
  # part from the main function
  
  in.Data$Spatial$Behav.mask[in.Data$Spatial$Behav.mask<behav.mask.low.value]<-behav.mask.low.value
  in.Data$Spatial$Phys.Mat<-in.Data$Spatial$Phys.Mat*in.Data$Spatial$Behav.mask
  #=========================================
  # get value for density for angles
  if (!is.na(k)) {
    D.kappa<-as.numeric(suppressWarnings(circular::dvonmises(0, mu=0, kappa=k)))
  }
  else {
    D.kappa<-as.numeric(suppressWarnings(circular::dvonmises(0, mu=0, kappa=0)))
  }
  # stage 1. create burn-in set
  #Last.Particles<-rep(as.integer(in.Data$Spatial$start.point), nParticles)
  #Last.Last.Particles<-Last.Particles
  
  # here I'll save results stack.
  #All.results
  
  Results.stack<-as.matrix(rep(as.integer(in.Data$Spatial$start.point), nParticles))
  
  Weights.stack<-as.matrix(rep(1/nParticles, nParticles))
  
  New.weights<-rep(1/nParticles, nParticles)
  #All.results<-NULL
  in.Data$outliers<-c()
	  #Trans<-vector(mode = "list", length = nrow(in.Data$Indices$Main.Index))
	  Trans<-vector(mode = "list")
  	if (check.outliers) {
	  in.Data$AB.distance<-c()
	  #in.Data$AC.distance<-c()
	  in.Data$AC.distance2<-c()
	  in.Data$Dif.ang-c()
	  #in.Data$BC.mean<-c()
  }
  
  propagate.particles<-function(Last.Particles, Current.Proposal, parallel=TRUE, Parameters, mycl) {
    Order.vector<-order(Last.Particles)
    Particles.rle<-bit::intrle(as.integer(Last.Particles[Order.vector]))
    if (is.null(Particles.rle)) 	Particles.rle<-rle(as.integer(Last.Particles[Order.vector]))
    Last.State<-cbind(Particles.rle$values, Particles.rle$length)
    Last.State<-cbind(Last.State,nMoving=stats::rbinom(dim(Last.State)[[1]], size=Last.State[,2], prob=Current.Proposal$Decision))
    Last.State.List<-split(Last.State, row(Last.State))
    nSeq<-nrow(Last.State)
    if (nSeq==1 | (!parallel)) {
      #cat("non.parallel\n")
      #New.Points<-t(lapply(Last.State, 1, pf.par.internal, Current.Proposal))			
      New.Points<-lapply(Last.State.List, FUN=function(x,Current.Proposal, Parameters) do.call(generate.points.dirs, c(x=list(x), Parameters, Current.Proposal=list(Current.Proposal))), Current.Proposal, Parameters)
      #New.Points<-lapply(Last.State.List, pf.par.internal, Current.Proposal)
    }
    else {
      cat(" initial diversity is ", nSeq)
      #New.Points<-clusterApplyLB(mycl, Last.State.List,  pf.par.internal, Current.Proposal)
      New.Points<-parallel::clusterApply(mycl, Last.State.List,  pf.par.internal, Current.Proposal)
	  cat("  Done\n")
    }
    #=======================================================
    New.Particles<-unlist(New.Points)[sort.list(Order.vector, na.last = NA, method =  "quick")]
    return(New.Particles)
  }
  
  
  ## for ver 1.6 
  #Prev.Weights<-rep(1, nParticles) 
  ResampleCount<-0
  steps.from.last<-2
  total_length<-nrow(in.Data$Indices$Main.Index)
  for (Time.Period in 1:total_length) {
    steps.from.last=steps.from.last+1

	cat("\n\n##########################\n     Time.Period", Time.Period, "of", total_length, "\n")
    #cat("prep. data:")
    Current.Proposal<-in.Data$Indices$Matrix.Index.Table[in.Data$Indices$Main.Index$Biol.Prev[Time.Period],]
    #=======================================
    cat("generating new particles")
    New.Particles<-propagate.particles(Last.Particles=Results.stack[,ncol(Results.stack)], Current.Proposal=Current.Proposal, parallel=parallel, Parameters=Parameters, mycl=mycl)
    
    #=====================================================
    # resampling step 
    # we need to estimate weights of resulting points and then resample them proportionally to weights..
    Current.Phys.Weights<-in.Data$Spatial$Phys.Mat[New.Particles,Time.Period]
    
    if (!is.na(k) & Time.Period >1) {
      get.directional.weights<-function(in.Data, Last.Last.Particles, Last.Particles, New.Particles, k, parallel, mycl, D.kappa) {
        #===================================
        #Prev.Dirs<-in.Data$Spatial$tmp$Azimuths[cbind(Last.Last.Particles, Last.Particles)]/180*pi
		
        Prev.Dirs<-apply(matrix(c(Last.Last.Particles, Last.Particles), ncol=2), 1, dir_fun, in.Data)/180*pi
		
        #New.Dirs<-in.Data$Spatial$tmp$Azimuths[cbind(Last.Particles,New.Particles)]/180*pi
        New.Dirs<-apply(matrix(c(Last.Particles, New.Particles), ncol=2), 1, dir_fun, in.Data)/180*pi
		
        FromTo=matrix(c(Prev.Dirs, New.Dirs), ncol=2)
        if (parallel) {Angles.probs<-parallel::parApply(mycl, FromTo, 1, flightr.dvonmises, mykap=k)
        } else {
          Angles.probs<-apply(FromTo,1, FUN=function(x, k) as.numeric(suppressWarnings(circular::dvonmises(x[[2]], mu=x[[1]], kappa=k))), k=k)
        }
        Angles.probs[is.na(Angles.probs)]<-D.kappa
        #cat("\n min.Kappa:", min(Angles.probs), ", max.Kappa", max(Angles.probs), "\n")
        return(Angles.probs)
      }
      
      Angles.probs<-get.directional.weights(in.Data, Results.stack[,ncol(Results.stack)-1], Results.stack[,ncol(Results.stack)], New.Particles, k, parallel, mycl, D.kappa)
      Current.Weights<- Current.Phys.Weights*Angles.probs
    }
    else {Current.Weights<-Current.Phys.Weights*D.kappa}

	if (smart.filter) {
	#=========================
	# here we should introduce clever mask..
	# this will be something like that
	Index.Stable<-which(New.Particles == Results.stack[,ncol(Results.stack)])
	Current.Weights[Index.Stable]<-Current.Weights[Index.Stable]*(in.Data$Spatial$Behav.mask[New.Particles][Index.Stable])
	}
	#=======================================================
	# now I want to add 3.3
	# the idea is following - we should compare distances from the last t-2 to t-1 and to t
	if (Time.Period>2 & check.outliers) {
	# here I want to try directiondl outliers..
	# so we need to pick up particles that moved..
	# 
	#=================================================================

	AB.distance<-stats::weighted.mean(	sp::spDists(in.Data$Spatial$Grid[Results.stack[,(ncol(Results.stack)-1)], c(1,2), drop=FALSE], in.Data$Spatial$Grid[Results.stack[,ncol(Results.stack)], c(1,2), drop=FALSE], longlat=TRUE, diagonal=TRUE), Weights.stack[,ncol(Weights.stack)])

	AC.distance2<-	stats::weighted.mean(sp::spDists(in.Data$Spatial$Grid[Results.stack[,(ncol(Results.stack)-1)], c(1,2), drop=FALSE], in.Data$Spatial$Grid[New.Particles, c(1,2), drop=FALSE], longlat=TRUE, diagonal=TRUE), Weights.stack[,ncol(Weights.stack)]*Current.Weights)	
	
	cat("AB.distance:", round(AB.distance, 2), "\n")
	cat("AC.distance2:", round(AC.distance2, 2), "\n")

	Dif.ang=180
	if (AB.distance>50) { # go for angles only if disctances are high!
	resample <- function(x, ...) x[sample.int(length(x), ...)]
	# the AB ones wil have folowing..

	BA.dir<-apply(Results.stack[,ncol(Results.stack):(ncol(Results.stack)-1), drop=FALSE], 1, dir_fun, in.Data)
	
	BA.moved<-which(!is.na(BA.dir))
	BA.mean<-circular::mean.circular(circular::circular(resample(BA.dir[BA.moved], replace=TRUE,prob=Weights.stack[,ncol(Weights.stack)][BA.moved]), units="degrees"), na.rm=TRUE)
	#cat(BA.mean)	#BC.dir<-in.Data$Spatial$tmp$Azimuths[cbind(Results.stack[,ncol(Results.stack)-1], New.Particles)]
	BC.dir<-apply(matrix(c(Results.stack[,ncol(Results.stack)], New.Particles), ncol=2), 1, dir_fun, in.Data)
	
	BC.moved<-which(!is.na(BC.dir))
	BC.mean<-circular::mean.circular(circular::circular(resample(BC.dir[BC.moved], replace=TRUE, prob=(Weights.stack[,ncol(Weights.stack)]*Current.Weights)[BC.moved]), units="degrees"), na.rm=TRUE)
	#cat(BC.mean)
	dif.ang<-function(x,y) {
	# this function provides the minimum angle between two angles..
	y=y*pi/180
	x=x*pi/180
	z<-c(y-x, y-x+2*pi, y-x-2*pi)
	z[which.min(abs(z))]*180/pi
	}
	if (!is.null(BA.mean) & !is.null(BC.mean)) {
	Dif.ang<-dif.ang(BA.mean, BC.mean)
	cat("anglular change", round(Dif.ang,2), "\n" )
	}
	}	
	
#=================================
# now I want to temporarily save this ..
#if (is.null(BA.mean)) BA.mean=NA
#if (is.null(BC.mean)) BC.mean=NA
  in.Data$AB.distance<-c(in.Data$AB.distance, AB.distance)
  #in.Data$AC.distance<-c(in.Data$AC.distance, AC.distance)
  in.Data$AC.distance2<-c(in.Data$AC.distance2, AC.distance2)
  in.Data$Dif.ang<-c(in.Data$Dif.ang, Dif.ang)
  #in.Data$BC.mean<-c(in.Data$BC.mean, BC.mean)
	# outlier detection 10
	# oulier detection
	#if (AB.distance>(AC.distance2*1.5) | (abs(Dif.ang)<90 & AB.distance>50)) {
	
	#if ( steps.from.last>2 & ((AB.distance>50 & AB.distance>(AC.distance2*1.3))| (abs(Dif.ang)<90))) {
	if ( steps.from.last>2 & ((AB.distance>50 & AB.distance>(AC.distance2*1.3))| (AB.distance>AC.distance2 & abs(Dif.ang)<100))) {
	steps.from.last=0
	cat("outlier number", length(in.Data$outliers)+1, "detected! removing observations from ", Time.Period-1, "!\n")
	in.Data$outliers<-c(in.Data$outliers, Time.Period-1)
	# and now the main question - how should we remove the weights??
	# probably the easiest way would be to add 1/nParticles into the last column..
	Weights.stack[,ncol(Weights.stack)]<-rep(1/nParticles, nParticles)
	# ok and now I have to change the resample to t-1
	}
	}
	
	
    #=====================================================
    # here is the change from 1.6 - Adding cumulative weights now..
	# and this came back at the version 3.2!
    Current.Weights.with.Prev.mat<-cbind(Weights.stack, Current.Weights)
    
	#rowProds <- function(a) exp(rowMeans(log(a)))
    rowProds <- function(a) exp(rowSums(log(a)))
	
	#Current.Weights.with.Prev<-	pmax(rowProds(Current.Weights.with.Prev.mat), 1e-323)
	
	# here is the place - I'll check weights without use of the last stage information! 
	Current.Weights.with.Prev<-	pmax(rowProds(Weights.stack), 1e-323)
	#
    #================================================================
    # ver 1.7. 
    # ADAPTIVE RESAMPLING
    if (adaptive.resampling!=1) {
#cat("min. current.weights.with.prev=", min(Current.Weights.with.Prev), "\n")
#cat("max. current.weights.with.prev=", max(Current.Weights.with.Prev), "\n")
#cat("sq.of.sums:",  sum(Current.Weights.with.Prev)^2, "\n")
#cat("sumofsq:", sum((Current.Weights.with.Prev)^2), "\n")
      ESS<-(sum(Current.Weights.with.Prev)^2)/sum((Current.Weights.with.Prev)^2)
      cat("ESS is ", ESS)
if (is.na(ESS)) {
	ESS=1
	save(Current.Weights.with.Prev, Current.Weights.with.Prev.mat, Current.Weights, file="tmp.RData")
	}		
    } else {
      ESS<-1
    }
    if ((ESS<(nParticles*adaptive.resampling) & Time.Period>1) | Time.Period == nrow(in.Data$Indices$Main.Index)) {	
	  if (length(unique(Current.Weights.with.Prev))==1) Current.Weights.with.Prev<-rep(1, nParticles)
	  Rows<-try(sample.int(nParticles, replace = TRUE, prob = Current.Weights.with.Prev))
			#if (class(Rows)=="try-error") {
			#	print(Current.Weights.with.Prev) 
			#	print(Rows)
			#	stop("weird..\n")														 
		#}#
		#}#
      ResampleCount<-ResampleCount+1
      cat(" - resampling", ResampleCount, "\n")
      Results.stack<-Results.stack[Rows,]
	  Weights.stack<-as.matrix(rep(1/nParticles, nParticles))
    } else {
      Rows<-1:nParticles # no resampling
      cat("\n")
    }
	New.weights<-Current.Weights[Rows]

    #=====================================================
    #               Smoothing
    #====================================================
    
    #======================================================
    # start OF SMOOTHING POINTS ESTIMATION
    
    #  ==================================
    # ver. 1.8 - now use one function for all particle propagations:
    #	cat("  smoothing new particles")
    #
    #	New.Particles.MH<-propagate.particles(Last.Particles= Results.stack[,ncol(Results.stack)], Current.Proposal=Current.Proposal, parallel=parallel, Parameters=Parameters, mycl=mycl)
    #
    #	Current.Phys.Weights.MH<-in.Data$Spatial$Phys.Mat[New.Particles.MH,Time.Period]
    #	
    #   if (!is.na(k) & Time.Period > 1) {
    #	  Angles.probs<-get.directional.weights(in.Data, Results.stack[,ncol(Results.stack)-1], Results.stack[,ncol(Results.stack)], New.Particles.MH, k, parallel, mycl, D.kappa)
    #      Current.Weights.MH<- Current.Phys.Weights.MH*Angles.probs
    #    }
    #    else {Current.Weights.MH<-Current.Phys.Weights.MH*D.kappa}
    
    #   end of smoothing estimation    
    #============================================
    #============================================
    
    # now I need to compare resampled results with new resuklts
    #according to MH
    #	v<-runif(nParticles)
    #	Accepted.Particles<-New.Particles.MH
    #	Ratio<-Current.Weights.MH/(Current.Weights[Rows])
    #	Ratio[is.na(Ratio)]<-1 # NA is caused by 0/0
    #	Eq<-!as.logical(floor(pmin(1, Ratio)-v)+1)
    #	Accepted.Particles[Eq]<-New.Particles[Rows][Eq]
    #	New.weights<-Current.Weights.MH
    #	New.weights[Eq]<-Current.Weights[Rows][Eq]
    #	cat("smoothed ", nParticles-sum(Eq), "Particles \n")
    
    #============================================
    # this line is needed instead of all MH stuff
    Accepted.Particles<-New.Particles[Rows]
    
    #All.results<-paste(All.results, Accepted.Particles, sep=".")
    ## adding points to the results matrix
    Results.stack<-cbind(Results.stack, Accepted.Particles)
    #==============================================================
    # 1.6
    Weights.stack<-cbind(Weights.stack, New.weights)
    
    #Prev.Weights<-(Prev.Weights[Rows]) * New.weights
    #Prev.Weights<-Prev.Weights/mean(Prev.Weights)
    # end of 1.6 addition
    #===================================================
    
    #Last.Last.Particles<-Last.Particles[Rows]
    #Last.Particles<-Accepted.Particles
    
    #	if (Time.Period==2) Particles.at.2<-Accepted.Particles
    #	if (Time.Period>2) {
    #		Particles.at.2<-Particles.at.2[Rows]
    cat("   unique P in point leaving stack", length(unique(Results.stack[,1])), "\n")
    #	}
    
    if (Time.Period<=L) cat("creating stack\n")
	if (Time.Period==L) 	Points<-vector(mode = "list")
    if (Time.Period>L) {
	# the new idea is that we could skip the saving all results and save just points and transitions - we are not outputting them anyways... THis will help avoiding the sort of All.results, that proved to be very slow..
      # save points
      #if (is.null(Points))  Points<-vector(mode = "list")
	  Rle<-bit::intrle(sort.int(Results.stack[,1], method="quick"))
      if (is.null(Rle)) Rle<-rle(sort.int(Results.stack[,1], method="quick"))
	  Points[length(Points)+1]<-list(Rle)
      #  All.results<-paste(Results.stack[,1], sep=".")
      #} else {
      #  All.results<-paste(All.results, Results.stack[,1], sep=".")
      #}
      # save transitions
      Trans[[Time.Period-L]]<-get.transition.rle(Results.stack[,1], Results.stack[,2])
      # clean Results.stack
      Results.stack<-Results.stack[,-1]
      # clean Weights.stack
      Weights.stack<-as.matrix(Weights.stack[,-1])
cat("******************\n")
    }
  }
  #####################################
  #save(All.results, file="All.results.usmoothed.RData")
  # now we need to add final point!
  
  if (known.last) {
    Results.stack<-pf.final.smoothing(in.Data, Results.stack, precision.sd=precision.sd, nParticles=nParticles, last.particles=Results.stack[,ncol(Results.stack)])
  }
  
  if (!is.list(Points)) Points<-vector(mode = "list")
  # and here we need to add a thing that will finish All.results and Trans from the points that are still in the stack
  cat("adding last points form the stack to the resutls\n")
  Length<-ncol(Results.stack)

  for (rest in 1:Length) {
    # save points
	  Rle<-bit::intrle(sort.int(Results.stack[,rest], method="quick"))
      if (is.null(Rle)) Rle<-rle(sort.int(Results.stack[,rest], method="quick"))
	  Points[length(Points)+1]<-list(Rle)

	#All.results<-paste(All.results, Results.stack[,rest], sep=".")
    if (rest<Length) {
      # save transitions
      #Trans[[Time.Period-L+rest]]<-get.transition.rle(Results.stack[,rest], Results.stack[,rest+1])
      Trans[[length(Trans)+1]]<-get.transition.rle(Results.stack[,rest], Results.stack[,rest+1])
    }
  }
  #if (parallel)   parallel::clusterEvalQ(mycl, rm(Parameters)) 
  #if (length(existing.cluster)==1) parallel::stopCluster(cl = mycl)
  if (sink2file) sink()
  tmp.results<-list(AB.distance=in.Data$AB.distance, AC.distance2=in.Data$AC.distance2, Dif.ang=in.Data$Dif.ang)

  return(list(Points=Points, Trans=Trans, Results=list(outliers=in.Data$outliers, tmp.results=tmp.results)))
}


get.coordinates.PF<-function(Points, in.Data, add.jitter=FALSE) {
  #library("aspace")
  # this function will extract point coordinates from the output matrix.. 
  # the question is do we need only mean and sd or also median and quantiles?
  # I will start from mean and SD
  #plot(c(min(in.Data$Spatial$Grid[,1]),max(in.Data$Spatial$Grid[,1])), c(min(in.Data$Spatial$Grid[,2]),max(in.Data$Spatial$Grid[,2])), type="n")
  #log <- capture.output({
  #  Means=aspace::calc_box(id=1,  points=in.Data$Spatial$Grid[inverse.rle(Points[[1]]),1:2])
  # 
  #for (i in 2:length(Points)) {
  #  Means[i,]=aspace::calc_box(id=i,  points=in.Data$Spatial$Grid[inverse.rle(Points[[i]]), 1:2])
    #plot_box(plotnew=FALSE, plotpoints=FALSE)
  #}
	#})
   if (length(in.Data$Results)==0) in.Data$Results<-list()
  #in.Data$Results$Final.Means<-Means
  
  #############
  # new part for medians
  cat("estimating quantiles for positions\n")
	
		# check whether Grid was over dateline:
	overdateline<-ifelse(attr(in.Data$Spatial$Grid, 'left')>	attr(in.Data$Spatial$Grid, 'right'), TRUE, FALSE)

	
	Quantiles<-c()
	CIntervals<-c()
	cur_Grid<-in.Data$Spatial$Grid
	cur_Grid[cur_Grid[,1]<0,1]<-cur_Grid[cur_Grid[,1]<0,1]+360
	for (i in 1:length(Points)) {
	

	#Mode_cur<-	cur_Grid[Points[[i]]$values[which.max(Points[[i]]$lengths)],1]
	#cur_Grid[,1]<-ifelse(cur_Grid[,1]<Mode_cur-180, cur_Grid[,1]+360, cur_Grid[,1])
	Quantiles<-rbind(Quantiles, c(summary(cur_Grid[inverse.rle(Points[[i]]),2]), Mode=cur_Grid[Points[[i]]$values[which.max(Points[[i]]$lengths)],2], summary(cur_Grid[inverse.rle(Points[[i]]),1]), Mode=cur_Grid[Points[[i]]$values[which.max(Points[[i]]$lengths)],1]))
	CIntervals<-rbind(CIntervals, c(stats::quantile(cur_Grid[inverse.rle(Points[[i]]),2], probs = c(0.025, 0.975)), stats::quantile(cur_Grid[inverse.rle(Points[[i]]),1], probs = c(0.025, 0.975))))
	}
	Quantiles<-as.data.frame(Quantiles)
	names(Quantiles)[1:6]<-paste(names(Quantiles)[1:6], "lat", sep="")
	names(Quantiles)[8:13]<-paste(names(Quantiles)[8:13], "lon", sep="")

	###########
	# doing jitter first
	if (add.jitter) {
	cat("adding jitter to medians\n")
	jitter_coords<-get.coords.jitter(in.Data)
	if (!is.null(jitter_coords)) {
	Quantiles$MedianlonJ<-jitter_coords[,1]
	Quantiles$MedianlatJ<-jitter_coords[,2]
	}
	} else {
	Quantiles$MedianlonJ<-Quantiles$Medianlon
	Quantiles$MedianlatJ<-Quantiles$Medianlat
	}
	names(Quantiles)<-gsub("\\s","", names(Quantiles))
	names(Quantiles)<-gsub("1","F", names(Quantiles))
	names(Quantiles)<-gsub("3","T", names(Quantiles))
	
	cat("adding 95% credibility intervals to medians\n")
	Quantiles$LCI.lat<-CIntervals[,1]
	Quantiles$UCI.lat<-CIntervals[,2]
	Quantiles$LCI.lon<-CIntervals[,3]
	Quantiles$UCI.lon<-CIntervals[,4]
	  
	if (overdateline) {
	   Columns<-c(8:15, 19, 20)
       for (i in Columns) {
           Quantiles[Quantiles[,i]>180,i]<-Quantiles[Quantiles[,i]>180,i]-360
       }	   
	}
	in.Data$Results$Quantiles<-Quantiles
	
    in.Data$Results$Points.rle<-Points
  return(in.Data)
}


estimate.movement.parameters<-function(Trans, in.Data, fixed.parameters=NA, a=45, b=500, parallel=FALSE, existing.cluster=NULL, estimatetruncnorm=FALSE, nParticles=1e6) {
  mycl=existing.cluster
  
  # old name update.proposal.PF
  # main function for proposal update
  # from the version 1.5 this function will do estimation of distance and it's SD..

  #=========================================
  # here is another place to change...
  # ver. 1.8
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #Trans<-vector(mode = "list", length = (dim(output.matrix)[2]-1))
  #cat("   extracting indexes\n")
  #for (i in 1:(dim(output.matrix)[2]-1)) {
  #  Trans[[i]]<-get.transition.rle(output.matrix[,i], output.matrix[,i+1])
  #}
  cat("   estimating distances\n")
  #####   let's try to get distance distribution:
  Distances<-Trans
  dist.fun<-function(x) {
   #in.Data$Spatial$tmp$Distance[x%/%1e5, x%%1e5]
    sp::spDists(in.Data$Spatial$Grid[x%/%1e5, c(1,2), drop=FALSE], in.Data$Spatial$Grid[x%%1e5, c(1,2), drop=FALSE],longlat=TRUE, diagonal=TRUE)
  }
  for (i in 1:length(Trans)) {
    Distances[[i]]$values<-sapply(Trans[[i]]$values, FUN=function(x) dist.fun(x))
  }
  cat("   estimating directions\n")
  #ok, now we want to get directions
  
  Directions<-Trans
  for (i in 1:length(Trans)) {
    Movement_Points<-matrix(c(Trans[[i]]$values%/%1e5, Trans[[i]]$values%%1e5), ncol=2)
    Directions[[i]]$values<-apply(Movement_Points, 1, FUN= dir_fun, in.Data)
  }
  cat("   estimating mean directions and kappas\n")
  
  Mean.Directions<-unlist(lapply(Directions, FUN=function(x) CircStats::circ.mean(circular::circular(inverse.rle(list(lengths=x$lengths[!is.na(x$values)], values=x$values[!is.na(x$values)]*pi/180))))*180/pi))
  #plot(Mean.Directions)
  cat("   estimating kappas\n") #CircStats::est.kappa
  Kappas<-unlist(lapply(Directions, FUN=function(x) CircStats::est.kappa(inverse.rle(list(lengths=x$lengths[!is.na(x$values)], values=x$values[!is.na(x$values)]*pi/180)))))
  
  
  # now we want to get mean ditance
  cat("   estimating mean dists\n")
  #
  #cat("   estimating mean and SD to report dists SD\n")
  # attempt to go just for simple distance estimation...
  Mean2report<-unlist(lapply(Distances, FUN=function(x) mean(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))
  SD2report<-unlist(lapply(Distances, FUN=function(x) stats::sd(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))
  cat("   estimating probs of migration\n")
  # ok, now we want to get parameters for distances
  #
  Probability.of.migration<-unlist(lapply(Distances, FUN=function(x) sum(x$lengths[x$values!=0])))/nParticles
  
  ## 
    if (estimatetruncnorm) {
  Mean.and.Sigma<-lapply(Distances, FUN=function(x)  mu.sigma.truncnorm(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])), a=a, b=b))
  #}
  Mean.Dists<-sapply(Mean.and.Sigma, "[[", i=1)
  cat("   estimating dists SD\n")
  Mean.SD<-sapply(Mean.and.Sigma, "[[", i=2)
  } else {
  Mean.Dists<-Mean2report
  Mean.SD<-SD2report
  }

  #unlist(lapply(Distances, FUN=function(x) mean(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))
  #plot(Mean.Dists)
  cat("   estimating median dists\n")
  Median.Dists<-unlist(lapply(Distances, FUN=function(x) stats::median(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))

  #unlist(lapply(Distances, FUN=function(x) sd(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))
  #plot(Kappas)
  
  
  
  
  cat("   creating output")
  ################
  # create new proposal from the posteriors
  #New.Matrix.Index.Table<-data.frame(Direction=Mean.Directions, M.mean=Mean.Dists, M.sd=Mean.SD, Decision=Probability.of.migration, Kappa=Kappas)
  ################
  
  Movement.results<-in.Data$Indices$Matrix.Index.Table
  #all.arrays.object<-in.Data
  
  #----------------
  # actually I do not want to update this - so I rather save it outside if possible..
  Movement.results$Decision<-Probability.of.migration
  Movement.results$Direction<-Mean.Directions
  Movement.results$Kappa<-Kappas
  Movement.results$M.mean<-Mean.Dists
  Movement.results$M.sd<-Mean.SD
  Movement.results$M.medians<-Median.Dists
  Movement.results$Mean2report<-Mean2report
  Movement.results$SD2report<-SD2report
  
  #all.arrays.object$Indices$Matrix.Index.Table<-New.Matrix.Index.Table
  
  #all.arrays.object$Indices$Main.Index$Biol.Prev=1:(nrow(all.arrays.object$Indices$Matrix.Index.Table))
    
  if (is.list(fixed.parameters)) {
    # this part can be moved upper not to estimate parametrs in case they are fixed.
    if ("M.sd" %in% names(fixed.parameters)) Movement.results$M.sd<-fixed.parameters$M.sd
    
    if ("M.mean" %in% names(fixed.parameters)) Movement.results$M.mean<-fixed.parameters$M.mean
    
    if ("Kappa" %in% names(fixed.parameters)) Movement.results$Kappa<-fixed.parameters$Kappa
    
    if ("Direction" %in% names(fixed.parameters)) Movement.results$Direction<-fixed.parameters$Direction
  }
  #all.arrays.object$Indices$Main.Index$Biol.Next=2:(nrow(all.arrays.object$Indices$Matrix.Index.Table))
  cat(" DONE!\n")
  Res<-list(Movement.results=Movement.results, Transitions.rle=Trans)
  return(Res)
}


get.transition.rle=function(From, To) {
  rle(sort.int(From*1e5+To, method = "quick"))
}


mu.sigma.truncnorm<-function(x, a=45, b=500) {
# this is used optionally
  if (length(unique(x))>1) {
    tr.norm<-function(prm) {
      sum(-log(truncnorm::dtruncnorm(as.numeric(x),a=a,b=b,mean=prm[1],sd=prm[2])))
    }
    Res=try(stats::optim(c(mean(x),stats::sd(x)), tr.norm, method="BFGS"))
    if (class(Res)=="try-error") {
      save(x, Res, file="x.RData")
      Res$par<-c(mean(x), stats::sd(x))
    }
    return(c(Res$par[1], Res$par[2]))
  } else {
    return(c(mean(x), stats::sd(x)))
  }
}


coords.aeqd.jitter <- function(coords, r, n)
{
# coords is a vector of leghth 2 c(lon, lat)
# r should be in meters
# made on te basis of this:
#"http://gis.stackexchange.com/questions/121489/1km-circles-around-lat-long-points-in-many-places-in-world"
 stopifnot(length(coords) == 2)

	p = sp::SpatialPoints(matrix(coords, ncol=2), proj4string=sp::CRS("+proj=longlat +datum=WGS84"))
    aeqd <- sprintf("+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                    p@coords[[2]], p@coords[[1]])
    projected <- sp::spTransform(p, sp::CRS(aeqd))
    buffered <- rgeos::gBuffer(projected, width=r, byid=TRUE)
    lambert <- sprintf("+proj=laea +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                    p@coords[[2]], p@coords[[1]])
	buffered_eqarea <- sp::spTransform(buffered, sp::CRS(lambert))
	random_points<-sp::spsample(buffered_eqarea,n=n,type="random")
	if (is.null(random_points)) random_points<-sp::spsample(buffered_eqarea,n=n,type="random", iter=40) 
	if (is.null(random_points)) random_points<-p
    sp::spTransform(random_points, p@proj4string)
}

# wrapper for jitter
get.coords.jitter<-function(in.Data) {
	Distance<-in.Data$Spatial$tmp$Distance
	if (is.null(Distance)) Distance=sp::spDists(in.Data$Spatial$Grid[,1:2], longlat=TRUE)
	Distance<-Distance[,1]
	JitRadius<-min(Distance[Distance>0])/2*1000 # in meters
	#now I want to generate random poitns in the radius of this
	coords=cbind(in.Data$Results$Quantiles$Medianlon, in.Data$Results$Quantiles$Medianlat)
	tmp<-try(apply(coords, 1,  coords.aeqd.jitter, r=JitRadius, n=1 ))
	jitter_coords<-NULL
	if (class(tmp)!="try-error") {
	jitter_coords<-t(sapply(tmp, sp::coordinates))
	}
	return(jitter_coords)
}


get.LL.PF<-function(in.Data, data) {
  # needed to estimate log likelihood of the optimization
  L=0
  if (is.list(data)) {
  for (i in 1:(length(data)-1)) { 
      L=L+log(mean(in.Data$Spatial$Phys.Mat[inverse.rle(data[[i]]),i]))
	  }
  } else{ 
    for (i in 1:(dim(data)[2]-1)) {
    L=L+log(mean(in.Data$Spatial$Phys.Mat[data[,i],i]))
	}
  }
  #L=L/(dim(All.results.mat)[2]-1)
  #LL=-sum(log(L))
  LL=-L
  return(LL)
}


pf.par.internal<-function(x, Current.Proposal) {
  # this simple function is needed to save I/0 during parallel run. Being executed at a slave it creates inoput for the next function taking Current proposal from the slave and Parameters from the master
  new.Parameters<-c(x=list(x), get("Parameters"), Current.Proposal=list(Current.Proposal))
  Res<-do.call( generate.points.dirs, new.Parameters)
  return(Res)
}


flightr.dvonmises<-function(x, mykap) {
  return(as.numeric(suppressWarnings(circular::dvonmises(x[[2]], mu=x[[1]], kappa=mykap))))}

  
pf.final.smoothing<-function(in.Data, results.stack, precision.sd=25, nParticles=1e6, last.particles=NA) {
  # this function simply resamples final points proportionally to the distance to known finish.
  Final.point.real<-in.Data$Spatial$stop.point
  # now we want to get distances.. I'll not index it as we will do this only once..
  Final.points.modeled=last.particles
  Weights<-stats::dnorm(  sp::spDists(in.Data$Spatial$Grid[Final.points.modeled, c(1,2), drop=FALSE], in.Data$Spatial$Grid[Final.point.real, c(1,2), drop=FALSE], longlat=TRUE), mean=0, sd=precision.sd)
  Rows<- try(suppressWarnings(sample.int(nParticles, replace = TRUE, prob = Weights/sum(Weights))))
  if (class(Rows) == 'try-error') {
    cat('final smoothing failed, error data saved to the working directory - smoothing.error.RData!\n')
	save(last.particles, Weights, file='smoothing.error.RData')
    return(results.stack)
	} else {
    return(results.stack[Rows,])
	}
}

dir_fun<-function(x, in.Data) {
	  maptools::gzAzimuth(in.Data$Spatial$Grid[x[[1]], c(1,2), drop=FALSE], in.Data$Spatial$Grid[x[[2]], c(1,2), drop=FALSE], type="abdali")
}
	
dist.fun<-function(x, Result) {
    sp::spDists(Result$Spatial$Grid[x%/%1e5, c(1,2), drop=FALSE], Result$Spatial$Grid[x%%1e5, c(1,2), drop=FALSE],longlat=TRUE, diagonal=TRUE)
}

lazy.result.plot<-function(Result) {
    graphics::plot(Meanlat~Meanlon, type="p", data=Result$Results$Quantiles, pch=3, col="blue", main="mean poistions")
    graphics::points(Meanlat~Meanlon, type="p", data=Result$Results$Quantiles, pch=3, col="blue")
    graphics::lines(Meanlat~Meanlon, data=Result$Results$Quantiles, col="blue")
	wrld_simpl<-NA
    load(system.file("data", "wrld_simpl.rda", package = "maptools"))
    sp::plot(wrld_simpl, add=TRUE)
}

