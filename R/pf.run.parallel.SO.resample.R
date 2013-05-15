

pf.run.parallel.SO.resample<-function(in.Data, cpus=2, nParticles=1e6, known.last=T, precision.sd=25, sea.value=0.01, k=1, parallel=T, plot=T, existing.cluster=NA, cluster.type="SOCK", a=45, b=500, sink2file=F, L=25, adaptive.resampling=F, RStudio=F) {
  # this dunction is doing main job. it works on the secondary master
  ### to make algorhytm work in a fast mode w/o directional proposal use k=NA
  if (sink2file & !RStudio)  sink(file=paste("pf.run.parallel.SO.resample", format(Sys.time(), "%H-%m"), "txt", sep="."))
   if (sink2file & RStudio) sink()
  #### if save.rle=T function will save a out.rle object in out.rle.RData file in working directory BUT to make it work There have to be out.rle column in Main.Index that will contain true for the twilights where results needed.
  ####
  require("parallel")
  require("circular")
  require("maptools")
  data("wrld_simpl", package="maptools")
  require("truncnorm")
  require("circular")
  require("bit")
  require("geosphere")

  # so, the idea here will be that we don't need to create these complicated Indexes..
  in.Data.short<-list(Main.Index=in.Data$Main.Index , Matrix.Index.Table=in.Data$Matrix.Index.Table,  Points.Land=in.Data$Points.Land, distance= in.Data$distance, Azimuths=in.Data$Azimuths, transitions=in.Data$transitions)
  Parameters<-list(in.Data=in.Data.short, a=a, b=b)
  if (parallel) {
    if (length(existing.cluster)==1) {
      mycl <- parallel:::makeCluster(cpus, type=cluster.type)
      parallel:::clusterSetRNGStream(mycl)
      ### we don' need to send all parameters to node. so keep it easy..
      # cleaning dataset
      parallel:::clusterExport(mycl,c("generate.points.dirs", "pf.par.internal", "my.dvonmises"))
      parallel:::clusterEvalQ(mycl, library("circular")) 
      parallel:::clusterEvalQ(mycl, library("truncnorm")) 
    }
    else {
      mycl<-existing.cluster
    }
    parallel:::clusterExport(mycl, "Parameters", envir=environment())
  }
  
  
  #########################################
  # part from the main function
  
  in.Data$Geogr.proposal[in.Data$Geogr.proposal==0]<-sea.value
  in.Data$Phys.Mat<-in.Data$Phys.Mat*in.Data$Geogr.proposal
  #=========================================
  # get value for density for angles
  if (!is.na(k)) {
    D.kappa<-as.numeric(suppressWarnings(circular:::dvonmises(0, mu=0, kappa=k)))
  }
  else {
    D.kappa<-as.numeric(suppressWarnings(circular:::dvonmises(0, mu=0, kappa=0)))
  }
  # stage 1. create burn-in set
  #Last.Particles<-rep(as.integer(in.Data$start.point), nParticles)
  #Last.Last.Particles<-Last.Particles

    # here I'll save results stack.
	#All.results
	
	Results.stack<-as.matrix(rep(as.integer(in.Data$start.point), nParticles))
	
	Weights.stack<-as.matrix(rep(1, nParticles))
	
	All.results<-NULL

	Trans<-vector(mode = "list", length = nrow(in.Data$Main.Index))

  	
	propagate.particles<-function(Last.Particles, Current.Proposal, parallel=T, Parameters, mycl) {
		Order.vector<-order(Last.Particles)
		Particles.rle<-bit:::intrle(as.integer(Last.Particles[Order.vector]))
		if (is.null(Particles.rle)) 	Particles.rle<-rle(as.integer(Last.Particles[Order.vector]))
		Last.State<-cbind(Particles.rle$values, Particles.rle$length)
		Last.State<-cbind(Last.State,nMoving=rbinom(dim(Last.State)[[1]], size=Last.State[,2], prob=Current.Proposal$Decision))
		Last.State.List<-split(Last.State, row(Last.State))
		nSeq<-nrow(Last.State)
		if (nSeq==1 | (!parallel)) {
		#cat("non.parallel\n")
		#New.Points<-t(lapply(Last.State, 1, pf.par.internal, Current.Proposal))			
		New.Points<-lapply(Last.State.List, FUN=function(x,Current.Proposal, Parameters) do.call(generate.points.dirs, c(x=list(x), Parameters, Current.Proposal=list(Current.Proposal))), Current.Proposal, Parameters)
		}
		else {
		  cat(" initial diversity is ", nSeq, "\n")
		  #New.Points<-clusterApplyLB(mycl, Last.State.List,  pf.par.internal, Current.Proposal)
		  New.Points<-parallel:::clusterApply(mycl, Last.State.List,  pf.par.internal, Current.Proposal)
		}
		#=======================================================
		New.Particles<-unlist(New.Points)[sort.list(Order.vector, na.last = NA, method =  "quick")]
		return(New.Particles)
	}

  
  ## for ver 1.6 
  Prev.Weights<-rep(1, nParticles) 
  ResampleCount<-0
for (Time.Period in 1:nrow(in.Data$Main.Index)) {
    cat("\n\n##########################\n     Time.Period", Time.Period, "\n")
    #cat("prep. data:")
    Current.Proposal<-in.Data$Matrix.Index.Table[in.Data$Main.Index$Biol.Prev[Time.Period],]
    #=======================================
	cat("generating new particles")
	New.Particles<-propagate.particles(Last.Particles=Results.stack[,ncol(Results.stack)], Current.Proposal=Current.Proposal, parallel=parallel, Parameters=Parameters, mycl=mycl)

    #=====================================================
    # resampling step 
    # we need to estimate weights of resulting points and then resample them proportionally to weights..
    Current.Phys.Weights<-in.Data$Phys.Mat[New.Particles,Time.Period]

    if (!is.na(k) & Time.Period >1) {
	get.directional.weights<-function(in.Data, Last.Last.Particles, Last.Particles, New.Particles, k, parallel, mycl, D.kappa) {
      #===================================
      Prev.Dirs<-in.Data$Azimuths[cbind(Last.Last.Particles, Last.Particles)]/180*pi
      New.Dirs<-in.Data$Azimuths[cbind(Last.Particles,New.Particles)]/180*pi
      FromTo=cbind(Prev.Dirs, New.Dirs)
      if (parallel) {Angles.probs<-parallel:::parApply(mycl, FromTo, 1, my.dvonmises, mykap=k)
      } else {
        Angles.probs<-apply(FromTo,1, FUN=function(x, k) as.numeric(suppressWarnings(circular:::dvonmises(x[[2]], mu=x[[1]], kappa=k))), k=k)
      }
      Angles.probs[is.na(Angles.probs)]<-D.kappa
      #cat("\n min.Kappa:", min(Angles.probs), ", max.Kappa", max(Angles.probs), "\n")
	  return(Angles.probs)
     }
 
	 Angles.probs<-get.directional.weights(in.Data, Results.stack[,ncol(Results.stack)-1], Results.stack[,ncol(Results.stack)], New.Particles, k, parallel, mycl, D.kappa)
      Current.Weights<- Current.Phys.Weights*Angles.probs
    }
    else {Current.Weights<-Current.Phys.Weights*D.kappa}
 	#=====================================================
	# here is the change from 1.6 - Adding cumulative weights now..
	Current.Weights.with.Prev.mat<-cbind(Weights.stack, Current.Weights)
	rowProds <- function(a) exp(rowMeans(log(a)))
	Current.Weights.with.Prev<-pmax(rowProds(Current.Weights.with.Prev.mat), 1e-323)	
	#print(head(Current.Weights.with.Prev))
	# if (any(Current.Weights.with.Prev==1e-323)) {
	# cat("precise\n")
	#rowProds <- function(a) exp(rowSums(log(a)))
	# Current.Weights.with.Prev<-pmax((Current.Weights*rowProds.precise(Weights.stack))^Degree, 1e-323)
	# }
	# I use ncol Weights stack to make surface a bit more flat..
	#================================================================
	# ver 1.7. 
	# ADAPTIVE RESAMPLING
	if (adaptive.resampling) {
	ESS<-(sum(Current.Weights.with.Prev)^2)/sum((Current.Weights.with.Prev)^2)
	cat("ESS is ", ESS)
	} else {
	ESS<-1
	}
	if (ESS<(nParticles/2) & Time.Period>1) {	
    Rows<-suppressWarnings(sample.int(nParticles, replace = TRUE, prob = Current.Weights.with.Prev))
	ResampleCount<-ResampleCount+1
	cat(" - resampling", ResampleCount, "\n")
    Results.stack<-Results.stack[Rows,]
	} else {
	Rows<-1:nParticles # no resampling
	cat("\n")
	}
	#=====================================================
	#               Smoothing
	#====================================================

#======================================================
# start OF SMOOTHING POINTS ESTIMATION

#  ==================================
# ver. 1.8 - now use one function for all particle propagations:
	cat("  smoothing new particles")

	New.Particles.MH<-propagate.particles(Last.Particles= Results.stack[,ncol(Results.stack)], Current.Proposal=Current.Proposal, parallel=parallel, Parameters=Parameters, mycl=mycl)

	Current.Phys.Weights.MH<-in.Data$Phys.Mat[New.Particles.MH,Time.Period]
	
    if (!is.na(k) & Time.Period > 1) {
	  Angles.probs<-get.directional.weights(in.Data, Results.stack[,ncol(Results.stack)-1], Results.stack[,ncol(Results.stack)], New.Particles.MH, k, parallel, mycl, D.kappa)
      Current.Weights.MH<- Current.Phys.Weights.MH*Angles.probs
    }
    else {Current.Weights.MH<-Current.Phys.Weights.MH*D.kappa}

#   end of smoothing estimation    
#============================================
#============================================

	# now I need to compare resampled results with new resuklts
	#according to MH
	v<-runif(nParticles)
	Accepted.Particles<-New.Particles.MH
	Ratio<-Current.Weights.MH/(Current.Weights[Rows])
	Ratio[is.na(Ratio)]<-1 # NA is caused by 0/0
	Eq<-!as.logical(floor(pmin(1, Ratio)-v)+1)
	Accepted.Particles[Eq]<-New.Particles[Rows][Eq]
	New.weights<-Current.Weights.MH
	New.weights[Eq]<-Current.Weights[Rows][Eq]
	cat("smoothed ", nParticles-sum(Eq), "Particles \n")

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

	if (Time.Period==2) Particles.at.2<-Accepted.Particles
#	if (Time.Period>2) {
#		Particles.at.2<-Particles.at.2[Rows]
		cat("   unique P in point leaving stack", length(unique(Results.stack[,1])), "\n")
#	}

	if (Time.Period<=L) cat("creating stack\n")
	if (Time.Period>L) {
		# save points
		if (is.null(All.results)) {
		All.results<-paste(Results.stack[,1], sep=".")
		} else {
		All.results<-paste(All.results, Results.stack[,1], sep=".")
		}
		# save transitions
		Trans[[Time.Period-L]]<-get.transition.rle(Results.stack[,1], Results.stack[,2])
		# clean Results.stack
		Results.stack<-Results.stack[,-1]
		# clean Weights.stack
		Weights.stack<-Weights.stack[,-1]
	}
 }
  #####################################
  #save(All.results, file="All.results.usmoothed.RData")
  # now we need to add final point!
  
  if (known.last) {
    Results.stack<-pf.final.smoothing(in.Data, Results.stack, precision.sd=precision.sd, nParticles=nParticles, save.memory=F, last.particles=Results.stack[,ncol(Results.stack)])
  }
  
  # and here we need to add a thing that will finish All.results and Trans from the points that are still in the stack
  Length<-ncol(Results.stack)
  cat("adding last points form the stack to the resutls\n")
  for (rest in 1:Length) {
		# save points
		All.results<-paste(All.results, Results.stack[,rest], sep=".")
		if (rest<Length) {
		# save transitions
		Trans[[Time.Period-L+rest]]<-get.transition.rle(Results.stack[,rest], Results.stack[,rest+1])
		}
   }
  #save(All.results, file="All.results.smoothed.RData")
  if (parallel)   parallel:::clusterEvalQ(mycl, rm(Parameters)) 
  if (length(existing.cluster)==1) parallel:::stopCluster(cl = mycl)
  if (sink2file) sink()
  return(list(All.results=All.results, Trans=Trans))
}
