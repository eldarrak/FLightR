# ver 3.7 - change to clever mask to filter already generated points or should we filter before the point generation???
# also added back the angle filters..
# ver 3.6 distance filter for AC.distance2*1.3 seems to be enough
# ver. 3.5 Feb 04 2014 added save.transitions option and check.outliers option..

# ver 3.4 decided to try directional outlier detection..

# ver 3.3. (22-10-2013) the idea is to add some stability to the outliers as it looks like better exclude such twilights...
# this will involve following parts:
# we will estimate distances from a to b and from a to c (that's very simple)
# and we will assume that the bird should not come back..
# it means that if AB is bigger than AC in more then say 95% of particles then we should skip stage B. How should we do it?
# the only problem is that at t we should resample only stage (t-1) but not t as we were doing.. 
# another problem is that we should estimate a weighted mean distance.
#
# 
#======================
# decided to add user defined degree of an adaptive resmapling..
# so now if it is 1 then we will resample every time if not then we will resample accordeing to ESS<adaptive.resampling
# ver.3.1 decided also to estimate weights as a sum of all the previous weights...
# I have had it somewhere before already..
# ver. 3.0 decided to change mask to allow bird to fly over but not allow it to stop where it is unlikely to do so..
#
# ver.2.1 going back for simple estimation of distances during update..
# ver 2.0 (08-02-2013)  added angle drift estimation...
# ver 1.9 (19-01-2013) decicided to go for a major revision of the whole thing
#   the idea is to make it block resampling.
# so no it is PF with MH smoothing, adaptive resamling and block sampling..

# ver 1.8 wrapped some stuff in functions
# ver 1.7 ESS added for adaptive resampling..

## ver 1.6 BUG found!!
## we needed to estimate weights accounting for all previous particles!!!!!
##

## ver 1.5 decided to add a new function into the M stgep of EM
## the idea is that it will be able to deal with growing distances..

## ver 1.4 decide to make sort - unsort faster!

## 1.3 MPI cluster support added through parameter cluster.type..
##     compilation added

## 1.2 faster version of get.matrix from char

## ver 1.1 addition of save.Res=F option when only first and last results will be saved. + patch not to plot when repeat..
## ver 1.0 MH smoothing
## ver 0.9 fixed movement
## ver 0.8 - likelihood evaluations and optimization added.
## ver 0.7 save.ram function added
## ver 0.6.. for the new input format..

# SO resample! ver 0.5
# I've got a feeling that way of making simulations is not good...
# we don't need to simulate by angle but !!! 
# WE NEED to check these angles...
# but first I would check on how results work!

# ver 0.4 Mat 16 - the prev version did not work as it was too slow... so the idea is that we need to create one number that will be ID for the transition - not for the point..

# ver 0.3 May 14 - attempt to incorporate Second order MC..
#
# ver 0.2 May 14 - the idea is to convert all result into character vectors as it will save nenory
# in this case we will need to pass All.Res together with last particles...
#
# version 0.1 May 10 2012
# stage 0!
# add geogr.proposal..
require(compiler)
enableJIT(3)

solar<-function (day) {
  tm <- as.POSIXlt(day, tz = "GMT")
  hh <- tm$hour
  mm <- tm$min
  ss <- tm$sec
  julday<-function (tm) {
    yr <- 1900 + tm$year
    mon <- tm$mon + 1
    day <- tm$mday
    yr[mon <= 2] <- yr[mon <= 2] - 1
    mon[mon <= 2] <- mon[mon <= 2] + 12
    a <- floor(yr/100)
    b <- 2 - a + floor(a/4)
    floor(365.25 * (yr + 4716) + floor(30.6001 * (mon + 1))) + 
      day + b - 1524.5
  }
  jday <- julday(tm) + (hh + (mm + ss/60)/60)/24
  t <- (jday - 2451545)/36525
  M <- 357.52911 + t * (35999.05029 - 0.0001537 * t)
  eqcent <- sin(pi/180 * M) * (1.914602 - t * (0.004817 + 1.4e-05 * 
    t)) + sin(pi/180 * 2 * M) * (0.019993 - 0.000101 * t) + 
    sin(pi/180 * 3 * M) * 0.000289
  L0 <- 280.46646 + t * (36000.76983 + 0.0003032 * t)
  L0 <- L0%%360
  lambda0 <- L0 + eqcent
  omega <- 125.04 - 1934.136 * t
  lambda <- lambda0 - 0.00569 - 0.00478 * sin(pi/180 * omega)
  seconds <- 21.448 - t * (46.815 + t * (0.00059 - t * (0.001813)))
  obliq0 <- 23 + (26 + (seconds/60))/60
  omega <- 125.04 - 1934.136 * t
  obliq <- obliq0 + 0.00256 * cos(pi/180 * omega)
  e <- 0.016708634 - t * (4.2037e-05 + 1.267e-07 * t)
  y <- tan(pi/180 * obliq/2)^2
  eqtime <- 180/pi * 4 * (y * sin(pi/180 * 2 * L0) - 2 * e * 
    sin(pi/180 * M) + 4 * e * y * sin(pi/180 * M) * cos(pi/180 * 
    2 * L0) - 0.5 * y^2 * sin(pi/180 * 4 * L0) - 1.25 * e^2 * 
    sin(pi/180 * 2 * M))
  solarDec <- asin(sin(pi/180 * obliq) * sin(pi/180 * lambda))
  sinSolarDec <- sin(solarDec)
  cosSolarDec <- cos(solarDec)
  solarTime <- (hh * 60 + mm + ss/60 + eqtime)/4
  list(solarTime = solarTime, sinSolarDec = sinSolarDec, cosSolarDec = cosSolarDec)
}
solar<-cmpfun(solar)

elevation<-function (lon, lat, sun) {
  hourAngle <- sun$solarTime + lon - 180
  cosZenith <- (sin(pi/180 * lat) * sun$sinSolarDec + cos(pi/180 * 
    lat) * sun$cosSolarDec * cos(pi/180 * hourAngle))
  cosZenith[cosZenith > 1] <- 1
  cosZenith[cosZenith < -1] <- -1
  90 - 180/pi * acos(cosZenith)
}
elevation<-cmpfun(elevation)

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
  
  if (x[[3]]>0) {
	# here is the addition of clever mask
	#if (in.Data$Points.Land[x[[1]],3]==0) x[[3]]=x[[2]]
	# end of addition
    Dists.distr<-in.Data$distance[x[[1]],]
    Dists.probs<-truncnorm:::dtruncnorm(Dists.distr, a=a, b=b, Current.Proposal$M.mean, Current.Proposal$M.sd)
    ###
    #  library(fields)
    # dists
    #image.plot(as.image(Dists.probs, x=in.Data$Points.Land, nrow=50, ncol=50))
    #
    if (Current.Proposal$Kappa>0) {
      Angles.dir<-in.Data$Azimuths[x[[1]],]
      Angles.probs<-as.numeric(suppressWarnings(circular:::dvonmises(Angles.dir/180*pi, mu=Current.Proposal$Direction/180*pi, kappa=Current.Proposal$Kappa)))
      Angles.probs[is.na(Angles.probs)]<-0
    }
    else {
      Angles.probs<-0.1591549
    }
    # this is neede to catch an error 
    if (which.max(c(min(c(Angles.probs, Dists.probs)), 0))==2) {
      cat("PF produced weird probs \n")
      tmp.out<-list(Angles.probs=Angles.probs, Dists.probs=Dists.probs, Current.Proposal=Current.Proposal, Data=x)
      save(tmp.out, file=paste("tmp.out", round(runif(n=1,min=1, max=1000)), "RData", sep="."))
      Biol.proposal<-pmax.int(Dists.probs*Angles.probs, 0)
    } else {
      Biol.proposal<-Dists.probs*Angles.probs
    }
    pos.biol<-suppressWarnings(sample.int(length(Biol.proposal), size=x[[3]], replace=T, prob=Biol.proposal))
    # ok, now we want to return more complicated stuff - final indexes!
    return(resample(as.integer(c(pos.biol, rep(x[[1]],(x[[2]]-x[[3]]))))))
  }
  else {
    return(as.integer(rep(x[[1]],(x[[2]]))))
  }
}
generate.points.dirs<-cmpfun(generate.points.dirs)

pf.final.smoothing<-function(in.Data, All.results, precision.sd=25, nParticles=1e6, save.memory=F, last.particles=NA) {
  # this function simply resamples final points proportionally to the distance to known finish.
  Final.point.real<-in.Data$stop.point
  # now we want to get distances.. I'll not index it as we will do this only once..
  Final.points.modeled=last.particles
  Weights<-dnorm(in.Data$distance[Final.points.modeled, Final.point.real], mean=0, sd=precision.sd)
  Rows<- suppressWarnings(sample.int(nParticles, replace = TRUE, prob = Weights/sum(Weights)))
  if (save.memory) {
    return(All.results[Rows])
  } else {		
    return(All.results[Rows,])
  }
}

pf.final.smoothing<-cmpfun(pf.final.smoothing)

#Res<-geologger.pf(all.out)
# now wen need to plot results..

plot.optimisation.results<-function(object, all.arrays.object, add=F, col="red", map.fill=T, type="mean", pch=3) {
  # simple plotting function
  if (as.character(type)=="mean") {
    Track.coords<-t(apply(object,2, FUN=function(x) cbind(mean( all.arrays.object$Points.Land[x,1]), mean(all.arrays.object$Points.Land[x,2]))))}
  else {
    cat("plotting medians\n")
    Track.coords<-t(apply(object,2, FUN=function(x) cbind(median( all.arrays.object$Points.Land[x,1]), median(all.arrays.object$Points.Land[x,2]))))
  }
  if (add) {lines(Track.coords, type="l", col=col)} 
  else {
    plot(Track.coords, type="l", col=col)
    data("wrld_simpl", package="maptools")
    
    plot(wrld_simpl, add=T, col=ifelse(map.fill, "lightgray" ,0))
  }
  points(Track.coords, pch=pch, col=col)
  lines(Track.coords, type="l", col=col)
}


return.matrix.from.char<-function(Res.txt) {
  #this function is needed to get matrix from character vector
  return(t(sapply(strsplit(Res.txt, "\\."), as.integer)))
}
return.matrix.from.char<-cmpfun(return.matrix.from.char)


# new type of plot.optim.results..
plot.final.track<-function(all.out, add=F, col="red", map.fill=T, pch=3) {  
  # add mean line:
  if (add) {
    lines(all.out$Final.Means$CENTRE.y~all.out$Final.Means$CENTRE.x, pch=3, col=col, lwd=2) } else {
      data("wrld_simpl", package="maptools")
      plot(all.out$Final.Means$CENTRE.y~all.out$Final.Means$CENTRE.x, pch=3, col=col, lwd=2, type="l")
      plot(wrld_simpl, add=T, col=ifelse(map.fill, "lightgray" ,0))
      lines(all.out$Final.Means$CENTRE.y~all.out$Final.Means$CENTRE.x, pch=3, col=col, lwd=2)
    }
  points(all.out$Final.Means$CENTRE.y~all.out$Final.Means$CENTRE.x, pch=pch, col=col)
}

# plot.final.track(all.out, add=F)
# plot.final.track(all.out, add=T, col="blue")
# plot.final.track(all.out, add=T, col="orange")


pf.par.internal<-function(x, Current.Proposal) {
  # this simple function is needed to save I/0 during parallel run. Being executed at a slave it creates inoput for the next function taking Current proposal from the slave and Parameters from the master
  new.Parameters<-c(x=list(x), get("Parameters"), Current.Proposal=list(Current.Proposal))
  Res<-do.call(generate.points.dirs, new.Parameters)
  return(Res)
}
pf.par.internal<-cmpfun(pf.par.internal)


my.dvonmises=function(x, mykap) {
  return(as.numeric(suppressWarnings(circular:::dvonmises(x[[2]], mu=x[[1]], kappa=mykap))))}

my.dvonmises<-cmpfun(my.dvonmises)




# in this version I'll try to make the function fast again..
# so I'll use only local weights
# so the things to eclude are
# 1. weights stack
# 2. MH step

# 1. so I want to use weights from just the last point.
# 

pf.run.parallel.SO.resample<-function(in.Data, cpus=2, nParticles=1e6, known.last=T, precision.sd=25, sea.value=0.01, k=1, parallel=T, plot=T, existing.cluster=NA, cluster.type="SOCK", a=45, b=500, sink2file=F, L=25, adaptive.resampling=0.5, RStudio=F, check.outliers=F) {
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
  
  
  if (any(in.Data$Points.Land[,3]!=1)) { 
	smart.filter=TRUE
	cat("smart filter is ON\n")
	} else {
	smart.filter=FALSE
	cat("smart filter is OFF\n")
	}
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
  
  Weights.stack<-as.matrix(rep(1/nParticles, nParticles))
  
  New.weights<-rep(1/nParticles, nParticles)
  All.results<-NULL
  in.Data$outliers<-c()
	  Trans<-vector(mode = "list", length = nrow(in.Data$Main.Index))
  	if (check.outliers) {
	  in.Data$AB.distance<-c()
	  #in.Data$AC.distance<-c()
	  in.Data$AC.distance2<-c()
	  in.Data$Dif.ang-c()
	  #in.Data$BC.mean<-c()
  }
  
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
      cat(" initial diversity is ", nSeq)
      #New.Points<-clusterApplyLB(mycl, Last.State.List,  pf.par.internal, Current.Proposal)
      New.Points<-parallel:::clusterApply(mycl, Last.State.List,  pf.par.internal, Current.Proposal)
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

  for (Time.Period in 1:nrow(in.Data$Main.Index)) {
    steps.from.last=steps.from.last+1

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

	if (smart.filter) {
	#=========================
	# here we should introduce clever mask..
	# this will be something like that
	Index.Stable<-which(New.Particles == Results.stack[,ncol(Results.stack)])
	Current.Weights[Index.Stable]<-Current.Weights[Index.Stable]*in.Data$Points.Land[New.Particles,3]
	}
	#=======================================================
	# now I want to add 3.3
	# the idea is following - we should compare distances from the last t-2 to t-1 and to t
	if (Time.Period>2 & check.outliers) {
	# here I want to try directiondl outliers..
	# so we need to pick up particles that moved..
	# 
	#=================================================================
	AB.distance<-weighted.mean(
	in.Data$distance[Results.stack[,(ncol(Results.stack)-1):ncol(Results.stack)]], Weights.stack[,ncol(Weights.stack)])
	#AC.distance<-	weighted.mean(in.Data$distance[cbind(Results.stack[,ncol(Results.stack)-1], New.Particles)], Current.Weights)
	AC.distance2<-	weighted.mean(in.Data$distance[cbind(Results.stack[,ncol(Results.stack)-1], New.Particles)], Weights.stack[,ncol(Weights.stack)]*Current.Weights)
	cat("AB.distance:", round(AB.distance, 2), "\n")
	#cat("AC.distance:", round(AC.distance, 2), "\n")
	cat("AC.distance2:", round(AC.distance2, 2), "\n")

	Dif.ang=180
	if (AB.distance>50) { # go for angles only if disctances are high!
	resample <- function(x, ...) x[sample.int(length(x), ...)]
	# the AB ones wil have folowing..
	BA.dir<-in.Data$Azimuths[Results.stack[,ncol(Results.stack):(ncol(Results.stack)-1)]]
	BA.moved<-which(!is.na(BA.dir))
	BA.mean<-mean.circular(circular(resample(BA.dir[BA.moved], replace=T,prob=Weights.stack[,ncol(Weights.stack)][BA.moved]), units="degrees"), na.rm=T)
	#cat(BA.mean)
	BC.dir<-in.Data$Azimuths[cbind(Results.stack[,ncol(Results.stack)-1], New.Particles)]
	BC.moved<-which(!is.na(BC.dir))
	BC.mean<-mean.circular(circular(resample(BC.dir[BC.moved], replace=T, prob=(Weights.stack[,ncol(Weights.stack)]*Current.Weights)[BC.moved]), units="degrees"), na.rm=T)
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
    if ((ESS<(nParticles*adaptive.resampling) & Time.Period>1) | Time.Period == nrow(in.Data$Main.Index)) {	
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
    #	Current.Phys.Weights.MH<-in.Data$Phys.Mat[New.Particles.MH,Time.Period]
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
      Weights.stack<-as.matrix(Weights.stack[,-1])
cat("******************\n")
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
  tmp<-list(AB.distance=in.Data$AB.distance, AC.distance2=in.Data$AC.distance2, Dif.ang=in.Data$Dif.ang)

  return(list(All.results=All.results, Trans=Trans, outliers=in.Data$outliers, tmp=tmp))
}

pf.run.parallel.SO.resample<-cmpfun(pf.run.parallel.SO.resample)


get.transition.rle=function(From, To) {
  rle(sort.int(From*1e5+To, method = "quick"))
}

get.transition.rle<-cmpfun(get.transition.rle)

# 
# mu.sigma.truncnorm<-function(points.current, dists, dists.real, azimuths, phys.proposal, Points.Land, point_ID, a=45, b=500, Current.mean, Current.sd, nonZero.Points) {
  # # this function is needed for the EM part for proposal update
  # points.current.nonZero<-points.current[nonZero.Points]
  # dists.current<-dists[points.current.nonZero]
  # dists.real.nz<-dists.real[nonZero.Points]
  
  # if (length(unique(dists.current))>1) {
    # azimuths.current<-azimuths[points.current.nonZero]
    # azimuths.current<-azimuths.current[!is.na(azimuths.current)]
    # Mean.Direction<- CircStats:::circ.mean(azimuths.current*pi/180)*180/pi
    # Angular.diff<-Mean.Direction-azimuths.current
    # # now I want to project point on the line..
    # #Projected.shift.X<-(dists.current*cos(Angular.diff*pi/180))
    # Second.Point<- destPoint(Points.Land[point_ID,],  b=Mean.Direction, d=1e6)
    # # ok now we have to points
    # # let's generate set of points for the lines
    # # so starting points will be at the border of our map
    # X_points_1<-seq(min(Points.Land[,1]), max(Points.Land[,1]), 0.25)
    # # Now Y
    # Y_points_1<-mean(Points.Land[,2])
    # # now I want to create second points for this line that will be perpendicular to our line..
    # Points_2<-destPoint(cbind( X_points_1,  Y_points_1),  b=Mean.Direction+90, d=1e6)
    # # ok, now I want to get intersections
    # Intersections<-gcIntersect(Points.Land[point_ID,], Second.Point, cbind( X_points_1,  Y_points_1), Points_2)
    # Intersections<-rbind(Intersections[,1:2], Intersections[,3:4])
    # # now I need to get points within our map:
    # Intersections<-Intersections[Intersections[,1]>min(Points.Land[,1]) & Intersections[,1]<max(Points.Land[,1]) & Intersections[,2]>min(Points.Land[,2])& Intersections[,2]<max(Points.Land[,2]),]
    # # now we want to go for distances
    # Closest.Points<-apply(Intersections, 1, FUN=function(x) which.min(spDists(Points.Land, cbind(x[[1]], x[[2]]), longlat=TRUE)))
    # Best.Point<- Closest.Points[which.max(phys.proposal[	Closest.Points])]
    # # and now I want to get simple distance from the current point to all other points and then to generate weights according to this probs..
    # Dists2best<-spDists(Points.Land, cbind(Points.Land[Best.Point,1],Points.Land[Best.Point,2]), longlat=TRUE)
    # # and to get the distribution we nee SD
    # SD.Dists2best<-sd(Dists2best[points.current.nonZero,1])
    # Weights<-dnorm(Dists2best[points.current.nonZero,1], mean=0, sd=SD.Dists2best)
    # tr.norm<-function(prm, dists.real.nz, Weights,a, b) {
      # return(sum(-log(truncnorm:::dtruncnorm(
        # dists.real.nz, a=a,b=b,mean=prm[1],sd=max(0.1,prm[2])))*(Weights^2)))
    # }
    # tmp.mean<-seq(a,b,by=50)
    # tmp.var <-seq(50,300,by=50)
    # tmp.prm <-expand.grid(tmp.mean,tmp.var)
    # tmp.LL  <-apply(tmp.prm, 1, FUN=tr.norm, dists.real.nz=dists.real.nz, a=a, b=b, Weights=Weights)
    # tmp.pos <-which.min(tmp.LL)
    # tmp.init<-tmp.prm[tmp.pos,]
    # Res=try(optim(tmp.init, tr.norm, gr=NULL, dists.real.nz, Weights, a, b, method="Nelder-Mead",control=list(trace=0)))
    
    # if (class(Res)=="try-error") Res=try(optim(tmp.init, tr.norm, gr=NULL, dists.real.nz, Weights, a, b, method="Nelder-Mead",control=list(trace=0)))
    # if (class(Res)=="try-error") save(dists.current, Res, file="x.RData")
    # return(c(Res$par[1], max(0.1,Res$par[2])))
  # } else {
    # return(c(mean(dists.real.nz), sd(dists.real.nz)))
  # }
# }

mu.sigma.truncnorm<-function(x, a=45, b=500) {
  if (length(unique(x))>1) {
    tr.norm<-function(prm) {
      sum(-log(truncnorm:::dtruncnorm(x,a=a,b=b,mean=prm[1],sd=prm[2])))
    }
    Res=try(optim(c(mean(x),sd(x)), tr.norm, method="BFGS"))
    #if (class(Res)=="try-error") Res=try(optim(c(mean(x),sd(x)), tr.norm,method="SANN"))
    if (class(Res)=="try-error") {
      save(x, Res, file="x.RData")
      Res$par<-c(mean(x), sd(x))
    }
    return(c(Res$par[1], Res$par[2]))
  } else {
    return(c(mean(x), sd(x)))
  }
}
mu.sigma.truncnorm<-cmpfun(mu.sigma.truncnorm)


update.proposal.PF<-function(output.matrix, Trans, in.Data, save.transitions=F, fixed.parameters=NA, a=45, b=500, parallel=F, existing.cluster=NULL, estimatetruncnorm=F) {
  mycl=existing.cluster
  # main function for proposal update
  # from the version 1.5 this function will do estimation of distance and it's SD..

  
  # fixed.parameters could be list(Kappa=1, M.sd=250)
  require("truncnorm")
  require("circular")
  require("CircStats")
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
    in.Data$distance[x%/%1e5, x%%1e5]
  }
  for (i in 1:length(Trans)) {
    Distances[[i]]$values<-sapply(Trans[[i]]$values, FUN=function(x) dist.fun(x))
  }
  cat("   estimating directions\n")
  #ok, now we want to get directions
  
  Directions<-Trans
  direct.fun<-function(x) {
    in.Data$Azimuths[x%/%1e5, x%%1e5]
  }
  for (i in 1:length(Trans)) {
    Directions[[i]]$values<-sapply(Trans[[i]]$values, FUN=function(x) direct.fun(x))
  }
  cat("   estimating mean directions\n")
  Mean.Directions<-unlist(lapply(Directions, FUN=function(x) CircStats:::circ.mean(inverse.rle(list(lengths=x$lengths[!is.na(x$values)], values=x$values[!is.na(x$values)]))*pi/180)*180/pi))
  #plot(Mean.Directions)
  cat("   estimating kappas\n")
  Kappas<-unlist(lapply(Directions, FUN=function(x) CircStats:::est.kappa(inverse.rle(list(lengths=x$lengths[!is.na(x$values)], values=x$values[!is.na(x$values)]))*pi/180)))
  
  
  # now we want to get mean ditance
  cat("   estimating mean dists\n")
  #
  #cat("   estimating mean and SD to report dists SD\n")
  # attempt to go just for simple distance estimation...
  Mean2report<-unlist(lapply(Distances, FUN=function(x) mean(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))
  SD2report<-unlist(lapply(Distances, FUN=function(x) sd(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))
  cat("   estimating probs of migration\n")
  # ok, now we want to get parameters for distances
  #
  Probability.of.migration<-unlist(lapply(Distances, FUN=function(x) sum(x$lengths[x$values!=0])))/(dim(output.matrix)[1])
  
  ## 
    if (estimatetruncnorm) {
  Mean.and.Sigma<-lapply(Distances, FUN=function(x) mu.sigma.truncnorm(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])), a=a, b=b))
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
  Median.Dists<-unlist(lapply(Distances, FUN=function(x) median(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))

  #unlist(lapply(Distances, FUN=function(x) sd(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))
  #plot(Kappas)
  
  
  
  
  cat("   creating output")
  ################
  # create new proposal from the posteriors
  #New.Matrix.Index.Table<-data.frame(Direction=Mean.Directions, M.mean=Mean.Dists, M.sd=Mean.SD, Decision=Probability.of.migration, Kappa=Kappas)
  ################
  
  all.arrays.object<-in.Data
  
  all.arrays.object$Matrix.Index.Table$Decision<-Probability.of.migration
  all.arrays.object$Matrix.Index.Table$Direction<-Mean.Directions
  all.arrays.object$Matrix.Index.Table$Kappa<-Kappas
  all.arrays.object$Matrix.Index.Table$M.mean<-Mean.Dists
  all.arrays.object$Matrix.Index.Table$M.sd<-Mean.SD
  all.arrays.object$Matrix.Index.Table$M.medians<-Median.Dists
  all.arrays.object$Matrix.Index.Table$Mean2report<-Mean2report
  all.arrays.object$Matrix.Index.Table$SD2report<-SD2report
  
  #all.arrays.object$Matrix.Index.Table<-New.Matrix.Index.Table
  
  all.arrays.object$Main.Index$Biol.Prev=1:(nrow(all.arrays.object$Matrix.Index.Table))
  if(save.transitions) all.arrays.object$Transitions.rle<-Trans
  
  
  if (is.list(fixed.parameters)) {
    # this part can be moved upper not to estimate parametrs in case they are fixed.
    if ("M.sd" %in% names(fixed.parameters)) all.arrays.object$Matrix.Index.Table$M.sd<-fixed.parameters$M.sd
    
    if ("M.mean" %in% names(fixed.parameters)) all.arrays.object$Matrix.Index.Table$M.mean<-fixed.parameters$M.mean
    
    if ("Kappa" %in% names(fixed.parameters)) all.arrays.object$Matrix.Index.Table$Kappa<-fixed.parameters$Kappa
    
    if ("Direction" %in% names(fixed.parameters)) all.arrays.object$Matrix.Index.Table$Kappa<-fixed.parameters$Direction
  }
  #all.arrays.object$Main.Index$Biol.Next=2:(nrow(all.arrays.object$Matrix.Index.Table))
  cat(" DONE!\n")
  return(all.arrays.object)
}

update.proposal.PF<-cmpfun(update.proposal.PF)

require(sp)
require(rgeos)
require(rgdal)

coords.aeqd.jitter <- function(coords, r, n)
{
# coords is a vector of leghth 2 c(lon, lat)
# r should be in meters
# made on te basis of this:
#"http://gis.stackexchange.com/questions/121489/1km-circles-around-lat-long-points-in-many-places-in-world"
 stopifnot(length(coords) == 2)

	p = SpatialPoints(matrix(coords, ncol=2), proj4string=CRS("+proj=longlat +datum=WGS84"))
    aeqd <- sprintf("+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                    p@coords[[2]], p@coords[[1]])
    projected <- spTransform(p, CRS(aeqd))
    buffered <- gBuffer(projected, width=r, byid=TRUE)
    lambert <- sprintf("+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                    p@coords[[2]], p@coords[[1]])
	buffered_eqarea <- spTransform(buffered, CRS(lambert))
	random_points<-spsample(buffered_eqarea,n=n,type="random")
	if (is.null(random_points)) random_points<-spsample(buffered_eqarea,n=n,type="random", iter=40) 
	if (is.null(random_points)) random_points<-p
    spTransform(random_points, p@proj4string)
}



get.coordinates.PF<-function(output.matrix, in.Data, save.points.distribution=F) {
  library("aspace")
  # this function will extract point coordinates from the output matrix.. 
  # the question is do we need only mean and sd or also median and quantiles?
  # I will start from mean and SD
  #plot(c(min(in.Data$Points.Land[,1]),max(in.Data$Points.Land[,1])), c(min(in.Data$Points.Land[,2]),max(in.Data$Points.Land[,2])), type="n")
  
  if (save.points.distribution) {
    require("bit")
    Points<-vector(mode = "list", length = (dim(output.matrix)[2]))
    cat("   extracting indexes for rle\n")
    for (i in 1:(dim(output.matrix)[2])) {
      Points[[i]]<-Rle<-bit:::intrle(sort.int(output.matrix[,i], method="quick"))
      if (is.null(Rle)) Points[[i]]<-rle(sort.int(output.matrix[,i], method="quick"))
    }
    in.Data$Points.rle<-Points
  }
  Means=aspace:::calc_box(id=1,  points=in.Data$Points.Land[output.matrix[,1],1:2])
  for (i in 2:(dim(output.matrix)[2])) {
    Means[i,]=aspace:::calc_box(id=i,  points=in.Data$Points.Land[output.matrix[,i], 1:2])
    #plot_box(plotnew=F, plotpoints=F)
  }
  in.Data$Final.Means<-Means
  
  #############
  # new part for medians
  cat("estimating quantiles for positions\n")
	Quantiles<-c()
	for (i in 1:length(Points)) {
	Quantiles<-rbind(Quantiles, c(summary(all.out$Points.Land[inverse.rle(Points[[i]]),2]), Mode=all.out$Points.Land[Points[[i]]$values[which.max(Points[[i]]$lengths)],2], summary(all.out$Points.Land[inverse.rle(Points[[i]]),1]), Mode=all.out$Points.Land[Points[[i]]$values[which.max(Points[[i]]$lengths)],1]))
	}
	Quantiles<-as.data.frame(Quantiles)
	names(Quantiles)[1:6]<-paste(names(Quantiles)[1:6], "lat", sep="")
	names(Quantiles)[8:13]<-paste(names(Quantiles)[8:13], "lon", sep="")

	###########
	# doing jitter first
	cat("adding jitter to medians\n")
	JitRadius<-min(all.out$distance[all.out$distance>0])/2*1000 # in meters
	#now I want to generate random poitns in the radius of this
	coords=cbind(Quantiles$Medianlon, Quantiles$Medianlat)
	tmp<-apply(coords, 1, coords.aeqd.jitter, r=JitRadius, n=1 )
	jitter_coords<-t(sapply(tmp, coordinates))

	Quantiles$MedianlonJ<-jitter_coords[,1]
	Quantiles$MedianlatJ<-jitter_coords[,2]

	names(Quantiles)<-gsub("\\s","", names(Quantiles))
	names(Quantiles)<-gsub("1","F", names(Quantiles))
	names(Quantiles)<-gsub("3","T", names(Quantiles))
	  
	in.Data$Quantiles<-Quantiles
  
  
  return(in.Data)
}
get.coordinates.PF<-cmpfun(get.coordinates.PF)

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

get.angle.drift<-cmpfun(get.angle.drift)
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

sun.matrix.internal<-cmpfun(sun.matrix.internal)
#Z<-sun.matrix.internal(in.Data$Matrix.Index.Table[1,], in.Data$Points.Land)
node.run<-function(x=x) {
  Points.Land<-get("Points.Land", globalenv())
  return(sun.matrix.internal(x, Points.Land))
}

node.run<-cmpfun(node.run)


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

create.spatial.sun.matrices<-cmpfun(create.spatial.sun.matrices)


# the next function will be one to estimate LL

get.LL.PF<-function(in.Data, All.results.mat) {
  # needed to estimate log likelihood of the optimization
  L=0
  for (i in 1:(dim(All.results.mat)[2]-1)) {
    L=L+log(mean(in.Data$Phys.Mat[All.results.mat[,i],i]))
  }
  #L=L/(dim(All.results.mat)[2]-1)
  #LL=-sum(log(L))
  LL=-L
  return(LL)
}
get.LL.PF<-cmpfun(get.LL.PF)

get.optimization.SD<-function(all.out.old, all.out, return.values=F) {
  require("sp")
  # this function will estimate SD between iterations.
  Coords<-cbind(all.out.old$Final.Means$CENTRE.x, all.out.old$Final.Means$CENTRE.y, all.out$Final.Means$CENTRE.x, all.out$Final.Means$CENTRE.y)
  Dists<-apply(Coords, 1, FUN=function(x) {sp:::spDistsN1(cbind(x[[1]], x[[2]]), cbind(x[[3]], x[[4]]), longlat=TRUE)})
  if(return.values) {return(Dists)
  } else {return(sd(Dists))}
}

get.optimization.SD<-cmpfun(get.optimization.SD)


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

EM.PF.wrapper<-cmpfun(EM.PF.wrapper)
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

EM.PF.safe.ram.wrapper<-cmpfun(EM.PF.safe.ram.wrapper)


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


# ok this function can work only wihtj final matrix..
# 
# this function is needed to estimate CI
#################################################################
####
create.credible.intervals.array<-function(in.Data, intervals=c(0.8, 0.6, 0.4, 0.2) ) {
  require(rgeos)
  require(maptools)
  # 
  # intervals should be a vector like c(0.8, 0.6, 0.4)
  intervals=sort(intervals, decreasing=T)
  Polygons<-vector(mode="list", length=length(intervals))
  names(Polygons)<-paste("ci.", intervals, sep="")
  nParticles=sum(in.Data$Points.rle[[1]]$lengths) #dim(All.results.mat)[1]
  Seq<-intervals 
  Polygons<-lapply(Polygons, list)
  cat("estimating rle..")
  for (Iteration in 1:length(in.Data$Points.rle)) {
    #cat("Iteration", Iteration, "\n")
    points2read<-inverse.rle(in.Data$Points.rle[[Iteration]])
    Density<-cumsum(sort.int(table(points2read), decreasing=T))/nParticles
    for (int in 1:length(Seq)) {
      Pointstopoly<-as.integer(names(Density)[Density<=Seq[int]])
      if (length(Pointstopoly)==0) {
        Pointstopoly<-as.integer(names(Density)[1])}
      Pol<-Polygons(list(Polygon(in.Data$Points.Land[c(Pointstopoly, Pointstopoly[1]), ])), ID=Iteration)
      Polygons[[int]][[Iteration]]<-Pol
    }
  }
  cat("   Done!\n")
  Convexed<-list()
  ###
  # now we need to combine points by 2..
  All.cr<-list()
  for (cred in 1:length(intervals)) {
    Convexed<-list()
    cat(intervals[cred], "\n")
    Curr.cr<-Polygons[[cred]]
    for (i in 1:(length(Curr.cr)-1)) {
      Spat.Polygons<-lapply(Curr.cr, FUN=function(x) SpatialPolygons(list(x), proj4string=CRS("+proj=longlat +datum=WGS84")))
      if(nrow(Spat.Polygons[[i]]@polygons[[1]]@Polygons[[1]]@coords)<=3) 	{ Spat.Polygons[[i]]@polygons[[1]]@Polygons[[1]]@coords<-rbind(Spat.Polygons[[i]]@polygons[[1]]@Polygons[[1]]@coords, Spat.Polygons[[i]]@polygons[[1]]@Polygons[[1]]@coords, Spat.Polygons[[i]]@polygons[[1]]@Polygons[[1]]@coords)
      }
      if(nrow(Spat.Polygons[[i+1]]@polygons[[1]]@Polygons[[1]]@coords)<=3) 	{ Spat.Polygons[[i+1]]@polygons[[1]]@Polygons[[1]]@coords<-rbind(Spat.Polygons[[i+1]]@polygons[[1]]@Polygons[[1]]@coords, Spat.Polygons[[i+1]]@polygons[[1]]@Polygons[[1]]@coords, Spat.Polygons[[i+1]]@polygons[[1]]@Polygons[[1]]@coords,Spat.Polygons[[i+1]]@polygons[[1]]@Polygons[[1]]@coords)
      }
      if (i==1) {
        Convexed<-gConvexHull(gUnion(gConvexHull(Spat.Polygons[[i]]), gBoundary(Spat.Polygons[[i+1]])))
      }
      else {Convexed<-gUnion(Convexed, gConvexHull(gUnion(gConvexHull(Spat.Polygons[[i]]), gConvexHull(Spat.Polygons[[i+1]]))))}
    }
    All.cr[[cred]]<-Convexed
  }
  names(All.cr)<-names(Polygons)
  return(All.cr)
}

plot.ci.array<-function(CI.array, xlim=c(-115, -94), ylim=c(16,38.5), all.out) {
  require(maptools)
  my.golden.colors <- colorRampPalette(
    c("white","#FF7100"),
    bias=2.0)	
  Col=my.golden.colors(1+length(CI.array))
  Col=Col[-1]
  #topo.colors(length(CI.array)+5)
  for (i in 1:length(CI.array)) {
    if(i==1) {
      add=F 
    } else { 
      add=T
    }
    if (class(CI.array[[i]])=="SpatialCollections") {
      
      suppressWarnings(plot(CI.array[[i]]@polyobj, col=Col[i], border=grey(0.5), xlim=xlim, ylim=ylim, add=add))
      suppressWarnings(plot(CI.array[[i]]@lineobj, col=Col[i], border=grey(0.5), xlim=xlim, ylim=ylim, add=add))
    } else  {
      suppressWarnings(plot(CI.array[[i]], col=Col[i], border=grey(0.5), xlim=xlim, ylim=ylim, add=add))
    }
  }
  data("wrld_simpl", package="maptools")
  plot(wrld_simpl, add=T, border=grey(0.8), lwd=2)
  
  # add mean line:
  points(all.out$Final.Means$CENTRE.y~all.out$Final.Means$CENTRE.x, pch=3, cex=0.5, col="red")
  lines(all.out$Final.Means$CENTRE.y~all.out$Final.Means$CENTRE.x, pch=3, col="red", lwd=2)
  
  # release position
  points(-98.7, 34.7, pch=1, cex=2, col="black", lwd=3)
  
  # add bars every month:
  Index<-which(as.POSIXlt(all.out$Matrix.Index.Table$Real.time)$mday==1 &  all.out$Main.Index$proposal.index=="Dusk")
  points(all.out$Final.Means$CENTRE.y[Index]~all.out$Final.Means$CENTRE.x[Index], pch=3, cex=3, col="black", lwd=2)
  
}

## movie function:
## main difference of this function is that it uses all.out..
make.movie2<-function(in.Data, Data, Result, video.name="test1.mp4", start=c(-98.7, 34.7), add.boundaries=T) {
  # add.boundaries makes function much slower!!!
  
  # ok, this function will simply use all.out, that we create anyway...
  op<-par(no.readonly=T)
  
  require(animation)
  require(maptools)
  require(fields)
  my.golden.colors <- colorRampPalette(
    c("white","#FF7100"),
    bias=2.0)	
  if (add.boundaries) data(wrld_simpl)
  # new addition - as far as spatial points may have 3 columns
  # 1 additional for water
  # we could check for this column and then exclude all the points which have it..
  Start=in.Data$Points.Land[in.Data$start[1],1:2]
  if (dim(in.Data$Points.Land)[2]>2) {
    Index<-  which(in.Data$Points.Land[,3]>0)
    in.Data$Points.Land=in.Data$Points.Land[Index,]
    #Dawn.matrix=Dawn.matrix[Index,]
    #Dusk.matrix=Dusk.matrix[Index,]
  }
  
  # adding rows ID
  Index.Dawn<-which(in.Data$Matrix.Index.Table$Dusk==F)
  
  Index.Dusk<-which(in.Data$Matrix.Index.Table$Dusk==T)
  # create pairs of twilights
  # it means that we want to get closest dusks and then closest dawns to this dusks
  Index.Dusk<-Index.Dusk[sapply(Index.Dawn, FUN=function(x) {Z<-Index.Dusk-x; which(Z>0)[1]})]
  Index.Dawn<-Index.Dawn[!is.na(Index.Dusk)]
  Index.Dusk<-Index.Dusk[!is.na(Index.Dusk)]
  Index.Dawn<-Index.Dawn[!duplicated(Index.Dusk, fromLast=T)]
  Index.Dusk<-Index.Dusk[!duplicated(Index.Dusk, fromLast=T)]
  #=======================================================
  saveVideo({
    close.screen(all=TRUE)	
    plot(c(0,1), c(0,1), type="n", axes=F, main=paste("video created by make.movie2 function", Sys.time()))
    for (i in 1:length(Index.Dawn)) {
      Dawn.Column<-in.Data$Phys.Mat[,Index.Dawn[i]]
      Dusk.Column<-in.Data$Phys.Mat[,Index.Dusk[i]]
      if (!length(Dusk.Column)==0) {
        par(mar=c(2,1.7,1.5,3), ps=9)
        split.screen(c(2,1))
        split.screen(c(2,2), screen = 1)
        split.screen(c(1,3), screen = 2)
        
        cat("producing image #", i, "\n")
        screen(7)
        
        Image<-as.image(Dawn.Column[Index], x= in.Data$Points.Land, nrow=50, ncol=50)
        image(Image , main="Dawn", col=my.golden.colors(64))
        box()
        if (add.boundaries) plot(wrld_simpl, add=T, border=grey(0.6), lwd=2)
        abline(v=Start[1], col=grey(0.5), lwd=2, lty=2)
        abline(h=Start[2], col=grey(0.5), lwd=2, lty=2)
        image.plot(Image, col=my.golden.colors(64), add=T,  legend.mar=3.1, legend.only=T)
        ###########################
        ## Multiplication graph
        #par( my.par)
        screen(8)
        Image<-as.image(Dusk.Column[Index]*Dawn.Column[Index], x= in.Data$Points.Land, nrow=50, ncol=50)
        image(Image, main="Multiplication", col=my.golden.colors(64))
        box()
        if (add.boundaries) plot(wrld_simpl, add=T, border=grey(0.6), lwd=2)
        abline(v=Start[1], col=grey(0.5), lwd=2, lty=2)
        abline(h=Start[2], col=grey(0.5), lwd=2, lty=2)
        image.plot(Image, col=my.golden.colors(64), add=T,  legend.mar=3.1, legend.only=T)
        ###########################
        ## Dusk graph			
        #par( my.par)
        screen(9)
        Image<-as.image(Dusk.Column[Index], x= in.Data$Points.Land, nrow=50, ncol=50)
        image(Image, main="Dusk", col=my.golden.colors(64))
        box()
        if (add.boundaries) plot(wrld_simpl, add=T, border=grey(0.6), lwd=2)
        abline(v=Start[1], col=grey(0.5), lwd=2, lty=2)
        abline(h=Start[2], col=grey(0.5), lwd=2, lty=2)
        image.plot(Image, col=my.golden.colors(64), add=T,  legend.mar=3.1, legend.only=T)					##########################################
        ### new set
        #########################################
        screen(5)
        par(mar=c(2,1.7,0.5,0.5))			
        plot(Result$Final.dawn$Data$Hour~Result$Final.dawn$Data$gmt, type="n", ylim=c(min(Result$Final.dawn$Loess.predict$fit, Result$Final.dawn$Data$Hour), max(Result$Final.dawn$Loess.predict$fit, Result$Final.dawn$Data$Hour))) 
        abline(v=Result$Final.dawn$Data$gmt[which(Result$Final.dawn$Loess.predict$Border==1)], col=colors()[432], lwd=2)
        points(Result$Final.dawn$Data$Hour~Result$Final.dawn$Data$gmt, col=grey(0.6), type="p", pch=3) 
        ## ok, now - loess
        lines(Result$Final.dawn$Loess.predict$fit~Result$Final.dawn$Data$gmt, col="darkgreen", lwd=2)
        box()
        #lines(Result$Final.dawn$Loess.predict$fit+Result$Final.dawn$Loess.predict$se.fit*1.96~Result$Final.dawn$Data$gmt, col="darkgreen", lwd=2, lty=2)
        #lines(Result$Final.dawn$Loess.predict$fit-Result$Final.dawn$Loess.predict$se.fit*1.96~Result$Final.dawn$Data$gmt, col="darkgreen", lwd=2, lty=2)
        Index.in.Result<-which(Result$Final.dawn$Data$gmt==in.Data$ Matrix.Index.Table$time[Index.Dawn[i]])			
        
        abline(v=Result$Final.dawn$Data$gmt[Index.in.Result], col="red", lwd=2)
        ##
        screen(3)
        par(mar=c(1.5,1.7,1.5,0.5))			
        Loess.gmt<-Result$Final.dawn$Data$gmt.adj[Index.in.Result]
        Point<-	which.min(abs(as.numeric(Data$d$gmt)-as.numeric(Result$Final.dawn$Data$gmt[Index.in.Result])))
        diap<-(Point-14):(Point+15)
        #diap<-c((Loess.gmt-14):(Loess.id.round+15))
        plot(light~gmt, data=Data$d[diap,], type="p", pch=3, ylim=c(1,127), main= paste(format(Result$Final.dawn$Data$gmt[Index.in.Result], tz="UTC"), " (UTC)"))
        lines(light~gmt, data=Data$d[diap,], lwd=2)
        abline(v=Result$Final.dawn$Data$gmt[Index.in.Result], col=grey(0.6), lwd=2)
        abline(v=Loess.gmt, col="darkgreen", lwd=3)
        screen(6)
        par(mar=c(2,1.7,0.5,0.5))	
        plot(Result$Final.dusk$Data$Hour~Result$Final.dusk$Data$gmt, type="n", 	ylim=c(min(Result$Final.dusk$Loess.predict$fit, Result$Final.dusk$Data$Hour), max(Result$Final.dusk$Loess.predict$fit, Result$Final.dusk$Data$Hour))) 
        abline(v=Result$Final.dusk$Data$gmt[which(Result$Final.dusk$Loess.predict$Border==1)], col=colors()[432], lwd=2)
        points(Result$Final.dusk$Data$Hour~Result$Final.dusk$Data$gmt, col=grey(0.6), type="p", pch=3) 
        ## ok, now - loess
        lines(Result$Final.dusk$Loess.predict$fit~Result$Final.dusk$Data$gmt, col="darkgreen", lwd=2)
        box()	#lines(Result$Final.dusk$Loess.predict$fit+Result$Final.dusk$Loess.predict$se.fit*1.96~Result$Final.dusk$Data$gmt, col="darkgreen", lwd=2, lty=2)
        #lines(Result$Final.dusk$Loess.predict$fit-Result$Final.dusk$Loess.predict$se.fit*1.96~Result$Final.dusk$Data$gmt, col="darkgreen", lwd=2, lty=2)
        
        ### dusk
        Index.in.Result.Dusk<-which(Result$Final.dusk$Data$gmt==in.Data$Matrix.Index.Table$time[Index.Dusk[i]])
        abline(v=Result$Final.dusk$Data$gmt[Index.in.Result.Dusk], col="red", lwd=2)
        
        Loess.gmt.Dusk<-Result$Final.dusk$Data$gmt.adj[Index.in.Result.Dusk]
        
        Point<-	which.min(abs(as.numeric(Data$d$gmt)-as.numeric(Result$Final.dusk$Data$gmt[Index.in.Result.Dusk])))
        
        #Loess.id.round<-floor(Loess.id)
        #Loess.gmt.round<-floor(Loess.gmt)
        diap<-(Point-14):(Point+15)
        #diap<-c((Loess.gmt-14):(Loess.id.round+15))
        screen(4)
        par(mar=c(1.5,1.7,1.5,0.5))			
        plot(light~gmt, data=Data$d[diap,], type="p", pch=3, ylim=c(1,127), main= paste(format(Result$Final.dusk$Data$gmt[Index.in.Result.Dusk], tz="UTC"), " (UTC)"))
        lines(light~gmt, data=Data$d[diap,], lwd=2)
        abline(v=Loess.gmt.Dusk, col="darkgreen", lwd=3)
        abline(v=Result$Final.dusk$Data$gmt[Index.in.Result.Dusk], col=grey(0.6), lwd=2)
        ######	
        close.screen(all=TRUE)	
      }
    }
  }, video.name = video.name, other.opts = "-b 300k", ffmpeg="C:/Program Files/ffmpeg/bin/ffmpeg.exe", outdir = getwd(), imgdir=getwd(), ani.width=1200, ani.height=800, clean = TRUE)
  par(op)
}

