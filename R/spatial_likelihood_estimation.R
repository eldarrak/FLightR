# spatial_likelihood_estiamtion.R


get.Phys.Mat.parallel<-function(all.out=NULL, Twilight.time.mat.dusk=NULL, Twilight.log.light.mat.dusk=NULL, Twilight.time.mat.dawn=NULL, Twilight.log.light.mat.dawn=NULL,  threads=2,  calibration=NULL, log.light.borders=NULL, log.irrad.borders=NULL, likelihood.correction=TRUE ) {

# let's say we have to submit all boundaries inside tha calibration object..

if (is.character(all.out)) all.out=get("all.out")
if (is.character(Twilight.time.mat.dusk)) Twilight.time.mat.dusk=get("Twilight.time.mat.dusk")
if (is.character(Twilight.time.mat.dawn)) Twilight.time.mat.dawn=get("Twilight.time.mat.dawn")
if (is.character(Twilight.log.light.mat.dusk)) Twilight.log.light.mat.dusk=get("Twilight.log.light.mat.dusk")
if (is.character(Twilight.log.light.mat.dawn)) Twilight.log.light.mat.dawn=get("Twilight.log.light.mat.dawn")
if (is.character(calibration)) calibration=get("calibration")


Grid<-all.out$Spatial$Grid

cat("making cluster\n")
mycl <- parallel:::makeCluster(threads)
    tmp<-parallel:::clusterSetRNGStream(mycl)
    tmp<-parallel:::clusterExport(mycl,c("Twilight.time.mat.dawn", "Twilight.time.mat.dusk", "Twilight.log.light.mat.dawn", "Twilight.log.light.mat.dusk", "Grid", "calibration"), envir=environment())
    tmp<-parallel:::clusterEvalQ(mycl, library("FLightR")) 

#====================
tryCatch( {
cat("estimating dusks\n")

	Twilight.vector<-1:(dim(Twilight.time.mat.dusk)[2])
	
	All.probs.dusk<-parSapplyLB(mycl, Twilight.vector, FUN=get.prob.surface, Twilight.log.light.mat=Twilight.log.light.mat.dusk, Twilight.time.mat=Twilight.time.mat.dusk, dusk=TRUE, Calib.param=calibration$Parameters$LogSlope, log.light.borders=calibration$Parameters$log.light.borders, log.irrad.borders=calibration$Parameters$log.irrad.borders, delta=NULL, Grid=Grid, calibration=calibration, impute.on.boundaries=calibration$Parameters$impute.on.boundaries, likelihood.correction=likelihood.correction)
		 	
cat("estimating dawns\n")
	 
	Twilight.vector<-1:(dim(Twilight.time.mat.dawn)[2])
	All.probs.dawn<-parSapplyLB(mycl, Twilight.vector, FUN=get.prob.surface, Twilight.log.light.mat=Twilight.log.light.mat.dawn, Twilight.time.mat=Twilight.time.mat.dawn, dusk=FALSE, Calib.param=calibration$Parameters$LogSlope, log.light.borders=calibration$Parameters$log.light.borders, log.irrad.borders=calibration$Parameters$log.irrad.borders, delta=NULL, Grid=Grid, calibration=calibration, impute.on.boundaries=calibration$Parameters$impute.on.boundaries, likelihood.correction=likelihood.correction)
	
	}, finally = stopCluster(mycl))
	
    #if(!is.null(mycl)) {
    #   parallel::stopCluster(mycl)
    #   mycl <- c()
    #}
cat("processing results\n")
All.probs.dusk.tmp<-All.probs.dusk
All.probs.dawn.tmp<-All.probs.dawn
	Phys.Mat<-c()
for (i in 1:nrow(all.out$Indices$Matrix.Index.Table)) {
	if (all.out$Indices$Matrix.Index.Table$Dusk[i]) {
		Phys.Mat<-cbind(Phys.Mat, All.probs.dusk.tmp[,1])
		All.probs.dusk.tmp<-as.matrix(All.probs.dusk.tmp[,-1])
		} else {
		Phys.Mat<-cbind(Phys.Mat, All.probs.dawn.tmp[,1])
		All.probs.dawn.tmp<-as.matrix(All.probs.dawn.tmp[,-1])
		}
}

Phys.Mat<-apply(Phys.Mat, 2, FUN=function(x) {x[x<=1e-70]=1e-70; return(x)})

return(Phys.Mat)
}


get.prob.surface<-function(Twilight.ID, dusk=TRUE, Twilight.time.mat, Twilight.log.light.mat, return.slopes=FALSE,  Calib.param, log.irrad.borders=c(-9, 3), delta=0, Grid, log.light.borders=log(c(2,64)), calibration=NULL, impute.on.boundaries=FALSE, likelihood.correction=TRUE) {
 
		if (Twilight.ID%%10== 1) cat("doing", Twilight.ID, "\n")	
		
		Twilight.solar.vector<-solar.FLightR(as.POSIXct(Twilight.time.mat[c(1:24, 26:49), Twilight.ID], tz="gmt", origin="1970-01-01"))
		Twilight.log.light.vector<-Twilight.log.light.mat[c(1:24, 26:49), Twilight.ID]
		Twilight.time.vector=Twilight.time.mat[c(1:24, 26:49), Twilight.ID]
		
			if (is.null(calibration)) {
			time_correction=Calib.param[1] # this is the only place where I use calib.param...
			} else {
			time_correction=calibration$time_correction_fun(get.declination(Twilight.time.vector[24]), as.numeric(dusk))
			}
		
		if (return.slopes) {
		Current.probs<-	apply(Grid, 1, get.current.slope.prob, Twilight.solar.vector=Twilight.solar.vector, Twilight.log.light.vector=Twilight.log.light.vector, plot=FALSE, verbose=FALSE,  log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, dusk=dusk, return.slopes=TRUE, Twilight.time.vector=Twilight.time.vector,  Calib.param= Calib.param, delta=delta, time_correction=time_correction, calibration=calibration, impute.on.boundaries=impute.on.boundaries, likelihood.correction=likelihood.correction)	
			} else {
		Current.probs<-	apply(Grid, 1, get.current.slope.prob, Twilight.solar.vector=Twilight.solar.vector, Twilight.log.light.vector=Twilight.log.light.vector, plot=FALSE, verbose=FALSE,  log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, dusk=dusk, return.slopes=FALSE, Twilight.time.vector=Twilight.time.vector,  Calib.param= Calib.param, delta=delta, time_correction=time_correction, calibration=calibration, impute.on.boundaries=impute.on.boundaries, likelihood.correction=likelihood.correction)
		}
		return(Current.probs)
	}

			
			
get.probs.lm<-function(Model, plot=FALSE, calibration=NULL, time_correction=NULL, Calib.param=NULL, return.slopes=FALSE, Twilight.time.vector=NULL,  delta=0, dusk=TRUE, x=NULL, likelihood.correction=TRUE) {
		# check for the intercept
		sum=0
			
		#require(mvtnorm)
		coef.coef<-Model$coefficients
		slope.sd<-Model$stderr[2]

		
		#if (coef(Model)[1] < calibration$calibration.bayesian.model$Intercept.boundary[1]) { 
		#if (is.na(coef.vcov[1])) coef.vcov[is.na(coef.vcov)]<-calibration$calibration.bayesian.model$Slope.integration$fast.tw.vcov # adding errors from MCMC
		if (!is.finite(slope.sd)) {
			if (is.null(calibration)) {
			slope.sd=0 
			} else slope.sd<- calibration$Parameters$mean.of.individual.slope.sigma
		}
		#if (length(resid(Model))== 3) coef.vcov<-coef.vcov*15
		# coef.vcov<-coef.vcov*100/(length(resid(Model))^2) # adding errors from MCMC
		#if (use.intercept) {
		#stop("use of intercept is not implemented!!!")
		#} else { 

		#test.Slope<-rnorm(1000, coef.coef[2], sqrt(coef.vcov[4]))
		if (slope.sd==0) {
      		test.Slope=coef.coef[2]
		} else {
     		test.Slope<-rnorm(1000, coef.coef[2], slope.sd)
		}
		if (is.null(time_correction)) {	

		if (is.null(calibration)) {
			time_correction=Calib.param[1] # this is the only place where I use calib.param...
			} else {
			#------------------------------
			# here is a change to sun declination from cosSolarDec..
			#time_correction=calibration$time_correction_fun(Twilight.solar.vector$cosSolarDec[1])
			time_correction=calibration$time_correction_fun(get.declination(Twilight.time.vector[24]), as.numeric(dusk))
			}
		}

		Expected.mean<-time_correction
		
	###########################################
        #LL correction fgoes in here	
		if (likelihood.correction) {
           if(is.null(calibration$c_fun)) {
		   warning('you must supply new calibration with calibration correction function!\n') 
		} else {
		test.Slope<-test.Slope-calibration$c_fun(slope.sd)
		}
		}
		
		if (is.null(delta)) {
		if (length(formals(calibration$lat_correction_fun))==1) {
		delta=calibration$lat_correction_fun(x[2])} else {
		delta=calibration$lat_correction_fun(x[2], Twilight.time.vector[13]) # have to check whether we always have time specified...
		}
		}
		#sum<-mean(dlnorm(test.Slope, Expected.mean+delta, Calib.param[2]))

		# correction for both parameters
		Probability<-mean(dlnorm(test.Slope, Expected.mean+delta, Calib.param[2]))
		#-----------------
		# new exp correction added 24 Apr 2015
		#sum<-mean(dlnorm(test.Slope, Expected.mean+delta, Calib.param[2])*test.Slope)
		#cat("time corr", time_correction ,"Expected.mean+delta", Expected.mean+delta, "\n")
		#-----------------------
				# correct for cos Lat
				# looks like cos is too much... should look at sqrt(cos)
				#sum=sum*(cos(x[2]/180*pi))
				#sum=sum*(cos(x[2]/180*pi)^0.5)
				#-----------------------
		
		#		}
		#if (plot) {
		#my.golden.colors <- colorRampPalette(
		#c("white","#FF7100"))
		#image(list(x=calibration$calibration.bayesian.model$Slope.integration$MCMC.All.densities$x, y=calibration$calibration.bayesian.model$Slope.integration$MCMC.All.densities$y,
		#z=matrix(dmvnorm(expand.grid(calibration$calibration.bayesian.model$Slope.integration$MCMC.All.densities$x, calibration$calibration.bayesian.model$Slope.integration$MCMC.All.densities$y), coef.coef, coef.vcov), nrow = length(calibration$calibration.bayesian.model$Slope.integration$MCMC.All.densities$x), ncol = length(calibration$calibration.bayesian.model$Slope.integration$MCMC.All.densities$y))), col=my.golden.colors(10))
		#contour(calibration$calibration.bayesian.model$Slope.integration$MCMC.All.densities, add=TRUE)
		#}
		#}
		if (!is.finite(Probability)) Probability<-0
		if (Probability<0) Probability<-0
		if (return.slopes) 	Probability<-c(Probability, coef.coef[2], Model$stderr[2])
		return(Probability)
		}

			
get.probs.nonparam.slope<-function(Slopes, plot=FALSE, calibration=NULL, time_correction=NULL, Calib.param=NULL, return.slopes=FALSE, Twilight.time.vector=NULL,  delta=0, dusk=TRUE) {

		if (is.null(time_correction)) {	

			if (is.null(calibration)) {
			time_correction=Calib.param[1] # this is the only place where I use calib.param...
			} else {
			time_correction=calibration$time_correction_fun(get.declination(Twilight.time.vector[24]), as.numeric(dusk))
			}
		}
		Expected.mean<-time_correction
	
		if (is.null(delta)) {
		if (length(formals(calibration$lat_correction_fun))==1) {
		delta=calibration$lat_correction_fun(x[2])} else {
		delta=calibration$lat_correction_fun(x[2], Twilight.time.vector[13]) # have to check whether we always have time specified...
		}
		}
		#Probability<-mean(dlnorm(Slopes, Expected.mean+delta, Calib.param[2]))
		Probability<-sum(dlnorm(Slopes, Expected.mean+delta, Calib.param[2]))
		if (!is.finite(Probability)) Probability<-0
		if (Probability<0) Probability<-0
		if (return.slopes) 	Probability<-c(Probability, mean(Slopes), sd(Slopes))
		return(Probability)
		}
		

get.current.slope.prob <-function (x, calibration = NULL, Twilight.solar.vector = NULL,   Twilight.log.light.vector, plot = FALSE, verbose = FALSE, log.light.borders = log(c(2, 64)), log.irrad.borders = c(-9, 1.5), dusk = TRUE, use.intercept = FALSE, return.slopes = FALSE, Twilight.time.vector = NULL, delta = 0, time_correction = NULL, Calib.param = NULL, impute.on.boundaries = FALSE, likelihood.correction=TRUE) {
    if (is.null(time_correction) & is.null(Twilight.solar.vector)) 
        stop("either time_correction or Twilight.solar.vector should be provided to get.current.slope.prob!")
		
	if (is.null(Twilight.solar.vector))  {
	Twilight.solar.vector<-solar.FLightR(as.POSIXct(Twilight.time.vector, tz="gmt", origin="1970-01-01"))
	}    
	Probability = 0
    Data <- check.boundaries(x, Twilight.solar.vector = Twilight.solar.vector, 
        Twilight.log.light.vector, plot = FALSE, verbose = verbose, 
        log.light.borders = log.light.borders, log.irrad.borders = log.irrad.borders, 
        dusk = dusk, Twilight.time.vector = Twilight.time.vector, 
        impute.on.boundaries = impute.on.boundaries)
    if (dim(Data)[1] > 1) {
        LogLight <- Data[, 1]
        LogIrrad <- Data[, 2]
        if (length(LogLight) >= 1) {
            if (calibration$Parameters$calibration.type == "parametric.slope") {
                #Model <- lm(LogLight ~ LogIrrad)
				Model<-fastLmPure(matrix(c(rep(1,length(LogIrrad)),LogIrrad), ncol=2), LogLight)

                #if (verbose) 
                  #print(summary(Model))
                Probability <- get.probs.lm(Model, plot = plot, 
                  calibration = calibration, time_correction = time_correction, 
                  Calib.param = Calib.param, return.slopes = return.slopes, 
                  delta = delta, Twilight.time.vector = Twilight.time.vector, 
                  dusk = dusk, x=x, likelihood.correction=likelihood.correction)
            }
            if (calibration$Parameters$calibration.type == "nonparametric.slope") {
                Slopes <- diff(LogLight)/diff(LogIrrad)
                Probability <- get.probs.nonparam.slope(Slopes, 
                  plot = plot, calibration = calibration, time_correction = time_correction, 
                  Calib.param = Calib.param, return.slopes = return.slopes, 
                  delta = delta, Twilight.time.vector = Twilight.time.vector, 
                  dusk = dusk)
            }
        }
        else {
            if (return.slopes) 
                Probability <- c(Probability, NA, NA)
            if (verbose) {
                print(str(calibration, max.level = 1))
                cat("calibration$Parameters:\n")
                print(calibration$Parameters)
            }
        }
    }
    else {
        Probability <- 0
        if (return.slopes) 
            Probability <- c(Probability, NA, NA)
    }
    return(Probability)
}

check.boundaries<-function(x, Twilight.solar.vector=NULL,  Twilight.log.light.vector, plot=FALSE, verbose=FALSE,  log.light.borders=log(c(2,64)), log.irrad.borders=c(-15, 50), dusk=TRUE, impute.on.boundaries=FALSE, Twilight.time.vector=NULL) {
# this function...
	if (is.null(Twilight.solar.vector))  {
	Twilight.solar.vector<-solar.FLightR(as.POSIXct(Twilight.time.vector, tz="gmt", origin="1970-01-01"))
	}
	Elevs<-elevation.FLightR(x[[1]], x[[2]], Twilight.solar.vector)
	LogIrrad<-log(get.Irradiance(Elevs*pi/180)+1e-20)
	LogLight<-Twilight.log.light.vector
	if (verbose) {
		cat("Elevs:\n")
		print(cbind(ID=1:length(LogLight), LogLight=LogLight, LogIrrad=LogIrrad))
			}
	#NotZero<-	which(LogLight>=log.light.borders[1] & LogLight<=log.light.borders[2] & LogIrrad<log.irrad.borders[2])
	#------------------------------------
	# version 0.3.6 after impute on.boundaries turned to FALSE figured out that there is not cut for low log irradiance so added it here..
	if (impute.on.boundaries) {
	NotZero<-	which(LogLight>=log.light.borders[1] & LogLight<=log.light.borders[2] & LogIrrad<log.irrad.borders[2])
	} else {
	NotZero<-	which(LogLight>=log.light.borders[1] & LogLight<=log.light.borders[2] & LogIrrad<log.irrad.borders[2] &  LogIrrad>log.irrad.borders[1])
	}
	# here we should make a new cut!
	# the idea is that we want to find at least two consequent points at the maximum
	Border.points<-which(LogLight>=log.light.borders[2])
	if (length(Border.points)==0) {
	if (dusk) {
	Border.points=1
	} else Border.points=48
	}
	if (dusk) {
	Border.points<-rev(Border.points)
	}
	if (verbose) print(Border.points)
	Diff.Border.points.lag1<- abs(diff(Border.points))
	Diff.Border.points.lag2<- abs(diff(Border.points, lag=2))
	Cut<-Border.points[intersect(which(Diff.Border.points.lag1==1) , which(Diff.Border.points.lag2==2))]
	if (length(Cut)==0) Cut=Border.points
	if (verbose) print(Cut)
	#NotZero<-NotZero[NotZero>=rev(Cut)[1]]
	if (dusk) {
	NotZero<-NotZero[NotZero>=Cut[1]]
	} else {
	NotZero<-NotZero[NotZero<=Cut[1]]
	}
	# and let's exclude points that are not in between crossing values..
	#if length(NotZero>0) {
	#Diap<-sort(c(LogLight[NotZero-1], LogLight[NotZero+1]))
	#NotZero<-NotZero[LogLight[NotZero]>Diap[1]&LogLight[NotZero]<Diap[2]]
	#}
	# here there is a trick - there is an artificial light, that does happen sometimes... so we can't skipt the whole thing in case we have one point that is lower than the boundary..
	
	#if (length(which(LogIrrad[NotZero]<log.irrad.borders[1] ))>0) {
		# here we have stop for the case when the sky is too clear
	#	return(matrix(nrow=0, ncol=0))
	#} else {
	# protection from the next twilight - extracting consequent measurements
	if (length(NotZero)>0) {
#print(LogIrrad)
#print(LogLight)
#print(NotZero)
	if (any(LogLight[NotZero[1]:NotZero[length(NotZero)]]>log.light.borders[2])) {
		if (Elevs[1]>rev(Elevs)[1]) {
			Cut.point<-rev(c(NotZero[1]:NotZero[length(NotZero)])[which(LogLight[NotZero[1]:NotZero[length(NotZero)]]>log.light.borders[2])])[1]
			NotZero<-NotZero[NotZero>Cut.point]
		} else {
			Cut.point<-c(NotZero[1]:NotZero[length(NotZero)])[which(LogLight[NotZero[1]:NotZero[length(NotZero)]]>log.light.borders[2])[1]]
			NotZero<-NotZero[NotZero<Cut.point]
		}
	}	
	}	
		
	#====================
	# here is the place to cut off the points that cross the maximum as in high latitudes there night can be very short...
	# 
	if (dusk) { 
		Index<-1:28
	} else {
		Index<-21:length(Elevs)
	}
	NotZero<-intersect(NotZero, Index)
	#=============================================
	# ok here I want to add light boundary estimates...
	# the idea is following - we assume that our boundary is exclusive
	# so we need to add a point that will be in half distance in time and will have the boundary value..
	# these points will always exist, so we don't want to check anything..
	# exept they could be crossing the Irrad boundary so this will need a check..
	First.LogIrrad<-Last.LogIrrad<-c()
	# LogIrrad
	if (verbose) print(NotZero)
	if (length(NotZero)>0) {
		#NotZero<-NotZero[1]:rev(NotZero)[1]
		First.LogIrrad<-NA
		if (dusk & !(LogLight[NotZero[1]] ==  log.light.borders[2])) First.LogIrrad<-LogIrrad[(NotZero[1]-1):NotZero[1]]
		if (!dusk & !(LogLight[NotZero[1]] ==  log.light.borders[1])) First.LogIrrad<-LogIrrad[(NotZero[1]-1):NotZero[1]]
		Last.LogIrrad<-NA
		if (dusk & !(LogLight[rev(NotZero)[1]] == log.light.borders[1])) Last.LogIrrad<-LogIrrad[(rev(NotZero)[1]):(rev(NotZero)[1]+1)]	
		if (!dusk & !(LogLight[rev(NotZero)[1]] ==  log.light.borders[2])) Last.LogIrrad<-LogIrrad[(rev(NotZero)[1]):(rev(NotZero)[1]+1)]	
		if (!is.null(Twilight.time.vector)) {
			First.Time<-Twilight.time.vector[(NotZero[1]-1):NotZero[1]]
			Last.Time<-Twilight.time.vector[(rev(NotZero)[1]):(rev(NotZero)[1]+1)]	
			}
		
	} else {
		# let's find the crossing
		Larger<-which(LogLight>= log.light.borders[1] & LogLight>= log.light.borders[2])
		Smaller<-which(LogLight<= log.light.borders[1] & LogLight<= log.light.borders[2])
		if (length(Larger)>0 & length(Smaller>0)) {
			#NotZero=c(24,25)
			if (dusk) {
				Interval<-c(rev(Larger)[1], Smaller[1])
				#LogLight[25]<-mean(LogLight[(Larger[1]-1):Larger[1]])
				#LogIrrad[25]<-mean(LogIrrad[(Larger[1]-1):Larger[1]])
			} else {
				Interval<-c(rev(Smaller)[1], Larger[1])
				#LogLight[25]<-mean(LogLight[(Smaller[1]-1):Smaller[1]])
				#LogIrrad[25]<-mean(LogIrrad[(Smaller[1]-1):Smaller[1]])
			}
			First.LogIrrad<-c(LogIrrad[Interval][1], mean(LogIrrad[Interval]))
			Last.LogIrrad<-c( mean(LogIrrad[Interval]), LogIrrad[Interval][2])
			if (!is.null(Twilight.time.vector)) {
			First.Time<-c(Twilight.time.vector[Interval][1], mean(Twilight.time.vector[Interval]))
			Last.Time<-c( mean(Twilight.time.vector[Interval]), Twilight.time.vector[Interval][2])
			}
			
		}
	}
	if (verbose) {
		cat("First.LogIrrad:", First.LogIrrad, "\n")
		cat("Last.LogIrrad:", Last.LogIrrad, "\n")
	}
	if (any(is.na(First.LogIrrad))) First.LogIrrad<-c()
	if (any(First.LogIrrad<=log.irrad.borders[1]) | any(First.LogIrrad>=log.irrad.borders[2])) First.LogIrrad<-c()
	
	if (any(is.na(Last.LogIrrad))) Last.LogIrrad<-c()

	if (any(Last.LogIrrad<log.irrad.borders[1]) | any(Last.LogIrrad>log.irrad.borders[2])) Last.LogIrrad<-c()
	if (!is.null(Twilight.time.vector)) {
			if(length(First.LogIrrad)==0) First.Time <-c()
			if(length(Last.LogIrrad)==0) Last.Time <-c()
			}
	if (impute.on.boundaries) {
	Imputing<-c() # Index for imputing
	if (length(First.LogIrrad)==2) Imputing<-c(Imputing, 1)
	if (length(Last.LogIrrad)==2) Imputing<-c(Imputing, 2)
	
	Left.Time<-Right.Time<-Right.LogIrrad<-Right.LogLight<-Left.LogIrrad<-Left.LogLight<-c()
	
	# light
	if (dusk) {
		if (1 %in% Imputing) {
			Left.LogLight<-max(log.light.borders)
			#Left.LogIrrad<-mean(First.LogIrrad)
			Left.LogIrrad<-log(mean(exp(First.LogIrrad)))
			#Left.LogIrrad<-max(First.LogIrrad)
			#if (!is.null(Twilight.time.vector)) Left.Time<-mean(First.Time)
			if (!is.null(Twilight.time.vector)) Left.Time<-min(First.Time)
		}
		if (2 %in% Imputing) {
			Right.LogLight<-min(log.light.borders)
			#Right.LogIrrad<-mean(Last.LogIrrad)
			Right.LogIrrad<-log(mean(exp(Last.LogIrrad)))
			#Right.LogIrrad<-min(Last.LogIrrad)
			#if (!is.null(Twilight.time.vector)) Right.Time<-mean(Last.Time)
			if (!is.null(Twilight.time.vector)) Right.Time<-max(Last.Time)

		}		
	} else {
		if (1 %in% Imputing) {
			Left.LogLight<-min(log.light.borders)
			#Left.LogIrrad<-mean(First.LogIrrad)
			Left.LogIrrad<-log(mean(exp(First.LogIrrad)))
			#Left.LogIrrad<-min(First.LogIrrad)
			#if (!is.null(Twilight.time.vector)) Left.Time<-mean(First.Time)
			if (!is.null(Twilight.time.vector)) Left.Time<-min(First.Time)


		}
		if (2 %in% Imputing) {
			Right.LogLight<-max(log.light.borders)
			#Right.LogIrrad<-mean(Last.LogIrrad)
			Right.LogIrrad<-log(mean(exp(Last.LogIrrad)))
			#Right.LogIrrad<-max(Last.LogIrrad)
			#if (!is.null(Twilight.time.vector)) Right.Time<-mean(Last.Time)
			if (!is.null(Twilight.time.vector)) Right.Time<-max(Last.Time)

		}		
	}
	New.Index<-which(c(Left.LogIrrad, LogIrrad[NotZero], Right.LogIrrad)>log.irrad.borders[1])

	if (verbose) {cat("nz\n"); print(NotZero)}

	Res<-cbind(c(Left.LogLight, LogLight[NotZero], Right.LogLight)[New.Index], c(Left.LogIrrad, LogIrrad[NotZero], Right.LogIrrad)[New.Index])
	if (!is.null(Twilight.time.vector)) Res<-cbind(Res, c(Left.Time, Twilight.time.vector[NotZero], Right.Time)[New.Index])
	
	# I want to output angle values now..
	Res<-cbind(Res, c(NA, Elevs[NotZero], NA)[New.Index])
	
	} else {
	Res<-cbind(LogLight[NotZero], LogIrrad[NotZero])
	if (!is.null(Twilight.time.vector)) Res<-cbind(Res,  Twilight.time.vector[NotZero])

	# I want to output angle values now..
	Res<-cbind(Res, Elevs[NotZero])

	}

	if (verbose) print(Res)
	# and now we want to return 
	if (plot) {
		Coef<-try(coef(lm(Res[,1]~Res[,2])))
		if (class(Coef)!='try-error') {
		    plot(LogLight~LogIrrad, main=paste(twilight=as.POSIXct(Twilight.time.vector[24], tz="gmt", origin="1970-01-01"), ifelse(dusk, "dusk", "dawn"), "intercept", Coef[1], "slope", Coef[2]), xlim=c(log.irrad.borders[1]-3, log.irrad.borders[2]+3), ylim=c(log.light.borders[1]-1, log.light.borders[2]+1))
		    print(Coef)
		} else {
		   plot(LogLight~LogIrrad, main=paste(twilight=as.POSIXct(Twilight.time.vector[24], tz= "gmt", origin="1970-01-01"), ifelse(dusk, "dusk", "dawn"), "intercept NA slope NA"), xlim=c(log.irrad.borders[1]-3, log.irrad.borders[2]+3))
		}
		   points(Res[,1]~Res[,2], col="red", lwd=2, pch="+")
		   par(ask = TRUE)
	    if (class(Coef)!='try-error') {
		plot(LogLight~LogIrrad, main=paste(twilight=as.POSIXct(Twilight.time.vector[24], tz="gmt", origin="1970-01-01"), ifelse(dusk, "dusk", "dawn"), "intercept", Coef[1], "slope", Coef[2]), xlim=c(log.irrad.borders[1]-3, log.irrad.borders[2]+3), ylim=c(log.light.borders[1]-1, log.light.borders[2]+1))
		} else {
		   plot(LogLight~LogIrrad, main=paste(twilight=as.POSIXct(Twilight.time.vector[24], tz= "gmt", origin="1970-01-01"), ifelse(dusk, "dusk", "dawn"), "intercept NA slope NA"), xlim=c(log.irrad.borders[1]-3, log.irrad.borders[2]+3))
		}
		par(ask=FALSE)
		if (Coef[1]<(-6)) warning("check twilight at around ", as.POSIXct(Twilight.time.vector[24], tz="gmt", origin="1970-01-01"), " it had strange shading\n")

	}
	return(Res)
	#}
}


