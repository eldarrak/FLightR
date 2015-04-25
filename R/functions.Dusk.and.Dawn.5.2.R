# ver 5.2 now I want to provide lat and time correction functions for every run and not to use the directory with saved files..
# so the calibration object now should have 
#	$time_correction_fun
#	$lat_correction_fun


# ver 5.1 will work with a new formalized calibration..
# there will be two functions one for time correction and one for Lat correction.
# time_correction_fun, file="time_correction_fun_600.RData")
# the correction shows what would be the estimated mean if you had cosSolarDec like that..
# so we could add a time_correction argument that will be jsut a number..





# 5.0 is based on the 4.0.d version. the addition is correction for time and latitude with GAM - I will first try to add gam without hte saving.period and will jsut assume that saving.period is 600..
# Also the parameters for the input should be changed to 0.67 and 0.4 (from 0.21 and 0.4)
 
# new 4.0.d will get cos(Lat) and then I'll go for a new gam - lat estimation
# 4.0.d now has got delta - this is a parameter we want to optimize..
# ver 4.0.b - just an attempt to use onoly one value - but not the distribution from LM
# ver 4.0 is version 3.4 with small changes taken from 4.0
#changed slope only prediction to density function.
# ver. 3.4 decided to make a better cutting off as we need to cut if there are 3 or more consequent values at the boundary.. at least two!

# ver 3.3 decided to try joint distribution of slope and Intercept

# ver.3.2 adding integration of intecept..
# ver.3.1 new model without any interactions.. based on general calibration


# 2.22  decided to try to delete interaction as it is probably not optimized welll..

# 2.21 - decide to try make boundary imputing better by making it at the log mean but not the log..
# 2.20 decided to multiply variances. to make sd of slope 0.15
# 

# 2.16 DECIDED TO ADD A TRUNCATION POINT..
# ver 15 - follows ver 13 22-10-2013 
# decided to combine the two models in one
# the idea is very simple - we just run the slope model then take the esitmated intercept velues and look at their distribution..
# the intercept model should be taken from that results..

# ver 13 - SAA interaction added.
# ver 12 - addition of SAA random effects integration.

## in this version I'll try to use an interaction of slope and intercept..

## in the version 10 I'll try to use both distribution - intercept and slope..

## in the new set of functions i'll focus on the correct mixed model and predictions from it.. 
## the basic idea is following:
## 1. during the calibration we will figure out how what is the slope of the base line..
## 2. the we will run a mixed model on the corrected for slope estimates
## 3. from this model we will get the distribution of slopes
## 4. for each day we will do a simple gamma regression and the will estimate the probability of the current slope..

# July 1 2013
# the new idea is that we could try to optimize the distribution of errors during the template fitting..
# sfor this we have to do more complicated calibration that will estimate separate thresholds but the joint shape of the distribution..
# I am not sure I know how can we do this..
# on of the ways will be 

## 24 June 2013
## Now I have decided to go for the lognormal distibution of errors...
## the idea is that error between the real twilight and clear sky twiligt should be 0 in case there is no error and should follow the lognormal distribution otherwise...
## so we will have to estimate parameters of the main distribution and then estimate values for the real points...
## how we will do that?
## 1. we need to calibrate logger again in a new way that will allow us to extract parameters of the distribution during calibration..
## 2. we also can add a refraction


require(compiler)

#get.Irradiance<-function(alpha, r=6378, s=6.9, refraction.correction=F) {

get.Irradiance<-function(alpha, r=6378, s=6.9) {
	# function from Ekstrom 2007
	erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
	## (see Abramowitz and Stegun 29.2.29)
	## and the so-called 'complementary error function'
	erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
	## and the inverses
	#erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
	#erfcinv <- function (x) qnorm(x/2, lower = FALSE)/sqrt(2)
	
	Res<-alpha
	u<-sqrt(r/(2*s))*sin(alpha)
	Res[(u<=0)]<-(exp(-u^2)/(1+erf(-u)))[(u<=0)]
	Res[(u>0)]<-(exp(-u^2)/(erfc(u)))[(u>0)]
	return(Res)
}

get.Irradiance<-cmpfun(get.Irradiance)
# the new idea is that we will give a vector of coordinates...

logger.template.calibrarion.internal<-function( Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk, positions=NA, plot=T, log.light.borders=NA,  log.irrad.borders=c(-8, 1.3), adjust.variance=T, impute.on.boundaries=T) {

	# =================
	# in this function I'll add a new lnorm calibration...
	#
	# now I'll try to do them both..
	#require(nlme)
	if (!is.list(positions)) {
		if (length (positions) !=2) stop("the positions should be a list with $dawn  and $dusk dataframes with pairs of coordiantes OR just one pair\n")
		dawn=matrix(ncol=2, nrow=dim(Twilight.time.mat.Calib.dawn)[2])
		dawn[,1]<-positions[1]
		dawn[,2]<-positions[2]
		#
		dusk=matrix(ncol=2, nrow=dim(Twilight.time.mat.Calib.dusk)[2])
		dusk[,1]<-positions[1]
		dusk[,2]<-positions[2]
		#
		positions=list(dawn=dawn, dusk=dusk)
	} 
	Twilight.time.mat.Calib.dawn<-Twilight.time.mat.Calib.dawn[-25,]
	Twilight.log.light.mat.Calib.dawn<-Twilight.log.light.mat.Calib.dawn[-25,]
	Twilight.time.mat.Calib.dusk<-Twilight.time.mat.Calib.dusk[-25,]
	Twilight.log.light.mat.Calib.dusk<-Twilight.log.light.mat.Calib.dusk[-25,]
	# let's try to create Calib.data.all first!!
	Calib.data.dawn<-data.frame()
	if (plot) par(ask=F)
	for (dawn in 1:dim(Twilight.time.mat.Calib.dawn)[2]) {
cat("checking dawn", dawn, "\n" )
		#Twilight.solar.vector<-solar(as.POSIXct(Twilight.time.mat.Calib.dawn[, dawn], tz="gmt", origin="1970-01-01"))
		Data<-check.boundaries(positions$dawn[dawn,], Twilight.solar.vector=NULL,  Twilight.log.light.vector = Twilight.log.light.mat.Calib.dawn[,dawn], plot=plot, verbose=F,  log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, dusk=F, Twilight.time.vector=Twilight.time.mat.Calib.dawn[, dawn], impute.on.boundaries=impute.on.boundaries)
		if (length(Data)==0) {
		cat ("dawn", dawn, "was excluded from the calibration\n")
		} else {
		Calib.data.dawn<-rbind(Calib.data.dawn, cbind(LogLight=Data[,1], LogIrrad=Data[,2], Day=dawn, Time=Data[,3], Elevs=Data[,4]))	
		}		
	}
	Calib.data.dawn$type<-"Dawn"

	Calib.data.dusk<-data.frame()
	for (dusk in 1:dim(Twilight.time.mat.Calib.dusk)[2]) {
cat("checking dusk", dusk, "\n" )

		#Twilight.solar.vector<-solar(as.POSIXct(Twilight.time.mat.Calib.dusk[, dusk], tz="gmt", origin="1970-01-01"))
		Data<-check.boundaries(positions$dusk[dusk,], Twilight.solar.vector=NULL,  Twilight.log.light.vector=Twilight.log.light.mat.Calib.dusk[,dusk], plot=plot, verbose=F,  log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, dusk=T,  Twilight.time.vector=Twilight.time.mat.Calib.dusk[, dusk], impute.on.boundaries=impute.on.boundaries)
#print(str(Data)	)
		if (length(Data)==0) {
		cat ("dusk", dusk, "was excluded from the calibration\n")
		} else {
		Calib.data.dusk<-rbind(Calib.data.dusk, cbind(LogLight=Data[,1], LogIrrad=Data[,2], Day=dim(Twilight.time.mat.Calib.dawn)[2]+dusk, Time=Data[,3], Elevs=Data[,4]))
		}
	}
	Calib.data.dusk$type<-"Dusk"
#print(str(Calib.data.dawn))
#print(str(Calib.data.dusk))
	Calib.data.all<-rbind(Calib.data.dawn, Calib.data.dusk)
	#print(str(Calib.data.all))	
	
	if (plot) {
		plot(LogLight~LogIrrad,data=Calib.data.all, type="n")
		for (Day in unique(Calib.data.all$Day)) {
		lines(LogLight~LogIrrad, data=Calib.data.all[Calib.data.all$Day==Day,], type="b", pch="+", col=rainbow(length(unique(Calib.data.all$Day)))[Day])
		}
	}
	
	
		# this works in case of rare points but with slope
	# this works in case of rare points but with slope
	Selected.days<-rle(Calib.data.all$Day)$values[which((rle(Calib.data.all$Day)$length)>1)]
	#==================================================
	# this is the part for the estimation of slope when we don't have enough data.
	# it looks like when we don't have enough points it is better just do estimate a slope value from all the slope values we have
			
	#m1<-lme(LogLight~LogIrrad,data=Calib.data.all[Calib.data.all$Day %in%Selected.days,],  random=~1|Day)

	#m2<-lme(LogLight~LogIrrad*type,data=Calib.data.all[Calib.data.all$Day %in%Selected.days,],  random=~1|Day)
	#Significance.of.dusk.dawn.diff=	prod(summary(m2)$tTable[3:4,5])

	####
	
	Calib.data.all$fDay<-as.factor(Calib.data.all$Day)

	#m2<-lmer(LogLight~LogIrrad+LogIrrad:as.factor(type)+fDay-1+(0+LogIrrad|fDay),data=Calib.data.all) 
	#Significance.of.dusk.dawn.diff<-rev(summary(m2)@coefs[,2])[1]
		
	#==================
	# here is the part for ver. 4.0
	# hm.. potentially we don't really want to calibrate inside the function..
	# so we will add this later but for now I would just do it outside..
		
	return(list(Calib.data.all=Calib.data.all
	#, Calib.data.all.lme=Calib.data.all.lme.no_interact, Significance.of.dusk.dawn.diff=Significance.of.dusk.dawn.diff, 
	#calibration.slope.par=list(
	#	Calibration.simple.slope=Calibration.simple.no_interact,
	#	tmp.x.slope=tmp.x.RE,
	#	calibration.slope.distr=tmp.dnorm.calib.no_interact.RE,
	#	RE.component=RE.component), 
	#calibration.intercept.par=list(
	#	Parameters=Parameters, CalibrationL.LL=Optim.result$value),
	#calibration.bayesian.model=calibration.bayesian.model # this is a stored list...
	))
	
}

# ok now we need template.calibration.function..
logger.template.calibration<-function(Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk, time.shift=0, positions, log.light.borders=NA,  log.irrad.borders=c(-15, 50), adjust.variance=F, plot=T, impute.on.boundaries=T) {
	
	Calibration.original<-logger.template.calibrarion.internal(Twilight.time.mat.Calib.dawn+time.shift, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk+time.shift, Twilight.log.light.mat.Calib.dusk, positions=positions, log.light.borders=log.light.borders,  log.irrad.borders= log.irrad.borders, adjust.variance=adjust.variance, plot=plot,impute.on.boundaries=impute.on.boundaries)
	
	# for lme
	#Calibration.simple<-summary(m1)$tTable[1,1:2]
	#Calibration.simple<- summary(model.lm$m2)$coefficients [1,1:2]
	#tmp.xmin<-(Calibration.simple[1]-3*Calibration.simple[2])
	#tmp.xmax<-(Calibration.simple[1]+3*Calibration.simple[2])
	#tmp.dx  <-(tmp.xmax-tmp.xmin)/100
	#tmp.x   <-seq(tmp.xmin,tmp.xmax,by=tmp.dx)
	#tmp.dnorm.calib<-dnorm(tmp.x,Calibration.simple[1],Calibration.simple[2])
	#Dusk.Dawn.difference<-summary(model.lm$m2)$coefficients[2,]
	#calibration<-list(calibration.simple=Calibration.simple, tmp.x=tmp.x, tmp.dnorm.calib=tmp.dnorm.calib, m1=model.lm$m1, m2=model.lm$m2 )
	return(Calibration.original)
	}

logger.template.calibration<-cmpfun(logger.template.calibration)
	
get.time.shift<-function(start, Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk,  log.light.borders=NA, diap=c(-600, 600), plot=T,  log.irrad.borders=c(-15, 50), verbose=F) {

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
	
# this will use the old idea about slopes..
get.current.slope.prob<-function(x, calibration=NULL, Twilight.solar.vector=NULL, Twilight.log.light.vector, plot=F, verbose=F,  log.light.borders=log(c(2,64)), log.irrad.borders=c(-9, 1.5), dusk=T, use.intercept=F, return.slopes=F, Twilight.time.vector=NULL,  delta=0, time_correction=NULL, Calib.param=NULL, impute.on.boundaries=T) {
	if (is.null(time_correction) & is.null(Twilight.solar.vector)) stop ("either time_correction or Twilight.solar.vector should be provided to get.current.slope.prob!")
	Probability=0
	Data<-check.boundaries(x, Twilight.solar.vector=Twilight.solar.vector,  Twilight.log.light.vector, plot=F, verbose=verbose,  log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, dusk=dusk, Twilight.time.vector=Twilight.time.vector, impute.on.boundaries=impute.on.boundaries)
	if (dim(Data)[1]>1) {
	LogLight<-Data[,1]
	LogIrrad<-Data[,2]
	
	#===========
	# here I'll check on how many points we have
		if (length(LogLight) >=2) { 
		Model<- lm(LogLight~LogIrrad)
		if (verbose) print(summary(Model))

	get.probs<-function(Model, plot=F, calibration=NULL, time_correction=NULL, Calib.param=NULL) {
		require(fields)
		# check for the intercept
		sum=0
			
		require(mvtnorm)
		coef.coef<-coef(Model)
		slope.sd<-sqrt(vcov(Model)[4])
		
		#if (coef(Model)[1] < calibration$calibration.bayesian.model$Intercept.boundary[1]) { 
		#if (is.na(coef.vcov[1])) coef.vcov[is.na(coef.vcov)]<-calibration$calibration.bayesian.model$Slope.integration$fast.tw.vcov # adding errors from MCMC
		if (is.na(slope.sd)) {
			if (is.null(calibration)) {
			slope.sd=0 
			} else slope.sd<- calibration$Parameters$mean.of.individual.slope.sigma
		}
		#if (length(resid(Model))== 3) coef.vcov<-coef.vcov*15
		# coef.vcov<-coef.vcov*100/(length(resid(Model))^2) # adding errors from MCMC
		if (use.intercept) {
		stop("use of intercept is not implemented!!!")
		} else { 

		#test.Slope<-rnorm(1000, coef.coef[2], sqrt(coef.vcov[4]))
		if (slope.sd==0) {
		test.Slope=coef.coef[2]
		
		} else {
		
	
		test.Slope<-rnorm(1000, coef.coef[2], slope.sd)
		}
		#-----------
		# here is an experimental correction to integrate out 
		# influence of the intercept
		#test.Slope=test.Slope+1.230660 -coef.coef[1]*0.405351 + coef.coef[1]^2*0.032388
		# this should always work is spectrum opacity of the tag does not change
		# end of experimental correction
		#-----------
				#lnorm.Slopes.fun<-function(x, Calib.param) {
				#	x=x[x>0]
					#return(dnorm(log(x), 0.2, 0.415))
				#	return(dnorm(log(x), Calib.param[1], Calib.param[2]))
				#}
				#sum<-sum(lnorm.Slopes.fun(test.Slope, Calib.param))/1000
				#sum<-mean(dlnorm(test.Slope, Calib.param[1], Calib.param[2]))
				#sum<-mean(dlnorm(test.Slope, Calib.param[1]+(Calib.param[2]^2)/2, Calib.param[2]))
		# ==========================
		# change for 5.0
		# addition of the GAM correction..
		# we will measure not Calib.param, but 
		# f(Calib.param[1], Lat, time)
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
		
		if (is.null(delta)) {
		if (length(formals(calibration$lat_correction_fun))==1) {
		delta=calibration$lat_correction_fun(x[2])} else {
		delta=calibration$lat_correction_fun(x[2], Twilight.time.vector[13]) # have to check whether we always have time specified...
		}
		}
		#sum<-mean(dlnorm(test.Slope, Expected.mean+delta, Calib.param[2]))

		# correction for both parameters
		#sum<-mean(dlnorm(test.Slope, Expected.mean+delta, Calib.param[2]))
		#-----------------
		# new exp correction added 24 Apr 2015
		sum<-mean(dlnorm(test.Slope, Expected.mean+delta, Calib.param[2])*test.Slope)
		#cat("time corr", time_correction ,"Expected.mean+delta", Expected.mean+delta, "\n")
		#-----------------------
				# correct for cos Lat
				# looks like cos is too much... should look at sqrt(cos)
				#sum=sum*(cos(x[2]/180*pi))
				sum=sum*(cos(x[2]/180*pi)^0.5)
				#-----------------------
		
				}
		#if (plot) {
		#my.golden.colors <- colorRampPalette(
		#c("white","#FF7100"))
		#image(list(x=calibration$calibration.bayesian.model$Slope.integration$MCMC.All.densities$x, y=calibration$calibration.bayesian.model$Slope.integration$MCMC.All.densities$y,
		#z=matrix(dmvnorm(expand.grid(calibration$calibration.bayesian.model$Slope.integration$MCMC.All.densities$x, calibration$calibration.bayesian.model$Slope.integration$MCMC.All.densities$y), coef.coef, coef.vcov), nrow = length(calibration$calibration.bayesian.model$Slope.integration$MCMC.All.densities$x), ncol = length(calibration$calibration.bayesian.model$Slope.integration$MCMC.All.densities$y))), col=my.golden.colors(10))
		#contour(calibration$calibration.bayesian.model$Slope.integration$MCMC.All.densities, add=T)
		#}
		#}
		return(sum)
		}
		
		Probability<-get.probs(Model, plot=plot, calibration=calibration, time_correction=time_correction, Calib.param=Calib.param)

		if (!is.finite(Probability)) Probability<-0
		if (Probability<0) Probability<-0
		if (return.slopes) 	Probability<-c(Probability, coef(Model)[2], sqrt(vcov(Model)[4]))

		} else {
		if (return.slopes) 	Probability<-c(Probability, NA, NA)

		if (verbose) {
				print(str(calibration, max.level=1))
				cat("calibration$Parameters:\n")
				print(calibration$Parameters)
				}
        }
		} else {Probability<-0} # this is done for testing purposes only...

	return(Probability)
}


get.current.slope.prob<-cmpfun(get.current.slope.prob)
# the new idea is that for the calibration and for the real run we want to have one function that will detect notZero..
# thi function should be able to work under the apply mode and return a dataframe with 2 columns:
# LogLight, LogIrrad, but already filtered from all incorrect points...


# the new idea of the 12 version is to add imputed points from the boundary..

check.boundaries<-function(x, Twilight.solar.vector=NULL,  Twilight.log.light.vector, plot=F, verbose=F,  log.light.borders=log(c(2,64)), log.irrad.borders=c(-15, 50), dusk=T, impute.on.boundaries=T, Twilight.time.vector=NULL) {
# this function...
	if (is.null(Twilight.solar.vector))  {
	Twilight.solar.vector<-solar(as.POSIXct(Twilight.time.vector, tz="gmt", origin="1970-01-01"))
	}
	Elevs<-elevation(x[[1]], x[[2]], Twilight.solar.vector)
	LogIrrad<-log(get.Irradiance(Elevs*pi/180)+1e-20)
	LogLight<-Twilight.log.light.vector
	if (verbose) {
		cat("Elevs:\n")
		print(cbind(ID=1:length(LogLight), LogLight=LogLight, LogIrrad=LogIrrad))
			}
	NotZero<-	which(LogLight>=log.light.borders[1] & LogLight<=log.light.borders[2] & LogIrrad<log.irrad.borders[2])
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
		NotZero<-NotZero[1]:rev(NotZero)[1]
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
		Coef<-coef(lm(Res[,1]~Res[,2]))
		plot(LogLight~LogIrrad, main=paste(twilight=as.POSIXct(Twilight.time.vector[24], tz="gmt", origin="1970-01-01"), ifelse(dusk, "dusk", "dawn"), "intercept", Coef[1], "slope", Coef[2]))
		points(Res[,1]~Res[,2], col="red", lwd=2, pch="+")
		print(Coef)
		par(ask = T)
		plot(LogLight~LogIrrad, main=paste(twilight=as.POSIXct(Twilight.time.vector[24], tz="gmt", origin="1970-01-01"), ifelse(dusk, "dusk", "dawn"), "intercept", Coef[1], "slope", Coef[2]))
		par(ask=F)
		if (Coef[1]<(-6)) warning("check twilight at around ", as.POSIXct(Twilight.time.vector[24], tz="gmt", origin="1970-01-01"), " it had strange shading\n")

	}
	return(Res)
	#}
}

check.boundaries<-cmpfun(check.boundaries)

get.prob.surface<-function(Twilight.ID, dusk=T, Twilight.time.mat, Twilight.log.light.mat, return.slopes=F,  Calib.param, log.irrad.borders=c(-9, 3), delta=0, Points.Land, log.light.borders=log(c(2,64)), calibration=NULL, impute.on.boundaries=T) {
 
		if (Twilight.ID%%10== 1) cat("doing", Twilight.ID, "\n")	
		
		Twilight.solar.vector<-solar(as.POSIXct(Twilight.time.mat[c(1:24, 26:49), Twilight.ID], tz="gmt", origin="1970-01-01"))
		Twilight.log.light.vector<-Twilight.log.light.mat[c(1:24, 26:49), Twilight.ID]
		Twilight.time.vector=Twilight.time.mat[c(1:24, 26:49), Twilight.ID]
		
			if (is.null(calibration)) {
			time_correction=Calib.param[1] # this is the only place where I use calib.param...
			} else {
			#------------------------------
			# here is a change to sun declination from cosSolarDec..
			#time_correction=calibration$time_correction_fun(Twilight.solar.vector$cosSolarDec[1])
			time_correction=calibration$time_correction_fun(get.declination(Twilight.time.vector[24]), as.numeric(dusk))
			}
		
		if (return.slopes) {
		Current.probs<-	apply(Points.Land, 1, get.current.slope.prob, Twilight.log.light.vector=Twilight.log.light.vector, plot=F, verbose=F,  log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, dusk=dusk, return.slopes=T, Twilight.time.vector=Twilight.time.vector,  Calib.param= Calib.param, delta=delta, time_correction=time_correction, calibration=calibration, impute.on.boundaries=impute.on.boundaries)	
			} else {
		Current.probs<-	apply(Points.Land, 1, get.current.slope.prob,  Twilight.log.light.vector=Twilight.log.light.vector, plot=F, verbose=F,  log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, dusk=dusk, return.slopes=F, Twilight.time.vector=Twilight.time.vector,  Calib.param= Calib.param, delta=delta, time_correction=time_correction, calibration=calibration, impute.on.boundaries=impute.on.boundaries)
		}
		return(Current.probs)
	}
	
get.prob.surface<-cmpfun(get.prob.surface)	

get.calib.param<-function(Calib.data.all, plot=F) {
Calib.data.all$fTwilight<-Calib.data.all$fDay
#Calib.data.all<-Calib.data.all[-which(Calib.data.all$fTwilight%in% c(191, 222, 249, 250)),]
# Calib.data.all<-Calib.data.all[-which(Calib.data.all$fTwilight%in% c( 222)),]
#============================
# we should do some outlier tests...

cur.data<-Calib.data.all
p.lm1<-lm(LogLight~fTwilight/LogIrrad,data=cur.data)
cur.slope<-matrix(coef(p.lm1),ncol=2) # foo[,2] - slope
cur.slope<-as.data.frame(cur.slope)
Diag<-diag(vcov(p.lm1))
Diag.slopes<-Diag[which(sapply(strsplit(names(Diag), ":"), length)==2)]

cur.slope$sd[!is.na(cur.slope[,2])]<-sqrt(Diag.slopes)
cur.slope$type<-aggregate(cur.data[,"type"],by=list(Day=cur.data$fTwilight),FUN=function(x) unique(x))[,2]
cur.slope$time<-aggregate(cur.data[,"Time"],by=list(Day=cur.data$fTwilight),FUN=function(x) x[1])[,2]

cur.slope$Duration<-aggregate(cur.data[,"Time"],by=list(Day=cur.data$fTwilight),FUN=function(x) max(x)-min(x))[,2]
names(cur.slope)[2]<-"slope"

#===========================
# ok I've decided to make separate lm
Slopes<-c()
Slopes.sd<-c()
Intercept=c()
Sigma=c()
Type=c()
Time=c()
Elevs<-c()
fDay<-c()
for (i in (unique(cur.data$fTwilight))) {
# lm
#plot(LogLight~LogIrrad, data=cur.data[cur.data$fTwilight==i,])
Data<-cur.data[cur.data$fTwilight==i,]
Data<-Data[Data$LogLight>0,]
Lm<-lm(LogLight~LogIrrad,data=Data)
Slopes<-c(Slopes, coef(Lm)[2])
Slopes.sd<-c(Slopes.sd, sqrt(vcov(Lm)[4]))
Intercept=c(Intercept, coef(Lm)[1])
Sigma<-c(Sigma, summary(Lm)$sigma)
Type=c(Type, cur.data$type[cur.data$fTwilight==i][1])
Time<-c(Time,ifelse(cur.data$type[cur.data$fTwilight==i][1]=="Dusk", max(cur.data$Time[cur.data$fTwilight==i]), min(cur.data$Time[cur.data$fTwilight==i])))
get.calib.param
Elevs=c(Elevs, mean(cur.data$Elevs[cur.data$fTwilight==i], na.rm=T))
fDay=c(fDay, cur.data$fDay[cur.data$fTwilight==i][1])
}
#plot(cur.slope$slope~Slopes)
#plot(cur.slope$sd~Slopes.sd)

# ok, we now want to replace old slopes with new ones.
# but first we want to add 0.5 instead of NA - at least this is how it works now inside the functions and how it was done in a simulations
cur.slope$slope<-Slopes
cur.slope$Intercept<-Intercept
cur.slope$Sigma<-Sigma
#Slopes.sd[is.na(Slopes.sd)]<-sd.fast
cur.slope$sd<-Slopes.sd

names(cur.slope)[2]<-"slope"
#=====================
if (plot) hist(log(cur.slope$slope))

Parameters<-list(Intercept=c(mean(cur.slope$Intercept, na.rm=T), sd(cur.slope$Intercept, na.rm=T)), LogSlope=c(mean(log(cur.slope$slope), na.rm=T), sd(log(cur.slope$slope), na.rm=T)), LogSigma=c(mean(log(cur.slope$Sigma[!is.na(cur.slope$sd)])), sd(log(cur.slope$Sigma[!is.na(cur.slope$sd)]))), mean.of.individual.slope.sigma=mean(cur.slope$sd, na.rm=T))

#cur.slope$time<-aggregate(cur.data[,"Time"],by=list(Day=cur.data$fTwilight),FUN=function(x) x[1])[,2]
#cur.slope$time<-aggregate(cur.data[,"Time"],by=list(Day=cur.data$fTwilight),FUN=mean)[,2]

Res<-list(Parameters=Parameters, Slopes=data.frame(Slope=cur.slope$slope, Time=Time, Intercept=cur.slope$Intercept, Sigma=cur.slope$Sigma, Slopes.sd=Slopes.sd, Type=Type, Elevs=Elevs, fDay=fDay))
return(Res) 
}



get.Phys.Mat.parallel<-function(all.out=NULL, Twilight.time.mat.dusk=NULL, Twilight.log.light.mat.dusk=NULL, Twilight.time.mat.dawn=NULL, Twilight.log.light.mat.dawn=NULL,  threads=2,  calibration=NULL, log.light.borders=NULL, log.irrad.borders=NULL ) {

# let's say we have to submit all boundaries inside tha calibration object..

if (is.character(all.out)) all.out=get("all.out")
if (is.character(Twilight.time.mat.dusk)) Twilight.time.mat.dusk=get("Twilight.time.mat.dusk")
if (is.character(Twilight.time.mat.dawn)) Twilight.time.mat.dawn=get("Twilight.time.mat.dawn")
if (is.character(Twilight.log.light.mat.dusk)) Twilight.log.light.mat.dusk=get("Twilight.log.light.mat.dusk")
if (is.character(Twilight.log.light.mat.dawn)) Twilight.log.light.mat.dawn=get("Twilight.log.light.mat.dawn")
if (is.character(calibration)) calibration=get("calibration")
#if (is.character(log.irrad.borders)) log.irrad.borders=get("log.irrad.borders")
#if (is.character(log.light.borders)) log.light.borders=get("log.light.borders")

#print(ls())
#print(str(Twilight.time.mat.dusk))

Points.Land<-all.out$Points.Land

cat("making cluster\n")
require(parallel)
mycl <- parallel:::makeCluster(Threads)
    tmp<-parallel:::clusterSetRNGStream(mycl)
    tmp<-parallel:::clusterExport(mycl,c("Twilight.time.mat.dawn", "Twilight.time.mat.dusk", "Twilight.log.light.mat.dawn", "Twilight.log.light.mat.dusk", "Points.Land", "calibration"), envir=environment())
    tmp<-parallel:::clusterEvalQ(mycl, library("circular")) 
    tmp<-parallel:::clusterEvalQ(mycl, library("truncnorm")) 
    tmp<-parallel:::clusterEvalQ(mycl, library("GeoLight")) 
    tmp<-parallel:::clusterEvalQ(mycl, library("FLightR")) 
    #tmp<-parallel:::clusterEvalQ(mycl, source("D:\\Geologgers\\LightR_development_code\\functions.Dusk.and.Dawn.5.1.r")) 

#====================
cat("estimating dusks\n")

	Twilight.vector<-1:(dim(Twilight.time.mat.dusk)[2])
	
	 All.probs.dusk<-parSapplyLB(mycl, Twilight.vector, FUN=get.prob.surface, Twilight.log.light.mat=Twilight.log.light.mat.dusk, Twilight.time.mat=Twilight.time.mat.dusk, dusk=T, Calib.param=calibration$Parameters$LogSlope, log.light.borders=calibration$Parameters$log.light.borders, log.irrad.borders=calibration$Parameters$log.irrad.borders, delta=NULL, Points.Land=Points.Land, calibration=calibration)
		 	
cat("estimating dawns\n")
	 
	Twilight.vector<-1:(dim(Twilight.time.mat.dawn)[2])
		 All.probs.dawn<-parSapplyLB(mycl, Twilight.vector, FUN=get.prob.surface, Twilight.log.light.mat=Twilight.log.light.mat.dawn, Twilight.time.mat=Twilight.time.mat.dawn, dusk=F, Calib.param=calibration$Parameters$LogSlope, log.light.borders=calibration$Parameters$log.light.borders, log.irrad.borders=calibration$Parameters$log.irrad.borders, delta=NULL, Points.Land=Points.Land, calibration=calibration)
stopCluster(mycl)

cat("processing results\n")
All.probs.dusk.tmp<-All.probs.dusk
All.probs.dawn.tmp<-All.probs.dawn
	Phys.Mat<-c()
for (i in 1:nrow(all.out$Matrix.Index.Table)) {
	if (all.out$Matrix.Index.Table$Dusk[i]) {
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


get.declination<-function(Dates) {

if (is.numeric(Dates[1])) Dates<-as.POSIXct(Dates, tz="UTC", origin="1970-01-01")

n=as.numeric(Dates-c(as.POSIXct("2000-01-01 12:00:00", tz="UTC")))
L=280.460+0.9856474*n
g=357.528+0.9856003*n
Lambda=(L+1.915*sin(g/180*pi)+0.020*sin(2*g/180*pi))%%360
epsilon = 23.439 - 0.0000004* n 
Dec = asin(sin(epsilon/180*pi)*sin(Lambda/180*pi))*180/pi
return(Dec)
}
