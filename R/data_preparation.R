# data_preparation functions

##############
# process.data
process.twilights<-function(All.p, Filtered_tw, measurement.period=60, saving.period=NULL) {
# this function just prepares data for the next steps..
##########
## Dusk
# processing Dusk
Dusk.all<-Filtered_tw$datetime[Filtered_tw$type==2]
Twilight.index.mat.dusk<-sapply(which(All.p$gmt %in% Dusk.all & All.p$type==2), FUN=function(x) (x-24):(x+24))
Twilight.index.mat.dusk<-apply(Twilight.index.mat.dusk, c(1,2), FUN=function(x) ifelse (x>0, x, NA))
Max.Index<-nrow(All.p)
Twilight.index.mat.dusk<-apply(Twilight.index.mat.dusk, c(1,2), FUN=function(x) ifelse (x>Max.Index, NA, x))

Twilight.time.mat.dusk<-apply(Twilight.index.mat.dusk, c(1,2), FUN=function(x) as.numeric(All.p$gmt[x]))
Twilight.time.mat.dusk<-apply(Twilight.time.mat.dusk, c(1,2), FUN=function(x) ifelse(is.finite(x), x, 0))
Twilight.log.light.mat.dusk<-apply(Twilight.index.mat.dusk, c(1,2), FUN=function(x) log(All.p$light[x]))
#Twilight.log.light.mat.dusk<-apply(Twilight.index.mat.dusk, c(1,2), FUN=function(x) All.p$light[x])
Twilight.log.light.mat.dusk<-apply(Twilight.log.light.mat.dusk, c(1,2), FUN=function(x) ifelse(is.finite(x), x, -1))


Twilight.time.mat.dusk<-Twilight.time.mat.dusk-(saving.period-measurement.period)

# processing Dawn
Dawn.all<-Filtered_tw$datetime[Filtered_tw$type==1]

Twilight.index.mat.dawn<-sapply(which(All.p$gmt %in% Dawn.all & All.p$type==1), FUN=function(x) (x-24):(x+24))
Twilight.index.mat.dawn<-apply(Twilight.index.mat.dawn, c(1,2), FUN=function(x) ifelse (x>0, x, NA))
Max.Index<-nrow(All.p)
Twilight.index.mat.dawn<-apply(Twilight.index.mat.dawn, c(1,2), FUN=function(x) ifelse (x>Max.Index, NA, x))

Twilight.time.mat.dawn<-apply(Twilight.index.mat.dawn, c(1,2), FUN=function(x) as.numeric(All.p$gmt[x]))
Twilight.time.mat.dawn<-apply(Twilight.time.mat.dawn, c(1,2), FUN=function(x) ifelse(is.finite(x), x, 0))
Twilight.log.light.mat.dawn<-apply(Twilight.index.mat.dawn, c(1,2), FUN=function(x) log(All.p$light[x]))
#Twilight.log.light.mat.dawn<-apply(Twilight.index.mat.dawn, c(1,2), FUN=function(x) All.p$light[x])
Twilight.log.light.mat.dawn<-apply(Twilight.log.light.mat.dawn, c(1,2), FUN=function(x) ifelse(is.finite(x), x, -1))

Res<-list(Twilight.time.mat.dusk=Twilight.time.mat.dusk, Twilight.log.light.mat.dusk=Twilight.log.light.mat.dusk, Twilight.time.mat.dawn=Twilight.time.mat.dawn,  Twilight.log.light.mat.dawn=Twilight.log.light.mat.dawn, measurement.period=measurement.period, saving.period=saving.period)
return(Res)
}




get.Irradiance<-function(alpha, r=6378, s=6.9, intigeo.template.correction=F) {
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
	
	if (intigeo.template.correction) Res[Res>0.001]=exp(log(Res[Res>0.001])*1.279+log(Res[Res>0.001])^2*0.091)
	
	return(Res)
}


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
		
	return(Calib.data.all)
	#, Calib.data.all.lme=Calib.data.all.lme.no_interact, Significance.of.dusk.dawn.diff=Significance.of.dusk.dawn.diff, 
	#calibration.slope.par=list(
	#	Calibration.simple.slope=Calibration.simple.no_interact,
	#	tmp.x.slope=tmp.x.RE,
	#	calibration.slope.distr=tmp.dnorm.calib.no_interact.RE,
	#	RE.component=RE.component), 
	#calibration.intercept.par=list(
	#	Parameters=Parameters, CalibrationL.LL=Optim.processed.light$value),
	#calibration.bayesian.model=calibration.bayesian.model # this is a stored list...
	
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
Day<-c()
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
Day=c(Day, cur.data$Day[cur.data$fTwilight==i][1])
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

Res<-list(Parameters=Parameters, Slopes=data.frame(Slope=cur.slope$slope, Time=Time, Intercept=cur.slope$Intercept, Sigma=cur.slope$Sigma, Slopes.sd=Slopes.sd, Type=Type, Elevs=Elevs, Day=Day))
return(Res) 
}


plot.slopes<-function(all.slopes, ylim=NULL) {
#old.par <- par(no.readonly = TRUE) 
#on.exit(par(old.par))
plot(log(all.slopes$Slopes$Slope)~as.POSIXct(all.slopes$Slopes$Time, tz="UTC", origin="1970-01-01"), type="n", main="red - dawn, black - dusk", xlab="time", ylab="log(Slope)", ylim=ylim)
lines(log(Slope)~as.POSIXct(Time, tz="UTC", origin="1970-01-01"), data=all.slopes$Slopes[all.slopes$Slopes$Type=="Dusk",])
points(log(Slope)~as.POSIXct(Time, tz="UTC", origin="1970-01-01"), data=all.slopes$Slopes[all.slopes$Slopes$Type=="Dusk",], pch="+")
points(log(Slope)~as.POSIXct(Time, tz="UTC", origin="1970-01-01"), data=all.slopes$Slopes[all.slopes$Slopes$Type=="Dawn",], pch="+", col="red")
lines(log(Slope)~as.POSIXct(Time, tz="UTC", origin="1970-01-01"), data=all.slopes$Slopes[all.slopes$Slopes$Type=="Dawn",], col="red")
#invisible()
}


get.calibration.parameters<-function(Calibration.periods, Proc.data, model.ageing=T, log.light.borders=log(c(2, 63)),  log.irrad.borders=c(-7,1.5)) {

# ageing.model this option should be used only in case there are two periods of calibration or there is one but very long.
if (nrow(Calibration.periods)==1 & model.ageing) warning("you have only one calibration period ageing estimation is unreliable and should be turned ot F!!!")

Calibration.periods.int<-Calibration.periods
Calibration.periods.int$calibration.start<-as.numeric(Calibration.periods.int$calibration.start)
Calibration.periods.int$calibration.stop<-as.numeric(Calibration.periods.int$calibration.stop)


Dusk.calib.days<-c()
Dawn.calib.days<-c()
Positions=list(dawn=c(), dusk=c())
#sum(apply(Calibration.periods.int, 1, FUN=function(x) x[2]-x[1]))/3600/24
for (i in 1:nrow(Calibration.periods)) {
Dusk.calib.days.cur<-which(Proc.data$Twilight.time.mat.dusk[1,]>Calibration.periods.int$calibration.start[i] & Proc.data$Twilight.time.mat.dusk[1,] < Calibration.periods.int$calibration.stop[i])
Dusk.calib.days<-unique(c(Dusk.calib.days, Dusk.calib.days.cur))

Dawn.calib.days.cur<-which(Proc.data$Twilight.time.mat.dawn[1,]>Calibration.periods.int$calibration.start[i] & Proc.data$Twilight.time.mat.dawn[1,] < Calibration.periods.int$calibration.stop[i] )
Dawn.calib.days<-unique(c(Dawn.calib.days,Dawn.calib.days.cur ))

	dawn.cur=matrix(ncol=2, nrow=length(Dawn.calib.days.cur))
		dawn.cur[,1]<-Calibration.periods.int$lon[i]
		dawn.cur[,2]<-Calibration.periods.int$lat[i]
		#
		dusk.cur=matrix(ncol=2, nrow=length(Dusk.calib.days.cur))
		dusk.cur[,1]<-Calibration.periods.int$lon[i]
		dusk.cur[,2]<-Calibration.periods.int$lat[i]
		Positions$dawn<-rbind(Positions$dawn, dawn.cur)
		Positions$dusk<-rbind(Positions$dusk, dusk.cur)

}
Twilight.time.mat.Calib.dusk<-Proc.data$Twilight.time.mat.dusk[,Dusk.calib.days]
Twilight.log.light.mat.Calib.dusk<-Proc.data$Twilight.log.light.mat.dusk[,Dusk.calib.days]

## Dawn

Twilight.time.mat.Calib.dawn<-Proc.data$Twilight.time.mat.dawn[,Dawn.calib.days]
Twilight.log.light.mat.Calib.dawn<-Proc.data$Twilight.log.light.mat.dawn[,Dawn.calib.days ]

#-------------------#

Calib.data.all<-logger.template.calibration(Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk, positions=Positions, log.light.borders=log.light.borders,  log.irrad.borders=log.irrad.borders, plot=F)

All.slopes<-get.calib.param(Calib.data.all, plot=F)

All.slopes$Slopes$logSlope<-log(All.slopes$Slopes$Slope)

if (model.ageing) {
All.slopes.int<-All.slopes
All.slopes.int$Slopes<-All.slopes.int$Slopes[is.finite(All.slopes.int$Slopes$logSlope),]
Time.start=min(c(Twilight.time.mat.Calib.dusk, Twilight.time.mat.Calib.dawn))
Model=lm(logSlope~I(Time-Time.start), data=All.slopes.int$Slopes)
Model$Time.start<-Time.start
calib_outliers<-All.slopes.int$Slopes$Time[which(abs(residuals(Model))>3*sd(residuals(Model)))]
Res<-list(calib_outliers=calib_outliers, ageing.model=Model, All.slopes=All.slopes)
} else {
calib_outliers<-All.slopes$Slopes$Time[which(abs(log(All.slopes$Slopes$Slope)-mean(log(All.slopes$Slopes$Slope), na.rm=T))>3*sd(log(All.slopes$Slopes$Slope)))]
Res<-list(calib_outliers=calib_outliers, All.slopes=All.slopes)
}
return(Res)
}


create.calibration<-function( All.slopes, Proc.data, FLightR.data, log.light.borders, log.irrad.borders, start, ageing.model=NULL) {
# Now we create 'parameters' object that will have all the details about the calibration
Parameters<-All.slopes$Parameters # LogSlope # 
Parameters$measurement.period<-Proc.data$measurement.period 
Parameters$saving.period<-Proc.data$saving.period 
Parameters$log.light.borders<-log.light.borders # these are the boundaries in which one should use the BAS tag.. for the other types of tags they will be different.
Parameters$min.max.values<-c(min(FLightR.data$Data$light), max(FLightR.data$Data$light))
Parameters$log.irrad.borders=log.irrad.borders
Parameters$start=start

Parameters$LogSlope_1_minute<-Parameters$LogSlope
 if (is.null(ageing.model)) {
time_correction_fun= eval(parse(text=paste("function (x,y) return(", Parameters$LogSlope[1], ")")))
lat_correction_fun<-function(x, y, z) return(0)
} else {
time_correction_fun= function(x, y) return(0)
lat_correction_fun<-eval(parse(text=paste("function (x,y) return(", Time.start, "+", coef(ageing.model)[1], "+", coef(ageing.model)[2], "* y)")))
}

Calibration<-list(Parameters=Parameters, time_correction_fun=time_correction_fun, lat_correction_fun=lat_correction_fun)

return(Calibration)
}


correct.hours<-function(datetime) {
	# this function is supposed to correct hours by adding some amount of them.
   hours <- as.numeric(format(datetime,"%H"))+as.numeric(format(datetime,"%M"))/60
   hours[as.POSIXlt(datetime)$isdst==1] <-hours[as.POSIXlt(datetime)$isdst==1]-1
	cor <- rep(NA, 24)
	for(i in 0:23){
		cor[i+1] <- max(abs((c(hours[1],hours)+i)%%24 - 
		            (c(hours,hours[length(hours)])+i)%%24),na.rm=T)
	}
	hours <- (hours + (which.min(round(cor,2)))-1)%%24
	return(hours)
}

make.result.list<-function(Data, raw.X, raw.Y) {
	Int<-c(min(raw.Y), max(raw.Y))
	Index<-sapply(raw.X,  FUN=function(x) which.min(abs(as.numeric(Data$d$gmt-x))))
	Res<-Data$d[Index,]
	# now I want to adjust time and light
	Res$light<-approx(x=Data$d$gmt, y=Data$d$light, xout=raw.X)$y
	Res$gmt<-raw.X
	Res$Hour<-raw.Y
	Result<-list(Data=Res)
	Result$Data$gmt.adj<-Result$Data$gmt
	Result$Data$gmt<-as.POSIXct(Result$Data$gmt, tz="UTC", origin="1970-01-01")
	Result$Data$gmt.adj<-as.POSIXct(Result$Data$gmt.adj, tz="UTC", origin="1970-01-01")
	return(Result)
}

make.processed.light.object<-function(FLightR.data) {
# dusk
raw.Y.dusk<-correct.hours(FLightR.data$twilights$datetime[FLightR.data$twilights$type==2 & FLightR.data$twilights$excluded==0])
raw.X.dusk<-as.numeric(FLightR.data$twilights$datetime[FLightR.data$twilights$type==2 & FLightR.data$twilights$excluded==0])
Data_tmp<-list(d=FLightR.data$Data)
Result.Dusk<-make.result.list(Data_tmp, raw.X.dusk, raw.Y.dusk)

# dusk
raw.Y.dawn<-correct.hours(FLightR.data$twilights$datetime[FLightR.data$twilights$type==1 & FLightR.data$twilights$excluded==0])
raw.X.dawn<-as.numeric(FLightR.data$twilights$datetime[FLightR.data$twilights$type==1 & FLightR.data$twilights$excluded==0])
Result.Dawn<-make.result.list(Data_tmp, raw.X.dawn, raw.Y.dawn)

# combine
processed.light<-list(Final.dusk=Result.Dusk, Final.dawn=Result.Dawn)

return(processed.light)
}

# this create proposal is corrected for the missing loess
create.proposal<-function(processed.light, start=c(-98.7, 34.7), end=NA, Grid) {
	if (is.na(end)) end=start
All.Days<-seq(min(min(processed.light$Final.dusk$Data$gmt), min(processed.light$Final.dawn$Data$gmt)),max(max(processed.light$Final.dusk$Data$gmt), max(processed.light$Final.dawn$Data$gmt)), by="days")

## Probabilty.of.migration
All.Days.extended<-seq(min(All.Days)-86400, max(All.Days)+86400, by="days")

# these are all potential twilights that we would have for the start point...
Potential.twilights<-sort(c(sunriset(matrix(start, nrow=1), All.Days.extended, direction="sunrise", POSIXct.out=TRUE)[,2], sunriset(matrix(start, nrow=1), All.Days.extended, direction="sunset", POSIXct.out=TRUE)[,2]))
# now we need to add dask or dawn to it..
Sunrise<-sunriset(matrix(start, nrow=1), Potential.twilights[1], direction="sunrise", POSIXct.out=TRUE)[,2]
Sunrises<-c(Sunrise-(3600*24), Sunrise, Sunrise+(3600*24))
Sunset<-sunriset(matrix(start, nrow=1), Potential.twilights[1], direction="sunset", POSIXct.out=TRUE)[,2]
Sunsets<-c(Sunset-(3600*24), Sunset, Sunset+(3600*24))
First.twilight<-ifelse(which.min(c(min(abs(difftime(Sunrises,Potential.twilights[1], units="mins"))), min(abs(difftime(Sunsets,Potential.twilights[1], units="mins")))))==1, "dawn", "dusk")

Index.tab<-data.frame(Date=Potential.twilights)

if (First.twilight=="dusk") {
	Index.tab$Dusk<-rep(c(T,F), times=length(All.Days.extended))
} else {
		Index.tab$Dusk<-rep(c(F,T), times=length(All.Days.extended))
}

Index.tab$Curr.mat<-NA # this will be NA if no data and row number if there are..
Index.tab$Real.time<-as.POSIXct(NA, tz="GMT")
Index.tab$time<-as.POSIXct(NA, tz="GMT")

######################################
# !!!! this will not work if bird will move for over 12 time zones!!!
# Dusk
for (i in 1:length(processed.light$Final.dusk$Data$gmt)) {
	Row2write<-which.min(abs(difftime(processed.light$Final.dusk$Data$gmt[i], Index.tab$Date[Index.tab$Dusk==T], units="mins")))
	Index.tab$time[Index.tab$Dusk==T][Row2write]<-processed.light$Final.dusk$Data$gmt[i]
	Index.tab$Real.time[Index.tab$Dusk==T][Row2write]<-processed.light$Final.dusk$Data$gmt.adj[i]
	Index.tab$Curr.mat[Index.tab$Dusk==T][Row2write]<-i
}

# Dawn
for (i in 1:length(processed.light$Final.dawn$Data$gmt)) {
	Row2write<-which.min(abs(difftime(processed.light$Final.dawn$Data$gmt[i], Index.tab$Date[Index.tab$Dusk==F], units="mins")))
	Index.tab$time[Index.tab$Dusk==F][Row2write]<-processed.light$Final.dawn$Data$gmt[i]
	Index.tab$Real.time[Index.tab$Dusk==F][Row2write]<-processed.light$Final.dawn$Data$gmt.adj[i]
	Index.tab$Curr.mat[Index.tab$Dusk==F][Row2write]<-i
}
# cutting empty ends..
while (is.na(Index.tab$Curr.mat[1])) Index.tab<-Index.tab[-1,]
while (is.na(Index.tab$Curr.mat[nrow(Index.tab)])) Index.tab<-Index.tab[-nrow(Index.tab),]
Index.tab$Point<-NA
First.Point<-which.min(spDistsN1(Grid[,1:2], start,  longlat=T))
Index.tab$Point[1]<-First.Point
#
# I decided that Curr.mat is not needed anymore
#Index.tab$Main.Index[,-which(names(Index.tab$Main.Index)=="Curr.mat")]
# this need to double checked
Index.tab$yday<-as.POSIXlt(Index.tab$Date, tz="GMT")$yday

 return(Index.tab)
}



geologger.sampler.create.arrays<-function(Index.tab, Grid, start, stop=start) {
	# this function wil take in Index.table with all the proposals and will create an array for case of missed data..
	
	# the main feature is that it can account for missing data..
	Index.tab.old<-Index.tab
	Index.tab$proposal.index<-NA
	Index.tab$proposal.index[Index.tab$Dusk==T]<-"Dusk"
	Index.tab$proposal.index[Index.tab$Dusk==F]<-"Dawn"
	
	Missed.twilights<-which(is.na(Index.tab$Curr.mat))
	
	if (length(Missed.twilights)>0) {
	Deleted.Rows.Count<-0
		for (i in Missed.twilights) {
			i=i-Deleted.Rows.Count
			Index.tab$proposal.index[i-1]<-"Comb"
			Index.tab$Decision[i-1]<-(Index.tab$Decision[i-1]+Index.tab$Decision[i])/2 # 
			if (Index.tab$Direction[i-1] != Index.tab$Direction[i])	{
				Index.tab$Direction[i-1]<-0
				Index.tab$Kappa[i-1]<-0
			}
			Index.tab$M.mean[i-1]<-Index.tab$M.mean[i-1]+Index.tab$M.mean[i]
			Index.tab$M.sd[i-1]<-sqrt((Index.tab$M.sd[i-1])^2+(Index.tab$M.sd[i])^2)
			Index.tab<-Index.tab[-i,]
			Deleted.Rows.Count=Deleted.Rows.Count+1
		}
	}
	output<-list()
	Matrix.Index.Table<-	Index.tab[,which(names(Index.tab) %in% c("Decision", "Direction", "Kappa", "M.mean", "M.sd", "yday", "Real.time", "time", "Dusk", "Loess.se.fit", "Loess.n")) ]
	# I  also remove last line as bird was not flying after last twilight
	Matrix.Index.Table<-Matrix.Index.Table[-nrow(Matrix.Index.Table),]

	
	# main index will have the same amount of rows we have in twilight matrices without first.. 
	Main.Index<-c()
	Main.Index$Biol.Prev<-1:nrow(Matrix.Index.Table)
	Main.Index$proposal.index<-Index.tab$proposal.index[-nrow(Index.tab)]
	Main.Index<-as.data.frame(Main.Index)
	
	Indices<-list(Matrix.Index.Table=Matrix.Index.Table, Main.Index=Main.Index)
	
	output$Indices=Indices
	
	# ok now we need to add to output all other parts..
	
	output$Spatial<-list(Grid=Grid) #[,1:2]

	output$Spatial$Behav.mask<-as.integer(Grid[,3])

	output$Spatial$start.point<-which.min(spDistsN1(Grid[,1:2], start,  longlat=T))

	if (!is.na(stop[1]) ) {
	output$Spatial$stop.point<-which.min(spDistsN1(Grid[,1:2], stop,  longlat=T))
	} else {
	output$Spatial$stop.point<-NA
	}
	
	#output$Spatial$tmp<-list(Distance=spDists(Grid[,1:2], longlat=T))

	#get.angles<-function(all.arrays.object) {
	#	return(apply(all.arrays.object$Spatial$Grid, 1, FUN=function(x) as.integer(round(gzAzimuth(from=all.arrays.object$Spatial$Grid, to=x)))))
	#}
	#output$Spatial$tmp$Azimuths<-get.angles(output)

	return(output)
	}
