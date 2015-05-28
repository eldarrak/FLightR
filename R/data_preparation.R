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
	#	Parameters=Parameters, CalibrationL.LL=Optim.result$value),
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


plot.slopes<-function(all.slopes) {
#old.par <- par(no.readonly = TRUE) 
#on.exit(par(old.par))
plot(log(all.slopes$Slopes$Slope)~as.POSIXct(all.slopes$Slopes$Time, tz="UTC", origin="1970-01-01"), type="n", main="red - dawn, black - dusk", xlab="time", ylab="log(Slope)")
lines(log(Slope)~as.POSIXct(Time, tz="UTC", origin="1970-01-01"), data=all.slopes$Slopes[all.slopes$Slopes$Type=="Dusk",])
points(log(Slope)~as.POSIXct(Time, tz="UTC", origin="1970-01-01"), data=all.slopes$Slopes[all.slopes$Slopes$Type=="Dusk",], pch="+")
points(log(Slope)~as.POSIXct(Time, tz="UTC", origin="1970-01-01"), data=all.slopes$Slopes[all.slopes$Slopes$Type=="Dawn",], pch="+", col="red")
lines(log(Slope)~as.POSIXct(Time, tz="UTC", origin="1970-01-01"), data=all.slopes$Slopes[all.slopes$Slopes$Type=="Dawn",], col="red")
#invisible()
}