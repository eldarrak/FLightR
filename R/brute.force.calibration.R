# brute.force.calibration.R
# ver 0.1 from 21.10.2014


get.1.minute.parameters<-function(parameters, measurement.period=60, saving.period=NULL, position, start.time, end.time, log.light.borders=c(log(2,64)), repeats=50, min.max.values=c(0,64)) {
# this stupid brute force function will just etimate what should be the parameter values on a 1 minute scale..

#========================================
# parameters used so far
# saving.period
# start
# min(Time)
# max(Time)
#log.light.borders
#========================================
Time.seq<-seq(from=start.time-7200
, to=end.time+7200, by=measurement.period)

# here the fixing period should be used.
#saving.period=600
Time.seq.saving<-seq(from=start.time-7200
, to=end.time+7200, by=saving.period)
# the complication here is that it looks like variance is changing now while switch time resolution..
# But Formally I would do this in two independent tasks..
# first I would optimize variance..

##===============================
## so here I simulate and solve the equations regarding SD
## this is not very much needed but we should probably go for that..
## the only problem I see here is that we will change the LogSlope, so the sd of
## the logSlope will also change
## how could we solve that?
## ok, let's for now leave it as -0.1 +0.4

To.run<-expand.grid(Slope.ideal= runif(1000,   parameters$LogSlope[1]-0.2,   parameters$LogSlope[1]+0.4), SD.ideal=runif(1000,  exp(log( parameters$LogSlope[2])-0.2),  exp(log( parameters$LogSlope[2])+0.2))) #
Res<-c()
Res<-as.data.frame(Res)
for (i in 1:repeats) {
To.run.cur<-To.run[sample( 1:nrow(To.run), 1),]
All.slope.runs=get.slopes(To.run=To.run.cur, Parameters=parameters, Lat=position[2], measurement.period=measurement.period, saving.period=saving.period, Time.seq=Time.seq, Time.seq.saving=Time.seq.saving, Lon=position[1], log.light.borders=log.light.borders, min.max.values=min.max.values)
Res<-rbind(Res, c(mean(All.slope.runs$Slope, na.rm=T), All.slope.runs$Slope.ideal[1], sd(All.slope.runs$Slope, na.rm=T), All.slope.runs$SD.ideal[1]))
names(Res)<-c("Slope", "Slope.ideal", "SD", "SD.ideal")
print(Res)
}

#Parameters=All.slopes$Parameters


#plot(Res$SD~Res$Slope.ideal) # no relationship
#plot(Res$SD~Res$Slope) # no relationship
# ok, this means that slope and SD are independent. GOOD.
#summary(lm(Res$Slope~Res$Slope.ideal))
#=ok I do see now absolute linear relationship
# with the same degree for slope and for SD
# interesting, but I would not pay too much attention to that

#plot(Res$SD~Res$SD.ideal)
Lm.sd<-lm(Res$SD~Res$SD.ideal)
summary(Lm.sd)
Slope.sd.ideal<-(parameters$LogSlope[2]-coef(Lm.sd)[1])/coef(Lm.sd)[2] #0.42 comparing to 0.4 from original calibration


plot(Res$Slope~Res$Slope.ideal)
Lm.slope<-lm(Res$Slope~Res$Slope.ideal)
summary(Lm.slope)
Slope.ideal<-(parameters$LogSlope[1]-coef(Lm.slope)[1])/coef(Lm.slope)[2] # 

# ok, now we have approximated slope also...
# but let's repeat that in a more detailed way

To.run<-expand.grid(Slope.ideal=runif(1000, Slope.ideal-0.3, Slope.ideal+0.3), SD.ideal=Slope.sd.ideal) #

All.slope.runs1=get.slopes(Repeats=repeats, To.run=To.run, Parameters=parameters, Lat=position[2], measurement.period=measurement.period, saving.period=saving.period, Time.seq=Time.seq, Time.seq.saving=Time.seq.saving, Lon=position[1], log.light.borders=log.light.borders, min.max.values=min.max.values)
#mean(All.slope.runs1$Slope) # works
#sd(All.slope.runs1$Slope)

# and now let's repeat the same story with Lm...
print(All.slope.runs1)
All.slope.runs1<-na.omit(All.slope.runs1)
Lm1<-lm(Slope~Slope.ideal, data=All.slope.runs1)
#summary(Lm1)

# this is stupid as what we do want to see are correction lms if they are independent then we should just use them to predict what is the 1 minute distribution..

Slope.parameter.final=(parameters$LogSlope[1]-coef(Lm1)[1])/coef(Lm1)[2] # 0.80 #

# Just realized that what we need to send our are the regression coeddicients.
RES<-rbind(coef(Lm1), coef(Lm.sd))
colnames(RES)<-c("a", "b")
rownames(RES)<-c("mu.slope.1.minute", "SD.slope.1.minute")
return(RES)
}

#
get.time.correction.function<-function(parameters, measurement.period=60, saving.period=NULL, position, log.light.borders=log(c(2,64)), Repeats=2, min.max.values=c(0,64), log.irrad.borders=c(-15, 50)) {
# as far as we do not provide time the slope function will runf for the whole year.

To.run<-expand.grid(Slope.ideal=Parameters$LogSlope_1_minute[1], SD.ideal=Parameters$LogSlope_1_minute[2]) #
All.slope.runs=get.slopes(Repeats=Repeats, To.run=To.run, Parameters=parameters, Lat=position[2], measurement.period=measurement.period, saving.period=saving.period, short.run=F, Lon=position[1], log.light.borders=log.light.borders, min.max.values=min.max.values, log.irrad.borders=log.irrad.borders)

Solar<-solar(as.POSIXct(All.slope.runs$gmt, tz="UTC", origin="1970-01-01"))
All.slope.runs$cosSolarDec<-Solar$cosSolarDec
#save(All.slope.runs, file="All.slope.runs_time_correction_600.RData")

# ok cos is linarly related!!!
# so let's make it linear
plot(Slope~cosSolarDec, data=All.slope.runs)
Lm2<-lm(Slope~cosSolarDec, data=All.slope.runs)
#summary(Lm2)
# ok, so this will be used as a time related pattern
# I think we should make a function of that and add it to the next step...

time_correction_fun<-approxfun(x=c(0.9, 1), y=predict(Lm2, newdata=data.frame(cosSolarDec=c(0.9, 1))))

Res<-list(time_correction_fun=time_correction_fun, coef=coef(Lm2), results=All.slope.runs)
return(Res)
}
#

get.lat.correction.function<-function(deltalim=c(-0.05, 0.35), measurement.period=60, saving.period=NULL, threads=2, mode="smart", Sigmas, LogSlope,  log.light.borders=log(c(2, 64)),log.irrad.borders=c(-50, 50), calibration=NULL, min.max.values=c(0,64)) {

cat("function will work in ", mode, "mode\n")

if (!mode %in% c("trial", "brute", "smart")) stop ("mode could be one of - trial, brute, smart")
require(parallel)

if (mode=="trial") {
	deltalim=c(0.15, 0.16)
	cat("deltalim set to ", deltalim, "\n")
	Res<-get.deltas.parallel(deltalim=deltalim, limits=c(-65,65), points=2, Sigmas=sigma, measurement.period=measurement.period, saving.period=saving.period, short.run=T, threads=threads, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders,  random.delta=T, calibration=calibration, min.max.values=min.max.values)
	Res<-as.data.frame(Res)
	print(str(Res))
	names(Res)<-c("Diff", "Sigma", "Delta", "Lat", "Diff.first", "Diff.second", "Sigma.init")
	Res$cosLat<-cos(Res$Lat/180*pi)
	RES<-list(simulations=Res)
}
if (mode=="brute") {
# real run
Res<-get.deltas.parallel(deltalim=deltalim, limits=c(-65,65), points=70, Sigmas=sigma, measurement.period=measurement.period, saving.period=saving.period, short.run=T, threads=threads, log.irrad.borders=log.irrad.borders, random.delta=T, calibration=calibration, min.max.values=min.max.values)
#=======
# !!! idealy there should be a check that will make sure that active variation occured at the deltalim specified
#=======
Res<-as.data.frame(Res)
names(Res)<-c("Diff", "Sigma", "Delta", "Lat", "Diff.first", "Diff.second", "Sigma.init")
Res$cosLat<-cos(Res$Lat/180*pi)
require(mgcv)
Model.all=try(gam(Delta~te(Diff, cosLat),  data=Res)) # this look the best
if (class(Model.all)=="try-error") {
RES<-list(simulations=Res)	
} else {
print(summary(Model.all))
#plot(predict(Model.all, newdata=data.frame(Lat=-90:90, Diff=0))~c(-90:90))
#lines(predict(Model.all, newdata=data.frame(Lat=-90:90, Diff=0))~c(-90:90))
#lines(predict(Model.all, newdata=data.frame(cosLat=cos(c(0:90)/180*pi), Diff=0)))

# look slike this si the simplest function something quadratic..

#plot(residuals(Model.all.2)~Res$cosLat)
#Res.limited<-Res[Res$Diff<2 & Res$Diff>-1.5,]

# let's make an approsfun now...
lat_correction_fun<-approxfun(y=predict(Model.all, newdata=data.frame(cosLat=cos(c(-90:90)/180*pi), Diff=0)), x=-90:90)

RES<-list(simulations=Res, lat_correction_fun=lat_correction_fun)
}
}
if (mode=="smart") {
 # in this mode we will try to make several attempts
 # so the first one is to make a very imprecise approximation
 # and it is better to make it at the maximum and minimum latitudes...
 # so we need to try at abs(min) and abs(max)
 deltalim_initial<-deltalim
 cat("...estimating compensation at latitude 0")


	Points<-ifelse(threads<7, 7, threads) # how many repeats to run..

	Res_min_lat<-get.deltas.parallel(deltalim=deltalim_initial, limits=c(0,0), points=Points, Sigmas=sigma, measurement.period=measurement.period, saving.period=saving.period, short.run=T, threads=threads, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, random.delta=T, calibration=calibration, min.max.values=min.max.values)
	Res_min_lat<-as.data.frame(Res_min_lat)
	names(Res_min_lat)<-c("Diff", "Sigma", "Delta", "Lat", "Diff.first", "Diff.second", "Sigma.init")
	Res_min_lat$cosLat<-cos(Res_min_lat$Lat/180*pi)
	# we do not need model here as the situation is very simples
	
	Res_min_lat<-Res_min_lat[Res_min_lat$Diff>-8 & Res_min_lat$Diff<8,]

	#plot(Delta~Diff, data=Res_min_lat)
	require(mgcv)
	Gam_min<-gam(Delta~s(Diff, k=3), data=Res_min_lat)
	Predict_min<-predict(Gam_min, se.fit=T, newdata=data.frame(Diff=0))
	
	
	cat("    Done!")	

	cat("estimated delta for 0 degrees is" , round(Predict_min$fit,3),  "+-"  , round(Predict_min$se.fit,3), "\n")
	plot(Delta~Diff, data=Res_min_lat, col="red")
	abline(v=0, col="red")
	if ((Predict_min$fit-3*Predict_min$se.fit) < deltalim_initial[1]) stop ("try to correct lower deltalim boundary to a smaller value\n")

	
 cat("...estimating compensation at latitude 65")
	Res_max_lat<-get.deltas.parallel(deltalim=deltalim_initial, limits=c(65,65), points=Points, Sigmas=sigma, measurement.period=measurement.period, saving.period=saving.period, short.run=T, threads=threads, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, random.delta=T, calibration=calibration, min.max.values=min.max.values)
	Res_max_lat<-as.data.frame(Res_max_lat)
	names(Res_max_lat)<-c("Diff", "Sigma", "Delta", "Lat", "Diff.first", "Diff.second", "Sigma.init")
	Res_max_lat<-Res_max_lat[Res_max_lat$Diff>-8 & Res_max_lat$Diff<8,]
	
	Res_max_lat$cosLat<-cos(Res_max_lat$Lat/180*pi)
	cat("    Done!")	

	Gam_max<-gam(Delta~s(Diff, k=3), data=Res_max_lat)
	Predict_max<-predict(Gam_max, se.fit=T, newdata=data.frame(Diff=0))

	
	cat("estimated delta for 65 degrees N is" , round(Predict_max$fit,3),  "+-"  , round(Predict_max$se.fit,3), "\n")
	
	plot(Delta~Diff, data=Res_max_lat)
	points(Delta~Diff, data=Res_min_lat, col="red")
	abline(v=0, col="red")

	if ((Predict_max$fit+3*Predict_max$se.fit) > deltalim_initial[2]) stop ("try to correct upper  deltalim boundary to a higher value\n")
	
 	plot(Delta~Diff, data=Res_max_lat)
	points(Delta~Diff, data=Res_min_lat, col="red")
	Gam_max<-gam(Delta~s(Diff, k=3), data=Res_max_lat)
	Predict_max<-predict(Gam_max, se.fit=T, newdata=data.frame(Diff=0))
	
	# ok and now we want at least  to take a diap from min to max and focus there..
	
	deltalim_corrected<-c(Predict_min$fit-0.1, 	Predict_max$fit+0.1)
	
	
	##################
	# checking for -65
	#Res_minus_max_lat<-get.deltas.parallel(deltalim=deltalim_initial, limits=c(-65,-65), points=Points, Sigmas=sigma, interval=saving.period, short.run=T, threads=threads, log.irrad.borders=c(-50, 50), random.delta=T, calibration=calibration)
	
	# ok now we have to go for the full run in a new boundaries...
	
	cat("doing the full run \n")
	
	Points<-(50%/%threads)*threads # how many repeats to run..
	
	Res<-get.deltas.parallel(deltalim=deltalim_corrected, limits=c(-65,65), points=Points, Sigmas=sigma, measurement.period=measurement.period, saving.period=saving.period, short.run=T, threads=threads, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, random.delta=T, calibration=calibration, min.max.values=min.max.values)

	cat(" ... done!")
	Res<-as.data.frame(Res)
	names(Res)<-c("Diff", "Sigma", "Delta", "Lat", "Diff.first", "Diff.second", "Sigma.init")
	
	Res<-Res[Res$Diff>-8 & Res$Diff<8,]

	
	Res$cosLat<-cos(Res$Lat/180*pi)
	
	#=== 
	# now we have to combine the three
	
	Res.all<-rbind(Res_min_lat,Res_max_lat, Res)
	
	require(mgcv)
	Model.all=try(gam(Delta~te(Diff, cosLat),  data=Res.all)) # this look the best
	if (class(Model.all)=="try-error") {
	RES<-list(simulations=Res)	
	} else {
	print(summary(Model.all))
	lat_correction_fun<-approxfun(y=predict(Model.all, newdata=data.frame(cosLat=cos(c(-90:90)/180*pi), Diff=0)), x=-90:90)

	RES<-list(simulations=Res, lat_correction_fun=lat_correction_fun)
	
 }

}
# ok this is it so far
return(RES)
}

#==================================
# here are the new functions for the latitude correction estimation...
# they are not wrapped yet though..


simulate.and.prepare.track<-function(measurement.period=60, saving.period=600, Parameters=Parameters, short.run=T, min.max.values=c(0, 64),log.light.borders=log(c(2,64)),Lat, first.date="2010-01-01 00:00:00") {

# for this we will have to create a To.run.object...
To.run<-data.frame(Slope.ideal=Parameters$LogSlope_1_minute[1], SD.ideal=Parameters$LogSlope_1_minute[2], Latitude=Lat) #

Track<-simulate.track(measurement.period=measurement.period, saving.period=saving.period, To.run=To.run, Parameters=Parameters, min.max.values=min.max.values, short.run=short.run, Time.seq=NULL, Time.seq.saving=NULL, Lon=0, first.date=first.date)
# ok, now we need to process it..
# to prepare for the run..

#==================================
# saving and reading track file
cat("   saving file\n")
lig.data<-cbind(format(as.POSIXct(Track$gmt, tz="gmt",origin="1970-01-01"), format="%d/%m/%Y"), format(as.POSIXct(Track$gmt, tz="gmt",origin="1970-01-01"), format="%H:%M:%S"), round(exp(Track$LogLight)))

File.name<-tempfile(pattern = "sim.no.move.", tmpdir = getwd(), fileext = ".csv")
write.table(lig.data, file=File.name, sep = ",", dec = ".", qmethod="double", quote = FALSE,row.names = FALSE, col.names=FALSE)

## no I wnat to work with the simulated data.

require(FLightR)
require(GeoLight)
require(maptools)
#data(wrld_simpl)
require(fields)
cat("   reading file\n")
# set wd
Data<-geologger.read.data(file=File.name)
try(unlink(File.name))

#========================================================
#========================================================

# new trick - let's try to load the real track
cat("   GeoLight step\n")

tw <- twilightCalc(Track$gmt, Track$light, allTwilights=T, ask=F, LightThreshold=log.light.borders[1]+1)
GLtab   <- tw[[2]] # Table to proceed with GeoLight
cat("automatically detected twilight:\n")
print(table(GLtab[,3]))

# as far as we have done an automatic detection we should check on whether the is no problems..
# the idea is that each dusk should be follwoed by dawn...
Index<-1:length(GLtab[,3])
GLtab[Index[Index%%2==1],3]<-round(mean(GLtab[Index[Index%%2==1],3]))
GLtab[Index[Index%%2==0],3]<-round(mean(GLtab[Index[Index%%2==0],3]))
# 
GLtab1<-GLtab
Filtered_tw <- data.frame(datetime=as.POSIXct(c(GLtab1$tFirst,GLtab1$tSecond),"UTC"),type=c(GLtab1$type,ifelse(GLtab1$type==1,2,1)))
Filtered_tw <- Filtered_tw[!duplicated(Filtered_tw$datetime),]
Filtered_tw <- Filtered_tw[order(Filtered_tw[,1]),]

# now I want to pair data and twilights..		  
Filtered_tw$light<-approx(x=Track$gmt, y=Track$light, xout=Filtered_tw$datetime)$y
Filtered_tw$id<-0
Track$type.real<-Track$type
Track$type<-0

Track.new<-Track[names(Track) %in% c("gmt", "light", "type")]
Track.new<-rbind(Track.new, data.frame(gmt= Filtered_tw$datetime, light=Filtered_tw$light, type=Filtered_tw$type))
All.p<-Track.new[order(Track.new$gmt),]
rownames(All.p)<-1:nrow(All.p)

raw.Y.dusk<-correct.hours(Filtered_tw$datetime[Filtered_tw$type==2])
raw.X.dusk<-as.numeric(Filtered_tw$datetime[Filtered_tw$type==2])
Result.Dusk<-make.result.list(Data, raw.X.dusk, raw.Y.dusk)

raw.Y.dawn<-correct.hours(Filtered_tw$datetime[Filtered_tw$type==1])
raw.X.dawn<-as.numeric(Filtered_tw$datetime[Filtered_tw$type==1])
Result.Dawn<-make.result.list(Data, raw.X.dawn, raw.Y.dawn)

Result.all<-list(Final.dusk=Result.Dusk, Final.dawn=Result.Dawn)
####

##START POINTS###
#all.out<-geologger.sampler.create.arrays(Index.tab, Points.Land, start=start)

#================= END== =====================================================

#=============================================================================

Proc.data<-process.twilights(All.p, Filtered_tw, measurement.period=measurement.period, saving.period=saving.period)

## Dusk
Twilight.time.mat.dusk<-Proc.data$Twilight.time.mat.dusk
Twilight.log.light.mat.dusk<-Proc.data$Twilight.log.light.mat.dusk

## Dawn
Twilight.time.mat.dawn<-Proc.data$Twilight.time.mat.dawn
Twilight.log.light.mat.dawn<-Proc.data$Twilight.log.light.mat.dawn


	Twilight.vector.dusk<-1:(dim(Twilight.time.mat.dusk)[2])
 	Twilight.vector.dawn<-1:(dim(Twilight.time.mat.dawn)[2])

##########################
# and now I want to save output..
	Res<-list(Dusk=list(Twilight.vector=Twilight.vector.dusk, Twilight.log.light.mat=Twilight.log.light.mat.dusk, Twilight.time.mat=Twilight.time.mat.dusk), Dawn=list(Twilight.vector=Twilight.vector.dawn, Twilight.log.light.mat=Twilight.log.light.mat.dawn, Twilight.time.mat=Twilight.time.mat.dawn), Lat=Lat, Filtered_tw=Filtered_tw)
	
	return(Res)
}


get.grid.of.tracks<-function(measurement.period=60, saving.period=600, Parameters=Parameters, short.run=T, min.max.values=c(0, 64), log.light.borders=log(c(2,64)),   Lats, cluster=NULL, first.date="2010-01-01 00:00:00") {
if (!is.null(cluster)) {
	require(parallel)
	#if (threads==-1) threads= detectCores()-1 # selecting all -1 available cores...
	#cluster<-makeCluster(threads)
		tmp<-parallel:::clusterSetRNGStream(cluster)
		### we don' need to send all parameters to node. so keep it easy..
		tmp<-parallel:::clusterExport(cluster, c( "simulate.and.prepare.track", "Parameters", "measurement.period","saving.period", "short.run", "min.max.values", "log.light.borders", "first.date", "simulate.track"),  envir=environment())
		tmp<-parallel:::clusterEvalQ(cluster, library("FLightR"))
		tmp<-parallel:::clusterEvalQ(cluster, library("GeoLight")) 
		tmp<-parallel:::clusterEvalQ(cluster, library("maptools")) 
		Tracks<-parLapply(cluster, Lats,fun=function(x) simulate.and.prepare.track(measurement.period=measurement.period, saving.period=saving.period,  Parameters=Parameters, min.max.values=min.max.values, short.run=short.run, log.light.borders=log.light.borders, Lat=x, first.date=first.date))
	#stopCluster(mycl)
} else {
Tracks<-lapply(Lats,FUN=function(x) simulate.and.prepare.track(measurement.period=measurement.period, saving.period=saving.period,  Parameters=Parameters, min.max.values=min.max.values, short.run=short.run, log.light.borders=log.light.borders, Lat=x, first.date=first.date))
}
return(Tracks)
}


get.diff<-function(prepared.data, calibration, min.max.values=c(1, 1150), log.light.borders=log(c(2, 1100)), log.irrad.borders=c(-15, 50)) {

	Points.Land<-as.matrix(expand.grid(0, seq(prepared.data$Lat-10, prepared.data$Lat+10, 0.5)))
	Points.Land<-cbind(Points.Land, 1)

	All.probs.dusk<-sapply(prepared.data$Dusk$Twilight.vector, FUN=get.prob.surface, Twilight.log.light.mat=prepared.data$Dusk$Twilight.log.light.mat, Twilight.time.mat=prepared.data$Dusk$Twilight.time.mat, dusk=T, Calib.param=calibration$Parameters$LogSlope, log.irrad.borders=log.irrad.borders, delta=prepared.data$delta, Points.Land=Points.Land,log.light.borders=log.light.borders, calibration=calibration)
	
	All.probs.dawn<-sapply(prepared.data$Dawn$Twilight.vector, FUN=get.prob.surface, Twilight.log.light.mat=prepared.data$Dawn$Twilight.log.light.mat, Twilight.time.mat=prepared.data$Dawn$Twilight.time.mat, dusk=F, Calib.param=calibration$Parameters$LogSlope, log.irrad.borders=log.irrad.borders, delta=prepared.data$delta, Points.Land=Points.Land,  log.light.borders=log.light.borders, calibration=calibration)

	Phys.Mat<-c()
	for (i in 1:nrow(prepared.data$Filtered_tw)) {
	if (prepared.data$Filtered_tw$type[i]==2) {
		Phys.Mat<-cbind(Phys.Mat, All.probs.dusk[,1])
		All.probs.dusk<-as.matrix(All.probs.dusk[,-1])
		} else {
		Phys.Mat<-cbind(Phys.Mat, All.probs.dawn[,1])
		All.probs.dawn<-as.matrix(All.probs.dawn[,-1])
		}
	}
	
	Diff<-try(Points.Land[which.max(apply(Phys.Mat[,1:(dim(Phys.Mat)[2])],1,  FUN=prod)),2]-prepared.data$Lat)
	#if (class(Diff)=="try-error") Diff=NA

	#Diff_1<-try(all.out$Points.Land[which.max(apply(all.out$Phys.Mat[,1:80],1,  FUN=prod)),2]-start[2])
	#if (class(Diff_1)=="try-error") Diff_1=NA

	return(Diff)
}


get.all.diffs<-function(Tracks, calibration, min.max.values=c(1, 1150), log.light.borders=log(c(2, 1100)), log.irrad.borders=c(-15, 50), cluster=NULL) {

if (!is.null(cluster)) {
	require(parallel)
	#if (threads==-1) threads= detectCores()-1 # selecting all -1 available cores...
	#mycl<-makeCluster(threads)
	#tmp<-parallel:::clusterSetRNGStream(mycl)
	### we don' need to send all parameters to node. so keep it easy..
	tmp<-parallel:::clusterExport(cluster, c( "get.diff", "log.irrad.borders", "log.light.borders","calibration"),  envir=environment())
	Diffs<-parSapply(cluster, Tracks, FUN=function(x) get.diff(prepared.data=x, calibration=calibration, min.max.values=min.max.values, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders))
	} else {
	Diffs<-sapply(Tracks, FUN=function(x) get.diff(prepared.data=x, calibration=calibration, min.max.values=min.max.values, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders))
	}
	return(Diffs)
	}

diffs.ll<-function(diffs, log=T) {
# least squares
ifelse(log, log(max(-1e70,mean(diffs^2, na.rm=T))), max(-1e70,mean(diffs^2, na.rm=T)))
}

test.deltas.old<-function(params, Tracks, Spline, calibration, min.max.values=c(1, 1150), log.light.borders=log(c(2, 1100)), log.irrad.borders=c(-15, 50), cluster=NULL) {
# params are parameters for delta...
print(params)
deltas=params[1] + Spline%*%params[2:4] # for params - first for the intercept.
Tracks <- mapply(append, Tracks, deltas, SIMPLIFY = FALSE)
Tracks <- lapply(Tracks, FUN=function(x) {names(x)[length(x)]="delta"; return(x)})
Diffs<-get.all.diffs(Tracks, calibration, min.max.values=min.max.values, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, cluster=cluster)
plot(Diffs~sapply(Tracks, "[[", i=3))
print(Diffs)
Res<-diffs.ll(Diffs, log=F)
cat("LL:", Res, "\n")
return(Res)
}


test.deltas<-function(params, Tracks, Spline, calibration, min.max.values=c(1, 1150), log.light.borders=log(c(2, 1100)), log.irrad.borders=c(-15, 50), cluster=NULL) {
# this function should send a spline as a calibration function..
# params are parameters for delta...
print(params)
#deltas=params[1] + Spline%*%params[2:4] # for params - first for the intercept.

#lat_correction_fun<-approxfun(y= bs(abs((-85:85)), knots=c(21, 42), degree=3, Boundary.knots=c(0,85), intercept=T)%*%params, x=-85:85)
#lat_correction_fun<-approxfun(y= cbind(1, c(-85:85)^2,c (-85:85)^4)%*%params, x=-85:85)
lat_correction_fun<-approxfun(y= cbind(1, cos(c(-85:85)/180*pi),cos(2*c(-85:85)/180*pi),cos(3*c(-85:85)/180*pi))%*%params, x=-85:85)

# trying to add knots..
calibration$lat_correction_fun=lat_correction_fun

Diffs<-get.all.diffs(Tracks, calibration, min.max.values=min.max.values, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, cluster=cluster)
par(mfrow=c(1,2))
plot(lat_correction_fun(-85:85)~ c(-85:85))

plot(Diffs~sapply(Tracks, "[[", i=3))
Res<-diffs.ll(Diffs, log=F)
cat("LL:", Res, "\n")
return(Res)
}





