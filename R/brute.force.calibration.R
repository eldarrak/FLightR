# brute.force.calibration.R
# ver 0.1 from 21.10.2014


get.1.minute.parameters<-function(parameters, saving.period=NULL, position, start.time, end.time, log.light.borders=c(log(2,64)), repeats=50) {
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
, to=end.time+7200, by=60)

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

To.run<-expand.grid(Slope.ideal= runif(1000,   parameters$LogSlope[1]-0.1,   parameters$LogSlope[1]+0.4), SD.ideal=runif(1000,  exp(log( parameters$LogSlope[2])-0.2),  exp(log( parameters$LogSlope[2])+0.2))) #
Res<-c()
Res<-as.data.frame(Res)
for (i in 1:repeats) {
To.run.cur<-To.run[sample( 1:nrow(To.run), 1),]
All.slope.runs=get.slopes(To.run=To.run.cur, Parameters=parameters, Lat=position[2], saving.period=saving.period, Time.seq=Time.seq, Time.seq.saving=Time.seq.saving, Lon=position[1], log.light.borders=log.light.borders)
Res<-rbind(Res, c(mean(All.slope.runs$Slope, na.rm=T), All.slope.runs$Slope.ideal[1], sd(All.slope.runs$Slope, na.rm=T), All.slope.runs$ SD.ideal[1]))
names(Res)<-c("Slope", "Slope.ideal", "SD", "SD.ideal")
print(Res)
}

#Parameters=All.slopes$Parameters


plot(Res$SD~Res$Slope.ideal) # no relationship
plot(Res$SD~Res$Slope) # no relationship
# ok, this means that slope and SD are independent. GOOD.
plot(Res$Slope~Res$Slope.ideal)
#summary(lm(Res$Slope~Res$Slope.ideal))
#=ok I do see now absolute linear relationship
# with the same degree for slope and for SD
# interesting, but I would not pay too much attention to that

plot(Res$SD~Res$SD.ideal)
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

All.slope.runs1=get.slopes(Repeats=25, To.run=To.run, Parameters=parameters, Lat=position[2], saving.period=600, Time.seq=Time.seq, Time.seq.saving=Time.seq.saving, Lon=position[1], log.light.borders=log.light.borders)
#mean(All.slope.runs1$Slope) # works
#sd(All.slope.runs1$Slope)

# and now let's repeat the same story with Lm...
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
get.time.correction.function<-function(parameters, saving.period=NULL, position, log.light.borders=log(c(2,64)), Repeats=2) {
# as far as we do not provide time the slope function will runf for the whole year.

To.run<-expand.grid(Slope.ideal=Parameters$LogSlope_1_minute[1], SD.ideal=Parameters$LogSlope_1_minute[2]) #
All.slope.runs=get.slopes(Repeats=Repeats, To.run=To.run, Parameters=parameters, Lat=position[2], saving.period=saving.period, short.run=F, Lon=position[1], log.light.borders=log.light.borders)

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

get.lat.correction.function<-function(deltalim=c(-0.05, 0.35), saving.period=NULL, threads=2, mode="smart", Sigmas, LogSlope,  log.irrad.borders=c(-50, 50), calibration=NULL) {

cat("function will work in ", mode, "mode\n")

if (!mode %in% c("trial", "brute", "smart")) stop ("mode could be one of - trial, brute, smart")
require(parallel)

if (mode=="rial") {
	deltalim=c(0.15, 0.16)
	cat("deltalim set to ", deltalim, "\n")
	Res<-get.deltas.parallel(deltalim=deltalim, limits=c(-65,65), points=2, Sigmas=sigma, interval=saving.period, short.run=T, threads=threads, log.irrad.borders=c(-50, 50),  random.delta=T, calibration=calibration)
	Res<-as.data.frame(Res)
	names(Res)<-c("Diff", "Sigma", "Delta", "Lat", "Diff.first", "Diff.second", "Sigma.init")
	Res$cosLat<-cos(Res$Lat/180*pi)
	RES<-list(simulations=Res)
}
if (mode=="brute") {
# real run
Res<-get.deltas.parallel(deltalim=deltalim, limits=c(-65,65), points=70, Sigmas=sigma, interval=saving.period, short.run=T, threads=threads, log.irrad.borders=c(-50, 50), random.delta=T, calibration=calibration)
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

	Res_min_lat<-get.deltas.parallel(deltalim=deltalim_initial, limits=c(0,0), points=Points, Sigmas=sigma, interval=saving.period, short.run=T, threads=threads, log.irrad.borders=c(-50, 50), random.delta=T, calibration=calibration)
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
	Res_max_lat<-get.deltas.parallel(deltalim=deltalim_initial, limits=c(65,65), points=Points, Sigmas=sigma, interval=saving.period, short.run=T, threads=threads, log.irrad.borders=c(-50, 50), random.delta=T, calibration=calibration)
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
	
	Res<-get.deltas.parallel(deltalim=deltalim_corrected, limits=c(-65,65), points=Points, Sigmas=sigma, interval=saving.period, short.run=T, threads=threads, log.irrad.borders=c(-50, 50), random.delta=T, calibration=calibration)

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
