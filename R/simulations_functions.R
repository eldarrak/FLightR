#simulations_functions.R
# these functions were used right befre the 0.2.1 version so they should work

get.shifts<-function(Track, Parameters, log.light.borders=log(c(2,64)), log.irrad.borders=c(-9,3), min.max.values=c(0,64), Grid, start, ask=FALSE, slopes.only=FALSE, delta=NULL, short.run=FALSE, measurement.period=60, saving.period=NULL, Time.seq=NULL, Time.seq.saving=NULL, calibration=NULL, GeoLight=FALSE) {
#========================================
if (length(unique(Track[,2])) !=1) stop("moving track is not implemented yet!")
# here is the lnorm distr that we currently use..
#Lnorm.param<-Parameters$LogSlope_1_minute

#===================== the new version will work with general simulation function \
#data(file.path(wd, "LightR_development_code\\get.slopes.5.0.r"))

# for this we will have to create a To.run.object...
To.run<-expand.grid(Slope.ideal=Parameters$LogSlope_1_minute[1], SD.ideal=Parameters$LogSlope_1_minute[2], Latitude=unique(Track[,2])) #
Track.initial<-Track
Track<-simulate.track(measurement.period=measurement.period, saving.period=saving.period, To.run=To.run, Parameters=Parameters, min.max.values=min.max.values, short.run=short.run, Time.seq=Time.seq, Time.seq.saving=Time.seq.saving)
#==================================
# saving and reading track file
cat("   saving file\n")

lig.data<-cbind(format(as.POSIXct(Track$gmt, tz="gmt",origin="1970-01-01"), format="%d/%m/%Y"), format(as.POSIXct(Track$gmt, tz="gmt",origin="1970-01-01"), format="%H:%M:%S"), round(exp(Track$LogLight)))

File.name<-tempfile(pattern = "sim.no.move.", tmpdir = getwd(), fileext = ".csv")
write.table(lig.data, file=File.name, sep = ",", dec = ".", qmethod="double", quote = FALSE,row.names = FALSE, col.names=FALSE)

## no I wnat to work with the simulated data.

cat("   reading file\n")

# set wd
Data<-geologger.read.data(file=File.name)
try(unlink(File.name))

#========================================================
#========================================================

# new trick - let's try to load the real track
cat("   GeoLight step\n")

tw <- twilightCalc(Track$gmt, Track$light, allTwilights=TRUE, ask=FALSE, LightThreshold=3)

#save(tw, file="tw.sim.no.move.RData")
#load("tw.sim.no.move.RData")


# so we will create a vector of points from Equator to Pole and will estimate the likelihood of out positions in that points..
#Grid<-cbind(start[1], seq(0, 67, 0.5))
#Grid<-cbind(-143, seq(50, 67, 0.5))

# now we want to solve the equation for every day..

GLtab   <- tw[[2]] # Table to proceed with GeoLight
cat("automatically detected twilight:\n")
print(table(GLtab[,3]))
 if (abs(diff(table(GLtab[,3])))>5 & ask) {
cat ("more than 5 twilights were detected incorrectly.. please do the selection by hand\n")
try(dev.off())
tw <- twilightCalc(Track$gmt, Track$light, allTwilights=TRUE, ask=TRUE, LightThreshold=3, nsee=1000)
save(tw, file=paste("tw.Lat", start[2], "attempt1.RData" , sep="."))
GLtab   <- tw[[2]]
} else {
# as far as we have done an automatic detection we should check on whether the is no problems..
# the idea is that each dusk should be follwoed by dawn...
Index<-1:length(GLtab[,3])
GLtab[Index[Index%%2==1],3]<-round(mean(GLtab[Index[Index%%2==1],3]))
GLtab[Index[Index%%2==0],3]<-round(mean(GLtab[Index[Index%%2==0],3]))

}
#==============================================================
# here is a brunch for geolight
if(GeoLight) {
GLtab_shifted <- twilightCalc(Track$gmt, Track$light, allTwilights=TRUE, ask=FALSE, LightThreshold=3, maxLight=saving.period/60)[[2]]
Elevation<-getElevation(GLtab_shifted$tFirst, GLtab_shifted$tSecond, GLtab_shifted$type, known.coord=start, plot=FALSE)
positionsGeoLight <- coord(GLtab_shifted$tFirst, GLtab_shifted$tSecond, GLtab_shifted$type, degElevation=Elevation)
positionsGeoLight<-as.data.frame(positionsGeoLight)
positionsGeoLight$tFirst<-GLtab_shifted$tFirst
#tripMap(positionsGeoLight)
# ok, good here are geolight positions.
}
#=================================================================
#filter <- loessFilter(GLtab[,1],GLtab[,2],GLtab[,3],k=10)

# 
#GLtab1<-GLtab[filter,]
GLtab1<-GLtab
Filtered_tw <- data.frame(datetime=as.POSIXct(c(GLtab1$tFirst,GLtab1$tSecond),"UTC"),type=c(GLtab1$type,ifelse(GLtab1$type==1,2,1)))


Filtered_tw <- Filtered_tw[!duplicated(Filtered_tw$datetime),]
Filtered_tw <- Filtered_tw[order(Filtered_tw[,1]),]



# ok, now the part for the template fit..
#

# now I want to pair data and twilights..		  
Filtered_tw$light<-approx(x=Track$gmt, y=Track$light, xout=Filtered_tw$datetime)$y
Filtered_tw$id<-0
Track$type.real<-Track$type
Track$type<-0

Track.new<-Track[names(Track) %in% c("gmt", "light", "type")]
Track.new<-rbind(Track.new, data.frame(gmt= Filtered_tw$datetime, light=Filtered_tw$light, type=Filtered_tw$type))

All.p<-Track.new[order(Track.new$gmt),]

rownames(All.p)<-1:nrow(All.p)


#=============================================================================
#================= START =====================================================
# NOW we need to add a stupid old part that will just create for the all out.object..
#processing dusk
raw.Y.dusk<-correct.hours(Filtered_tw$datetime[Filtered_tw$type==2])

#######################################################
raw.X.dusk<-as.numeric(Filtered_tw$datetime[Filtered_tw$type==2])
Result.Dusk<-make.result.list(Data, raw.X.dusk, raw.Y.dusk)


raw.Y.dawn<-correct.hours(Filtered_tw$datetime[Filtered_tw$type==1])

#######################################################################

raw.X.dawn<-as.numeric(Filtered_tw$datetime[Filtered_tw$type==1])

Result.Dawn<-make.result.list(Data, raw.X.dawn, raw.Y.dawn)

Result.all<-list(Final.dusk=Result.Dusk, Final.dawn=Result.Dawn)
####

# ok now I want to create an output..

Index.tab<-create.proposal(Result.all, start=start, Grid=Grid)
#save(Index.tab, file="Index.tab.RData")

Index.tab$yday<-as.POSIXlt(Index.tab$Date, tz="GMT")$yday

Index.tab$Decision<-0.1 # prob of migration
Index.tab$Direction<- 0 # direction 0 - North
Index.tab$Kappa<-0 # distr concentration 0 means even
Index.tab$M.mean<- 300 # distance
Index.tab$M.sd<- 500 # distance sd


##START POINTS###
all.out<-geologger.sampler.create.arrays(Index.tab, Grid, start=start)

#================= END== =====================================================

#=============================================================================

Proc.data<-process.twilights(All.p, Filtered_tw, measurement.period=measurement.period, saving.period=saving.period)

## Dusk
Twilight.time.mat.dusk<-Proc.data$Twilight.time.mat.dusk

Twilight.log.light.mat.dusk<-Proc.data$Twilight.log.light.mat.dusk

## Dawn

Twilight.time.mat.dawn<-Proc.data$Twilight.time.mat.dawn

Twilight.log.light.mat.dawn<-Proc.data$Twilight.log.light.mat.dawn


###############
# ok and now we should be ready for the estimates...
# 
#============== ver 4.0.c
#Lnorm.param.estimation[1]<-Lnorm.param.estimation[1]+(Lnorm.param.estimation[1]^2+exp(Parameters$LogSigma)^2)/2
#============================
# before going gor the simulation I would probably like to go for a simple estimation in just one point as this could help us a lot!
do.linear.regresion<-function(Twilight.ID, start, dusk=TRUE, Twilight.time.mat, Twilight.log.light.mat, return.slopes=FALSE,  Calib.param, log.irrad.borders, verbose=FALSE, log.light.borders=log(c(2,64))) {
#=========================================================================

		Twilight.log.light.vector<-Twilight.log.light.mat[c(1:24, 26:49), Twilight.ID]
		Twilight.time.vector=Twilight.time.mat[c(1:24, 26:49), Twilight.ID]
		Data<-check.boundaries(start, Twilight.solar.vector=NULL,  Twilight.log.light.vector=Twilight.log.light.vector, plot=FALSE, verbose=verbose,  log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, dusk=dusk, Twilight.time.vector=Twilight.time.vector)
		Res<-c(NA, NA)
		if (dim(Data)[1]>1) {
		LogLight<-Data[,1]
		LogIrrad<-Data[,2]
		
		#===========
		# here I'll check on how many points we have
			if (length(LogLight) >=2) { 
			Model<- stats::lm(LogLight~LogIrrad)
			Model<- stats::lm(LogLight~LogIrrad)
			if (verbose) print(summary(Model))
		# now we need to get the results:
		Res[1]<-stats::coef(Model)[2] 
		Res[2]<-summary(Model)$sigma

	}}
	return(Res)
	}


	Twilight.vector<-1:(dim(Twilight.time.mat.dusk)[2])

		Slopes.dusk<-sapply(Twilight.vector, FUN=do.linear.regresion, start=start, Twilight.log.light.mat=Twilight.log.light.mat.dusk, Twilight.time.mat=Twilight.time.mat.dusk, dusk=TRUE, log.irrad.borders=log.irrad.borders, log.light.borders=log.light.borders)
		
	Twilight.vector<-1:(dim(Twilight.time.mat.dawn)[2])

		Slopes.dawn<-sapply(Twilight.vector, FUN=do.linear.regresion, start=start, Twilight.log.light.mat=Twilight.log.light.mat.dawn, Twilight.time.mat=Twilight.time.mat.dawn, dusk=FALSE,log.irrad.borders=log.irrad.borders, log.light.borders=log.light.borders)

	All.probs.dusk<-c()
	All.probs.dawn<-c()
		cat("\n detected dusk", dim(Twilight.time.mat.dusk)[2], "\n")
		cat("\n detected dawn",dim(Twilight.time.mat.dawn)[2] ,"\n")

	if (!slopes.only) {
		#cat("detected dusk", dim(Twilight.time.mat.dusk)[2], "\n")
		#cat("detected dawn",dim(Twilight.time.mat.dawn)[2] ,"\n")
	Twilight.vector<-1:(dim(Twilight.time.mat.dusk)[2])
 
		 All.probs.dusk<-sapply(Twilight.vector, FUN=get.prob.surface, Twilight.log.light.mat=Twilight.log.light.mat.dusk, Twilight.time.mat=Twilight.time.mat.dusk, dusk=TRUE, Calib.param=Parameters$LogSlope, log.irrad.borders=log.irrad.borders, delta=delta, Grid=Grid,log.light.borders=log.light.borders, calibration=calibration)
	
	Twilight.vector<-1:(dim(Twilight.time.mat.dawn)[2])
		All.probs.dawn<-sapply(Twilight.vector, FUN=get.prob.surface, Twilight.log.light.mat=Twilight.log.light.mat.dawn, Twilight.time.mat=Twilight.time.mat.dawn, dusk=FALSE, Calib.param=Parameters$LogSlope,log.irrad.borders=log.irrad.borders, delta=delta, Grid=Grid,  log.light.borders=log.light.borders, calibration=calibration)
# ok, what should we do now?
# first I'd found a maximum of the probabilities for each day and comapre it with the means 

All.probs.dusk.tmp<-All.probs.dusk
All.probs.dawn.tmp<-All.probs.dawn
}
	Phys.Mat<-c()
	Slopes<-c()
for (i in 1:nrow(Filtered_tw)) {
	if (Filtered_tw$type[i]==2) {
		if (!slopes.only) {
		Phys.Mat<-cbind(Phys.Mat, All.probs.dusk.tmp[,1])
		All.probs.dusk.tmp<-as.matrix(All.probs.dusk.tmp[,-1])
		}
		Slopes<-rbind(Slopes, Slopes.dusk[,1])
		Slopes.dusk<-as.matrix(Slopes.dusk[,-1])
		} else {
		if (!slopes.only) {
		Phys.Mat<-cbind(Phys.Mat, All.probs.dawn.tmp[,1])
		All.probs.dawn.tmp<-as.matrix(All.probs.dawn.tmp[,-1])
		}
		Slopes<-rbind(Slopes, Slopes.dawn[,1])
		Slopes.dawn<-as.matrix(Slopes.dawn[,-1])
		}
}
all.out$Phys.Mat<-Phys.Mat
all.out$Phys.Mat.time.label<-Filtered_tw$datetime
all.out$Slopes<-Slopes
if (GeoLight) all.out$positionsGeoLight<-positionsGeoLight
		if (!slopes.only) {
		return(all.out)
		} else {
		return(Slopes)
		}
}


get.slopes<-function(Repeats=1, file.head="tmp", Lon=0, Lat=NULL, measurement.period=60, saving.period=NULL, To.run, Parameters=NULL, short.run=FALSE, Time.seq=NULL, Time.seq.saving=NULL, log.light.borders=log(c(2,64)), min.max.values=c(0, 64), log.irrad.borders=c(-15, 50) , plot=TRUE) {
To.run.initial<-To.run
Lat.initial<-Lat
All.slope.runs<-c()
 for (i in 1:Repeats) {
 cat("  repeat", i, "from", Repeats, "\n")
if (is.null(Lat.initial)) {
	
	if (ncol(To.run.initial)==2) { 
		To.run$Latitude<-round(runif(nrow(To.run), -55, 55))
		} else {
		cat("latitudes were provided with To.run object \n")
		}
	} else {
	To.run$Latitude<-Lat.initial
	}

Track<-simulate.track(measurement.period=measurement.period, saving.period=saving.period, To.run=To.run, Parameters=Parameters, short.run=short.run, Time.seq=Time.seq, Time.seq.saving=Time.seq.saving, Lon=Lon,  min.max.values= min.max.values, plot=plot)
	
#===========
# input
#SD.ideal<-To.run$SD.ideal[i]
#Slope.ideal<-To.run$Slope.ideal[i]
#Latitude<-To.run$Latitude[i]
#==============
# function

#start=c(0,Latitude) # let's make it random also.. but this should already come from TO.run

#Grid<-cbind(start[1], start[2], 1)

# now we will be measuring evey 60 seconds and save max inside a period for the future.




#================================================================
# estimation
#====================================
# this is the estimation part, so it should be after the generation part..
tw <- twilightCalc(Track$gmt, Track$light, allTwilights=TRUE, ask=FALSE, LightThreshold=3)

# now we want to solve the equation for every day..
GLtab   <- tw[[2]] # Table to proceed with GeoLight

# as far as we have done an automatic detection we should check on whether the is no problems..
# the idea is that each dusk should be follwoed by dawn...
Index<-1:length(GLtab[,3])
GLtab[Index[Index%%2==1],3]<-round(mean(GLtab[Index[Index%%2==1],3]))
GLtab[Index[Index%%2==0],3]<-round(mean(GLtab[Index[Index%%2==0],3]))

GLtab1<-GLtab

Filtered_tw <- data.frame(datetime=as.POSIXct(c(GLtab1$tFirst,GLtab1$tSecond),"UTC"),type=c(GLtab1$type,ifelse(GLtab1$type==1,2,1)))

Filtered_tw <- Filtered_tw[!duplicated(Filtered_tw$datetime),]
Filtered_tw <- Filtered_tw[order(Filtered_tw[,1]),]



# ok, now the part for the template fit..
#

# now I want to pair data and twilights..		  
Filtered_tw$light<-approx(x=Track$gmt, y=Track$light, xout=Filtered_tw$datetime)$y
Filtered_tw$id<-0
Track$type.real<-Track$type
Track$type<-0

Track.new<-Track[names(Track) %in% c("gmt", "light", "type")]
Track.new<-rbind(Track.new, data.frame(gmt= Filtered_tw$datetime, light=Filtered_tw$light, type=Filtered_tw$type))

All.p<-Track.new[order(Track.new$gmt),]

rownames(All.p)<-1:nrow(All.p)

Proc.data<-process.twilights(All.p, Filtered_tw, measurement.period=measurement.period, saving.period=saving.period)

## Dusk
Twilight.time.mat.dusk<-Proc.data$Twilight.time.mat.dusk

Twilight.log.light.mat.dusk<-Proc.data$Twilight.log.light.mat.dusk

## Dawn

Twilight.time.mat.dawn<-Proc.data$Twilight.time.mat.dawn

Twilight.log.light.mat.dawn<-Proc.data$Twilight.log.light.mat.dawn

#log.light.borders<-log(c(2,64))

#		cat("doing", Twilight.ID, "\n")	

# now I need to create Grid
Row<-sapply(as.numeric(Filtered_tw$datetime), FUN=function(x) which.min(abs(x-Track$Time.seq)))


# ok and here we need to get around somehow...
	Twilight.vector<-1:(dim(Twilight.time.mat.dusk)[2])
#Grid<-cbind(start[1],start[2],1)

All.probs.dusk<-c()
for (Twilight.ID in Twilight.vector) {
	Lon<-Track$Lon[Row[Filtered_tw$type==2][Twilight.ID]]
	Lat<-Track$Lat[Row[Filtered_tw$type==2][Twilight.ID]]
	Grid<-cbind(Lon, Lat)
	
	Prob.surf<-try(get.prob.surface(Twilight.ID=Twilight.ID, Twilight.log.light.mat=Twilight.log.light.mat.dusk, Twilight.time.mat=Twilight.time.mat.dusk, dusk=TRUE, return.slopes=TRUE, log.irrad.borders=log.irrad.borders, Calib.param=Parameters$LogSlope, Grid=Grid, delta=0, log.light.borders=log.light.borders))
	#print(str(Prob.surf))
	All.probs.dusk<-cbind(All.probs.dusk, Prob.surf)
 }
	Twilight.vector<-1:(dim(Twilight.time.mat.dawn)[2])

	All.probs.dawn<-c()
	for (Twilight.ID in Twilight.vector) {
	Lon<-Track$Lon[Row[Filtered_tw$type==1][Twilight.ID]]
	Lat<-Track$Lat[Row[Filtered_tw$type==1][Twilight.ID]]
	#print(Grid)
	Grid<-cbind(Lon, Lat)
	Prob.surf<-try(get.prob.surface(Twilight.ID=Twilight.ID, Twilight.log.light.mat=Twilight.log.light.mat.dawn, Twilight.time.mat=Twilight.time.mat.dawn, dusk=FALSE, return.slopes=TRUE, log.irrad.borders=log.irrad.borders, Calib.param=Parameters$LogSlope, Grid=Grid, delta=0, log.light.borders=log.light.borders))
	#print(str(Prob.surf))
	All.probs.dawn<-cbind(All.probs.dawn, Prob.surf)
	}

	# ok now we need to bring the slopes back to scale...

All.slopes.dusk<-t(All.probs.dusk[2:3,])
All.slopes.dawn<-t(All.probs.dawn[2:3,])

All.slopes<-c()
for (k in 1:nrow(Filtered_tw)) {
	if (Filtered_tw$type[k]==2) {
		All.slopes<-rbind(All.slopes, All.slopes.dusk[1,])
		All.slopes.dusk<-as.matrix(All.slopes.dusk[-1,])
		} else {
		All.slopes<-rbind(All.slopes, All.slopes.dawn[1,])
		All.slopes.dawn<-as.matrix(All.slopes.dawn[-1,])
		}
}

All.slopes<-as.data.frame(All.slopes)
names(All.slopes)[1:2]<-c("Slope", "Slope.sd")
All.slopes$gmt<-as.numeric(Filtered_tw$datetime)
All.slopes$type<-as.numeric(Filtered_tw$type)


All.slopes<-cbind(All.slopes, Track[Row, c("Slope.ideal", "SD.ideal", "Lat")])
All.slopes$Slope<- log(All.slopes$Slope)
if (plot) graphics::plot(All.slopes$Slope~All.slopes$Slope.ideal)

All.slope.runs<-rbind(All.slope.runs, All.slopes)
if (plot)  {
graphics::par(mfrow=c(1,2))
graphics::plot((All.slope.runs$Slope)~All.slope.runs$Slope.ideal)
mean((All.slope.runs$Slope), na.rm=TRUE)
graphics::plot(Slope~Lat, data=All.slope.runs)
}
}
return(All.slope.runs)
}
# and now we want to save that... 


simulate.track<-function(measurement.period=60, saving.period=600, To.run, Parameters=Parameters, short.run=FALSE, Time.seq=NULL, Time.seq.saving=NULL, Lon=0, min.max.values=c(0, 64), first.date="2010-01-01 00:00:00", last.date="2010-03-20 23:59:59", plot=TRUE) {
# important here is that min and max values may be different from light.borders.
# and it is actually better to make them different if there is enough point to make an estimation...

if (saving.period%%measurement.period !=0) stop("saving period / measurement.period has to be integer!")
time.shift<-sample(1:saving.period, 1)
if (is.null(Time.seq) | is.null(Time.seq.saving)) {
	if (!short.run) {
	Time.seq<-seq(from=as.numeric(as.POSIXct(first.date, tz="UTC")), to=as.numeric(as.POSIXct("2010-12-31 23:59:59", tz="UTC")), by=measurement.period)+time.shift

	Time.seq.saving<-seq(from=as.numeric(as.POSIXct(first.date, tz="UTC")), to=as.numeric(as.POSIXct("2010-12-31 23:59:59", tz="UTC")), by=saving.period)+time.shift
	} else {
	Time.seq<-seq(from=as.numeric(as.POSIXct(first.date, tz="UTC")), to=as.numeric(as.POSIXct(last.date, tz="UTC")), by=measurement.period)+time.shift

	Time.seq.saving<-seq(from=as.numeric(as.POSIXct(first.date, tz="UTC")), to=as.numeric(as.POSIXct(last.date, tz="UTC")), by=saving.period)+time.shift
	}
}

Track<-cbind(Lon, NA, Time.seq)
Track<-as.data.frame(Track)
names(Track)<-c("Lon","Lat", "Time.seq")
Track$type<-ifelse(as.numeric(format(as.POSIXct(Time.seq, tz="utc", origin="1970-01-01"), "%H"))<12, "Dawn", "Dusk")

Track$type<-as.numeric(as.factor(Track$type))

Rle<-rle(Track$type)
Rle$values<-1:length(Rle$values)
Track$Day<-inverse.rle(Rle)

#=============================
# ok now we have everything we need to add corrdinates, slopes and sd
#=============================
print(To.run)
Index<-sample.int(nrow(To.run), max(Track$Day), replace=TRUE)


tmpRle<-Rle
tmpRle$values<-To.run$Latitude[Index]
Track$Lat<-inverse.rle(tmpRle)

tmpRle<-Rle
tmpRle$values<-To.run$Slope.ideal[Index]
Track$Slope.ideal<-inverse.rle(tmpRle)

tmpRle<-Rle
tmpRle$values<-To.run$SD.ideal[Index]
Track$SD.ideal<-inverse.rle(tmpRle)


#======================================
# ok now we could go for the estimation

# step 2. generating light curve.
# Angles
cat("   Estimating solar angles\n")
cat("Time...")
Time<-as.POSIXct(Track[,3], tz="gmt", origin="1970-01-01")
cat("Solar...")
S<-solar.FLightR(Time)
Track.row<-1:dim(Track)[1]
cat("Angles...")
Angles<-sapply(Track.row,  FUN=function(x) { elevation(Track[x,1], Track[x, 2], lapply(S, "[", i=x))})
cat("Irradiance\n")
# Irradiance
Irradiance<-sapply(Angles, FUN=function(x) get.Irradiance(x*pi/180))
# ok, now I want to use the output from model as it has everything we need...
Track<-cbind(Track, Irradiance, Angles)

Track<-as.data.frame(Track)
cat("defining dusks and dawns\n")
# we need to add Day now..
# let's do it by taking making switch att the maxuimum - minimum of a daily curve..

cat("   generating Light levels from the model\n")

Track$LogLight<-NA

	for (j in unique(Track$Day)) {
		
		Current.sample.size<-length(which((Track$Day==j)))
		Current.values<- log(Track$Irradiance[Track$Day==j])*exp(rnorm(1, (Track$Slope.ideal[Track$Day==j][1]) , Track$SD.ideal[Track$Day==j][1]))  + rnorm(1, Parameters$Intercept[1], Parameters$Intercept[2]) + rnorm(Current.sample.size, 0, exp(rnorm(1,  Parameters$LogSigma[1], Parameters$LogSigma[2])))

		# plot(Current.values)
		# and now we need to add these valeus to the lines
		Track$LogLight[Track$Day==j]<-Current.values
		}
		
#======================================
# now we just want to estimate values for each line in Track

Track$LogLight[Track$LogLight>log(min.max.values[2])] <-log(min.max.values[2])

Track$LogLight[Track$LogLight<max(log(min.max.values[1]), -1)] <-max(log(min.max.values[1]), -1)

Track$light<-exp(Track$LogLight)

Track$light[Track$light<min.max.values[1]] <-min.max.values[1]


 if (!short.run & plot) graphics::plot(Track$light[5000:6000], type="b", pch=".")

 Track$gmt<-as.POSIXct(Track$Time.seq, tz="gmt", origin="1970-01-01")

 # now we need to get the estimates without saving file I'd say 
Track.new<-Track[Track$Time.seq %in% Time.seq.saving,] # creating new track

if (saving.period!=measurement.period) {
Track.new<-Track.new[-1,]
New.light<-Track.new$light
for (i in 2:(saving.period/measurement.period)) {
	#New.light<-pmax(New.light, Track$light[Track$Time.seq %in% (Time.seq.saving-(60*(i-1)))] )
	New.light<-pmax(New.light, Track$light[Track$Time.seq %in% ((Time.seq.saving[-1])-(measurement.period*(i-1)))] )
}
} else  {
New.light<-Track.new$light
}
Track.new$light<-New.light
# ok this worked..
Track.old<-Track
#Track<-Track.old
Track<-Track.new

return(Track)
}


## the version is copied from
# ver.4.0.-4.4.estimating.error.fixed.sd.fast.1-10.minutes.r

# decided to add a function that will complete the task at different latitudes with the parallel package..
# the idea is very simple:
# we just make a general wrapper that runs the whole thing at each latitude and then combines results...


get.deltas.one.basic<-function(delta=0, start=c(0,0), Sigma=0.5, return.all.out=FALSE, measurement.period=60, saving.period=600, short.run=TRUE, calibration=NULL, log.light.borders=log(c(2,64)), min.max.values=c(0,64), log.irrad.borders=c(-15, 50)) {
# so this function will have to return 1 0 or -1

Grid<-as.matrix(expand.grid(start[1], seq(start[2]-8, start[2]+8, 0.5)))
Grid<-cbind(Grid, 1)
Time.seq<-seq(from=as.numeric(as.POSIXct("2010-01-01 00:00:00", tz="UTC")), to=as.numeric(as.POSIXct("2010-03-22 23:59:59", tz="UTC")), by=saving.period)
Track<-cbind(start[1], start[2], Time.seq)

Parameters=calibration$Parameters
Parameters$LogSigma=c(log(Sigma), 0.0)
all.out<-get.shifts(Track=Track, Grid=Grid, start=start, Parameters=Parameters, log.light.borders=log.light.borders, min.max.values=min.max.values, log.irrad.borders=log.irrad.borders, slopes.only=FALSE, delta=delta, short.run=short.run, measurement.period=measurement.period, saving.period=saving.period, calibration=calibration)
Real.Sigma<-mean(all.out$Slopes[,2], na.rm=TRUE)
# and now we need to get resulting value

Diff<-try(all.out$Spatial$Grid[which.max(apply(all.out$Phys.Mat[,1:(dim(all.out$Phys.Mat)[2])],1,  FUN=prod)),2]-start[2])
if (class(Diff)=="try-error") Diff=NA

Diff_1<-try(all.out$Spatial$Grid[which.max(apply(all.out$Phys.Mat[,1:80],1,  FUN=prod)),2]-start[2])
if (class(Diff_1)=="try-error") Diff_1=NA

Diff_2<-try(all.out$Spatial$Grid[which.max(apply(all.out$Phys.Mat[,80:(dim(all.out$Phys.Mat)[2])],1,  FUN=prod)),2]-start[2])
if (class(Diff_2)=="try-error") Diff_2=NA

Res<-c(Diff, Real.Sigma, delta, start[2], Diff_1=Diff_1, Diff_2=Diff_2)
if (return.all.out) {return(all.out)
} else {return(Res)}
}

get.deltas.intermediate<-function(deltalim=c(-0.2, 0.2), start=c(0,0), Sigma=0.5, measurement.period=60, saving.period=600, short.run=TRUE, repeats=3, random.delta=TRUE, fast=FALSE, calibration=NULL, log.light.borders=log(c(2,64)), log.irrad.borders, min.max.values=c(0,64)) {
	if (random.delta) {
		if (fast) {
			Deltas<-rep(runif(5, deltalim[1], deltalim[2]), repeats)
		} else {
			Deltas<-rep(runif(ceiling((deltalim[2]-deltalim[1])/0.02), deltalim[1], deltalim[2]), repeats)
		}
	} else {
	Deltas<-rep(seq(deltalim[1], deltalim[2], 0.01), repeats)
	}
Res<-c()
for (i in Deltas) {
Res<-rbind(Res, get.deltas.one.basic(delta=i, start=start, Sigma=Sigma, measurement.period=measurement.period, saving.period=saving.period, short.run=short.run, calibration=calibration, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, min.max.values=min.max.values))
#try(print(tail(Res, 20)))
}
return(Res)
}

# and now the next on that will iterate Sigma
get.deltas.main<-function(deltalim=c(-0.2, 0.2), start=c(0,0), Sigmas=seq(0, 0.8, 0.1), measurement.period=60, saving.period=600, short.run=TRUE, repeats=3, random.delta=TRUE, fast=FALSE, calibration=NULL, log.light.borders=log(c(2,64)),  log.irrad.borders=c(-15, 50), min.max.values=c(0,60)) {

Res<-c()
for (i in Sigmas) {
cat(Sigmas, "\n")
Res_local<-try(get.deltas.intermediate(deltalim=deltalim, start=start, Sigma=i, measurement.period=measurement.period, saving.period=saving.period, short.run=short.run, repeats=repeats, random.delta=random.delta, fast=fast, calibration=calibration, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, min.max.values=min.max.values))
if (class(Res_local) == "try-error") {
		save(list = ls(all.names = TRUE), file=paste("Res", start[2],"tmp.RData", sep="."), envir=environment())
		
Res<-rbind(Res,cbind(rep(NA, 6), i))

} else {
Res<-rbind(Res,cbind(Res_local, i))
}
#print(Res)
#save(Res, file=paste("Res", start[2],"tmp.RData", sep="."))
}

return(Res)
}


get.deltas.parallel<-function(deltalim=c(-0.2, 0.2), limits=c(-65,65), points=20, Sigmas=seq(0, 0.8, 0.1), measurement.period=60, saving.period=600, short.run=TRUE, threads=2, log.light.borders=log(c(2,64)), log.irrad.borders=c(-15, 50), repeats=1, random.delta=TRUE, calibration=NULL, fast=FALSE, min.max.values=c(0, 64)) {

if (is.character(measurement.period)) measurement.period=get("measurement.period")
if (is.character(saving.period)) saving.period=get("saving.period")
if (is.character(Sigmas)) Sigmas=get("Sigmas")
if (is.character(calibration)) calibration=get("calibration")
if (is.character(deltalim)) deltalim=get("deltalim")
#Parameters=list(Intercept=c(3.71, 1.25), LogSlope=c(0.72, 0.4)),

# points means number of latitudes that should be used for the run..
# the question is what we have to download before we can run this on cluster
mycl<-parallel::makeCluster(threads)
Lats<-runif(points, min(limits), max(limits))
Coords<-cbind(0, Lats)
    tmp<-parallel::clusterSetRNGStream(mycl)
    ### we don' need to send all parameters to node. so keep it easy..
    tmp<-parallel::clusterExport(mycl, c("calibration","log.irrad.borders", "log.light.borders","min.max.values", "deltalim"), envir=environment())
    tmp<-parallel::clusterEvalQ(mycl, library("FLightR"))
	tryCatch(Res<-parallel::parApply(mycl, Coords, 1, FUN=function(x) as.data.frame(get.deltas.main(start=x,  deltalim=deltalim, Sigmas=Sigmas, measurement.period=measurement.period, saving.period=saving.period, short.run=short.run, repeats=1, random.delta=random.delta, fast=fast, calibration=calibration, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, min.max.values=min.max.values))), finally = parallel::stopCluster(mycl))
	#Res1<-apply(Coords, 1, FUN=function(x) as.data.frame(get.deltas.main(start=x,  deltalim=deltalim, Sigmas=Sigmas, interval=interval, short.run=short.run, LogSlope=LogSlope, Parameters=Parameters, repeats=1, random.delta=random.delta)))
	parallel::stopCluster(cl = mycl)
	Res<-do.call(rbind.data.frame, Res)
	return(Res)
}


# brute.force.calibration.R
# ver 0.3 from 04.08.2015


get.1.minute.parameters<-function(initial=NULL, parameters, position, start.time, end.time, print=TRUE) {
# this stupid brute force function will just etimate what should be the parameter values on a 1 minute scale..
if (is.null(initial))initial=parameters$LogSlope
if (print) print(initial)

Time.seq<-seq(from=start.time-7200
, to=end.time+7200, by=parameters$measurement.period)

# here the fixing period should be used.
#saving.period=600
Time.seq.saving<-seq(from=start.time-7200
, to=end.time+7200, by=parameters$saving.period)
To.run.cur=data.frame(Slope.ideal=initial[1],SD.ideal=initial[2])
All.slope.runs=get.slopes(To.run=To.run.cur, Parameters=parameters, Lat=position[2], measurement.period=parameters$measurement.period, saving.period=parameters$saving.period, Time.seq=Time.seq, Time.seq.saving=Time.seq.saving, Lon=position[1], log.light.borders=parameters$log.light.borders, min.max.values=parameters$min.max.values, log.irrad.borders=parameters$log.irrad.borders, plot=FALSE)

#########
# LL
LL=log((mean(All.slope.runs$Slope, na.rm=TRUE)-parameters$LogSlope[1])^2+(stats::sd(All.slope.runs$Slope, na.rm=TRUE)-parameters$LogSlope[2])^2)
cat("\n\n   ############## LL:", LL, "\n\n")
return(LL)
}


#
get.time.correction.function<-function(parameters, measurement.period=60, saving.period=NULL, position, log.light.borders=log(c(2,64)), Repeats=2, min.max.values=c(0,64), log.irrad.borders=c(-15, 50)) {
# as far as we do not provide time the slope function will runf for the whole year.

To.run<-expand.grid(Slope.ideal=Parameters$LogSlope_1_minute[1], SD.ideal=Parameters$LogSlope_1_minute[2]) #
All.slope.runs=get.slopes(Repeats=Repeats, To.run=To.run, Parameters=parameters, Lat=position[2], measurement.period=measurement.period, saving.period=saving.period, short.run=FALSE, Lon=position[1], log.light.borders=log.light.borders, min.max.values=min.max.values, log.irrad.borders=log.irrad.borders)

Solar<-solar.FLightR(as.POSIXct(All.slope.runs$gmt, tz="UTC", origin="1970-01-01"))
All.slope.runs$cosSolarDec<-Solar$cosSolarDec
#save(All.slope.runs, file="All.slope.runs_time_correction_600.RData")

# ok cos is linarly related!!!
# so let's make it linear
graphics::plot(Slope~cosSolarDec, data=All.slope.runs)
All.slope.runs<-All.slope.runs[is.finite(All.slope.runs$Slope),]
Lm2<-stats::lm(Slope~cosSolarDec, data=All.slope.runs)
#summary(Lm2)
# ok, so this will be used as a time related pattern
# I think we should make a function of that and add it to the next step...

time_correction_fun<-stats::approxfun(x=c(0.9, 1), y=predict(Lm2, newdata=data.frame(cosSolarDec=c(0.9, 1))))

Res<-list(time_correction_fun=time_correction_fun, coef=stats::coef(Lm2), results=All.slope.runs)
return(Res)
}
#

get.lat.correction.function<-function(deltalim=c(-0.05, 0.35), measurement.period=60, saving.period=NULL, threads=2, mode="smart", Sigmas, LogSlope,  log.light.borders=log(c(2, 64)),log.irrad.borders=c(-50, 50), calibration=NULL, min.max.values=c(0,64)) {

cat("function will work in ", mode, "mode\n")

if (!mode %in% c("trial", "brute", "smart")) stop ("mode could be one of - trial, brute, smart")
if (mode=="trial") {
	deltalim=c(0.15, 0.16)
	cat("deltalim set to ", deltalim, "\n")
	Res<-get.deltas.parallel(deltalim=deltalim, limits=c(-65,65), points=2, Sigmas=sigma, measurement.period=measurement.period, saving.period=saving.period, short.run=TRUE, threads=threads, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders,  random.delta=TRUE, calibration=calibration, min.max.values=min.max.values)
	Res<-as.data.frame(Res)
	print(str(Res))
	names(Res)<-c("Diff", "Sigma", "Delta", "Lat", "Diff.first", "Diff.second", "Sigma.init")
	Res$cosLat<-cos(Res$Lat/180*pi)
	RES<-list(simulations=Res)
}
if (mode=="brute") {
# real run
Res<-get.deltas.parallel(deltalim=deltalim, limits=c(-65,65), points=70, Sigmas=sigma, measurement.period=measurement.period, saving.period=saving.period, short.run=TRUE, threads=threads, log.irrad.borders=log.irrad.borders, random.delta=TRUE, calibration=calibration, min.max.values=min.max.values)
#=======
# !!! idealy there should be a check that will make sure that active variation occured at the deltalim specified
#=======
Res<-as.data.frame(Res)
names(Res)<-c("Diff", "Sigma", "Delta", "Lat", "Diff.first", "Diff.second", "Sigma.init")
Res$cosLat<-cos(Res$Lat/180*pi)
Model.all=try(mgcv::gam(Delta~te(Diff, cosLat),  data=Res)) # this look the best
if (class(Model.all)=="try-error") {
RES<-list(simulations=Res)	
} else {
print(summary(Model.all))
lat_correction_fun<-stats::approxfun(y=mgcv::predict.gam(Model.all, newdata=data.frame(cosLat=cos(c(-90:90)/180*pi), Diff=0)), x=-90:90)

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

	Res_min_lat<-get.deltas.parallel(deltalim=deltalim_initial, limits=c(0,0), points=Points, Sigmas=sigma, measurement.period=measurement.period, saving.period=saving.period, short.run=TRUE, threads=threads, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, random.delta=TRUE, calibration=calibration, min.max.values=min.max.values)
	Res_min_lat<-as.data.frame(Res_min_lat)
	names(Res_min_lat)<-c("Diff", "Sigma", "Delta", "Lat", "Diff.first", "Diff.second", "Sigma.init")
	Res_min_lat$cosLat<-cos(Res_min_lat$Lat/180*pi)
	# we do not need model here as the situation is very simples
	
	Res_min_lat<-Res_min_lat[Res_min_lat$Diff>-8 & Res_min_lat$Diff<8,]

	Gam_min<-mgcv::gam(Delta~s(Diff, k=3), data=Res_min_lat)
	Predict_min<-mgcv::predict.gam(Gam_min, se.fit=TRUE, newdata=data.frame(Diff=0))
	
	
	cat("    Done!")	

	cat("estimated delta for 0 degrees is" , round(Predict_min$fit,3),  "+-"  , round(Predict_min$se.fit,3), "\n")
	graphics::plot(Delta~Diff, data=Res_min_lat, col="red")
	abline(v=0, col="red")
	if ((Predict_min$fit-3*Predict_min$se.fit) < deltalim_initial[1]) stop ("try to correct lower deltalim boundary to a smaller value\n")

	
 cat("...estimating compensation at latitude 65")
	Res_max_lat<-get.deltas.parallel(deltalim=deltalim_initial, limits=c(65,65), points=Points, Sigmas=sigma, measurement.period=measurement.period, saving.period=saving.period, short.run=TRUE, threads=threads, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, random.delta=TRUE, calibration=calibration, min.max.values=min.max.values)
	Res_max_lat<-as.data.frame(Res_max_lat)
	names(Res_max_lat)<-c("Diff", "Sigma", "Delta", "Lat", "Diff.first", "Diff.second", "Sigma.init")
	Res_max_lat<-Res_max_lat[Res_max_lat$Diff>-8 & Res_max_lat$Diff<8,]
	
	Res_max_lat$cosLat<-cos(Res_max_lat$Lat/180*pi)
	cat("    Done!")	

	Gam_max<-mgcv::gam(Delta~s(Diff, k=3), data=Res_max_lat)
	Predict_max<-mgcv::predict.gam(Gam_max, se.fit=TRUE, newdata=data.frame(Diff=0))

	
	cat("estimated delta for 65 degrees N is" , round(Predict_max$fit,3),  "+-"  , round(Predict_max$se.fit,3), "\n")
	
	graphics::plot(Delta~Diff, data=Res_max_lat)
	graphics::points(Delta~Diff, data=Res_min_lat, col="red")
	graphics::abline(v=0, col="red")

	if ((Predict_max$fit+3*Predict_max$se.fit) > deltalim_initial[2]) stop ("try to correct upper  deltalim boundary to a higher value\n")
	
 	graphics::plot(Delta~Diff, data=Res_max_lat)
	graphics::points(Delta~Diff, data=Res_min_lat, col="red")
	Gam_max<-mgcv::gam(Delta~s(Diff, k=3), data=Res_max_lat)
	Predict_max<-mgcv::predict.gam(Gam_max, se.fit=TRUE, newdata=data.frame(Diff=0))
	
	# ok and now we want at least  to take a diap from min to max and focus there..
	
	deltalim_corrected<-c(Predict_min$fit-0.1, 	Predict_max$fit+0.1)
	
	
	##################
	# checking for -65	
	# ok now we have to go for the full run in a new boundaries...
	
	cat("doing the full run \n")
	
	Points<-(50%/%threads)*threads # how many repeats to run..
	
	Res<-get.deltas.parallel(deltalim=deltalim_corrected, limits=c(-65,65), points=Points, Sigmas=sigma, measurement.period=measurement.period, saving.period=saving.period, short.run=TRUE, threads=threads, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, random.delta=TRUE, calibration=calibration, min.max.values=min.max.values)

	cat(" ... done!")
	Res<-as.data.frame(Res)
	names(Res)<-c("Diff", "Sigma", "Delta", "Lat", "Diff.first", "Diff.second", "Sigma.init")
	
	Res<-Res[Res$Diff>-8 & Res$Diff<8,]

	
	Res$cosLat<-cos(Res$Lat/180*pi)
	
	#=== 
	# now we have to combine the three
	
	Res.all<-rbind(Res_min_lat,Res_max_lat, Res)
	
	Model.all=try(mgcv::gam(Delta~te(Diff, cosLat),  data=Res.all)) # this look the best
	if (class(Model.all)=="try-error") {
	RES<-list(simulations=Res)	
	} else {
	print(summary(Model.all))
	lat_correction_fun<-approxfun(y=mgcv::predict.gam(Model.all, newdata=data.frame(cosLat=cos(c(-90:90)/180*pi), Diff=0)), x=-90:90)

	RES<-list(simulations=Res, lat_correction_fun=lat_correction_fun)
	
 }

}
# ok this is it so far
return(RES)
}

#==================================
# here are the new functions for the latitude correction estimation...
# they are not wrapped yet though..


simulate.and.prepare.track<-function(measurement.period=60, saving.period=600, Parameters=Parameters, short.run=TRUE, min.max.values=c(0, 64),log.light.borders=log(c(2,64)),Lat, first.date="2010-01-01 00:00:00", last.date="2010-03-20 23:59:59") {

# for this we will have to create a To.run.object...
To.run<-data.frame(Slope.ideal=Parameters$LogSlope_1_minute[1], SD.ideal=Parameters$LogSlope_1_minute[2], Latitude=Lat) #

Track<-simulate.track(measurement.period=measurement.period, saving.period=saving.period, To.run=To.run, Parameters=Parameters, min.max.values=min.max.values, short.run=short.run, Time.seq=NULL, Time.seq.saving=NULL, Lon=0, first.date=first.date, last.date=last.date)
# ok, now we need to process it..
# to prepare for the run..

#==================================
# saving and reading track file
cat("   saving file\n")
lig.data<-cbind(format(as.POSIXct(Track$gmt, tz="gmt",origin="1970-01-01"), format="%d/%m/%Y"), format(as.POSIXct(Track$gmt, tz="gmt",origin="1970-01-01"), format="%H:%M:%S"), round(exp(Track$LogLight)))

File.name<-tempfile(pattern = "sim.no.move.", tmpdir = getwd(), fileext = ".csv")
write.table(lig.data, file=File.name, sep = ",", dec = ".", qmethod="double", quote = FALSE,row.names = FALSE, col.names=FALSE)

## no I wnat to work with the simulated data.
cat("   reading file\n")
# set wd
Data<-geologger.read.data(file=File.name)
try(unlink(File.name))

#========================================================
#========================================================

# new trick - let's try to load the real track
cat("   GeoLight step\n")

tw <- twilightCalc(Track$gmt, Track$light, allTwilights=TRUE, ask=FALSE, LightThreshold=log.light.borders[1]+1)
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
#all.out<-geologger.sampler.create.arrays(Index.tab, Grid, start=start)

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


get.grid.of.tracks<-function(measurement.period=60, saving.period=600, Parameters=Parameters, short.run=TRUE, min.max.values=c(0, 64), log.light.borders=log(c(2,64)),   Lats, cluster=NULL, first.date="2010-01-01 00:00:00",last.date="2010-03-20 23:59:59") {
if (!is.null(cluster)) {
	#if (threads==-1) threads= detectCores()-1 # selecting all -1 available cores...
	#cluster<-makeCluster(threads)
		tmp<-parallel::clusterSetRNGStream(cluster)
		### we don' need to send all parameters to node. so keep it easy..
		tmp<-parallel::clusterExport(cluster, c( "simulate.and.prepare.track", "Parameters", "measurement.period","saving.period", "short.run", "min.max.values", "log.light.borders", "first.date", "simulate.track", "last.date"),  envir=environment())
		tmp<-parallel::clusterEvalQ(cluster, library("FLightR"))
		tryCatch(Tracks<-parallel::parLapply(cluster, Lats,fun=function(x) simulate.and.prepare.track(measurement.period=measurement.period, saving.period=saving.period,  Parameters=Parameters, min.max.values=min.max.values, short.run=short.run, log.light.borders=log.light.borders, Lat=x, first.date=first.date, last.date=last.date)), finally = parallel::stopCluster(mycl))
	#stopCluster(mycl)
} else {
Tracks<-lapply(Lats,FUN=function(x) simulate.and.prepare.track(measurement.period=measurement.period, saving.period=saving.period,  Parameters=Parameters, min.max.values=min.max.values, short.run=short.run, log.light.borders=log.light.borders, Lat=x, first.date=first.date, last.date=last.date))
}
return(Tracks)
}


get.diff<-function(prepared.data, calibration, min.max.values=c(1, 1150), log.light.borders=log(c(2, 1100)), log.irrad.borders=c(-15, 50)) {

	Grid<-as.matrix(expand.grid(0, seq(prepared.data$Lat-10, prepared.data$Lat+10, 0.5)))
	Grid<-cbind(Grid, 1)

	All.probs.dusk<-sapply(prepared.data$Dusk$Twilight.vector, FUN=get.prob.surface, Twilight.log.light.mat=prepared.data$Dusk$Twilight.log.light.mat, Twilight.time.mat=prepared.data$Dusk$Twilight.time.mat, dusk=TRUE, Calib.param=calibration$Parameters$LogSlope, log.irrad.borders=log.irrad.borders, delta=prepared.data$delta, Grid=Grid,log.light.borders=log.light.borders, calibration=calibration)
	
	All.probs.dawn<-sapply(prepared.data$Dawn$Twilight.vector, FUN=get.prob.surface, Twilight.log.light.mat=prepared.data$Dawn$Twilight.log.light.mat, Twilight.time.mat=prepared.data$Dawn$Twilight.time.mat, dusk=FALSE, Calib.param=calibration$Parameters$LogSlope, log.irrad.borders=log.irrad.borders, delta=prepared.data$delta, Grid=Grid,  log.light.borders=log.light.borders, calibration=calibration)

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
	
	Diff<-try(Grid[which.max(apply(Phys.Mat[,1:(dim(Phys.Mat)[2])],1,  FUN=prod)),2]-prepared.data$Lat)
	#if (class(Diff)=="try-error") Diff=NA

	#Diff_1<-try(all.out$Grid[which.max(apply(all.out$Phys.Mat[,1:80],1,  FUN=prod)),2]-start[2])
	#if (class(Diff_1)=="try-error") Diff_1=NA

	return(Diff)
}


get.all.diffs<-function(Tracks, calibration, min.max.values=c(1, 1150), log.light.borders=log(c(2, 1100)), log.irrad.borders=c(-15, 50), cluster=NULL) {

if (!is.null(cluster)) {
	tmp<-parallel::clusterExport(cluster, c( "get.diff", "log.irrad.borders", "log.light.borders","calibration"),  envir=environment())
	tryCatch(Diffs<-parallel::parSapply(cluster, Tracks, FUN=function(x) get.diff(prepared.data=x, calibration=calibration, min.max.values=min.max.values, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders)))
	} else {
	Diffs<-sapply(Tracks, FUN=function(x) get.diff(prepared.data=x, calibration=calibration, min.max.values=min.max.values, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders))
	}
	return(Diffs)
	}

diffs.ll<-function(diffs, log=TRUE) {
# least squares
ifelse(log, log(max(-1e70,mean(diffs^2, na.rm=TRUE))), max(-1e70,mean(diffs^2, na.rm=TRUE)))
}

test.deltas.old<-function(params, Tracks, Spline, calibration, min.max.values=c(1, 1150), log.light.borders=log(c(2, 1100)), log.irrad.borders=c(-15, 50), cluster=NULL) {
# params are parameters for delta...
print(params)
deltas=params[1] + Spline%*%params[2:4] # for params - first for the intercept.
Tracks <- mapply(append, Tracks, deltas, SIMPLIFY = FALSE)
Tracks <- lapply(Tracks, FUN=function(x) {names(x)[length(x)]="delta"; return(x)})
Diffs<-get.all.diffs(Tracks, calibration, min.max.values=min.max.values, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, cluster=cluster)
graphics::plot(Diffs~sapply(Tracks, "[[", i=3))
print(Diffs)
Res<-diffs.ll(Diffs, log=FALSE)
cat("LL:", Res, "\n")
return(Res)
}


test.deltas<-function(params, Tracks, Spline, calibration, min.max.values=c(1, 1150), log.light.borders=log(c(2, 1100)), log.irrad.borders=c(-15, 50), cluster=NULL) {
# this function should send a spline as a calibration function..
# params are parameters for delta...
print(params)

lat_correction_fun<-approxfun(y= cbind(1, cos(c(-85:85)/180*pi),cos(2*c(-85:85)/180*pi),cos(3*c(-85:85)/180*pi))%*%params, x=-85:85)

# trying to add knots..
calibration$lat_correction_fun=lat_correction_fun

Diffs<-get.all.diffs(Tracks, calibration, min.max.values=min.max.values, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, cluster=cluster)
graphics::par(mfrow=c(1,2))
graphics::plot(lat_correction_fun(-85:85)~ c(-85:85))
graphics::plot(Diffs~sapply(Tracks, "[[", i=3))
Res<-diffs.ll(Diffs, log=FALSE)
cat("LL:", Res, "\n")
return(Res)
}

test.deltas3<-function(params, Tracks, calibration, min.max.values=c(1, 1150), log.light.borders=log(c(2, 1100)), log.irrad.borders=c(-15, 50), cluster=NULL) {
# this function should send a spline as a calibration function..
# params are parameters for delta...
print(params)
#deltas=params[1] + Spline%*%params[2:4] # for params - first for the intercept.

lat_correction_fun<-stats::approxfun(y= splines::bs((-85:85), degree=5, Boundary.knots=c(-85,85), intercept=TRUE)%*%params, x=-85:85)
# trying to add knots..
calibration$lat_correction_fun=lat_correction_fun

Diffs<-get.all.diffs(Tracks, calibration, min.max.values=min.max.values, log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, cluster=cluster)

Diff_index<-which(abs(mean(Diffs)-Diffs)<(5*stats::sd(Diffs))) # exclude outliers
Diffs_cor<-Diffs[Diff_index]

graphics::par(mfrow=c(1,2))
graphics::plot(lat_correction_fun(-85:85)~ c(-85:85))
graphics::plot(Diffs~sapply(Tracks, "[[", i=3))
graphics::points(Diffs_cor~sapply(Tracks, "[[", i=3)[Diff_index], col="red", pch="+")
Res<-diffs.ll(Diffs_cor, log=FALSE)
cat("LL:", Res, "\n")
return(Res)
}



