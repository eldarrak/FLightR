
##########################
# This is an example of use of FLightR
# with the BAS tag...
# ver 0.2 from 04 Apr 2015


############################################
## PART 1. install.package and get the data#
############################################
##Load the necessary packages##
require(GeoLight)
require(maptools)
require(rgeos)
require(geosphere)
require(raster)
require(fields)

require(tsoutliers)
require(forecast)

require(circular)
require(truncnorm)
require(parallel)
require(bit)
require(aspace)
require(rgdal)
require(CircStats)

require(devtools)
install_github("eldarrak/FLightR@current_stable_version")
require(FLightR)

# for now I will assume that one use the TAGS service (http://tags.animalmigration.org) and saved light-twilight data from there...

require(RCurl)

# read script lines from website
text <- getURL("https://raw.githubusercontent.com/eldarrak/FLightr/master/examples/MEE_manuscript_example/749.csv", ssl.verifypeer = FALSE, followlocation = TRUE)
lig.raw<-read.csv(text=text, stringsAsFactors =F)
	
FLightR.data<-read.tags.light.twilight(lig.raw, end.date="2012-05-12")

Proc.data<-process.twilights(FLightR.data$Data, FLightR.data$twilights, measurement.period=60, saving.period=120)

# tag measures data every minute (measurement.period=60)
# and saves maximum over 2 minutes (saving.period=120)

start=c(-80.46,	42.62) # where the bird was captured...

################################################
## PART 2. Calibration
################################################

##----------------------------------------------
##   Search for a proper calibration period 
##   and manual check of outliers
##----------------------------------------------

## we need to select days when bird was in a known location
## these are typically days in the beginning or in the end
## of migration. To do we first will plot all sun slopes over 
## the whole period and then will decide when is our calibration period

## Dusk
Twilight.time.mat.Calib.dusk<-Proc.data$Twilight.time.mat.dusk

Twilight.log.light.mat.Calib.dusk<-Proc.data$Twilight.log.light.mat.dusk

## Dawn

Twilight.time.mat.Calib.dawn<-Proc.data$Twilight.time.mat.dawn

Twilight.log.light.mat.Calib.dawn<-Proc.data$Twilight.log.light.mat.dawn

# in the next screens you will see twilights..
# it is very important to go through them and check on whether there are some weird ones that should be excluded..
# weird in this case means seriously nonlinear - such that linear regression over them will be seriously biased...
# in the scrip window you will number of dawn or dusk..
# write down numbers that you want to exclude and add them at the next step..

# this is especially important for the nest box or cavity breeders.
# the abnormally fast twilights corresponding to twilight missed inside a cavity should be excluded..

Calib.data.all<-logger.template.calibration(Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk, positions=start, log.light.borders=log(c(2, 64)),  log.irrad.borders=c(-1000, 1000), plot=F) 

# plot=T - will plot al the twilights and you might exclude weird ones.
# plot=F will just go to the end without manual check possibility

Calib.data.all<-Calib.data.all[[1]]
All.slopes<-get.calib.param(Calib.data.all, plot=F)

All.slopes$Slopes<-as.data.frame(All.slopes$Slopes)
plot(log(All.slopes$Slopes$Slope)~All.slopes$Slopes$Time, type="n")

lines(log(Slope)~Time, data=All.slopes$Slopes[All.slopes$Slopes$Type=="Dusk",])
points(log(Slope)~Time, data=All.slopes$Slopes[All.slopes$Slopes$Type=="Dusk",], pch="+")
points(log(Slope)~Time, data=All.slopes$Slopes[All.slopes$Slopes$Type=="Dawn",], pch="+", col="red")
lines(log(Slope)~Time, data=All.slopes$Slopes[All.slopes$Slopes$Type=="Dawn",], col="red")

## we conclude that in the end of the trip our bird was around for about a month. - that's a nice period for the calibration..

as.POSIXct(All.slopes$Slopes$Time[620], tx="UTC", origin="1970-01-01")
# we assume it to be 04 Apr - as it have started to sit in a next box
abline(v=All.slopes$Slopes$Time[620])

# also we want to exclude following twilights

Dusk_to_exclude<-c(1, 72, 122, 267)
Dawn_to_exclude<-c(6, 7,8, 10, 11, 12, 298, 301, 322,329, 330, 331)
Proc.data$Twilight.time.mat.dusk<-Proc.data$Twilight.time.mat.dusk[,-Dusk_to_exclude]
Proc.data$Twilight.log.light.mat.dusk<-Proc.data$Twilight.log.light.mat.dusk[,-Dusk_to_exclude]
Proc.data$Twilight.time.mat.dawn<-Proc.data$Twilight.time.mat.dawn[,-Dawn_to_exclude]
Proc.data$Twilight.log.light.mat.dawn<-Proc.data$Twilight.log.light.mat.dawn[,-Dawn_to_exclude]

Joint_ind<-sort(c(which(FLightR.data$twilights$type==1)[-Dawn_to_exclude], which(FLightR.data$twilights$type==2)[-Dusk_to_exclude]))
FLightR.data$twilights$excluded=1
FLightR.data$twilights$excluded[Joint_ind]<-0


##----------------------------------------------
##   Calibration for the selected period
##----------------------------------------------

Dusk.calib.days<-which(as.POSIXct(Proc.data$Twilight.time.mat.dusk[1,], tz="UTC", origin="1970-01-01") >= as.POSIXct("2012-04-10", tz="UTC"))
Twilight.time.mat.Calib.dusk<-Proc.data$Twilight.time.mat.dusk[,Dusk.calib.days]
Twilight.log.light.mat.Calib.dusk<-Proc.data$Twilight.log.light.mat.dusk[,Dusk.calib.days]

Dawn.calib.days<-which(as.POSIXct(Proc.data$Twilight.time.mat.dawn[1,], tz="UTC", origin="1970-01-01") >=as.POSIXct("2012-04-10", tz="UTC") )
Twilight.time.mat.Calib.dawn<-Proc.data$Twilight.time.mat.dawn[,Dawn.calib.days]
Twilight.log.light.mat.Calib.dawn<-Proc.data$Twilight.log.light.mat.dawn[,Dawn.calib.days ]

Calib.data.all<-logger.template.calibration(Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk, positions=start,log.light.borders=log(c(2, 64)),  log.irrad.borders=c(-50, 50), plot=F) 
Calib.data.all<-Calib.data.all[[1]]

All.slopes<-get.calib.param(Calib.data.all, plot=T)

plot(log(All.slopes$Slopes$Slope)~All.slopes$Slopes$Time)

# Now we create 'parameters' object that will have all the details about the calibration
Parameters<-All.slopes$Parameters # LogSlope 0.24 0.17
Parameters$measurement.period<-Proc.data$measurement.period 
Parameters$saving.period<-Proc.data$saving.period 
Parameters$log.light.borders<-log(c(2, 64)) # these are the boundaries in which one should use the BAS tag.. for the other types of tags they will be different.
Parameters$min.max.values<-c(min(FLightR.data$Data$light), max(FLightR.data$Data$light))
Parameters$log.irrad.borders=c(-50, 50)
Parameters$start=start

# now we have to get the parameters on the 1 minute scale where we assume the process is really happening (do not forget tag saves maximum value over 2 minute period)

tmp<-optim(par=Parameters$LogSlope, get.1.minute.parameters, gr=NULL, parameters=Parameters, method="L-BFGS-B", lower=c(-3, 0), upper=c(3, Inf), position=Parameters$start, start.time=min(All.slopes$Slopes$Time), end.time=max(All.slopes$Slopes$Time), print=T)

Slope.1.minute<-tmp$par[1]
SD.slope.1.minute<-tmp$par[2]

Parameters$LogSlope_1_minute<-c(Slope.1.minute, SD.slope.1.minute)

# and now we could estimate time and latitude correction functions but I am still not sure how much they are actually needed
# as far functions after expect to have compensation functios ready we will just create a simple fake correction functions

time_correction_fun= eval(parse(text=paste("function (x) return(", Parameters$LogSlope[1], ")")))
Calibration<-list(Parameters=Parameters, time_correction_fun=time_correction_fun)

lat_correction_fun<-function(x, y, z) return(0)
Calibration$lat_correction_fun<-lat_correction_fun

#----------------------------------------------------------
# automated outlier detection:
#----------------------------------------------------------
require(parallel)
Threads=detectCores()-1
Outliers<-detect.tsoutliers(Calibration, Proc.data, plot=T, Threads=Threads)

Proc.data<-Outliers$Proc.data

FLightR.data$twilights$excluded[which(!as.numeric(FLightR.data$twilights$datetime) %in% c(Proc.data$Twilight.time.mat.dusk[25,]+Calibration$Parameters$saving.period-Calibration$Parameters$measurement.period,  Proc.data$Twilight.time.mat.dawn[25,]) & FLightR.data$twilights$excluded!=1 )]<-2

#--------------------------------------------------------
# recalibration with outliers excluded..
#--------------------------------------------------------
Dusk.calib.days<-which(as.POSIXct(Proc.data$Twilight.time.mat.dusk[1,], tz="UTC", origin="1970-01-01") >= as.POSIXct("2012-04-10", tz="UTC"))
Twilight.time.mat.Calib.dusk<-Proc.data$Twilight.time.mat.dusk[,Dusk.calib.days]
Twilight.log.light.mat.Calib.dusk<-Proc.data$Twilight.log.light.mat.dusk[,Dusk.calib.days]

Dawn.calib.days<-which(as.POSIXct(Proc.data$Twilight.time.mat.dawn[1,], tz="UTC", origin="1970-01-01") >=as.POSIXct("2012-04-10", tz="UTC") )
Twilight.time.mat.Calib.dawn<-Proc.data$Twilight.time.mat.dawn[,Dawn.calib.days]
Twilight.log.light.mat.Calib.dawn<-Proc.data$Twilight.log.light.mat.dawn[,Dawn.calib.days ]

Calib.data.all<-logger.template.calibration(Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk, positions=start,log.light.borders=log(c(2, 64)),  log.irrad.borders=c(-50, 50), plot=F) 
Calib.data.all<-Calib.data.all[[1]]

All.slopes<-get.calib.param(Calib.data.all, plot=T)

plot(log(All.slopes$Slopes$Slope)~All.slopes$Slopes$Time)

# Now we create 'parameters' object that will have all the details about the calibration
Parameters<-All.slopes$Parameters # LogSlope 0.23 0.16
Parameters$measurement.period<-Proc.data$measurement.period 
Parameters$saving.period<-Proc.data$saving.period 
Parameters$log.light.borders<-log(c(2, 64)) # these are the boundaries in which one should use the BAS tag.. for the other types of tags they will be different.
Parameters$min.max.values<-c(min(FLightR.data$Data$light), max(FLightR.data$Data$light))
Parameters$log.irrad.borders=c(-50, 50)
Parameters$start=start

# now we have to get the parameters on the 1 minute scale where we assume the process is really happening (do not forget tag saves maximum value over 2 minute period)

tmp<-optim(par=Parameters$LogSlope, get.1.minute.parameters, gr=NULL, parameters=Parameters, method="L-BFGS-B", lower=c(-3, 0), upper=c(3, Inf), position=Parameters$start, start.time=min(All.slopes$Slopes$Time), end.time=max(All.slopes$Slopes$Time), print=T)

Slope.1.minute<-tmp$par[1]
SD.slope.1.minute<-tmp$par[2]

Parameters$LogSlope_1_minute<-c(Slope.1.minute, SD.slope.1.minute)

# and now we could estimate time and latitude correction functions but I am still not sure how much they are actually needed
# as far functions after expect to have compensation functios ready we will just create a simple fake correction functions

time_correction_fun= eval(parse(text=paste("function (x) return(", Parameters$LogSlope[1], ")")))
Calibration<-list(Parameters=Parameters, time_correction_fun=time_correction_fun)

lat_correction_fun<-function(x, y, z) return(0)
Calibration$lat_correction_fun<-lat_correction_fun

###########################################################
##  Part 3. Define spatial extent for the optimisation   ##
###########################################################

# we will define just  a box, 
# but one could use any other more complicated shape

ylim = c(15, 45)
xlim = c(-92, -70)

library(geosphere)
Globe.Points<-regularCoordinates(200) # 50 km between each point

All.Points.Focus<-Globe.Points[Globe.Points[,1]>xlim[1] & Globe.Points[,1]<xlim[2] & Globe.Points[,2]>ylim[1] & Globe.Points[,2]<ylim[2],]

# here we could cut by the sea but we will not do it now

library(raster)
plot(All.Points.Focus, type="n")
map('state',add=TRUE, lwd=1,  col=grey(0.5))
map('world',add=TRUE, lwd=1.5,  col=grey(0.8))
abline(v=start[1])
abline(h=start[2])

Points.Land<-cbind(All.Points.Focus, Land=1)
# the masks work in a follwing way:
#
# Spatial constrain: if one wants to prevent animal being at this point at twilight 
# the point should just be excluded
# 
# Spatiobehavioural constrain: if you do allow to be flying at the point but not be stationary
# then set Land=0 for this point

##########################################################
## Part 4. Preestimation of all the matrices           ##
##########################################################
 
raw.Y.dusk<-correct.hours(FLightR.data$twilights$datetime[FLightR.data$twilights$type==2 & FLightR.data$twilights$excluded==0])
raw.X.dusk<-as.numeric(FLightR.data$twilights$datetime[FLightR.data$twilights$type==2 & FLightR.data$twilights$excluded==0])
Data_tmp<-list(d=FLightR.data$Data)
Result.Dusk<-make.result.list(Data_tmp, raw.X.dusk, raw.Y.dusk)

raw.Y.dawn<-correct.hours(FLightR.data$twilights$datetime[FLightR.data$twilights$type==1 & FLightR.data$twilights$excluded==0])
raw.X.dawn<-as.numeric(FLightR.data$twilights$datetime[FLightR.data$twilights$type==1 & FLightR.data$twilights$excluded==0])
Result.Dawn<-make.result.list(Data_tmp, raw.X.dawn, raw.Y.dawn)
Result.all<-list(Final.dusk=Result.Dusk, Final.dawn=Result.Dawn)
####

Index.tab<-create.proposal(Result.all, start=start, Points.Land=Points.Land)
Index.tab$yday<-as.POSIXlt(Index.tab$Date, tz="GMT")$yday
Index.tab$Decision<-0.1 # prob of migration
Index.tab$Direction<- 0 # direction 0 - North
Index.tab$Kappa<-0 # distr concentration 0 means even
Index.tab$M.mean<- 300 # distance mu
Index.tab$M.sd<- 500 # distance sd

all.in<-geologger.sampler.create.arrays(Index.tab, Points.Land, start=start)

all.in$Matrix.Index.Table$Decision<-0.1
all.in$M.mean<-300

# we also want to restrict irradiance values to c(-12, 5)
Calibration$Parameters$log.irrad.borders<-c(-12, 5)

# the next step might have some time
# with the current example it takes about 5 min at 24 core workstation

library(parallel)
Threads= detectCores()-1
Phys.Mat<-get.Phys.Mat.parallel(all.in, Proc.data$Twilight.time.mat.dusk, Proc.data$Twilight.log.light.mat.dusk, Proc.data$Twilight.time.mat.dawn, Proc.data$Twilight.log.light.mat.dawn,  threads=Threads, calibration=Calibration)

all.in$Phys.Mat<-Phys.Mat

#---------------------------
# we can also check the pre-estimated matrices:

par(mfrow=c(3,3))
par(ask=T)
my.golden.colors <- colorRampPalette(
			c("white","#FF7100"))
			
for (i in 1:dim(all.in$Phys.Mat)[2]) {
image(as.image(all.in$Phys.Mat[,i], x=all.in$Points.Land[,1:2], nrow=30, ncol=30), col=my.golden.colors(64), main=paste( ifelse((all.in$Matrix.Index.Table$Dusk[i]), "Dusk","Dawn"), all.in$Matrix.Index.Table$Real.time[i], "i=", i))
	map('state',add=TRUE, lwd=1,  col=grey(0.5))
	map('world',add=TRUE, lwd=1.5,  col=grey(0.8))
	abline(v=all.in$Points.Land[all.in$start,1], col="grey", lwd=2)
	abline(h=all.in$Points.Land[all.in$start,2], col="grey", lwd=2)
}

# if there are any outliers it will be better to exclude them now
#setting all their values to 0 e.g all.in$Phys.Mat[,c(255,303, 306, 458, 476)]<-1

# saving object...
save(all.in, file="TRES.all.in.RData")

###############################################
## Part 5. Main estimation
###############################################

# 1e6 particles take a while to run, so we highly recommend to run it first with 1e3 particles.. the result will be bad but it will check whether the PF works on yuor workstation.

# do not forget about RAM - each node will eat 2 - 2.5 Gb of it!!!

Result<-run.particle.filter(all.in, save.Res=F, cpus=min(Threads,6), nParticles=1e6, known.last=T,
 precision.sd=25, sea.value=1, save.memory=T, k=NA, parallel=T, 
 plot=T, prefix="pf", extend.prefix=T, max.kappa=100, 
 min.SD=25, min.Prob=0.01, max.Prob=0.99, 
 fixed.parameters=list(M.mean=300, M.sd=500, Kappa=0), 
 cluster.type="SOCK", a=45, b=1000, L=90, update.angle.drift=F, adaptive.resampling=0.99, save.transitions=T, check.outliers=T)
gc()

save(Result, file="result.TRES.RData")


###############################################
## Part 6. Plotting
###############################################

# 1. basic plotting could be done as follows:

pdf("FLightR_map.pdf",width=5,height=5)
par(mfrow=c(1,1))
par(mar=c(4,4,3,1),las=1,mgp=c(2.25,1,0))

#--------------------------------------------------------
# we could plot either mean or median.

#Mean_coords<-cbind(Result$Quantiles$Meanlon, Result$Quantiles$Meanlat)
Median_coords<-cbind(Result$Quantiles$MedianlonJ, Result$Quantiles$MedianlatJ)
plot(Median_coords, type = "n",ylab="Latitude",xlab="Longitude")
data(wrld_simpl)
plot(wrld_simpl, add = T, col = "grey95", border="grey70")
lines(Median_coords, col = "darkgray", cex = 0.1)
points(Median_coords, pch = 16, cex = 0.75, col = "darkgray")
#lines(Mean_coords, col = "blue", cex = 0.1)
#points(Mean_coords, pch = 16, cex = 0.75, col = "blue")
title("FLightR analysis", line = 0.3)
box("plot")
dev.off()

#######################################
##Plot FLgihtR and known path by lat and long##
#######################################
Quantiles<-Result$Quantiles[1:length(Result$Matrix.Index.Table$Real.time),]
Quantiles$Time<-Result$Matrix.Index.Table$Real.time

pdf("FLightR_lat_lon.pdf",width=5,height=5)

par(mfrow=c(2,1))
par(mar=c(2,4,3,1),cex=1)
 Sys.setlocale("LC_ALL", "English")  

 #Longitude
plot(Quantiles$Medianlon~Quantiles$Time, las=1,col=grey(0.1),pch=16,ylab="Longitude",xlab="",lwd=2, ylim=range(c( Quantiles$LCI.lon, Quantiles$UCI.lon )), type="n")


polygon(x=c(Quantiles$Time, rev(Quantiles$Time)), y=c(Quantiles$LCI.lon, rev(Quantiles$UCI.lon)), col=grey(0.9), border=grey(0.5))

polygon(x=c(Quantiles$Time, rev(Quantiles$Time)), y=c(Quantiles$TrdQu.lon, rev(Quantiles$FstQu.lon)), col=grey(0.7), border=grey(0.5))

lines(Quantiles$Medianlon~Quantiles$Time, col=grey(0.1),lwd=2)

abline(v=as.POSIXct("2013-09-22 21:34:30 EDT"), col="red", lwd=1)
abline(v=as.POSIXct("2014-03-22 21:34:30 EDT"), col="red", lwd=1)
title("FLightR analysis", line = 0.3)

#Latitude
par(mar=c(3,4,1,1))

plot(Quantiles$Medianlat~Quantiles$Time, las=1,col=grey(0.1),pch=16,ylab="Latitude",xlab="",lwd=2, ylim=range(c( Quantiles$UCI.lat, Quantiles$LCI.lat )), type="n")

polygon(x=c(Quantiles$Time, rev(Quantiles$Time)), y=c(Quantiles$LCI.lat, rev(Quantiles$UCI.lat)), col=grey(0.9), border=grey(0.5))

polygon(x=c(Quantiles$Time, rev(Quantiles$Time)), y=c(Quantiles$TrdQu.lat, rev(Quantiles$FstQu.lat)), col=grey(0.7), border=grey(0.5))

lines(Quantiles$Medianlat~Quantiles$Time, col=grey(0.1),lwd=2)

abline(v=as.POSIXct("2013-09-22 21:34:30 EDT"), col="red", lwd=1)
abline(v=as.POSIXct("2014-03-22 21:34:30 EDT"), col="red", lwd=1)

dev.off()

# there is also possibility to plot prob of migration and  migratory ditance..

#####################################










