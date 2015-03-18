# in this script will show basic usage of a current version of a package

# and I will do it with an example of tree swallow data.

# first we need to install the package
# via devtools


library(devtools)
install_github("eldarrak/FLightR@towards_0.4") # this should be changed to the master at some point

library(FLightR)
library(GeoLight)

#===================================================
# for no I will assume that one use the TAGS service (http://tags.animalmigration.org) and saved light-twilight data from there...

 require(RCurl)

  # read script lines from website
  text <- getURL("https://raw.githubusercontent.com/eldarrak/FLightr/towards_0.4/examples/749.csv", ssl.verifypeer = FALSE, followlocation = TRUE)
  lig.raw<-read.csv(text=text, stringsAsFactors =F)

  #lig.raw<-read.csv("749.csv", stringsAsFactors =F)
	
	FLightR.data<-read.tags.light.twilight(lig.raw, end.date="2012-05-12")

	Proc.data<-process.twilights(FLightR.data$Data, FLightR.data$twilights, measurement.period=60, saving.period=120)

str(Proc.data)


start=c(-80.46,	42.62) # where the bird was captured...

# calibration
## PART 2 NEW Calibration

#####################
## first we need to select days in the beginning or in the end of the calibration..
## let's first estimates slopes for the whole year assuming 
## that tag was in the starting point and look at the graph


## Dusk
Twilight.time.mat.Calib.dusk<-Proc.data$Twilight.time.mat.dusk

Twilight.log.light.mat.Calib.dusk<-Proc.data$Twilight.log.light.mat.dusk

## Dawn

Twilight.time.mat.Calib.dawn<-Proc.data$Twilight.time.mat.dawn

Twilight.log.light.mat.Calib.dawn<-Proc.data$Twilight.log.light.mat.dawn


#####################################
# in the next screens you will see twilights..
# it is very important to go through them and check on whther there are some weird ones that should be excluded..
# weird in this case means seriously nonlinear - such that linear regression over them will be seriously biased...
# in the scrip window you will number of dawn or dusk..
# write down numbers that you want to exclude and add them at the next step..
Calib.data.all<-logger.template.calibration(Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk, positions=start, log.light.borders=log(c(2, 64)),  log.irrad.borders=c(-100, 1000), plot=T) 


Calib.data.all<-Calib.data.all[[1]]
All.slopes<-get.calib.param(Calib.data.all, plot=F)

All.slopes$Slopes<-as.data.frame(All.slopes$Slopes)

plot(log(All.slopes$Slopes$Slope)~All.slopes$Slopes$Time)
points(log(All.slopes$Slopes$Slope)~All.slopes$Slopes$Time)
abline(v=All.slopes$Slopes$Time[27])
as.POSIXct(All.slopes$Slopes$Time[27], tx="UTC", origin="1970-01-01")

# we want  to exclude following twilights

Dusk_to_exclude<-c(1, 72, 122, 267)
Dawn_to_exclude<-c(6, 7,8, 10, 11, 12, 298, 301, 322,329, 330, 331)
Proc.data$Twilight.time.mat.dusk<-Proc.data$Twilight.time.mat.dusk[,-Dusk_to_exclude]
Proc.data$Twilight.log.light.mat.dusk<-Proc.data$Twilight.log.light.mat.dusk[,-Dusk_to_exclude]
Proc.data$Twilight.time.mat.dawn<-Proc.data$Twilight.time.mat.dawn[,-Dawn_to_exclude]
Proc.data$Twilight.log.light.mat.dawn<-Proc.data$Twilight.log.light.mat.dawn[,-Dawn_to_exclude]


Joint_ind<-sort(c(which(FLightR.data$twilights$type==1)[-Dawn_to_exclude], which(FLightR.data$twilights$type==2)[-Dusk_to_exclude]))
FLightR.data$twilights<-FLightR.data$twilights[Joint_ind,]

# the start end the end a very different!
# the question is what to use for the calibration
# I think we should try the first values.
# though during the initial calibration we used the last two weeks..

# so let's try to wotrk with the firs tweeks.

# let's look at the values before the twilight 27
#"2011-07-12 11:13:49 CEST"
Dusk.calib.days<-which(as.POSIXct(Proc.data$Twilight.time.mat.dusk[1,], tz="UTC", origin="1970-01-01") < as.POSIXct("2011-07-12", tz="UTC"))

Twilight.time.mat.Calib.dusk<-Proc.data$Twilight.time.mat.dusk[,Dusk.calib.days]

Twilight.log.light.mat.Calib.dusk<-Proc.data$Twilight.log.light.mat.dusk[,Dusk.calib.days]

## Dawn

Dawn.calib.days<-which(as.POSIXct(Proc.data$Twilight.time.mat.dawn[1,], tz="UTC", origin="1970-01-01") < as.POSIXct("2011-07-12", tz="UTC") )

Twilight.time.mat.Calib.dawn<-Proc.data$Twilight.time.mat.dawn[,Dawn.calib.days]
Twilight.log.light.mat.Calib.dawn<-Proc.data$Twilight.log.light.mat.dawn[,Dawn.calib.days ]



Calib.data.all<-logger.template.calibration(Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk, positions=start,log.light.borders=log(c(2, 64)),  log.irrad.borders=c(-50, 50), plot=F) 
Calib.data.all<-Calib.data.all[[1]]

All.slopes<-get.calib.param(Calib.data.all, plot=T)

plot(log(All.slopes$Slopes$Slope)~All.slopes$Slopes$Time)

# it looks liek there is an outlier but I'll leave it in for now, as it looks absolutely correct..


#=======================================================
# step 2
#=======================================================
Parameters<-All.slopes$Parameters
Parameters$measurement.period<-Proc.data$measurement.period 
Parameters$saving.period<-Proc.data$saving.period 
Parameters$log.light.borders<-log(c(2, 64))
Parameters$min.max.values<-c(min(FLightR.data$Data$light), max(FLightR.data$Data$light))
Parameters$log.irrad.borders=c(-50, 50)
Parameters$start=start

############################
# we assume that process happens on the measurement scale of 1 minute, so we want to estimate what will be the parameters of the process on this 1 minute scale

require(GeoLight)
Res<-get.1.minute.parameters(parameters=Parameters, measurement.period=Parameters$measurement.period, saving.period=Parameters$saving.period, position=Parameters$start, start.time=min(All.slopes$Slopes$Time), end.time=max(All.slopes$Slopes$Time), log.light.borders=Parameters$log.light.borders, min.max.values=Parameters$min.max.values, repeats=25, log.irrad.borders=Parameters$log.irrad.borders)

Slope.1.minute<-(Parameters$LogSlope[1]-Res[1,1])/Res[1,2] # 
SD.slope.1.minute<- (Parameters$LogSlope[2]-Res[2,1])/Res[2,2] #

Parameters$LogSlope_1_minute<-c(Slope.1.minute, SD.slope.1.minute)

# and now we want to estimate time correction function

Res_time_correction<-get.time.correction.function(parameters=Parameters, measurement.period=Parameters$measurement.period, saving.period=Parameters$saving.period, position=Parameters$start, log.light.borders=Parameters$log.light.borders, min.max.values=Parameters$min.max.values ,log.irrad.borders=Parameters$log.irrad.borders, Repeats=3)

#========
#Res_time_correction$time_correction_fun(1)
#Res_time_correction$time_correction_fun(0.9)

Calibration<-list(Parameters=Parameters, time_correction_fun=Res_time_correction$time_correction_fun)

#############
save(Calibration, file="Calibration.TRES.749.FLightR.0.3.first2weeks.RData")
save(FLightR.data, Proc.data,file="Proc.data.TRES.749.FLightR.0.3.RData")
load("Proc.data.Filtered_tw.Data.TRES.749.FLightR.0.3.RData")

#======================================
# Latitude correction function.
# and now we are coming to currently the least optimized part.. 
# though I am intensively working on it, so there easy solution will appear soon.

load(file="XXX.RData")
Calibration$lat_correction_fun<-lat_correction_fun
Calibration$Gamlatcalib<-Gamlatcalib

Calibration$lat_correction_fun(75, "2010-03-20", Calibration$Gamlatcalib)


#
#===============================================================
# Caibration is done and we start with the real data now..
########## 
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

#mypolygon <- drawPoly()  # click on the map to draw a polygon and press ESC when finished
#summary(mypolygon)
#save(mypolygon, file="mypolygon.SouthAmerica.RData")
#All.Points.Focus.sp<-SpatialPoints(All.Points.Focus)
#All.Points.Focus<-All.Points.Focus[ !is.na(over(All.Points.Focus.sp,mypolygon )),] # 10485
#rm(All.Points.Focus.sp)

Points.Land<-cbind(All.Points.Focus, Land=1)


#==================================
# making some checks - having a look at some twilight to see on whether the whole thing make any sense


Twilight.ID=1
my.golden.colors <- colorRampPalette(
			c("white","#FF7100"))
Twilight.solar.vector<-solar(as.POSIXct(Proc.data$Twilight.time.mat.dawn[c(1:24, 26:49), Twilight.ID], tz="gmt", origin="1970-01-01"))
Twilight.log.light.vector<-Proc.data$Twilight.log.light.mat.dawn[c(1:24, 26:49), Twilight.ID]
par(mfrow=c(1,2))
		Current.probs1<-	apply(Points.Land, 1, get.current.slope.prob, calibration=Calibration,  Twilight.time.vector=Proc.data$Twilight.time.mat.dawn[c(1:24, 26:49), Twilight.ID], Twilight.log.light.vector=Twilight.log.light.vector, plot=F, verbose=F,  log.light.borders=Calibration$Parameters$log.light.borders, log.irrad.borders=Calibration$Parameters$log.irrad.borders, dusk=F, Calib.param=Calibration$Parameters$LogSlope, Twilight.solar.vector=Twilight.solar.vector, delta=NULL)
		image.plot(as.image(Current.probs1, x=Points.Land[,1:2], nrow=30, ncol=30), col=my.golden.colors(64), main=paste("Dawn", Twilight.ID))
map('state',add=TRUE, lwd=1,  col=grey(0.5))
		map('world',add=TRUE, lwd=1.5,  col=grey(0.8))
abline(v=start[1])
abline(h=start[2])


# let's test dusk

# Dusk
Twilight.ID=1
Twilight.log.light.vector<-Proc.data$Twilight.log.light.mat.dusk[c(1:24, 26:49), Twilight.ID]

Twilight.solar.vector<-solar(as.POSIXct(Proc.data$Twilight.time.mat.dusk[c(1:24, 26:49), Twilight.ID], tz="gmt", origin="1970-01-01"))

Twilight.time.vector<-Proc.data$Twilight.time.mat.dusk[c(1:24, 26:49), Twilight.ID]

		Current.probs1<-	apply(Points.Land, 1, get.current.slope.prob, calibration=Calibration,  Twilight.time.vector=Proc.data$Twilight.time.mat.dusk[c(1:24, 26:49), Twilight.ID], Twilight.log.light.vector=Twilight.log.light.vector, plot=F, verbose=F,   log.light.borders=Calibration$Parameters$log.light.borders, log.irrad.borders=Calibration$Parameters$log.irrad.borders, dusk=T, Calib.param=Calibration$Parameters$LogSlope, Twilight.solar.vector=Twilight.solar.vector, delta=NULL)

		image.plot(as.image(Current.probs1, x=Points.Land[,1:2], nrow=30, ncol=30), col=my.golden.colors(64), main=paste("Dusk", Twilight.ID))
map('state',add=TRUE, lwd=1,  col=grey(0.5))
		map('world',add=TRUE, lwd=1.5,  col=grey(0.8))
abline(v=start[1])
abline(h=start[2])


################################
## ok now we want to make the rest but for this I feel I should better do on the other machine.
 
raw.Y.dusk<-correct.hours(Filtered_tw$datetime[Filtered_tw$type==2])
raw.X.dusk<-as.numeric(Filtered_tw$datetime[Filtered_tw$type==2])

Result.Dusk<-make.result.list(Data, raw.X.dusk, raw.Y.dusk)

raw.Y.dawn<-correct.hours(Filtered_tw$datetime[Filtered_tw$type==1])
raw.X.dawn<-as.numeric(Filtered_tw$datetime[Filtered_tw$type==1])
Result.Dawn<-make.result.list(Data, raw.X.dawn, raw.Y.dawn)
Result.all<-list(Final.dusk=Result.Dusk, Final.dawn=Result.Dawn)
####
library(maptools)
Index.tab<-create.proposal(Result.all, start=start, Points.Land=Points.Land)

Index.tab$yday<-as.POSIXlt(Index.tab$Date, tz="GMT")$yday
Index.tab$Decision<-0.05 # prob of migration
Index.tab$Direction<- 0 # direction 0 - North
Index.tab$Kappa<-0 # distr concentration 0 means even
Index.tab$M.mean<- 300 # distance
Index.tab$M.sd<- 500 # distance sd

all.in<-geologger.sampler.create.arrays(Index.tab, Points.Land, start=start)

#===========================================
# and now we want to add Phys.Mat

#===============================
library(parallel)
Threads= detectCores()-1
Phys.Mat<-get.Phys.Mat.parallel(all.in, Proc.data$Twilight.time.mat.dusk, Proc.data$Twilight.log.light.mat.dusk, Proc.data$Twilight.time.mat.dawn, Proc.data$Twilight.log.light.mat.dawn,  threads=Threads, calibration=Calibration)


all.in$Phys.Mat<-Phys.Mat

par(mfrow=c(3,3), ask=T)
for (t in seq(1,dim(all.in$Phys.Mat)[2], by=30)) {
# ok now I want to see how stable my estimates are.
	image.plot(as.image(apply(all.in$Phys.Mat[,t:(t+30)],1,  FUN=prod), x=all.in$Points.Land[,1:2], nrow=60, ncol=60), col=my.golden.colors(64), main=paste("twilight number", t))
library(maps)
map('world', add=T)
map('state', add=T)
abline(v=start[1])
abline(h=start[2])
abline(h=start[2]+1)

}

#=====================
# plotting to exclude outliers

par(mfrow=c(3,3))
par(ask=T)
for (i in 1:dim(all.in$Phys.Mat)[2]) {
image(as.image(all.in$Phys.Mat[,i], x=all.in$Points.Land[,1:2], nrow=30, ncol=30), col=my.golden.colors(64), main=paste( ifelse((all.in$Matrix.Index.Table$Dusk[i]), "Dusk","Dawn"), all.in$Matrix.Index.Table$Real.time[i], "i=", i))
	map('state',add=TRUE, lwd=1,  col=grey(0.5))
	map('world',add=TRUE, lwd=1.5,  col=grey(0.8))
	abline(v=all.in$Points.Land[all.in$start,1], col="grey", lwd=2)
	abline(h=all.in$Points.Land[all.in$start,2], col="grey", lwd=2)
}

# if there were any outliers it will be better to exclude them now setting all their values to 0









