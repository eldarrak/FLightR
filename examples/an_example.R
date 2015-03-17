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

	Proc.data<-process.twilights(FLightR.data, Filtered_tw, measurement.period=60, saving.period=120)

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
Calib.data.all<-logger.template.calibration(Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk, positions=start, log.light.borders=log(c(2, 64)),  log.irrad.borders=c(-100, 1000)) 


Calib.data.all<-Calib.data.all[[1]]
All.slopes<-get.calib.param(Calib.data.all, plot=F)

All.slopes$Slopes<-as.data.frame(All.slopes$Slopes)

plot(log(All.slopes$Slopes$Slope)~All.slopes$Slopes$Time)
points(log(All.slopes$Slopes$Slope)~All.slopes$Slopes$Time)
abline(v=All.slopes$Slopes$Time[27])
as.POSIXct(All.slopes$Slopes$Time[27], tx="UTC", origin="1970-01-01")
# we had to exclude dawns
# 6, 7, 70, 11, 12, 298, 322,329, 330, 331
# and dusk 1, 72
Proc.data$Twilight.time.mat.dusk<-Proc.data$Twilight.time.mat.dusk[,-c(1, 72, 122, 267)]
Proc.data$Twilight.log.light.mat.dusk<-Proc.data$Twilight.log.light.mat.dusk[,-c(1, 72, 122, 267)]
Proc.data$Twilight.time.mat.dawn<-Proc.data$Twilight.time.mat.dawn[,-c(6, 7,8, 10, 11, 12, 298, 301, 322,329, 330, 331)]
Proc.data$Twilight.log.light.mat.dawn<-Proc.data$Twilight.log.light.mat.dawn[,-c(6, 7,8, 10, 11, 12, 298, 301, 322,329, 330, 331)]

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

#install_github("eldarrak/FLightR@towards_0.4")

#time_shift=-64
time_shift=0
Calib.data.all<-logger.template.calibration(Twilight.time.mat.Calib.dawn+time_shift, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk+time_shift, Twilight.log.light.mat.Calib.dusk, positions=start,log.light.borders=log(c(2, 64)),  log.irrad.borders=c(-50, 50)) 
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
Parameters$min.max.values<-c(1, 64)
Parameters$log.irrad.borders=c(-50, 50)
Parameters$start=start

Times.start=Sys.time()
Res<-get.1.minute.parameters(parameters=Parameters, measurement.period=Parameters$measurement.period, saving.period=Parameters$saving.period, position=Parameters$start, start.time=min(All.slopes$Slopes$Time), end.time=max(All.slopes$Slopes$Time), log.light.borders=Parameters$log.light.borders, min.max.values=Parameters$min.max.values, repeats=25, log.irrad.borders=Parameters$log.irrad.borders)

Times.end=Sys.time()
print(Times.end-Times.start) #6 minutes for not optimized function not too bad


Slope.1.minute<-(Parameters$LogSlope[1]-Res[1,1])/Res[1,2] # 

SD.slope.1.minute<- (Parameters$LogSlope[2]-Res[2,1])/Res[2,2] #
print(Slope.1.minute)
print(SD.slope.1.minute)


Parameters$LogSlope_1_minute<-c(Slope.1.minute, SD.slope.1.minute)

Times.start=Sys.time()


# and now instead of this part we want to go for a new part optimization part
# we assume that it will work best..


Res_time_correction<-get.time.correction.function(parameters=Parameters, measurement.period=Parameters$measurement.period, saving.period=Parameters$saving.period, position=Parameters$start, log.light.borders=Parameters$ log.light.borders, min.max.values=Parameters$min.max.values ,log.irrad.borders=Parameters$log.irrad.borders, Repeats=3)

Times.end=Sys.time()
print(Times.end-Times.start) #6 minutes fr not optimized function not too bad

#========



Res_time_correction$time_correction_fun(1)
Res_time_correction$time_correction_fun(0.9)

hist(Res_time_correction$results$Slope)

sigma=exp(Parameters$LogSigma[1])
#LogSlope_1_minute=c(Slope.1.minute, SD.slope.1.minute)
#Parameters$LogSlope_1_minute<-LogSlope_1_minute



Calibration<-list(Parameters=Parameters, time_correction_fun=Res_time_correction$time_correction_fun)

#############
save(Calibration, file="Calibration.TRES.749.FLightR.0.3.first2weeks.RData")
save(Proc.data, Filtered_tw, Data, file="Proc.data.Filtered_tw.TRES.749.FLightR.0.3.RData")
load("Proc.data.Filtered_tw.Data.TRES.749.FLightR.0.3.RData")



