
get.slopes<-function(Repeats=1, file.head="tmp", Lon=0, Lat=NULL, measurement.period=60, saving.period=NULL, To.run, Parameters=NULL, short.run=F, Time.seq=NULL, Time.seq.saving=NULL, log.light.borders=log(c(2,64)), min.max.values=c(0, 64), log.irrad.borders=c(-15, 50) , plot=T) {
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

#Points.Land<-cbind(start[1], start[2], 1)

# now we will be measuring evey 60 seconds and save max inside a period for the future.




#================================================================
# estimation
#====================================
# this is the estimation part, so it should be after the generation part..
tw <- twilightCalc(Track$gmt, Track$light, allTwilights=T, ask=F, LightThreshold=3)

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

#All.p<-All.p[!duplicated(All.p[,2:3], fromLast=T),]
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

# now I need to create Points.Land
Row<-sapply(as.numeric(Filtered_tw$datetime), FUN=function(x) which.min(abs(x-Track$Time.seq)))


# ok and here we need to get around somehow...
	Twilight.vector<-1:(dim(Twilight.time.mat.dusk)[2])
#Points.Land<-cbind(start[1],start[2],1)

All.probs.dusk<-c()
for (Twilight.ID in Twilight.vector) {
	Lon<-Track$Lon[Row[Filtered_tw$type==2][Twilight.ID]]
	Lat<-Track$Lat[Row[Filtered_tw$type==2][Twilight.ID]]
	Points.Land<-cbind(Lon, Lat)
	
	Prob.surf<-try(get.prob.surface(Twilight.ID=Twilight.ID, Twilight.log.light.mat=Twilight.log.light.mat.dusk, Twilight.time.mat=Twilight.time.mat.dusk, dusk=T, return.slopes=T, log.irrad.borders=log.irrad.borders, Calib.param=Parameters$LogSlope, Points.Land=Points.Land, delta=0, log.light.borders=log.light.borders))
	#print(str(Prob.surf))
	All.probs.dusk<-cbind(All.probs.dusk, Prob.surf)
 }

	#	 All.probs.dusk<-sapply(Twilight.vector, FUN=get.prob.surface, cor.model=Gam2, Twilight.log.light.mat=Twilight.log.light.mat.dusk, Twilight.time.mat=Twilight.time.mat.dusk, dusk=T, correct.points=correct.points, return.slopes=T, log.irrad.borders=c(-100, 100))
#	cat("minimum dusk duration:", min(All.probs.dusk[4,]), "\n" )

	Twilight.vector<-1:(dim(Twilight.time.mat.dawn)[2])

	All.probs.dawn<-c()
	for (Twilight.ID in Twilight.vector) {
	Lon<-Track$Lon[Row[Filtered_tw$type==1][Twilight.ID]]
	Lat<-Track$Lat[Row[Filtered_tw$type==1][Twilight.ID]]
	#print(Points.Land)
	Points.Land<-cbind(Lon, Lat)
	Prob.surf<-try(get.prob.surface(Twilight.ID=Twilight.ID, Twilight.log.light.mat=Twilight.log.light.mat.dawn, Twilight.time.mat=Twilight.time.mat.dawn, dusk=F, return.slopes=T, log.irrad.borders=log.irrad.borders, Calib.param=Parameters$LogSlope, Points.Land=Points.Land, delta=0, log.light.borders=log.light.borders))
	#print(str(Prob.surf))
	All.probs.dawn<-cbind(All.probs.dawn, Prob.surf)
	}
	#All.probs.dawn<-sapply(Twilight.vector, FUN=get.prob.surface, cor.model=Gam2, Twilight.log.light.mat=Twilight.log.light.mat.dawn, Twilight.time.mat=Twilight.time.mat.dawn, dusk=F, correct.points=correct.points, return.slopes=T, log.irrad.borders=c(-100, 100))
#cat("minimum dawn duration:", min(All.probs.dawn[4,]), "\n" )
# ok, what should we do now?
# first I'd found a maximum of the probabilities for each day and comapre it with the means 

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

#str(All.slopes)
All.slopes<-as.data.frame(All.slopes)
names(All.slopes)[1:2]<-c("Slope", "Slope.sd")
All.slopes$gmt<-as.numeric(Filtered_tw$datetime)
All.slopes$type<-as.numeric(Filtered_tw$type)


All.slopes<-cbind(All.slopes, Track[Row, c("Slope.ideal", "SD.ideal", "Lat")])
#All.slopes<-All.slopes[-which(All.slopes$Duration<saving.period),]
#plot(All.slopes$Slope~All.slopes$Duration)
All.slopes$Slope<- log(All.slopes$Slope)
if (plot) plot(All.slopes$Slope~All.slopes$Slope.ideal)
#lm(All.slopes$Slope~All.slopes$Slope.ideal)
#Gam1<-gam(Slope.ideal~s(Slope, Slope.sd), data=All.slopes)

#predict(Gam1, newdata=data.frame(Slope=0.2, Slope.sd=0.4))

#Gam2<-gam(SD.ideal~s(Slope, Slope.sd), data=All.slopes)

#plot(Gam2)
#predict(Gam2, newdata=data.frame(Slope=0.2, Slope.sd=0.4))



#plot(All.slopes$Slope.sd~All.slopes$SD.ideal)
#lm(All.slopes$Slope.SD~All.slopes$SD.ideal)

All.slope.runs<-rbind(All.slope.runs, All.slopes)
#save(All.slope.runs, file=paste(file.head, "All.slope.runs.RData", sep="."))
if (plot)  {
par(mfrow=c(1,2))
plot((All.slope.runs$Slope)~All.slope.runs$Slope.ideal)
mean((All.slope.runs$Slope), na.rm=T)

plot(Slope~Lat, data=All.slope.runs)
}
}
return(All.slope.runs)
}
# and now we want to save that... 


simulate.track<-function(measurement.period=60, saving.period=600, To.run, Parameters=Parameters, short.run=F, Time.seq=NULL, Time.seq.saving=NULL, Lon=0, min.max.values=c(0, 64), first.date="2010-01-01 00:00:00", last.date="2010-03-20 23:59:59", plot=T) {
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
Index<-sample.int(nrow(To.run), max(Track$Day), replace=T)


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

require(FLightR)
# step 2. generating light curve.
# Angles
cat("   Estimating solar angles\n")
cat("Time...")
Time<-as.POSIXct(Track[,3], tz="gmt", origin="1970-01-01")
cat("Solar...")
S<-solar(Time)
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


 if (!short.run & plot) plot(Track$light[5000:6000], type="b", pch=".")

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


require(compiler)
simulate.track<-cmpfun(simulate.track)
get.slopes<-cmpfun(get.slopes)
