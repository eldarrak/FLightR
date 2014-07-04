# ver 5.0 will work with functions 5.0 and with the 1 minute simulation..
# ver 0.4 works without mcmc model
# ver.0.3 corrected for the shift!!
# ver 0.2 point correction possibility
# ver 0.1 13.02.2014

get.shifts<-function(Track, Parameters, log.light.borders=log(c(2,64)), log.irrad.borders=c(-9,3), Points.Land, start, ask=F, slopes.only=F, delta=NULL, short.run=F, saving.period=600, Time.seq=NULL, Time.seq.saving=NULL) {
#========================================
if (length(unique(Track[,2])) !=1) stop("moving track is not implemented yet!")
# here is the lnorm distr that we currently use..
Lnorm.param<-Parameters$LogSlope



require(FLightR)

#===================== the new version will work with general simulation function \
#data(file.path(wd, "LightR_development_code\\get.slopes.5.0.r"))

# for this we will have to create a To.run.object...
To.run<-expand.grid(Slope.ideal=Lnorm.param[1], SD.ideal=Lnorm.param[2], Latitude=unique(Track[,2])) #

Track.initial<-Track
Track<-simulate.track(saving.period=saving.period, To.run=To.run, Parameters=Parameters, short.run=short.run, Time.seq=Time.seq, Time.seq.saving=Time.seq.saving, log.light.borders=log.light.borders)
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
data(wrld_simpl)
require(fields)

cat("   reading file\n")

# set wd
Data<-geologger.read.data(file=File.name)
try(unlink(File.name))

#========================================================
#========================================================

# new trick - let's try to load the real track
cat("   GeoLight step\n")

tw <- twilightCalc(Track$gmt, Track$light, allTwilights=T, ask=F, LightThreshold=3)

#save(tw, file="tw.sim.no.move.RData")
#load("tw.sim.no.move.RData")


# so we will create a vector of points from Equator to Pole and will estimate the likelihood of out positions in that points..
#Points.Land<-cbind(start[1], seq(0, 67, 0.5))
#Points.Land<-cbind(-143, seq(50, 67, 0.5))

# now we want to solve the equation for every day..

GLtab   <- tw[[2]] # Table to proceed with GeoLight
cat("automatically detected twilight:\n")
print(table(GLtab[,3]))
 if (abs(diff(table(GLtab[,3])))>5 & ask) {
cat ("more than 5 twilights were detected incorrectly.. please do the selection by hand\n")
try(dev.off())
tw <- twilightCalc(Track$gmt, Track$light, allTwilights=T, ask=T, LightThreshold=3, nsee=1000)
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
Elevation<-getElevation(GLtab$tFirst, GLtab$tSecond, GLtab$type, known.coord=start, plot=F)
positionsGeoLight <- coord(GLtab$tFirst, GLtab$tSecond, GLtab$type, degElevation=Elevation)
positionsGeoLight<-as.data.frame(positionsGeoLight)
positionsGeoLight$tFirst<-GLtab$tFirst
#tripMap(positionsGeoLight)
# ok, good here are geolight positions.

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

#All.p<-All.p[!duplicated(All.p[,2:3], fromLast=T),]
rownames(All.p)<-1:nrow(All.p)

#=============================================================================
#================= START =====================================================
# NOW we need to add a stupid old part that will just create for the all out.object..
#processing dusk
raw.Y.dusk<-correct.hours(Filtered_tw$datetime[Filtered_tw$type==2])
#Dusk.segments<-get.discontinuities(as.numeric(Filtered_tw$datetime[Filtered_tw$type==2]), raw.Y.dusk , plot=T, p.level=0.001)

#######################################################
raw.X.dusk<-as.numeric(Filtered_tw$datetime[Filtered_tw$type==2])
Result.Dusk<-make.result.list(Data, raw.X.dusk, raw.Y.dusk)


raw.Y.dawn<-correct.hours(Filtered_tw$datetime[Filtered_tw$type==1])

#Dawn.segments<-get.discontinuities(Filtered_tw$datetime[Filtered_tw$type==1], raw.Y.dawn, p.level=0.005, plot=T)
#######################################################################

raw.X.dawn<-as.numeric(Filtered_tw$datetime[Filtered_tw$type==1])

Result.Dawn<-make.result.list(Data, raw.X.dawn, raw.Y.dawn)
#Result.Dawn<-run.segmented.lnorm.loess(Data, raw.X.dawn, raw.Y.dawn,  Segments=Dawn.segments, span.correction=15, dusk=F, window.size=9, cpus=1, maxiter=100, use.first.interval=T, plot=T)

Result.all<-list(Final.dusk=Result.Dusk, Final.dawn=Result.Dawn)
####

# ok now I want to create an output..

Index.tab<-create.proposal(Result.all, start=start, Points.Land=Points.Land)
#save(Index.tab, file="Index.tab.RData")

Index.tab$yday<-as.POSIXlt(Index.tab$Date, tz="GMT")$yday

Index.tab$Decision<-0.1 # prob of migration
Index.tab$Direction<- 0 # direction 0 - North
Index.tab$Kappa<-0 # distr concentration 0 means even
Index.tab$M.mean<- 300 # distance
Index.tab$M.sd<- 500 # distance sd


##START POINTS###
all.out<-geologger.sampler.create.arrays(Index.tab, Points.Land, start=start)

#================= END== =====================================================

#=============================================================================
#########
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

Twilight.time.mat.dusk<-Twilight.time.mat.dusk-(saving.period-60)

# now we need to load the calibration that we already have...
#load("Calibration.3.3.RData")


###############
# ok and now we should be ready for the estimates...
# 
#============== ver 4.0.c
#Lnorm.param.estimation[1]<-Lnorm.param.estimation[1]+(Lnorm.param.estimation[1]^2+exp(Parameters$LogSigma)^2)/2
#============================
# before going gor the simulation I would probably like to go for a simple estimation in just one point as this could help us a lot!
do.linear.regresion<-function(Twilight.ID, start, dusk=T, Twilight.time.mat, Twilight.log.light.mat, return.slopes=F,  Calib.param, log.irrad.borders, verbose=F, log.light.borders=log(c(2,64))) {
#=========================================================================
		#cat("doing", Twilight.ID, "\n")	
		#Twilight.solar.vector<-solar(as.POSIXct(Twilight.time.mat[c(1:24, 26:49), Twilight.ID], tz="gmt", origin="1970-01-01"))
		Twilight.log.light.vector<-Twilight.log.light.mat[c(1:24, 26:49), Twilight.ID]
		Twilight.time.vector=Twilight.time.mat[c(1:24, 26:49), Twilight.ID]
		Data<-check.boundaries(start, Twilight.solar.vector=NULL,  Twilight.log.light.vector=Twilight.log.light.vector, plot=F, verbose=verbose,  log.light.borders=log.light.borders, log.irrad.borders=log.irrad.borders, dusk=dusk, Twilight.time.vector=Twilight.time.vector)
		Res<-c(NA, NA)
		if (dim(Data)[1]>1) {
		LogLight<-Data[,1]
		LogIrrad<-Data[,2]
		
		#===========
		# here I'll check on how many points we have
			if (length(LogLight) >=2) { 
			Model<- lm(LogLight~LogIrrad)
			if (verbose) print(summary(Model))
		# now we need to get the results:
		Res[1]<-coef(Model)[2] 
		Res[2]<-summary(Model)$sigma

	}}
	return(Res)
	}


	Twilight.vector<-1:(dim(Twilight.time.mat.dusk)[2])

		Slopes.dusk<-sapply(Twilight.vector, FUN=do.linear.regresion, start=start, Twilight.log.light.mat=Twilight.log.light.mat.dusk, Twilight.time.mat=Twilight.time.mat.dusk, dusk=T, log.irrad.borders=log.irrad.borders, log.light.borders=log.light.borders)
		
	Twilight.vector<-1:(dim(Twilight.time.mat.dawn)[2])

		Slopes.dawn<-sapply(Twilight.vector, FUN=do.linear.regresion, start=start, Twilight.log.light.mat=Twilight.log.light.mat.dawn, Twilight.time.mat=Twilight.time.mat.dawn, dusk=F,log.irrad.borders=log.irrad.borders, log.light.borders=log.light.borders)

	All.probs.dusk<-c()
	All.probs.dawn<-c()
		cat("\n detected dusk", dim(Twilight.time.mat.dusk)[2], "\n")
		cat("\n detected dawn",dim(Twilight.time.mat.dawn)[2] ,"\n")

	if (!slopes.only) {
		cat("detected dusk", dim(Twilight.time.mat.dusk)[2], "\n")
		cat("detected dawn",dim(Twilight.time.mat.dawn)[2] ,"\n")
	Twilight.vector<-1:(dim(Twilight.time.mat.dusk)[2])
 
		 All.probs.dusk<-sapply(Twilight.vector, FUN=get.prob.surface, Twilight.log.light.mat=Twilight.log.light.mat.dusk, Twilight.time.mat=Twilight.time.mat.dusk, dusk=T, Calib.param=Lnorm.param, log.irrad.borders=log.irrad.borders, delta=delta, Points.Land=Points.Land, interval=saving.period, log.light.borders=log.light.borders)
		 #All.probs.dusk<-sapply(Twilight.vector, FUN=get.prob.surface, Twilight.log.light.mat=Twilight.log.light.mat.dusk, Twilight.time.mat=Twilight.time.mat.dusk, dusk=T, Calib.param=Lnorm.param, log.irrad.borders=log.irrad.borders, delta=delta, Points.Land=cbind(0,38), return.slopes=T)
	
	Twilight.vector<-1:(dim(Twilight.time.mat.dawn)[2])
		All.probs.dawn<-sapply(Twilight.vector, FUN=get.prob.surface, Twilight.log.light.mat=Twilight.log.light.mat.dawn, Twilight.time.mat=Twilight.time.mat.dawn, dusk=F,Calib.param=Lnorm.param,log.irrad.borders=log.irrad.borders, delta=delta, Points.Land=Points.Land, interval=saving.period, log.light.borders=log.light.borders)
		#All.probs.dawn<-sapply(Twilight.vector, FUN=get.prob.surface, Twilight.log.light.mat=Twilight.log.light.mat.dawn, Twilight.time.mat=Twilight.time.mat.dawn, dusk=F,Calib.param=Lnorm.param,log.irrad.borders=log.irrad.borders, delta=delta, Points.Land=Points.Land)
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
all.out$positionsGeoLight<-positionsGeoLight
		if (!slopes.only) {
		return(all.out)
		} else {
		return(Slopes)
		}
}
