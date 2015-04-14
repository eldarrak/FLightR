############
# detect ts outliers

get.equatorial.max<-function(Proc.data, calibration, dusk=T, Twilight.ID, center=NULL) {
if (dusk) {
Twilight.solar.vector<-solar(as.POSIXct(Proc.data$Twilight.time.mat.dusk[c(1:24, 26:49), Twilight.ID], tz="gmt", origin="1970-01-01"))
Twilight.log.light.vector<-Proc.data$Twilight.log.light.mat.dusk[c(1:24, 26:49), Twilight.ID]
Twilight.time.vector=Proc.data$Twilight.time.mat.dusk[c(1:24, 26:49), Twilight.ID]
} else {
Twilight.solar.vector<-solar(as.POSIXct(Proc.data$Twilight.time.mat.dawn[c(1:24, 26:49), Twilight.ID], tz="gmt", origin="1970-01-01"))
Twilight.log.light.vector<-Proc.data$Twilight.log.light.mat.dawn[c(1:24, 26:49), Twilight.ID]
Twilight.time.vector=Proc.data$Twilight.time.mat.dawn[c(1:24, 26:49), Twilight.ID]
}

#ok let's now create a line at equator
Points.Land<-cbind(seq(-180, 180, length.out=360*2+1), 0, 1)
 if (is.null(center)) {
	Current.probs1<-	apply(Points.Land, 1, get.current.slope.prob, calibration=calibration,  Twilight.time.vector=Twilight.time.vector, Twilight.log.light.vector=Twilight.log.light.vector, plot=F, verbose=F,  log.light.borders=calibration$Parameters$log.light.borders, log.irrad.borders=calibration$Parameters$log.irrad.borders, dusk=dusk, Calib.param=calibration$Parameters$LogSlope, Twilight.solar.vector=Twilight.solar.vector, delta=NULL)

	Points.Land<-cbind(seq(Points.Land[which.max(Current.probs1),1]
	-5, Points.Land[which.max(Current.probs1),1]
	+5, length.out=100), 0, 1)
	Points.Land[,1]<-Points.Land[,1]%%360
	} else {
	Points.Land<-cbind(seq(center-30, center+30, by=0.1), 0, 1)
	Points.Land[,1]<-Points.Land[,1]%%360
	}
Current.probs1<-	apply(Points.Land, 1, get.current.slope.prob, calibration=calibration,  Twilight.time.vector=Twilight.time.vector, Twilight.log.light.vector=Twilight.log.light.vector, plot=F, verbose=F,  log.light.borders=calibration$Parameters$log.light.borders, log.irrad.borders=calibration$Parameters$log.irrad.borders, dusk=dusk, Calib.param=calibration$Parameters$LogSlope, Twilight.solar.vector=Twilight.solar.vector, delta=NULL)

Final.max=Points.Land[which.max(Current.probs1),1]
return(Final.max)
}


detect.outliers<-function(Lons, plot=T) {
require('tsoutliers')

.Lons.ts<-ts(Lons)
.Lons.ts<<-.Lons.ts
#fit <- auto.arima(Lons.ts, max.p=1, max.d=1, max.q=0)
fit <- arima(.Lons.ts, order=c(1,1,0))

otypes <- c("AO", "TC", "LS")

# I want to limit cval to allow it no more than 10% of outliers.. this could become a parameter later on..

mo2 <- locate.outliers.oloop(.Lons.ts, fit, types = otypes, maxit=15)

# if 


 if (nrow(mo2$outliers)<length(Lons)*0.075) {
 cat("adjusting cval down\n")
 Cval=3.5
 while(nrow(mo2$outliers)<length(Lons)*0.075) {
 Cval=Cval*0.95
 mo2 <- locate.outliers.oloop(.Lons.ts, fit, types = otypes, maxit=15, cval=Cval)
 } 
 cat("cval adjusted to", Cval, "\n")
 }
 
 if (nrow(mo2$outliers)>length(Lons)*0.1) {
 cat("adjusting cval up\n")
 Cval=3.5
 while(nrow(mo2$outliers)>length(Lons)*0.1) {
 Cval=Cval*1.05
 mo2 <- locate.outliers.oloop(.Lons.ts, fit, types = otypes, maxit=15, cval=Cval)
 } 
 cat("cval adjusted to", Cval, "\n")
 }

Outliers1=remove.outliers(mo2, .Lons.ts, method = "en-masse",  tsmethod.call = fit$call, cval=1)$outliers
#Outliers2=remove.outliers(mo2, .Lons.ts, method = "bottom-up",  tsmethod.call = fit$call, cval=1)$outliers
rm(".Lons.ts", envir=globalenv())

Outliers1_c<-Outliers1$ind[Outliers1$type %in% c("AO", "TC")]
if (plot) {
plot(.Lons.ts)
abline(v=Outliers1$ind, col="black")
abline(v=Outliers1$ind[Outliers1$type=="AO"], col="blue")
abline(v=Outliers1$ind[Outliers1$type=="TC"], col="brown")
}
return(Outliers1_c)
}


detect.tsoutliers<-function(calibration, Proc.data, plot=T, Threads=NULL) {

if (is.character(Proc.data)) Proc.data=get("Proc.data")
if (is.character(calibration)) calibration=get("calibration")

calibration$Parameters$log.irrad.borders<-c(-15, 5)

if (!is.null(Threads)) {
	cat("making cluster\n")
	require(parallel)
	mycl <- parallel:::makeCluster(Threads)
	tmp<-parallel:::clusterSetRNGStream(mycl)
    tmp<-parallel:::clusterExport(mycl,c("calibration", "Proc.data"), envir=environment())
    #tmp<-parallel:::clusterEvalQ(mycl, library("circular")) 
    #tmp<-parallel:::clusterEvalQ(mycl, library("truncnorm")) 
    tmp<-parallel:::clusterEvalQ(mycl, library("GeoLight")) 
    tmp<-parallel:::clusterEvalQ(mycl, library("FLightR")) 
	cat("estimating dusk errors projection on equator\n")
	
	Dusks<-1:(dim(Proc.data$Twilight.time.mat.dusk)[2])
	Lons.dusk<-parSapply(mycl, Dusks, FUN=function(x) get.equatorial.max(Proc.data, calibration, dusk=T, x))
	
	cat("estimating dawn errors projection on equator\n")
	
	Dawns<-1:(dim(Proc.data$Twilight.time.mat.dawn)[2])
	Lons.dawn<-parSapply(mycl, Dawns, FUN=function(x) get.equatorial.max(Proc.data, calibration, dusk=F, x))
	stopCluster(mycl)

} else {
	###########
	cat("estimating dusk errors projection on equator\n")
	Lons=c()
	for (i in 1:(dim(Proc.data$Twilight.time.mat.dusk)[2])) {
	cat(i, "\n")
	Lons=c(Lons, get.equatorial.max(Proc.data, calibration, dusk=T, i))
	plot(Lons)
	}

	Lons.dusk<-Lons

	cat("estimating dawn errors projection on equator\n")

	Lons=c()
	for (i in 1:(dim(Proc.data$Twilight.time.mat.dawn)[2])) {
	cat(i, "\n")
	Lons=c(Lons, get.equatorial.max(Proc.data, calibration, dusk=F, i))
	plot(Lons)
	}
	Lons.dawn<-Lons
}
# now I have to optimize Lons in order to find the best cut point
ll.cut.point<-function(delta, Lons) {
Res<-log(sd((Lons+delta)%%360))
return(Res)
}

Deltas<--180:180
Delta_cur_dusk<-Deltas[which.min(sapply(Deltas, ll.cut.point, Lons=Lons.dusk))]
Lons.dusk.orig<-Lons.dusk
Lons.dusk<-((Lons.dusk+Delta_cur_dusk)%%360)-Delta_cur_dusk

Delta_cur_dawn<-Deltas[which.min(sapply(Deltas, ll.cut.point, Lons=Lons.dawn))]
Lons.dawn.orig<-Lons.dawn

Lons.dawn<-((Lons.dawn+Delta_cur_dawn)%%360)-Delta_cur_dawn

#-------------------------------------------------
# now we need to speciall focus on the outliers - maybe they were just estimated at the wrong maximum..
loess.filter<-function(x, y, k=3) {
Loess<-loess(y~x)
Resid<-residuals(Loess)
outliers<-which(abs(Resid)>(sd(Resid)*k))
Res=list(outliers=outliers, center=predict(Loess)[outliers])
return(Res)
}

Problematic.Dusks<-loess.filter(x=Proc.data$Twilight.time.mat.dusk[25,], y=Lons.dusk, k=3)
Problematic.Dawns<-loess.filter(x=Proc.data$Twilight.time.mat.dawn[25,], y=Lons.dawn, k=3)

par(mfrow=c(1,2))
plot(x=Proc.data$Twilight.time.mat.dusk[25,], y=Lons.dusk, pch="+")
if (length(Problematic.Dusks$outliers) >0) abline(v=Proc.data$Twilight.time.mat.dusk[25,Problematic.Dusks$outliers])
plot(x=Proc.data$Twilight.time.mat.dawn[25,], y=Lons.dawn, pch="+")
if (length(Problematic.Dawns$outliers) >0) abline(v=Proc.data$Twilight.time.mat.dawn[25,Problematic.Dawns$outliers])

# for these points I'd like to rerun the estimation..
if (length(Problematic.Dusks$outliers) >0| length(Problematic.Dawns$outliers)>0) {
if (!is.null(Threads)) {
	cat("making cluster\n")
	require(parallel)
	mycl <- parallel:::makeCluster(Threads)
	tmp<-parallel:::clusterSetRNGStream(mycl)
    tmp<-parallel:::clusterExport(mycl,c("calibration", "Proc.data"), envir=environment())
    #tmp<-parallel:::clusterEvalQ(mycl, library("circular")) 
    #tmp<-parallel:::clusterEvalQ(mycl, library("truncnorm")) 
    tmp<-parallel:::clusterEvalQ(mycl, library("GeoLight")) 
    tmp<-parallel:::clusterEvalQ(mycl, library("FLightR")) 
	
	if (length(Problematic.Dusks$outliers) >0) {
	cat("reestimating dusk errors projection on equator\n")
	Dusks<-cbind(Problematic.Dusks$outliers,  ((Problematic.Dusks$center+Delta_cur_dusk)%%360)-Delta_cur_dusk)
	Lons.dusk_short<-parApply(mycl, Dusks, 1, FUN=function(x) get.equatorial.max(Proc.data, calibration, dusk=T, x[1], center=x[2]))
	Lons.dusk_short<-((Lons.dusk_short+Delta_cur_dusk)%%360)-Delta_cur_dusk
	Lons.dusk[Problematic.Dusks$outliers]<-Lons.dusk_short

	}
	
	if (length(Problematic.Dawns$outliers) >0) {
	cat("reestimating dawn errors projection on equator\n")
	Dawns<-cbind(Problematic.Dawns$outliers, Problematic.Dawns$center-Delta_cur_dawn)
	Lons.dawn_short<-parApply(mycl, Dawns, 1, FUN=function(x) get.equatorial.max(Proc.data, calibration, dusk=F, x[1], center=x[2]))
	Lons.dawn_short<-((Lons.dawn_short+Delta_cur_dawn)%%360)-Delta_cur_dawn
	Lons.dawn[Problematic.Dawns$outliers]<-Lons.dawn_short

	}
	stopCluster(mycl)
} else {
	if (length(Problematic.Dusks$outliers) >0) {
	cat("reestimating dusk errors projection on equator\n")
	Dusks<-cbind(Problematic.Dusks$outliers,  ((Problematic.Dusks$center+Delta_cur_dusk)%%360)-Delta_cur_dusk)
	Lons.dusk_short<-apply(Dusks, 1, FUN=function(x) get.equatorial.max(Proc.data, calibration, dusk=T, x[1], center=x[2]))
	Lons.dusk_short<-((Lons.dusk_short+Delta_cur_dusk)%%360)-Delta_cur_dusk
	Lons.dusk[Problematic.Dusks$outliers]<-Lons.dusk_short
	}
	
	if (length(Problematic.Dawns$outliers) >0) {
	cat("reestimating dawn errors projection on equator\n")
	Dawns<-cbind(Problematic.Dawns$outliers, Problematic.Dawns$center-Delta_cur_dawn)
	Lons.dawn_short<-apply(Dawns, 1, FUN=function(x) get.equatorial.max(Proc.data, calibration, dusk=F, x[1], center=x[2]))
	Lons.dawn_short<-((Lons.dawn_short+Delta_cur_dawn)%%360)-Delta_cur_dawn
	}
}
}
cat("detecting outliers\n")

#par(mfrow=c(2,1))
#plot(Lons.dusk~Proc.data$Twilight.time.mat.dusk[1,])
#plot(Lons.dawn~Proc.data$Twilight.time.mat.dawn[1,])
Dusk.outliers=detect.outliers(Lons=Lons.dusk, plot=F)
Dawn.outliers=detect.outliers(Lons=Lons.dawn, plot=F)
cat(length(Dusk.outliers), "detected for Dusks and", length(Dawn.outliers), "for Dawns\n" )
if (plot) {
par(mfrow=c(2,1))
plot(Lons.dusk~Proc.data$Twilight.time.mat.dusk[1,], main="Dusk")
abline(v=Proc.data$Twilight.time.mat.dusk[1,][Dusk.outliers])
plot(Lons.dawn~Proc.data$Twilight.time.mat.dawn[1,], main="Dawn")
abline(v=Proc.data$Twilight.time.mat.dawn[1,][Dawn.outliers])
}
Proc.data$Twilight.time.mat.dusk<-Proc.data$Twilight.time.mat.dusk[,-Dusk.outliers]
Proc.data$Twilight.log.light.mat.dusk<-Proc.data$Twilight.log.light.mat.dusk[,-Dusk.outliers]
Proc.data$Twilight.time.mat.dawn<-Proc.data$Twilight.time.mat.dawn[,-Dawn.outliers]
Proc.data$Twilight.log.light.mat.dawn<-Proc.data$Twilight.log.light.mat.dawn[,-Dawn.outliers]
Res<-list(Proc.data=Proc.data, Lons.dusk=Lons.dusk, Lons.dawn=Lons.dawn)
return(Res)
}
