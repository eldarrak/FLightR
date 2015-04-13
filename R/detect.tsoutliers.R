############
# detect ts outliers

get.equatorial.max<-function(Proc.data, calibration, dusk=T, Twilight.ID) {
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

Current.probs1<-	apply(Points.Land, 1, get.current.slope.prob, calibration=calibration,  Twilight.time.vector=Twilight.time.vector, Twilight.log.light.vector=Twilight.log.light.vector, plot=F, verbose=F,  log.light.borders=calibration$Parameters$log.light.borders, log.irrad.borders=calibration$Parameters$log.irrad.borders, dusk=dusk, Calib.param=calibration$Parameters$LogSlope, Twilight.solar.vector=Twilight.solar.vector, delta=NULL)

Points.Land<-cbind(seq(Points.Land[which.max(Current.probs1),1]
-5, Points.Land[which.max(Current.probs1),1]
+5, length.out=100), 0, 1)

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
 cat("adjusting cval\n")
 Cval=3.5
 while(nrow(mo2$outliers)<length(Lons)*0.075) {
 Cval=Cval*0.95
 mo2 <- locate.outliers.oloop(.Lons.ts, fit, types = otypes, maxit=15, cval=Cval)
 } 
 cat("cval adjusted to", Cval, "\n")
 }
 
 if (nrow(mo2$outliers)>length(Lons)*0.1) {
 cat("adjusting cval\n")
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

calibration$Parameters$log.irrad.borders<-c(-15, 50)

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
Delta_cur<-Deltas[which.min(sapply(Deltas, ll.cut.point, Lons=Lons.dusk))]
Lons.dusk_cor<-Lons.dusk+Delta_cur

Delta_cur<-Deltas[which.min(sapply(Deltas, ll.cut.point, Lons=Lons.dawn))]
Lons.dawn_cor<-Lons.dawn+Delta_cur

cat("detecting outliers\n")

#par(mfrow=c(2,1))
#plot(Lons.dusk~Proc.data$Twilight.time.mat.dusk[1,])
#plot(Lons.dawn~Proc.data$Twilight.time.mat.dawn[1,])
Dusk.outliers=detect.outliers(Lons=Lons.dusk_cor, plot=F)
Dawn.outliers=detect.outliers(Lons=Lons.dawn_cor, plot=F)
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
