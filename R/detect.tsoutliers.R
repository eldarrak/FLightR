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

Current.probs1<-	apply(Points.Land, 1, get.current.slope.prob, calibration=Calibration,  Twilight.time.vector=Twilight.time.vector, Twilight.log.light.vector=Twilight.log.light.vector, plot=F, verbose=F,  log.light.borders=Calibration$Parameters$log.light.borders, log.irrad.borders=Calibration$Parameters$log.irrad.borders, dusk=dusk, Calib.param=Calibration$Parameters$LogSlope, Twilight.solar.vector=Twilight.solar.vector, delta=NULL)

Points.Land<-cbind(seq(Points.Land[which.max(Current.probs1),1]
-5, Points.Land[which.max(Current.probs1),1]
+5, length.out=100), 0, 1)

Current.probs1<-	apply(Points.Land, 1, get.current.slope.prob, calibration=Calibration,  Twilight.time.vector=Twilight.time.vector, Twilight.log.light.vector=Twilight.log.light.vector, plot=F, verbose=F,  log.light.borders=Calibration$Parameters$log.light.borders, log.irrad.borders=Calibration$Parameters$log.irrad.borders, dusk=dusk, Calib.param=Calibration$Parameters$LogSlope, Twilight.solar.vector=Twilight.solar.vector, delta=NULL)

Final.max=Points.Land[which.max(Current.probs1),1]
return(Final.max)
}


detect.outliers<-function(Lons, plot=T) {
require('tsoutliers')
Lons.ts<-ts(Lons)
fit <- arima(Lons.ts, order = c(1, 1, 0))
resid <- residuals(fit)
pars <- coefs2poly(fit)
otypes <- c("AO", "TC", "LS")
mo2 <- locate.outliers.oloop(Lons.ts, fit, types = otypes, maxit =10)
if (plot) {
plot(Lons.ts)
abline(v=mo2$outliers$ind, col="black")
abline(v=mo2$outliers$ind[mo2$outliers$type=="AO"], col="blue")
abline(v=mo2$outliers$ind[mo2$outliers$type=="TC"], col="brown")
}
return(mo2$outliers$ind[mo2$outliers$type %in% c("AO", "TC")])
}


detect.tsoutliers<-function(Calibration, Proc.data, plot=T) {

Calibration$Parameters$log.irrad.borders<-c(-1000, 1000)
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

cat("detecting outliers\n")

#par(mfrow=c(2,1))
#plot(Lons.dusk~Proc.data$Twilight.time.mat.dusk[1,])
#plot(Lons.dawn~Proc.data$Twilight.time.mat.dawn[1,])
Dusk.outliers=detect.outliers(Lons.dusk, plot=F)
Dawn.outliers=detect.outliers(Lons.dawn, plot=F)
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
return(Proc.data)
}
