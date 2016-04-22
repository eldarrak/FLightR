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
Grid<-cbind(seq(-180, 180, length.out=360*2+1), 0, 1)
 if (is.null(center)) {
	Current.probs1<-	apply(Grid, 1, get.current.slope.prob, calibration=calibration,  Twilight.time.vector=Twilight.time.vector, Twilight.log.light.vector=Twilight.log.light.vector, plot=F, verbose=F,  log.light.borders=calibration$Parameters$log.light.borders, log.irrad.borders=calibration$Parameters$log.irrad.borders, dusk=dusk, Calib.param=calibration$Parameters$LogSlope, Twilight.solar.vector=Twilight.solar.vector, delta=NULL, impute.on.boundaries=Proc.data$impute.on.boundaries)

	Grid<-cbind(seq(Grid[which.max(Current.probs1),1]
	-5, Grid[which.max(Current.probs1),1]
	+5, length.out=100), 0, 1)
	Grid[,1]<-Grid[,1]%%360
	} else {
	Grid<-cbind(seq(center-45, center+45, by=0.1), 0, 1)
	Grid[,1]<-Grid[,1]%%360
	}
Current.probs1<-	apply(Grid, 1, get.current.slope.prob, calibration=calibration,  Twilight.time.vector=Twilight.time.vector, Twilight.log.light.vector=Twilight.log.light.vector, plot=F, verbose=F,  log.light.borders=calibration$Parameters$log.light.borders, log.irrad.borders=calibration$Parameters$log.irrad.borders, dusk=dusk, Calib.param=calibration$Parameters$LogSlope, Twilight.solar.vector=Twilight.solar.vector, delta=NULL, impute.on.boundaries=Proc.data$impute.on.boundaries)

Final.max=Grid[which.max(Current.probs1),1]
return(Final.max)
}


tsoutliers.test <- function(x,plot=FALSE)
{
# this function is taken from here:
# http://stats.stackexchange.com/questions/1142/simple-algorithm-for-online-outlier-detection-of-a-generic-time-series
    x <- as.ts(x)
    if(frequency(x)>1)
        resid <- stl(x,s.window="periodic",robust=TRUE)$time.series[,3]
    else
    {
        tt <- 1:length(x)
        resid <- residuals(loess(x ~ tt))
    }
    resid.q <- quantile(resid,prob=c(0.25,0.75))
    iqr <- diff(resid.q)
    limits <- resid.q + 1.5*iqr*c(-1,1)
    score <- abs(pmin((resid-limits[1])/iqr,0) + pmax((resid - limits[2])/iqr,0))
    if(plot)
    {
        plot(x)
        x2 <- ts(rep(NA,length(x)))
        x2[score>0] <- x[score>0]
        tsp(x2) <- tsp(x)
        points(x2,pch=19,col="red")
        return(invisible(score))
    }
    else
        return(score)
}

detect.outliers<-function(Lons, plot=T, max.outlier.proportion=0.2) {
require('tsoutliers')

.Lons.ts<-ts(Lons)
.Lons.ts<<-.Lons.ts
#fit <- auto.arima(Lons.ts, max.p=1, max.d=1, max.q=0)
fit <- arima(.Lons.ts, order=c(1,1,0))

otypes <- c("AO", "TC", "LS")

# I want to limit cval to allow it no more than 10% of outliers.. this could become a parameter later on..

mo2 <- locate.outliers.oloop(.Lons.ts, fit, types = otypes, maxit=100)

# if 

 Cval=3.5
 if (nrow(mo2$outliers)<length(Lons)*0.075) {
 cat("adjusting cval down\n")
 while(nrow(mo2$outliers)<length(Lons)*0.075) {
 Cval=Cval*0.95
 mo2 <- locate.outliers.oloop(.Lons.ts, fit, types = otypes, maxit=100, cval=Cval)
 } 
 cat("cval adjusted to", Cval, "\n")
 }
 
 if (nrow(mo2$outliers)>length(Lons)*max.outlier.proportion) {
 cat("adjusting cval up\n")
 while(nrow(mo2$outliers)>length(Lons)*max.outlier.proportion) {
 Cval=Cval*1.05
 mo2 <- locate.outliers.oloop(.Lons.ts, fit, types = otypes, maxit=100, cval=Cval)
 } 
 cat("cval adjusted to", Cval, "\n")
 }

Outliers1=try(remove.outliers(mo2, .Lons.ts, method = "en-masse",  tsmethod.call = fit$call, cval=1)$outliers)

if (class(Outliers1)=="try-error") {
cat("error detected, switching detection method to bottom-up")
Outliers1=remove.outliers(mo2, .Lons.ts, method = "bottom-up",  tsmethod.call = fit$call, cval=1)$outliers
}
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


detect.tsoutliers<-function(calibration, Proc.data, plot=T, Threads=NULL, max.outlier.proportion=0.2, simple.version=F) {
if (!requireNamespace("tsoutliers", quietly = TRUE)) {
    stop("Pkg tsoutliers needed for this function to work. Please install it.",
      call. = FALSE)
  }
if (is.character(Proc.data)) Proc.data=get("Proc.data")
if (is.character(calibration)) calibration=get("calibration")

#calibration$Parameters$log.irrad.borders<-c(-4.5, 4.5)

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

if (!simple.version) {
#-------------------------------------------------
# now we need to specially focus on the outliers - maybe they were just estimated at the wrong maximum..
loess.filter<-function(x, y, k=3, exclude=NULL) {
if (is.null(exclude)) {
Loess<-loess(y~x)
} else {
y_new<-y
y_new[exclude]<-NA
x_new<-x
x_new[exclude]<-NA
Loess<-loess(y_new~x_new)
}
Predict<-predict(Loess)
Predict<-approx(y=Predict[!is.na(x)], x=x[!is.na(x)], xout=x, rule=2)$y
#Resid<-y-Predict
Resid<-residuals(Loess)
outliers<-sort(c(as.numeric(names(Resid))[(which(abs(Resid)>(sd(Resid)*k)))], exclude))
Res=list(outliers=outliers, center=Predict[outliers])
return(Res)
}

# first we want to eclude obvious outliers 
# wwe will do less then

Dusk.obvious.outliers<-which(tsoutliers.test(Lons.dusk)>1.5)
Problematic.Dusks<-loess.filter(x=Proc.data$Twilight.time.mat.dusk[25,], y=Lons.dusk, k=3, exclude=Dusk.obvious.outliers)

Dawn.obvious.outliers<-which(tsoutliers.test(Lons.dawn)>1.5)
Problematic.Dawns<-loess.filter(x=Proc.data$Twilight.time.mat.dawn[25,], y=Lons.dawn, k=3, exclude=Dawn.obvious.outliers)

par(mfrow=c(2,1))
plot(x=as.POSIXct(Proc.data$Twilight.time.mat.dusk[25,], tz="UTC", origin="1970-01-01"), y=Lons.dusk, pch="+")
if (length(Problematic.Dusks$outliers) >0) abline(v=Proc.data$Twilight.time.mat.dusk[25,Problematic.Dusks$outliers])
plot(x=as.POSIXct(Proc.data$Twilight.time.mat.dawn[25,], tz="UTC", origin="1970-01-01"), y=Lons.dawn, pch="+")
if (length(Problematic.Dawns$outliers) >0) abline(v=Proc.data$Twilight.time.mat.dawn[25,Problematic.Dawns$outliers])

# for these points I'd like to rerun the estimation..
if (length(Problematic.Dusks$outliers) >0| length(Problematic.Dawns$outliers)>0) {
	if (length(Problematic.Dusks$outliers) >0) {
	cat("reestimating dusk errors projection on equator\n")
	Dusks<-cbind(Problematic.Dusks$outliers,  ((Problematic.Dusks$center+Delta_cur_dusk)%%360)-Delta_cur_dusk)
	Lons.dusk_short<-apply(Dusks, 1, FUN=function(x) get.equatorial.max(Proc.data, calibration, dusk=T, x[1], center=x[2]))
	Lons.dusk_short<-((Lons.dusk_short+Delta_cur_dusk)%%360)-Delta_cur_dusk
	Lons.dusk[Problematic.Dusks$outliers]<-Lons.dusk_short
	}
	
	if (length(Problematic.Dawns$outliers) >0) {
	cat("reestimating dawn errors projection on equator\n")
	Dawns<-cbind(Problematic.Dawns$outliers, ((Problematic.Dawns$center+Delta_cur_dawn)%%360)-Delta_cur_dawn)
	Lons.dawn_short<-apply(Dawns, 1, FUN=function(x) get.equatorial.max(Proc.data, calibration, dusk=F, x[1], center=x[2]))
	Lons.dawn_short<-((Lons.dawn_short+Delta_cur_dawn)%%360)-Delta_cur_dawn
	Lons.dawn[Problematic.Dawns$outliers]<-Lons.dawn_short

	}

}
}
cat("detecting outliers\n")

#par(mfrow=c(2,1))
#plot(Lons.dusk~Proc.data$Twilight.time.mat.dusk[1,])
#plot(Lons.dawn~Proc.data$Twilight.time.mat.dawn[1,])
Dusk.outliers=detect.outliers(Lons=Lons.dusk, plot=F, max.outlier.proportion=max.outlier.proportion)
Dawn.outliers=detect.outliers(Lons=Lons.dawn, plot=F, max.outlier.proportion=max.outlier.proportion)
cat(length(Dusk.outliers), "detected for Dusks and", length(Dawn.outliers), "for Dawns\n" )
if (plot) {
par(mfrow=c(2,1))
plot(Lons.dusk~as.POSIXct(Proc.data$Twilight.time.mat.dusk[1,], tz="UTC", origin="1970-01-01"), main="Dusk")
if (length(Dusk.outliers)>0) abline(v=Proc.data$Twilight.time.mat.dusk[1,][Dusk.outliers])
plot(Lons.dawn~as.POSIXct(Proc.data$Twilight.time.mat.dawn[1,], tz="UTC", origin="1970-01-01"), main="Dawn")
if (length(Dawn.outliers)>0) abline(v=Proc.data$Twilight.time.mat.dawn[1,][Dawn.outliers])
}
Proc.data$Twilight.time.mat.dusk<-Proc.data$Twilight.time.mat.dusk
if (length(Dusk.outliers)>0)  Proc.data$Twilight.time.mat.dusk<-Proc.data$Twilight.time.mat.dusk[,-Dusk.outliers]
Proc.data$Twilight.log.light.mat.dusk<-Proc.data$Twilight.log.light.mat.dusk
if (length(Dusk.outliers)>0) Proc.data$Twilight.log.light.mat.dusk<-Proc.data$Twilight.log.light.mat.dusk[,-Dusk.outliers]
Proc.data$Twilight.time.mat.dawn<-Proc.data$Twilight.time.mat.dawn
if (length(Dawn.outliers)>0) Proc.data$Twilight.time.mat.dawn<-Proc.data$Twilight.time.mat.dawn[,-Dawn.outliers]
Proc.data$Twilight.log.light.mat.dawn<-Proc.data$Twilight.log.light.mat.dawn

if (length(Dawn.outliers)>0) Proc.data$Twilight.log.light.mat.dawn<-Proc.data$Twilight.log.light.mat.dawn[,-Dawn.outliers]
Res<-list(Proc.data=Proc.data, Lons.dusk=Lons.dusk, Lons.dawn=Lons.dawn)
return(Res)
}
