# solar_functions.R
solar.tripEstimation<-function (day) {
  tm <- as.POSIXlt(day, tz = "GMT")
  hh <- tm$hour
  mm <- tm$min
  ss <- tm$sec
  julday<-function (tm) {
    yr <- 1900 + tm$year
    mon <- tm$mon + 1
    day <- tm$mday
    yr[mon <= 2] <- yr[mon <= 2] - 1
    mon[mon <= 2] <- mon[mon <= 2] + 12
    a <- floor(yr/100)
    b <- 2 - a + floor(a/4)
    floor(365.25 * (yr + 4716) + floor(30.6001 * (mon + 1))) + 
      day + b - 1524.5
  }
  jday <- julday(tm) + (hh + (mm + ss/60)/60)/24
  t <- (jday - 2451545)/36525
  M <- 357.52911 + t * (35999.05029 - 0.0001537 * t)
  eqcent <- sin(pi/180 * M) * (1.914602 - t * (0.004817 + 1.4e-05 * 
    t)) + sin(pi/180 * 2 * M) * (0.019993 - 0.000101 * t) + 
    sin(pi/180 * 3 * M) * 0.000289
  L0 <- 280.46646 + t * (36000.76983 + 0.0003032 * t)
  L0 <- L0%%360
  lambda0 <- L0 + eqcent
  omega <- 125.04 - 1934.136 * t
  lambda <- lambda0 - 0.00569 - 0.00478 * sin(pi/180 * omega)
  seconds <- 21.448 - t * (46.815 + t * (0.00059 - t * (0.001813)))
  obliq0 <- 23 + (26 + (seconds/60))/60
  omega <- 125.04 - 1934.136 * t
  obliq <- obliq0 + 0.00256 * cos(pi/180 * omega)
  e <- 0.016708634 - t * (4.2037e-05 + 1.267e-07 * t)
  y <- tan(pi/180 * obliq/2)^2
  eqtime <- 180/pi * 4 * (y * sin(pi/180 * 2 * L0) - 2 * e * 
    sin(pi/180 * M) + 4 * e * y * sin(pi/180 * M) * cos(pi/180 * 
    2 * L0) - 0.5 * y^2 * sin(pi/180 * 4 * L0) - 1.25 * e^2 * 
    sin(pi/180 * 2 * M))
  solarDec <- asin(sin(pi/180 * obliq) * sin(pi/180 * lambda))
  sinSolarDec <- sin(solarDec)
  cosSolarDec <- cos(solarDec)
  solarTime <- (hh * 60 + mm + ss/60 + eqtime)/4
  list(solarTime = solarTime, sinSolarDec = sinSolarDec, cosSolarDec = cosSolarDec)
}
solar.tripEstimation<-cmpfun(solar.tripEstimation)

elevation.tripEstimation<-function (lon, lat, sun) {
  hourAngle <- sun$solarTime + lon - 180
  cosZenith <- (sin(pi/180 * lat) * sun$sinSolarDec + cos(pi/180 * 
    lat) * sun$cosSolarDec * cos(pi/180 * hourAngle))
  cosZenith[cosZenith > 1] <- 1
  cosZenith[cosZenith < -1] <- -1
  90 - 180/pi * acos(cosZenith)
}
elevation.tripEstimation<-cmpfun(elevation.tripEstimation)


#====================
# these are a bit different solar functions based on the GeoLight's formulas


solar.GeoLight<-function(gmt) {
    n <- gmt - as.POSIXct(strptime("2000-01-01 12:00:00", "%Y-%m-%d %H:%M:%S"), 
        "UTC")
    L <- 280.46 + 0.9856474 * n
    L <- as.numeric(L)
    g <- 357.528 + 0.9856003 * n
    g <- as.numeric(g)
    t.v <- floor(g/360)
    g <- g - 360 * t.v
    g.rad <- g * pi/180
    t.l <- floor(L/360)
    L <- L - 360 * t.l
    L.rad <- L * pi/180
    LAMBDA <- L + 1.915 * sin(g.rad) + 0.02 * sin(2 * g.rad)
    LAMBDA.rad <- LAMBDA * pi/180
    epsilon <- 23.439 - 4e-07 * n
    epsilon.rad <- as.numeric(epsilon) * pi/180
    alpha.rad <- atan(cos(epsilon.rad) * sin(LAMBDA.rad)/cos(LAMBDA.rad))
    alpha.rad <- ifelse(cos(LAMBDA.rad) < 0, alpha.rad + pi, 
        alpha.rad)
    alpha <- alpha.rad * 180/pi
    deklination.rad <- asin(sin(epsilon.rad) * sin(LAMBDA.rad))
    deklination <- deklination.rad * 180/pi
    tag <- paste(format(gmt, format="%Y"), "-", format(gmt, format="%m"), "-", format(gmt, format="%d"), " 00:00:00", sep = "")
    JD0 <- as.POSIXct(strptime(tag, "%Y-%m-%d %H:%M:%S"), "UTC")
    JD0 <- JD0 - as.POSIXct(strptime("2000-01-01 12:00:00", "%Y-%m-%d %H:%M:%S"), 
        "UTC")
    T0 <- JD0/36525
	
    Time <- as.numeric(format(gmt, format="%H")) + as.numeric(format(gmt, format="%M"))/60 + as.numeric(format(gmt, format="%S"))/60/100
    theta.Gh <- 6.697376 + 2400.05134 * T0 + 1.002738 * Time
    theta.Gh <- as.numeric(theta.Gh)
    t.d <- floor(theta.Gh/24)
    theta.Gh <- theta.Gh - t.d * 24
    theta.G <- theta.Gh * 15
	Res<-list(
		theta.G = theta.G,
		alpha=alpha,
		sinSolarDec = sin(deklination.rad),
		cosSolarDec = cos(deklination.rad)
		)
	return(Res)	
	}

solar.GeoLight<-cmpfun(solar.GeoLight)

elevation.GeoLight<-function(lon, lat, solarFLightR) {
    theta <- solarFLightR$theta.G + lon
    tau <- theta - solarFLightR$alpha
    tau.rad <- tau/180 * pi
    h <- asin(solarFLightR$cosSolarDec * cos(tau.rad) * cos(lat/180 * 
        pi) + solarFLightR$sinSolarDec * sin(lat/180 * pi))
    h.grad <- h/pi * 180
    R <- 1.02/(tan((h.grad + 10.3/(h.grad + 5.11))/180 * pi))
    hR.grad <- h.grad + R/60
    return(hR.grad)
}

elevation.GeoLight<-cmpfun(elevation.GeoLight)


solar<-function(gmt, mode=c("tripEstimation", "GeoLight")) {
if (mode[1]=="tripEstimation") {
	solar.tripEstimation(gmt)
	} else {
	solar.GeoLight(gmt)
	}
}

elevation<-function(lon, lat, sun, mode=c("tripEstimation", "GeoLight")) {
if (mode[1]=="tripEstimation") {
	elevation.tripEstimation(lon, lat, sun)
	} else {
	elevation.GeoLight(lon, lat, sun)
	}
}


get.declination<-function(Dates) {

if (is.numeric(Dates[1])) Dates<-as.POSIXct(Dates, tz="UTC", origin="1970-01-01")

n=as.numeric(Dates-c(as.POSIXct("2000-01-01 12:00:00", tz="UTC")))
L=280.460+0.9856474*n
g=357.528+0.9856003*n
Lambda=(L+1.915*sin(g/180*pi)+0.020*sin(2*g/180*pi))%%360
epsilon = 23.439 - 0.0000004* n 
Dec = asin(sin(epsilon/180*pi)*sin(Lambda/180*pi))*180/pi
return(Dec)
}
