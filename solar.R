## this is the function from tripEstimation
## look there for all details
solar<-function (day) {
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