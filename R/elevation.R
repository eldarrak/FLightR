## this is the function from tripEstimation
## look there for all details
elevation<-function (lon, lat, sun) {
    hourAngle <- sun$solarTime + lon - 180
    cosZenith <- (sin(pi/180 * lat) * sun$sinSolarDec + cos(pi/180 * 
        lat) * sun$cosSolarDec * cos(pi/180 * hourAngle))
    cosZenith[cosZenith > 1] <- 1
    cosZenith[cosZenith < -1] <- -1
    90 - 180/pi * acos(cosZenith)
}
