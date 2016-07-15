#' @export
BAStags2TAGS <- function(raw, twl, threshold) {
  # function by Simeon
  names(raw) <- c("Twilight", "Light")
  twl$Light <- threshold
  
  tmp01 <- merge(raw, twl, all.y = T, all.x = T)
  out <- data.frame(datetime = tmp01[,1], light = tmp01[,2],
                    twilight = ifelse(tmp01$Rise==TRUE & !is.na(tmp01$Rise), 1,
                                      ifelse(tmp01$Rise==FALSE & !is.na(tmp01$Rise), 2, 0)),
                    interp = FALSE, excluded = ifelse(tmp01$Deleted & !is.na(tmp01$Deleted), TRUE, FALSE))
  out$interp[out$twilight>0] <- TRUE
  out<-out[order(out[,1]),]
  out[,1]<-format(out[,1], format="%Y-%m-%dT%H:%M:%S.000Z")
  return(out)
 }

#' @export
BAStag2TAGS <- BAStags2TAGS

#' @export
GeoLight2TAGS<-function (raw, gl_twl) {
   names(raw) <- c("datetime", "light")
   raw$twilight<-0
   twl <- data.frame(datetime = as.POSIXct(c(gl_twl$tFirst, 
      gl_twl$tSecond), "UTC"), twilight = c(gl_twl$type, ifelse(gl_twl$type == 
      1, 2, 1)))
   twl <- twl[!duplicated(twl$datetime),]
   twl <- twl[order(twl[, 1]), ]
   twl$light <- approx(x = raw$datetime, y = raw$light, xout = twl$datetime)$y

   tmp01 <- merge(raw, twl, all.y = T, all.x = T)

   out <- data.frame(datetime = tmp01[,1], light = tmp01[,2],
                    twilight = tmp01[,3],
                    interp = FALSE, excluded = FALSE)
  out$interp[out$twilight>0] <- TRUE
  out<-out[order(out[,1]),]
  out[,1]<-format(out[,1], format="%Y-%m-%dT%H:%M:%S.000Z")
  return(out)
   }
