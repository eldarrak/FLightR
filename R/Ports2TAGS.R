#' Function to write down twilights annotated in BAStag package data in so-called TAGS format
#'
#' this function converts combines twilights detected in BAStag with raw data and writes them down in TAGS format that can be easily read by \code{\link{get.tags.data}}
#' @param raw original data - dataframe with two columns first column must contain time and second measured light levels
#' @param twl twilights object from \code{preprocess.light} function
#' @param threshold threshold value used for twilight definition in \code{preprocess.light}
#' @param filename if NULL data.frame in TAGS format will be returned otherwise .csv file in TAGS format will be written
#' @return \code{NULL} if \code{filename} is provided or TAGS formatted dataframe.
#' @details TAGS format returned or written as .csv by this function is a dataframe with columns
#'  \itemize{
#'    \item{\code{datetime}}{ date and time in ISO 8601 format e.g. 2013-06-16T00:00:11.000Z}
#'   \item{\code{light}}{ light value measured by tag}
#'   \item{\code{twilight}}{ assigned by the software numeric indication of whether the record belongs to sunrise (1), sunset (2) or none of those (0)}
#'   \item{\code{excluded}}{ indication of whether a twilight was excluded during manual inspection (logical, \code{TRUE | FALSE})}
#'   \item{\code{interp}}{ indication of whether the light value at twilight was interpolated (logical, \code{TRUE | FALSE})}
#'  }
#' The fields \code{excluded} and \code{interp} may have values of \code{TRUE} only for \code{twilight > 0}. 
#' @seealso \code{\link{twGeos2TAGS}} and \code{\link{GeoLight2TAGS}}
#' 
#' @author Eldar Rakhimberdiev & Simeon Lisovski
#' @export
BAStag2TAGS <- function(raw, twl, threshold, filename=NULL) {
  if (any(is.na(c(twl[,1], twl[,2])))) stop('NA detected in twl, check what went wrong there!')
  names(raw) <- c("Twilight", "Light")
  twl$Light <- threshold
  if (!'POSIXct' %in% class(twl$Twilight)) stop ('Twilight column in twl object should have POSIX format!')

  tmp01 <- merge(raw, twl, all.y = TRUE, all.x = TRUE)
  out <- data.frame(datetime = tmp01[,1], light = tmp01[,2],
                    twilight = ifelse(tmp01$Rise==TRUE & !is.na(tmp01$Rise), 1,
                                      ifelse(tmp01$Rise==FALSE & !is.na(tmp01$Rise), 2, 0)),
                    interp = FALSE, excluded = ifelse(tmp01$Deleted & !is.na(tmp01$Deleted), TRUE, FALSE))
  out$interp[out$twilight>0] <- TRUE
  out<-out[order(out[,1]),]
  out[,1]<-format(out[,1], format="%Y-%m-%dT%H:%M:%S.000Z")
  if (!is.null(filename)) {
  utils::write.csv(out, file=paste(strsplit(filename, '.csv')[[1]][1], ".csv", sep=''), quote=FALSE, row.names=FALSE)
  return(NULL)
  } else { return(out)}
 }

#' Function to write down twilights annotated in twGeos package data in so-called TAGS format
#'
#' this function converts combines twilights detected in twGeos with raw data and writes them down in TAGS format that can be easily read by \code{\link{get.tags.data}}
#' @param raw original data - dataframe with two columns first column must contain time and second measured light levels
#' @param twl twilights object from \code{preprocess.light} function
#' @param threshold threshold value used for twilight definition in \code{preprocess.light}
#' @param filename if NULL data.frame in TAGS format will be returned otherwise .csv file in TAGS format will be written
#' @return \code{NULL} if \code{filename} is provided or TAGS formatted dataframe.
#' @details TAGS format returned or written as .csv by this function is a dataframe with columns
#'  \itemize{
#'    \item{\code{datetime}}{ date and time in ISO 8601 format e.g. 2013-06-16T00:00:11.000Z}
#'   \item{\code{light}}{ light value measured by tag}
#'   \item{\code{twilight}}{ assigned by the software numeric indication of whether the record belongs to sunrise (1), sunset (2) or none of those (0)}
#'   \item{\code{excluded}}{ indication of whether a twilight was excluded during manual inspection (logical, \code{TRUE | FALSE})}
#'   \item{\code{interp}}{ indication of whether the light value at twilight was interpolated (logical, \code{TRUE | FALSE})}
#'  }
#' The fields \code{excluded} and \code{interp} may have values of \code{TRUE} only for \code{twilight > 0}. 
#' @seealso \code{\link{BAStag2TAGS}} and \code{\link{GeoLight2TAGS}}
#'
#' @author Eldar Rakhimberdiev & Simeon Lisovski
#' @export
twGeos2TAGS <- function(raw, twl, threshold, filename=NULL) {
  if (any(is.na(c(twl[,1], twl[,2])))) stop('NA detected in twl, check what went wrong there!')
  if (!'POSIXct' %in% class(twl$Twilight)) stop ('Twilight column in twl object should have POSIX format!')

  names(raw) <- c("Twilight", "Light")
  twl$Light <- threshold
  
  tmp01 <- merge(raw, twl, all.y = TRUE, all.x = TRUE)
  out <- data.frame(datetime = tmp01[,1], light = tmp01[,2],
                    twilight = ifelse(tmp01$Rise==TRUE & !is.na(tmp01$Rise), 1,
                                      ifelse(tmp01$Rise==FALSE & !is.na(tmp01$Rise), 2, 0)),
                    interp = FALSE, excluded = ifelse(tmp01$Deleted & !is.na(tmp01$Deleted), TRUE, FALSE))
  out$interp[out$twilight>0] <- TRUE
  out<-out[order(out[,1]),]
  out[,1]<-format(out[,1], format="%Y-%m-%dT%H:%M:%S.000Z")
  if (!is.null(filename)) {
  utils::write.csv(out, file=paste(strsplit(filename, '.csv')[[1]][1], ".csv", sep=''), quote=FALSE, row.names=FALSE)
  return(NULL)
  } else { return(out)}
 }
 

#' Function to write down twilights annotated in GeoLight package data in so-called TAGS format
#'
#' this function converts combines twilights detected in BAStag ot twGeos with raw data and writes them down in TAGS format that can be easily read by \code{\link{get.tags.data}}
#' @param raw original data - dataframe with two columns first column must contain time and second measured light levels
#' @param gl_twl twilights object from GeoLight
#' @param threshold threshold value used for twilight definition in GeoLight
#' @param filename if NULL data.frame in TAGS format will be returned otherwise .csv file in TAGS format will be written
#' @return \code{NULL} if \code{filename} is provided or TAGS formatted dataframe.
#' @details TAGS format returned or written as .csv by this function is a dataframe with columns
#'  \itemize{
#'    \item{\code{datetime}}{ date and time in ISO 8601 format e.g. 2013-06-16T00:00:11.000Z}
#'   \item{\code{light}}{ light value measured by tag}
#'   \item{\code{twilight}}{ assigned by the software numeric indication of whether the record belongs to sunrise (1), sunset (2) or none of those (0)}
#'   \item{\code{excluded}}{ indication of whether a twilight was excluded during manual inspection (logical, \code{TRUE | FALSE})}
#'   \item{\code{interp}}{ indication of whether the light value at twilight was interpolated (logical, \code{TRUE | FALSE})}
#'  }
#' The fields \code{excluded} and \code{interp} may have values of \code{TRUE} only for \code{twilight > 0}. 
#' @seealso \code{\link{twGeos2TAGS}} and \code{\link{BAStag2TAGS}}
#' 
#' @author Eldar Rakhimberdiev & Simeon Lisovski
#' @export
GeoLight2TAGS<-function (raw, gl_twl, threshold, filename=NULL) {
   if (any(is.na(c(gl_twl[,1], gl_twl[,2])))) stop('NA detected in gl_twl, check what went wrong there!')
   names(raw) <- c("datetime", "light")
   raw$twilight<-0
   twl <- data.frame(datetime = as.POSIXct(c(gl_twl$tFirst, 
      gl_twl$tSecond), "UTC"), twilight = c(gl_twl$type, ifelse(gl_twl$type == 
      1, 2, 1)))
   twl <- twl[!duplicated(twl$datetime),]
   twl <- twl[order(twl[, 1]), ]
   twl$light <- mean(stats::approx(x = raw$datetime, y = raw$light, xout = twl$datetime)$y, na.rm=T)

   tmp01 <- merge(raw, twl, all.y = TRUE, all.x = TRUE)

   out <- data.frame(datetime = tmp01[,1], light = tmp01[,2],
                    twilight = tmp01[,3],
                    interp = FALSE, excluded = FALSE)
  out$interp[out$twilight>0] <- TRUE
  out<-out[order(out[,1]),]
  out[,1]<-format(out[,1], format="%Y-%m-%dT%H:%M:%S.000Z")
    if (!is.null(filename)) {
	  utils::write.csv(out, file=paste(strsplit(filename, '.csv')[[1]][1], ".csv", sep=''), quote=FALSE, row.names=FALSE)
  return(NULL)
  } else { return(out)}
}
