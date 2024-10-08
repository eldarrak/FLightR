% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_preparation.R
\name{plot_slopes_by_location}
\alias{plot_slopes_by_location}
\title{plots log of observed versus expected slope by time for a known location}
\usage{
plot_slopes_by_location(
  Proc.data,
  location,
  log.light.borders = "auto",
  log.irrad.borders = "auto",
  ylim = NULL,
  xlim = NULL
)
}
\arguments{
\item{Proc.data}{processed data object generated by \code{\link{get.tags.data}}}

\item{location}{vector with longitude and latitude of calibration location (degrees).}

\item{log.light.borders}{numeric vector with length of 2 for minimum and maximum log(light) levels to use. Default value 'auto', will take these values from the Proc.data object.}

\item{log.irrad.borders}{numeric vector with length of 2 for minimum and maximum log(irradiance) values to use. Default value 'auto', will take these values from the Proc.data object.}

\item{ylim}{the y limits of the plot. The default value, NULL, indicates that the range of the finite values to be plotted should be used.}

\item{xlim}{the x limits of the plot. The default value, NULL, otherwise can be POSIXct or character in a form readable by \code{\link{as.POSIXct}}.}
}
\value{
'NULL'
}
\description{
The function calculates and plots calibration slopes for sunsets and sunrises for every day of the tracking period, based on the assumption that the tag remained in the same (calibration) location all the time.
}
\details{
The plot of calibration slopes is used for finding start and end dates of a calibration period (the time period, during which the tag remained in the calibration location with coordinates (x,y)). During the calibration period, the calibration slopes vary little both, between the twilight events (sunrises and sunsets) and in time. When the tag changes location, the slopes for sunrises and sunsets start to deviate. There may potentially be several calibration periods for the same location (if the bird returned to the same location several times). The boundaries (start and end dates) of each of these periods are captured visually. If there were more than one calibration location, the procedure is repeated, once for each location. 
All the obtained calibration periods can be entered in a data frame 'Calibration.periods', for further analysis. Each line of the data frame contains start and end dates (if applicable) of the calibration period and geographic coordinates of the location.
}
\examples{
File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
Proc.data<-get.tags.data(File)
plot_slopes_by_location(Proc.data=Proc.data, location=c(5.43, 52.93))
abline(v=as.POSIXct("2013-08-20", tz='GMT')) # end of first calibration period
abline(v=as.POSIXct("2014-05-05", tz='GMT')) # start of the second calibration period

}
\author{
Eldar Rakhimberdiev
}
