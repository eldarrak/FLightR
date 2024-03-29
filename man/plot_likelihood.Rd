% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/new_plotting_functions.R
\name{plot_likelihood}
\alias{plot_likelihood}
\title{plot likelihood surface over map}
\usage{
plot_likelihood(object, date = NULL, twilight.index = NULL)
}
\arguments{
\item{object}{either output from \code{\link{make.prerun.object}} or \code{\link{run.particle.filter}}}

\item{date}{either NULL or a date (possibly with time) closest to the twilight you wan to be plotted}

\item{twilight.index}{number of likelihood surface to be plotted}
}
\value{
'NULL'
}
\description{
plots specific likelihood surface over map
}
\details{
function plots likelihoods before particle filter run, so these are pure results of calibrations without any movement model
}
\examples{
File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
# to run example fast we will cut the real data file by 2013 Aug 20
Proc.data<-get.tags.data(File, end.date=as.POSIXct('2013-07-02', tz='GMT'))
Calibration.periods<-data.frame(
       calibration.start=as.POSIXct(c(NA, "2014-05-05"), tz='GMT'),
       calibration.stop=as.POSIXct(c("2013-08-20", NA), tz='GMT'),
       lon=5.43, lat=52.93) 
       #use c() also for the geographic coordinates, if you have more than one calibration location
       # (e. g.,  lon=c(5.43, 6.00), lat=c(52.93,52.94))

# NB Below likelihood.correction is set to FALSE for fast run! 
# Leave it as default TRUE for real examples
Calibration<-make.calibration(Proc.data, Calibration.periods, likelihood.correction=FALSE)

Grid<-make.grid(left=0, bottom=50, right=10, top=56,
  distance.from.land.allowed.to.use=c(-Inf, Inf),
  distance.from.land.allowed.to.stay=c(-Inf, Inf))

all.in<-make.prerun.object(Proc.data, Grid, start=c(5.43, 52.93),
                             Calibration=Calibration, threads=2)
plot_likelihood(all.in, twilight.index=10)

}
\author{
Eldar Rakhimberdiev
}
