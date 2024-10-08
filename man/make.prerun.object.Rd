% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_preparation.R
\name{make.prerun.object}
\alias{make.prerun.object}
\title{combines data, calibration and sets up priors}
\usage{
make.prerun.object(
  Proc.data,
  Grid,
  start,
  end = start,
  Calibration,
  threads = -1,
  Decision = 0.05,
  Direction = 0,
  Kappa = 0,
  M.mean = 300,
  M.sd = 500,
  likelihood.correction = TRUE
)
}
\arguments{
\item{Proc.data}{Processed data object created by \code{\link{get.tags.data}}.}

\item{Grid}{Spatial grid created by \code{\link{make.grid}}.}

\item{start}{release location (lat, lon).}

\item{end}{end of the track location. Will use \code{start} by default. Use NA in case of unknown end point.}

\item{Calibration}{Calibration object created by \code{\link{make.calibration}}.}

\item{threads}{number of parallel threads to use. default is -1, which means FLightR will use all available threads except 1. Value 1 will force sequential evaluation}

\item{Decision}{prior for migration probability values from 0 to 1 are allowed}

\item{Direction}{Direction prior for direction of migration (in degrees) with 0 pointing to the North}

\item{Kappa}{concentration parameter for vonMises distribution, 0 means uniform or even distribution. Will set some prior for direction for all the track, so is not recommended to be changed}

\item{M.mean}{Prior for mean distance travelled between consecutive twilights, km}

\item{M.sd}{Prior for sd of distance travelled between consecutive twilights, the higher the value is the wider is the the distribution}

\item{likelihood.correction}{Should likelihood correction estimated during \code{\link{make.calibration}} run be used?}
}
\value{
Object to be uses in the \code{\link{run.particle.filter}}
}
\description{
This function is one step before \code{\link{run.particle.filter}}. It combines data, calibration, spatial extent and movement priors and estimates spatial likelihoods that used later in the particle filter.
}
\examples{
File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
# to run example fast we will cut the real data file by 2013 Aug 20
Proc.data<-get.tags.data(File, end.date=as.POSIXct('2013-07-02', tz='GMT'))
Calibration.periods<-data.frame(
       calibration.start=NA,
       calibration.stop=as.POSIXct("2013-08-20", tz='GMT'),
       lon=5.43, lat=52.93) 
       #use c() also for the geographic coordinates, if you have more than one calibration location
       # (e. g.,  lon=c(5.43, 6.00), lat=c(52.93,52.94))
print(Calibration.periods)

# NB Below likelihood.correction is set to FALSE for fast run! 
# Leave it as default TRUE for real examples
Calibration<-make.calibration(Proc.data, Calibration.periods, likelihood.correction=FALSE)

Grid<-make.grid(left=0, bottom=50, right=10, top=56,
  distance.from.land.allowed.to.use=c(-Inf, Inf),
  distance.from.land.allowed.to.stay=c(-Inf, Inf))

all.in<-make.prerun.object(Proc.data, Grid, start=c(5.43, 52.93),
                             Calibration=Calibration, threads=2)

}
\author{
Eldar Rakhimberdiev
}
