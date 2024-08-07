% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial.mask.R
\name{make.grid}
\alias{make.grid}
\title{makes spatial grid}
\usage{
make.grid(
  left = -180,
  bottom = -90,
  right = 180,
  top = 90,
  distance.from.land.allowed.to.use = c(-Inf, Inf),
  distance.from.land.allowed.to.stay = c(-Inf, Inf),
  plot = TRUE,
  return.distances = FALSE,
  probability.of.staying = 0.5
)
}
\arguments{
\item{left}{- left boundary in degrees (-180 <= left <= 180)}

\item{bottom}{- lower boundary in degrees (-90 <= bottom <= 90)}

\item{right}{- right boundary in degrees (-180 <= right <= 180)}

\item{top}{- top boundary in degrees (-90 <= right <= 90)}

\item{distance.from.land.allowed.to.use}{- define how far from the shore animal could occur. Unit - km, negative values are for inland and positive for offshore directions. \code{Inf} stays for infinity}

\item{distance.from.land.allowed.to.stay}{- define how far from the shore animal could stay stationary between twilights. Unit - km, negative values are for inland and positive for offshore directions. \code{Inf} stays for infinity}

\item{plot}{show a plot of final grid.}

\item{return.distances}{- return distances to the shoreline}

\item{probability.of.staying}{- assigned probability value for grid cells that do not satisfy \code{distance.from.water.allowed.to.stay}}
}
\value{
dataframe with coordinates(lon and lat) and \code{probability.of.staying}
}
\description{
This function makes a rectangular grid with use defined boundaries and probabilities of being stationary.
}
\examples{
Grid<-make.grid(left=-14, bottom=30, right=13, top=57,
  distance.from.land.allowed.to.use=c(-Inf, Inf),
  distance.from.land.allowed.to.stay=c(-Inf, Inf))

}
\author{
Eldar Rakhimberdiev
}
