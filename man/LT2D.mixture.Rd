% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2DLTfunctions mixture.r
\name{LT2D.mixture}
\alias{LT2D.mixture}
\title{Fitting function for mixture model LT2D}
\usage{
LT2D.mixture(
  DataFrameInput,
  hr,
  b,
  ystart,
  pi.x,
  logphi1,
  logphi2,
  w,
  lambda,
  formulas = NULL,
  ipars = NULL,
  xpars = NULL,
  ypars = NULL,
  control = list(),
  hessian = FALSE,
  corrFlag = 0.7,
  debug = FALSE
)
}
\arguments{
\item{DataFrameInput, }{data.frame ; data frame with required columns
 stratum (stratum number of observation), transect (transect number of
observation, object(object number of detection, NA if transect with no
detections), size (group size of detection), area (stratum area of observation)
and L (transect length of transect).}

\item{hr}{character; name of hazard rate to fit}

\item{b}{numeric; vector of start parameters for hr}

\item{ystart}{numeric; furthest possible forwards distance at which we can
detect animals.}

\item{pi.x}{character; name of perpendicular density to fit}

\item{w}{numeric; perpendicular truncation distance}

\item{hessian}{boolean}

\item{corrFlag}{numeric; value above which correlation flag is raised}

\item{logphi}{numeric; vector of start parameters for pi.x}

\item{rmin}{numeric; radial distance below which all values were rounded down to 0}
}
\description{
The equivalent function to \code{\link{LT2D.fit}}, when a 2
component mixture model is desired for the perpendicular density.
}
