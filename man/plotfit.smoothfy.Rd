% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2DLTfunctions mixture.r
\name{plotfit.smoothfy}
\alias{plotfit.smoothfy}
\title{Plot smooth fitted f(y) and forward distance distribution for small xs}
\usage{
plotfit.smoothfy(
  fit,
  nclass = 12,
  nfys = 200,
  xmax = max(fit$unrounded.points.with.betas[[1]]$x),
  main = "",
  plot = TRUE
)
}
\arguments{
\item{fit}{fitted object output by \code{\link{fityx}}}

\item{nclass}{number of histogram bins to use}

\item{nfys}{number of points to use in plotting smooth}

\item{xmax}{maxumum perp. dist. to use}
}
\value{
Invisibly returns a list with these elements
hst=hst,y=ys,smfy=smfy
\code{$hst} = histogram object from call to \link{\code{hist}}, containing data for
histogram of the detections within perp. dist \code{xmax}.
\code{$y} = y values for plot of smooth f(y)
\code{$smfy} = smooth f(y)
}
\description{
Plot spline smooth of f(y) and forward distance distribution
resulting from a call of \code{\link{plotfit.y}}.
}
\details{
Plot f(y) and forward distance distribution resulting from a call of
\code{\link{fityx}}. This is a post-hoc fix of \link{\code{plotfit.y}}, which
produces f(y) with some sharp and implausible bends.
}
\seealso{
\code{\link{fityx}}
}
