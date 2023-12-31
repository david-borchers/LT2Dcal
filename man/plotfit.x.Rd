% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2DLTfunctions mixture.r
\name{plotfit.x}
\alias{plotfit.x}
\title{Plot fitted hazard and perpendicular denisty distribution}
\usage{
plotfit.x(
  est,
  nclass = 10,
  nint = 100,
  plot = TRUE,
  dotitle = "FALSE",
  addTruth = FALSE,
  true.pi.x = NULL,
  true.logphi = NULL,
  true.hr = NULL,
  true.b = NULL,
  N = NULL,
  pi.x.additional = NULL,
  ...
)
}
\arguments{
\item{nclass}{number of histogram classes}

\item{nint}{number of intervals in numerical integration}

\item{plot}{boolean, plot results}

\item{addTruth}{boolean add true hazard and perp. density when simulating}

\item{true.pi.x}{true perpendicular density distribution function used when simulating}

\item{true.logphi}{true perpendicular density distribution function
parameters used when simulating}

\item{true.hr}{true hazard rate function used when simulating.}

\item{true.b}{true hazard rate function parameters used when simulating.}

\item{N}{true number of animals in simulated distribution, used to calculate bias (see details).}

\item{x}{perpendicular distance observations}

\item{fit}{return from a call of \code{\link{fityx}}}

\item{true.legend}{If true (and addTruth) plots legend for true functions in bottom left.}
}
\value{
list with:
\code{$gridx} = x values used in plotting
\code{$p.xpifit} = product of perpendicular distance det probability p(x)
and perpendicular animal  distribution \eqn{\pi(x)}.
\code{$mufit} = effective strip width \eqn{\hat p}
\code{$f.xfit} =
\code{$p.xfit} =
\code{$ptot} =
\code{$p.xfit.std} =
\code{$adbn} =
\code{$N} = true number of animals in population
\code{$n} = number of seen animals
\code{$Nhat} = estimated number of animals.
\code{$bias} = \eqn{\hat N} bias.
}
\description{
Plot fitted hazard and perpendicular denisty distribution
resulting from a call of \code{\link{fityx}}.  When simulating, true hazard functions
and perpendicular density distribution can also be added to the plot.
}
\details{
When \code{N} is specified, bias in estimated number of animals ,\eqn{\hat N},
is calculated.
}
\examples{
\dontrun{
ystart=4;w=1
hr=h2; b=log(c(0.75,1))
pi.x=pi.norm; logphi=c(0.5,log(0.2))
N=50 #true number of animals
#generate some observations
simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
hr=hr,b=b,w=w,ystart=ystart)
x=simDat$locs$x; y=simDat$locs$y
est.yx=fityx(y,x,b,hr,ystart,pi.x,logphi,w)
plotdat.yx=plotfit.x(x,est.yx,addTruth=TRUE,true.logphi=logphi,true.b=b,N=N)
}
}
\seealso{
\code{\link{fityx}}
}
