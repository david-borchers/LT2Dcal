% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2DLTfunctions mixture.r
\name{px}
\alias{px}
\title{Numerical calculation of perpendicular detection function from a hazard}
\usage{
px(x, b, hrname, ystart, notInf = 10000, almostzero = 1/1000)
}
\arguments{
\item{x}{scale or vector; perp. distance}

\item{b}{two-element vector of hazard rate parameters, some of whihc may be logged}

\item{hrname}{hazard rate function}

\item{ystart}{max forward distance at which could possibly detect animal (see details).}

\item{notInf}{Value of hrname function evaluation, above which the hazard is treated as
infinite when integrating, and therefore detection probabiity is assumed to be 1.}

\item{almostzero}{Smallest value of forward distance used in integration, expressed 
as a fracton of ystart This is used as an ad-hoc way of avoiding divergent integrals 
due to the hazard being infinite at forward distance zero.}
}
\value{
probability of detection at x
}
\description{
Calculates the perpendicular detection function, \eqn{p(x)}, for a 
given hazard.
}
\examples{
gridx=seq(0,1,length=50)
p.x=px(gridx,b=c(-0.2876821, -2.3025851),
 hr=h2,ystart=4,nint=100)
plot(gridx,p.x,type="l",ylim=c(0,max(p.x)),
xlab="prep. distance, x",ylab="p(x)")
}
