% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2DLTfunctions mixture.r
\name{h1.2b}
\alias{h1.2b}
\title{Detection hazard function \code{h1b2} prob(detect | available at x,y)}
\usage{
h1.2b(y, x, b)
}
\arguments{
\item{y}{Forward distance}

\item{x}{perpendicular distance
theta[1]*(y^2+x^2)^(-theta[2]/2).}
}
\value{
probability of detection given that an animal is availabe at location x,y
}
\description{
This hazard function has the form k(r,y)=a*r^(-b) from Hayes and Buckland (1983)
p36, with b fixed equal to 2. Note: This function uses x for perp. dist., they use y.
}
\examples{
h1(0.5,0.5,b=log(c(0.001,1)))
}
\references{
Hayes, R. J., and S. T. Buckland. "Radial-distance models for the line-transect method." Biometrics (1983): 29-42.
}
\seealso{
\code{\link{h2}}, \code{\link{ghy}}, \code{\link{ghy2}}
}
