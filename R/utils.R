#'@title Calculates perp. dist. detection function approximation
#'
#'@description  Jitters the perpendicular distances (\code{x}s) within \code{xcut} of 
#'zero and the forward distances (\code{y}s) within \code{ycut} of 0. Does this by 
#'adding a random angle from a uniform random variable on (0,\code{anground.degree}) 
#'degrees to the angle of the \code{x}s within \code{xcut} of zero and subtracting 
#'a uniform random variable on (0,\code{anground.degree}) degrees from the angle of 
#'the \code{y}s within \code{ycut} of zero.
#'
#'@param x Observed perpendicular distance data
#'@param y Observed forward distance data
#'@param xcut Perpendicular distance within which rounding is assumed to have occurred
#'@param ycut Forward distance within which rounding is assumed to have occurred
#'@param anground.degree Max angle to jitter x- and y-data
#'@param seed Random number seed for reproducibility. If NULL, random seed is generated
#'
#'@return Data frame with elements \code{$x} and \code{$y} being the new (jittered) 
#'perpendicular and forward distance observations, respectively.
#'@examples
#'xs = seq(0,w,length=100)
#'ys = seq(0,ystart,length=100)
#'h.fun.name = "h1"
#'b=c(-7.3287948, 0.9945317)
#'h.fun = match.fun(h.fun.name) # make h.fun the function specified via a character variable
#'p.vals = p.approx(ys,xs,h.fun,b) # detection function values to plot
#'plot(ys,p.vals,type='l',xlab='Perp. distance, x',ylab=expression(p(x)))
#'
#'@export
jitterzeros = function(x,y,xcut,ycut,anground.degree,seed=NULL) {
  anground = anground.degree*pi/180 # angle rounding range in radians
  xadj0 = xobs
  yadj0 = yobs
  robs = sqrt(xobs^2 + yobs^2)
  angobs = atan(xobs/yobs)
  if(is.null(seed)) seed = round(runif(1,1,10000))
  set.seed(seed)
  xrounds = which(xobs<10) # xs with angles within 10m
  yrounds = which(yobs<10) # ys with angles within 10m
  nxadj = length(xrounds)
  nyadj = length(yrounds)
  xadj0[xrounds] = robs[xrounds] * sin(angobs[xrounds] + runif(nxadj,0,anground))
  yadj0[yrounds] = robs[yrounds] * cos(angobs[yrounds] - runif(nyadj,0,anground))
  return(data.frame(x=xadj0, y=yadj0))
}
