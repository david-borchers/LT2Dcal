#'@title 2D line transect survey simulator
#'
#'@description  Generates detections and associated forward distances (\code{y}) for a population
#'of objects with perpendicular distances given by \code{x}, using the hazard function 
#'specified by \code{h.fun}, with parameters specified by \code{b}, and maximum forward 
#'detection distance \code{ystart}.
#'
#'@param x perpendicular distances at which detection probability is to be calculated
#'@param hr The hazard function (either a function or a character variable of the function name)
#'@param b parameter vector for hazard function \code{hr}.
#'@param ystart Max forward detectin distance, beyond which the detection hazard is 
#'assumed to be zero.
#'@param ny Number of equally-spaced points in the interval (0,\code{ystart}) to use to 
#'specify the hazard function. Defaults to 1000.
#'@param getIDs If TRUE, returns the row index of \code{x} of the detected objects.
#'
#'@return A data frame with elements \code{$x} and  \code{$y} being the perpendicular and 
#'forward distances to detections. If \code{getIDs}=TRUE, also returns the row index of \code{x} 
#'of the detected objects in an additional element \code{$seenID}.
#'
#'@examples
#'ystart = 1300
#'w = 1500
#'hr = "h1"
#'b = c(5,0.7)
#'pi.x = "pi.hnorm"
#'logphi = 6.6
#'L=1000
#'En = 100 # expected population size
#'dat = simpop2DLT(L,w,pi.x,logphi,En,fixed.n=TRUE) # fixed.n=TRUE fixes popn size a En
#'par(mfrow=c(1,2))
#'plot(dat,pch="+",main=paste("Population (N=",nrow(dat),")"))
#'abline(v=0,lty=2)
#'hist(dat$x,breaks=seq(-w,w,length=21),xlab="Perp. distance",main="")
#'# now do LT2D survey of these points (only need their x-values for that)
#'seen = detect2DLT(dat$x,hr,b,ystart,ny=1000)
#'layout(matrix(c(1,1,2,3),nrow=2))
#'plot(seen,xlim=c(-w,w),pch="+",main=paste(nrow(seen),"Detections"))
#'abline(v=0,lty=2)
#'hist(seen$x,breaks=seq(-w,w,length=21),xlab="Perp. distance",main="")
#'hist(seen$y,breaks=seq(0,ystart,length=21),xlab="Forward distance",main="")
#'
#'@export
detect2DLT = function(x,hr,b,ystart,ny=1000, getIDs=FALSE) {
  require(spatstat)
  ys = seq(0,ystart,length=ny) # possible y values (discrete approx of continuous possibles)
  # calculate detection probabilities for all y's for every x:
  pstuff = p.approx(ys,x,hr,b,what="all") # get p(x) and f(y|x,seen)
  p = pstuff$p
  fmat = pstuff$fmat
  nx = length(x)
  seen = which(p>=runif(nx,0,1)) # these are indices of detected animals
  n = length(seen)
  if(length(seen)>0) {
    fmat = fmat[seen,,drop=FALSE] # only keep f(y|x,seen) for detected animals
    x = x[seen] # only keep x for detected animals
    y = rep(NA,n) # to hold the sampled forward distances
    for(i in 1:n) {
      goty = FALSE
      while(!goty) { # keep sampling until have at least one y
        f = fmat[i,]*100/(mean(fmat[i,]*ystart)) # scale to get expected number points = 100 
        fxfun = approxfun(ys,f) # continuous approximation to perp dist thinning function
        intensityfun = function(x,y) return(fxfun(y))
        # now use thinned poisson process function to generate points
        window=owin(c(0,1),c(0,ystart)) # create strip window
        pp=rpoispp(intensityfun,win=window,nsim=1)
        if(pp$n>0) {
          y[i] = sample(rep(pp$y,2),1) # randomly chosen y is draw from f(y|x,seen)
          goty = TRUE
        }
      }
    }
    if(getIDs) df = data.frame(x,y,seenID=seen)
    else df = data.frame(x,y)
    return(df)
  } else {
    return(NA)
  }
}



#'@title Simulate population for 2D line transect survey
#'
#'@description  Generates locations of a population with perpendicular distance
#'distribution specified by \code{pi.x}, with parameters specified by \code{logphi}, 
#'in a strip of length \code{L} and width \code{2w}. The number of individuals in 
#'the population is a Poisson random variable with mean \code{En}, unless \code{fixed.n}
#'is TRUE, in which case exactly \code{En} individuals are simulated.
#'
#'@param L Line length (length of strip)
#'@param w Strip half-width
#'@param pi.x Perpendicular distance distribution function (character or function)
#'@param logphi Parameters of \code{pi.x}.
#'@param En Expected number of objects to simulate (or exact number if \code{fixed.n} is TRUE)
#'@param fixed.n If TRUE exactly \code{En} objects are simulated, else Poisson number
#'of objects is simulated, with expected valuue \code{En}.
#'
#'@return A data frame with elements \code{$x} and  \code{$y} being the perpendicular and 
#'along-transect distances of objects.
#'@examples
#'ystart = 1300
#'w = 1500
#'hr = "h1"
#'b = c(5,0.7)
#'pi.x = "pi.hnorm"
#'logphi = 6.6
#'L=1000
#'En = 100 # expected population size
#'dat = simpop2DLT(L,w,pi.x,logphi,En,fixed.n=TRUE) # fixed.n=TRUE fixes popn size a En
#'par(mfrow=c(1,2))
#'plot(dat,pch="+",main=paste("Population (N=",nrow(dat),")"))
#'abline(v=0,lty=2)
#'hist(dat$x,breaks=seq(-w,w,length=21),xlab="Perp. distance",main="")
#'
#'@export
simpop2DLT=function(L,w,pi.x,logphi,En,fixed.n=FALSE){
  require(spatstat)
  pix=match.fun(pi.x)
  xs = seq(-w,w,length=)
  pis = pix(xs,logphi,w)
  meanpi= mean(pis)
  px = pis*En/(meanpi*2*w*L) # scale to get expected number points = En 
  pxfun = approxfun(xs,px) # continuous approximation to perp dist thinning function
  intensityfun = function(x,y) return(pxfun(x))
  # now use thinned poisson process function to generate points
  window=owin(c(-w,w),c(0,L)) # create strip window
  pp=rpoispp(intensityfun,win=window,nsim=1)
  n = pp$n
  dat = data.frame(x=pp$x,y=pp$y)
  while(n<En & fixed.n) {
    pp=rpoispp(intensityfun,win=window,nsim=1)
    dat = rbind(dat,data.frame(x=pp$x,y=pp$y))
    n = n + pp$n
  }
  if(fixed.n) dat = dat[1:En,]
  return(dat)
}
