sim.data <- function(n=400, density, move, match){
  ## PART I. Simulated data
  
  df <- data.frame(id=rep(1:n, each = 2), obs=rep(1:2, n),
                   x=NA, y=NA, forw.dist=NA, detect=NA, 
                   angle=NA)
  
  # (a) initial animal distribution
  if (density==0){pi.fun.name <- pi.const; logphi <- NA}  # uniform
  else if(density==1){pi.fun.name <- pi.chnorm; logphi <- c(0, 6)}  # avoid
  else{pi.fun.name <- pi.hnorm; logphi <- 6.5}  # attracted
  
  # simulated positions
  pos <- simpop2DLT(L=600,w = 1600, pi.x = pi.fun.name, logphi = logphi, En = n, fixed.n = T)
  pos$obs <- 0
  df$x[df$obs==1] <- pos$x  # original coordinates
  df$y[df$obs==1] <- pos$y
  
  detect.2d <- detect2DLT(pos$x, hr = ip0, b=c(4.9, 0.036), ystart = 1300, ny = 1000)
  pos$obs[pos$x %in% detect.2d$x] <- 1  
  df$detect[df$obs==1] <- pos$obs  # first observer detection
  df$forw.dist[df$obs==1 & df$detect==1] <- detect.2d$y
  
  
  # (b) animal movement
  if (move==0){  # avoidance
    nleft <- nrow(df[df$x <= 0 & df$obs==1, ])
    
    angleleft <- rwrappedcauchy(nleft, mu=circular(pi), rho=0.9, control.circular=list(units="radian"))
    df$angle[df$x <= 0 & df$obs==1] <- as.numeric(angleleft)
    angleright <- rwrappedcauchy((n-nleft), mu=circular(0), rho=0.9, control.circular=list(units="radian"))
    df$angle[df$x > 0 & df$obs==1] <- as.numeric(angleright)
    
    dist <- (max(abs(df$x[df$obs==1]))-abs(df$x[df$obs==1]))*runif(n, 0, 0.1)*5
    df$x[df$obs==2] <- df$x[df$obs==1] + dist * cos(df$angle[df$obs==1])
    df$y[df$obs==2] <- df$y[df$obs==1] + dist * sin(df$angle[df$obs==1])
    }
  else{  # attracted 
    df$x[df$obs==2] <- df$x[df$obs==1]
    df$y[df$obs==2] <- df$y[df$obs==1]  # PLEASE MODIFY THIS!
  }
  
  ys <- seq(0, 1300, length.out=100)  
  obs2.probs <- p.approx(ys, df$x[df$obs==2], ip0, b=c(4.9, 0.036), what = "px")
  df$detect[df$obs==2] <- rbinom(n, 1, obs2.probs)  # second observer detection
  df[ , c("angle")] <- list(NULL)
  
  return(df)
}


plot(df$x[df$obs==1], df$y[df$obs==1], xlim=c(-1600, 1600))
abline(v=0, col="red", lty=2)
plot(df$x[df$obs==2], df$y[df$obs==2], xlim=c(-1600, 1600), ylim = c(0, 600))

par(mfrow=c(1,2))
plot(pos,pch="+",main=paste("Population (N=",nrow(pos),")"))
abline(v=0,lty=2)
hist(pos$x,breaks=seq(-1600,1600,length=21),xlab="Perp. distance",main="")
