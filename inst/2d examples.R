load("./inst/avoiddf.Rdata")
load("./inst/attractdf.Rdata")
pdf(file="./inst/ydbn.pdf",h=5,w=10)
par(mfrow=c(1,2))
hist(abs(attractdf$y))
hist(abs(avoiddf$y))
dev.off()

fit.2d <- function(df, density){
  if (density==0){pi.fun.name <- "pi.const"; logphi <- NA}  # uniform
  else if(density==1){pi.fun.name <- "pi.chnorm"; logphi <- c(0, 6)}  # avoid
  else if (density==2){pi.fun.name <- "pi.hnorm"; logphi <- 6.5}  # attracted
  simDat <- df[df$obs == 1 & df$detect == 1,]
  all.1s <- rep(1,length(simDat$x))
  obj <- 1:length(simDat$x)
  sim.df <- data.frame(x = simDat$x,
                       y = simDat$y,
                       stratum = all.1s,
                       transect = all.1s,
                       L = 600,
                       area = 2*1600*600,
                       object = obj,
                       size = all.1s)
  fit <- LT2D.fit(DataFrameInput = sim.df,
                  hr = 'ip0',
                  b = c(4.9, 0.036),
                  ystart = 1300,
                  pi.x = pi.fun.name,
                  logphi = logphi,
                  w = 1600,
                  hessian = TRUE,
                  control = list(trace=5))
  est <- fit$ests[nrow(fit$ests),ncol(fit$ests)]
  lci<-uci<- NA
  try({boot <- LT2D.bootstrap(fit)
  lci <- boot$ci[1]
  uci <- boot$ci[2]
  }, silent = TRUE)
  output <- c(est, lci, uci)
  names(output) = NULL
  return(list(fit,output))
}
#attractdf <- sim.data(800,2,0)

density=2
df = attractdf

atfit = fit.2d(attractdf, 2) #works but underestimate like all the others
avfit = fit.2d(avoiddf, 2) 

# Function used to generate data ------------------------------------------
sim.data <- function(n=400, density, move){
  #browser()
  ## PART I. Simulated data
  
  df <- data.frame(id=rep(1:n, each = 2), obs=rep(1:2, n),
                   x=NA, y=NA, forw.dist=NA, detect=NA, 
                   angle=NA)
  
  # (a) initial animal distribution
  if(density==0){pi.fun.name <- 'pi.const'; logphi <- NA}  # uniform
  else if(density==1){pi.fun.name <- 'pi.chnorm'; logphi <- c(0, 6)}  # avoid
  else if(density==2){pi.fun.name <- 'pi.hnorm'; logphi <- 6.6}  # attracted
  
  # simulated positions
  pos <- simpop2DLT(L=600,w = 2000, pi.x = pi.fun.name,logphi = logphi, En = n, fixed.n = T)
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
    
    dist <- (max(abs(df$x[df$obs==1]))-abs(df$x[df$obs==1]))*runif(n, 0, 0.1)*2.5
    df$x[df$obs==2] <- df$x[df$obs==1] + dist * cos(df$angle[df$obs==1])
    df$y[df$obs==2] <- df$y[df$obs==1] + dist * sin(df$angle[df$obs==1])/1000
  }else if(move==1){ #random
    df$angle <- rwrappedcauchy(n, mu = circular(0), rho = 0, control.circular=list(units="radian"))
    dist <- rlnorm(n,log(12) + 1.5,1.5)
    df$x[df$obs==2] <- df$x[df$obs==1] + dist * cos(df$angle[df$obs==1])
    df$y[df$obs==2] <- df$y[df$obs==1] + dist * sin(df$angle[df$obs==1])/1000
  }
  
  ys <- seq(0, 1300, length.out=100)  
  obs2.probs <- p.approx(ys, df$x[df$obs==2], ip0, b=c(4.9, 0.036), what = "px")
  df$detect[df$obs==2] <- rbinom(n, 1, obs2.probs)  # second observer detection
  
  df$detect[abs(df$x)>1600] <- 0
  # return dataset
  df$keep <- rep(!(df$detect[df$obs==1] == 0 & df$detect[df$obs==2] == 0), each = 2)
  df <- subset(df, df$keep==TRUE)
  df[ , c("angle", "keep")] <- list(NULL)
  return(df)
}


