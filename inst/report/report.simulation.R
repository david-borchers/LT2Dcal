library(LT2D)
library(circular)
library(Distance)
library(fields)

###-----------------------------------------------------------------------------
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
  
  detect.2d <- detect2DLT(pos$x, hr = ip0, b=c(4.9, 0.036), ystart = 1700, ny = 1000)
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
  
  ys <- seq(0, 1700, length.out=100)  
  obs2.probs <- p.approx(ys, df$x[df$obs==2], ip0, b=c(4.9, 0.036), what = "px")
  df$detect[df$obs==2] <- rbinom(n, 1, obs2.probs)  # second observer detection
  
  df$detect[abs(df$x)>1600] <- 0
  # return dataset
  df$keep <- rep(!(df$detect[df$obs==1] == 0 & df$detect[df$obs==2] == 0), each = 2)
  df <- subset(df, df$keep==TRUE)
  df[ , c("angle", "keep")] <- list(NULL)
  return(df)
}


###-----------------------------------------------------------------------------
sim.mismatch <- function(df){
  #browser()
  
  ## 1. get distance between every pair of detected objects
  df1 <- subset(df, df$obs==1 & df$detect==1)[ , c("x", "y")]  
  df2 <- subset(df, df$obs==2 & df$detect==1)[ , c("x", "y")]  
  dist.pair <- as.data.frame(rdist(df1, df2))
  
  dist.pair$unique <- 1:nrow(dist.pair)
  df1$id <- 1:nrow(df1); df1$detect <- df1$obs <- rep(1, nrow(df1))
  
  ## 2. using min distance to decide mismatching
  for (i in 1:(ncol(dist.pair)-1)) {
    min.index <- which.min(dist.pair[, i])
    if (length(min.index)>0){  
      if (dist.pair[min.index, i] < 300) {detect2 <- 1}
      else if (dist.pair[min.index, i] > 1000) {detect2 <- 0}
      else{detect2 <- rbinom(1, 1, p.approx(ys <- seq(0, 1700, length.out=100), dist.pair[min.index, i], ip0, b=c(6, 0.000005), what = "px"))}
    }else{detect2 <- 0}  # if no obs1 detection to match
    
    if (detect2==1){  # if matched
      df1[nrow(df1) + 1, ] <- c(df2$x[i], df2$y[i], dist.pair$unique[min.index], 2, 1)
      dist.pair <- dist.pair[-min.index, ]
    }else{  # if no match
      id <- max(unique(df1$id))+1
      df1[nrow(df1) + 1, ] <- c(df2$x[i], df2$y[i], id, 1, 0)
      df1[nrow(df1) + 1, ] <- c(df2$x[i], df2$y[i], id, 2, 1)}}
  
  for (i in dist.pair$unique){ 
    df1[nrow(df1) + 1, ] <- c(df1$x[df1$id==i], df1$y[df1$id==i], i, 2, 0)}
  
  ## 3. return new dataset
  df1 <- df1[order(df1$id), ]
  return(df1)
}


###-----------------------------------------------------------------------------
chapman.mr <- function(df, mismatch){
  if (mismatch==TRUE){df <- sim.mismatch(df)}
  
  S1 <- nrow(df[df$obs==1 & df$detect==1, ])  # first occasion
  S2 <- nrow(df[df$obs==2 & df$detect==1, ])  # second occasion
  B <- df$detect[df$obs==1]==1 & df$detect[df$obs==2]==1  
  B <- length(B[B==TRUE])  # caught by both occasions
  N.hat <- (S1+1)*(S2+1)/(B+1)-1  # abundance estimate
  var.N <- (S1+1)*(S2+1)*(S1-B)*(S2-B)/(((B+1)^2)*(B+2))
  d <- exp(1.96*sqrt(log(1+(var.N/(N.hat^2)))))
  lcl <- N.hat/d; ucl <- N.hat*d
  return(c(N.hat, lcl, ucl))
}


###-----------------------------------------------------------------------------
ds.analysis <- function(df){
  #browser()
  df1 <- subset(df, df$obs==1 & df$detect==1)
  n <- nrow(df1)
  new_df1 <- data.frame(Region.Label = rep(1, n), Area = rep(1.92e+09, n), 
                        Sample.Label = rep(1, n), Effort = rep(600000, n),
                        distance = abs(df1$x))
  df.ds <- ds(new_df1, truncation=1600, transect="line", key="hr", order=0, monotonicity = "none")

  return(c(df.ds$ddf$Nhat, df.ds$dht$individuals$N$lcl, df.ds$dht$individuals$N$ucl))
}


# -------------------------------------------------------------------------
fit.mrds <- function(df, mismatch){
  if(mismatch){
    df <- sim.mismatch(df)
    names(df) <- c( "x","y","object","observer","detected")
  }else{
    names(df) <- c("object","observer", "x","y","forw.dist","detected")
    }
  df$distance <- abs(df$x)
  df$distance[df$distance > 1600] <- 1600
  df$Region.Label = rep(1,dim(df)[1])
  df$Sample.Label = rep(1,dim(df)[1])
  model <- ddf(method = "io", dsmodel =~cds(key ="hr"),
               mrmodel =~glm(link = "logit", formula = ~distance),
               data = df, meta.data = list(width = 1600), control = list(refit = T, nrefit = 5, debug = T))
  ests <- dht(model, region.table = data.frame(Region.Label = 1, Area = 600000*1600*2),
              sample.table = data.frame(Region.Label = 1, Sample.Label = 1,Effort = 600000))

  N <- ests$individuals$N$Estimate
  lci <- ests$individuals$N$lcl
  uci <- ests$individuals$N$ucl
  return(c(N, lci, uci))
}


# -------------------------------------------------------------------------
fit.2d <- function(df, density){
  if (density==0){pi.fun.name <- "pi.const"; logphi <- NA}  # uniform
  else if(density==1){pi.fun.name <- "pi.chnorm"; logphi <- c(0, 6)}  # avoid
  else if (density==2){pi.fun.name <- "pi.hnorm"; logphi <- 6.5}  # attracted
  simDat <- df[df$obs == 1 & df$detect == 1,]
  all.1s <- rep(1,length(simDat$x))
  obj <- 1:length(simDat$x)
  sim.df <- data.frame(x = simDat$x,
                       y = simDat$forw.dist,
                       stratum = all.1s,
                       transect = all.1s,
                       L = 600,
                       area = 2*1600*600,
                       object = obj,
                       size = all.1s)
  fit <- LT2D.fit(DataFrameInput = sim.df,
                  hr = 'ip0',
                  b = c(4.9, 0.036),
                  ystart = 1700,
                  pi.x = pi.fun.name,
                  logphi = logphi,
                  w = 1600,
                  hessian = TRUE)
  est <- fit$ests[nrow(fit$ests),ncol(fit$ests)]
  lci<-uci<- NA
  try({boot <- LT2D.bootstrap(fit)
  lci <- boot$ci[1]
  uci <- boot$ci[2]
  }, silent = TRUE)
  output <- c(est, lci, uci)
  names(output) = NULL
  return(output)
}


###-----------------------------------------------------------------------------
simulation <- function(n=400, b=99, density, move, mismatch){
  #browser()
  input <- rep(n, b)
  df <- lapply(input, sim.data, density, move)
  
  ests.mr <- lapply(df, chapman.mr, mismatch)
  ests.ds <- lapply(df, ds.analysis)
  ests.mrds <- lapply(df, fit.mrds, mismatch)
  ests.2d <- lapply(df, fit.2d, density)
  
  df.ests.mr <- data.frame(t(sapply(ests.mr,c)))
  colnames(df.ests.mr) <- c("N.hat", "lcl", "ucl")
  df.ests.ds <- data.frame(t(sapply(ests.ds,c)))
  colnames(df.ests.ds) <- c("N.hat", "lcl", "ucl")
  df.ests.mrds <- data.frame(t(sapply(ests.mrds,c)))
  colnames(df.ests.mrds) <- c("N.hat", "lcl", "ucl")
  df.ests.2d <- data.frame(t(sapply(ests.2d,c)))
  colnames(df.ests.2d) <- c("N.hat", "lcl", "ucl")
  
  result <- function(method, n){
    sd <- sqrt(var(method$N.hat))
    bias <- mean((method$N.hat[method$N.hat < n+2*sd]-n)/n)  # mean relative bias
    
    check <- n > method$lcl[!is.na(method$lcl)] & n < method$ucl[!is.na(method$ucl)]
    cover.p <- length(check[check==TRUE])/length(check)  # coverage probability
    
    return(c(bias, cover.p))
    }
  
  list.method <- list(df.ests.mr, df.ests.ds, df.ests.mrds, df.ests.2d)
  out <- lapply(list.method, result, n)
  output <- data.frame(t(sapply(out,c)), row.names = c("MR", "DS", "MRDS", "2D"))
  colnames(output) <- c("mean relative bias", "coverage probability")

  return(output)
}


