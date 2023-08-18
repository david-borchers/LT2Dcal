library(LT2D)

sim.2d <- function(n, density){
  #browser()

  # (a) initial animal distribution
  if(density==0){pi.fun.name <- 'pi.const'; logphi <- NA}  # uniform
  else if(density==1){pi.fun.name <- 'pi.chnorm'; logphi <- c(0, 6)}  # avoid
  else if(density==2){pi.fun.name <- 'pi.hnorm'; logphi <- 6.6}  # attracted
  
  # simulated positions
  pos <- simpop2DLT(L=600,w = 2000, pi.x = pi.fun.name,logphi = logphi, En = n, fixed.n = T)
  
  detect.2d <- detect2DLT(pos$x, hr = ip0, b=c(4.9, 0.036),ystart = 1700, ny = 1000)
  
  detect.2d <- detect.2d[abs(detect.2d$x)<=1600, ]
  
  return(data.frame(x=detect.2d$x, y=detect.2d$y))
}
df <- sim.2d(1200, 1)


fit.sim.2d <- function(df, density){
  #browser()
  if (density==0){pi.fun.name <- "pi.const"; logphi <- NULL}  # uniform
  else if(density==1){pi.fun.name <- "pi.chnorm"; logphi <- c(0, 6)}  # avoid
  else if (density==2){pi.fun.name <- "pi.hnorm"; logphi <- 6.6}  # attracted
  
  simDat <- df
  all.1s <- rep(1,length(simDat$x))
  obj <- 1:length(simDat$x)
  sim.df <- data.frame(x = abs(simDat$x),
                       y = simDat$y,
                       stratum = all.1s,
                       transect = all.1s,
                       L = 600,
                       area = 2*2000*600,
                       object = obj,
                       size = all.1s)
  
  try({
    fit <- LT2D.fit(DataFrameInput = sim.df,
                    hr = 'ip0',
                    b = c(4.9, 0.036),
                    ystart = 1700,
                    pi.x = pi.fun.name,
                    logphi = logphi,
                    w = 2000, 
                    hessian = TRUE)
  })
  
  return(fit$fit$par)
}
out <- fit.2d(df, 1)


sim.fit.2d <- function(b, n, density){
  browser()
  for (i in 1:b){
    df <- sim.2d(n, density)
    par <- sim.fit.2d(n, density)
    write.table(df, file = "df.csv", row.names = FALSE, 
                append = TRUE, col.names = FALSE, sep = ", ")
    write.table(as.data.frame(par), file = "par.csv", row.names = FALSE, 
                append = TRUE, col.names = FALSE, sep = ", ")
  }
}
result <- sim.fit.2d(3, 1200, 1)
