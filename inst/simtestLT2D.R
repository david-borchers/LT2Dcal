library(tcltk2) # for progress bar function
library(LT2D)

pi.fun.name <- 'pi.chnorm'; logphi <- c(0, 6)
hr.name = "h1.2"; hr.pars=c(15, 0.25)
ystart = 2000
w=1500
par(mfrow=c(1,2))
plotdetmodel(hr.name,hr.pars,ystart,xlim=c(0,w),what="px")
plotpimodel(pi.fun.name,logphi,w)

En = 1200

En = 400


set.seed(123)
Nsim = 100
ests.ip0 = fits.ip0 = ests.h1 = fits.h1 = ests.h1.2b = fits.h1.2b = vector("list", Nsim)

# Set up progress bar, before looping:
pb <- tkProgressBar(title=paste("Simulation Progress (Nsim=",Nsim,")",sep=""), min=0, max=Nsim, width=400)

for(i in 1:Nsim) {
  pos <- simpop2DLT(L=600, w=w, pi.x=pi.fun.name, logphi=logphi, En=En, fixed.n=TRUE)
#  plot(pos,pch="+")
  detect.2d <- detect2DLT(pos$x, hr=hr.name, b=hr.pars, ystart=ystart, ny=1000, getIDs=TRUE)
#  points(detect.2d$x,pos$y[detect.2d$seen], col="red")
#  hist(detect.2d$x)
  
  simDat = detect.2d
  simDat$forw.dist = detect.2d$y
  all.1s <- rep(1,length(simDat$x))
  obj <- 1:length(simDat$x)
  sim.df <- data.frame(x = abs(simDat$x),
                       y = simDat$forw.dist,
                       stratum = all.1s,
                       transect = all.1s,
                       L = 600,
                       area = 2*2000*600,
                       object = obj,
                       size = all.1s)
  
#  fit.ip0 <- LT2D.fit(DataFrameInput = sim.df,
#                  hr = 'ip0',
#                  b = c(4.9, 0.036),
#                  ystart = 1700,
#                  pi.x = pi.fun.name,
#                  logphi = logphi,
#                  w = 2000, 
#                  hessian = TRUE)
#  ests.ip0[[i]] = fit.ip0$est
#  fits.ip0[[i]] = fit.ip0$fit
#  
#  fit.h1 <- LT2D.fit(DataFrameInput = sim.df,
#                    hr = 'h1',
#                    b = c(3, 0.5),
#                    ystart = 1700,
#                    pi.x = pi.fun.name,
#                    logphi = logphi,
#                    w = 2000, 
#                    hessian = TRUE)
#  ests.h1[[i]] = fit.h1$est
#  fits.h1[[i]] = fit.h1$fit

  # fit fixed-b model
  fit.h1.2 <- LT2D.fit(DataFrameInput = sim.df,
                       hr = 'h1.2',
                       b = hr.pars,
                       ystart = ystart,
                       pi.x = pi.fun.name,
                       logphi = logphi,
                       w = w, 
                       hessian = TRUE,
                       control=list(trace=5))
  ests.h1.2[[i]] = fit.h1.2$est
  fits.h1.2[[i]] = fit.h1.2$fit  
    # Progress, inside loop
  setTkProgressBar(pb, i, label=paste( round(i/Nsim*100, 0),"% done"))
  
}
# Close progress bar after looping
close(pb)

Nest = rep(NA,Nsim)
for(i in 1:Nsim)  Nest[i] = ests.h1.2b[[i]][2,11]

saveRDS(list(ests.h1.2b=ests.h1.2b, fits.h1.2b=fits.h1.2b),file="simests.h1.2b.Rds")


# Fit hazard rate with b>=2 AND SIMUATE FROM SAME MODEL
ests.h1.2 = fits.h1.2 = vector("list", Nsim)

# Set up progress bar, before looping:
pb <- tkProgressBar(title=paste("Simulation Progress (Nsim=",Nsim,")",sep=""), min=0, max=Nsim, width=400)

for(i in 1:Nsim) {
  pos <- simpop2DLT(L=600,w = 2000, pi.x = pi.fun.name,logphi = logphi, En = n, fixed.n = T)
  detect.2d <- detect2DLT(pos$x, hr = h1.2, b=c(4.9, -2),ystart = 1700, ny = 1000, getIDs=TRUE)

  simDat = detect.2d
  simDat$forw.dist = detect.2d$y
  all.1s <- rep(1,length(simDat$x))
  obj <- 1:length(simDat$x)
  sim.df <- data.frame(x = abs(simDat$x),
                       y = simDat$forw.dist,
                       stratum = all.1s,
                       transect = all.1s,
                       L = 600,
                       area = 2*2000*600,
                       object = obj,
                       size = all.1s)
  
  fit.h1.2 <- LT2D.fit(DataFrameInput = sim.df,
                        hr = 'h1.2',
                        b = c(6,-10),
                        ystart = 1700,
                        pi.x = pi.fun.name,
                        logphi = logphi,
                        w = 2000, 
                        hessian = TRUE)
  ests.h1.2[[i]] = fit.h1.2$est
  fits.h1.2[[i]] = fit.h1.2$fit  
  # Progress, inside loop
  setTkProgressBar(pb, i, label=paste( round(i/Nsim*100, 0),"% done"))
  
}
# Close progress bar after looping
close(pb)

plot(fit.h1.2)
gof.LT2D(fit.h1.2)

Nest = rep(NA,Nsim)
for(i in 1:Nsim)  Nest[i] = ests.h1.2b[[i]][2,11]




tt=plot(fits.h1[[3]])
p = approxfun(tt$x$gridx,tt$x$p.xfit)
x = seq(0.01,2000,length=500)
plot(x,p(x),type="l")
x = seq(0.001,50,length=500)
dx = diff(x)[1]
diff(p(x))[1]/dx

h1.to.HB(fits.h1[[16]]$par[1:2])
h1.to.HB(fits.h1[[17]]$par[1:2])


df17 = sim.df
fit.h1.2 <- LT2D.fit(DataFrameInput = df17,
                   hr = 'h1.2',
                   b = c(3, 0.5),
                   ystart = 1700,
                   pi.x = pi.fun.name,
                   logphi = logphi,
                   w = 2000, 
                   hessian = TRUE)
plot(fit.h1.2)

df17trunc = df17
# left-truncate so smallest x is at zero
df17trunc$x = df17trunc$x - min(df17trunc$x)+1
fit.h1.2tr <- LT2D.fit(DataFrameInput = df17trunc,
                     hr = 'h1.2',
                     b = c(6, -90),
                     ystart = 1700,
                     pi.x = pi.fun.name,
                     logphi = logphi,
                     w = 2000- min(df17trunc$x)+1, 
                     hessian = TRUE)
plot(fit.h1.2tr)


# fit fixed-b model
fit.h1.2b <- LT2D.fit(DataFrameInput = df17,
                       hr = 'h1.2b',
                       b = c(6),
                       ystart = 1700,
                       pi.x = pi.fun.name,
                       logphi = logphi,
                       w = 2000, 
                       hessian = TRUE)
plot(fit.h1.2b)

