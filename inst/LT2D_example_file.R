library(devtools)

install_github("david-borchers/LT2Dcal",force=TRUE)

library('LT2D')

# set perp trunc, forward trunc
# total transect length and survey
# area:
w = 0.15 ; ystart = 0.55
L = 10 ; A = 2*w*L

# set value of 'true' parameters
# for simulated data:
b=c(-7.3287948, 0.9945317)
logphi <- c(0.02,-4.42)

# produce simulated data:
set.seed(3)
simDat = simXY(50, 'pi.norm',
               logphi, 'h1', 
               b, w, 
               ystart)$locs

# create the data.frame:
all.1s <- rep(1,length(simDat$x))
obj <- 1:length(simDat$x)
sim.df <- data.frame(x = simDat$x,
                     y = simDat$y,
                     stratum = all.1s,
                     transect = all.1s,
                     L = L,
                     area = A,
                     object = obj,
                     size = all.1s)

# fit an LT2D model
fit <- LT2D.fit(DataFrameInput = sim.df,
                hr = 'h1',
                # start values for b:
                b = b,
                ystart = ystart,
                pi.x = 'pi.norm',
                # start values for logphi:
                logphi = logphi,
                w = w,
                hessian = TRUE)

# plot fitted functions
par(mfrow=c(1,2))
plot(fit)

# print table of point estimates
fit$ests

# look at goodness of fit
par(mfrow=c(1,2))
gof.LT2D(fit, plot=TRUE)

# Look at AIC
fit$fit$AIC

# bootstrap for variance and interval estimates of total abundance
boot <- LT2D.bootstrap(fit,r=999,alpha = 0.05)
boot$ci
# plot bootstrap histogram:
hist(boot$Ns,main='',xlab='Estimate of Abundance')
# add point estimate and CI
mle.N = fit$ests[nrow(fit$ests),ncol(fit$ests)]
points(mle.N,0,pch=19,col="red")
arrows(boot$ci[1],0,boot$ci[2],0,angle=90,code=3,col="red")
