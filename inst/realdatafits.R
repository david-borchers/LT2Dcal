dfs = readRDS("./inst/gobidf.Rds")
names(dfs)
argdf = dfs$argalijitter
ibdf = dfs$ibexjitter
combdf = dfs$jitterdf


# Exploratory plots
xbreaks = seq(0,max(na.omit(combdf$x)),length=11)
ybreaks = seq(0,max(na.omit(combdf$y)),length=11)
pdf("./inst/argibplots.pdf")
par(mfrow=c(3,3))
# Argali
plot(argdf$x,argdf$y,pch="+",col="blue")
hist(argdf$x,col="blue",breaks=xbreaks)
hist(argdf$y,col="blue",breaks=ybreaks)
# Ibex
plot(ibdf$x,ibdf$y,pch="+",col="red")
hist(ibdf$x,col="red",breaks=xbreaks)
hist(ibdf$y,col="red",breaks=ybreaks)
# combied
plot(ibdf$x,ibdf$y,pch="+",col="blue",xlim=range(na.omit(combdf$x)),ylim=range(na.omit(combdf$y)))
points(argdf$x,argdf$y,pch="+",col="red")
hist(combdf$x,col="purple",breaks=xbreaks)
hist(combdf$y,col="purple",breaks=ybreaks)
dev.off()



# ======
# Argali
# ======
pi.fun.name <- 'pi.hnorm'; logphi <- log(c(700))

# Fit with hazard rate, constrained to have a shoulder
hr.name = "h1.2"; hr.pars=c(12, 0.25)
ystart = 2000
w=1000
par(mfrow=c(1,2))
plotdetmodel(hr.name,hr.pars,ystart,xlim=c(0,w),what="px")
plotpimodel(pi.fun.name,logphi,w)
#dev.off()

arg.h1.2_chnorm <- LT2D.fit(DataFrameInput = argdf,
                     hr = hr.name,
                     b = hr.pars,
                     ystart = ystart,
                     pi.x = pi.fun.name,
                     logphi = logphi,
                     w = w, 
                     hessian = TRUE,
                     control=list(trace=5))

# Fit with model ip0
ystart = 2000
w=1000
hr.name = "ip0"; hr.pars=c(7, 0.25)
nx=101
deriv.x = dpdx(hr.name,hr.pars,ystart,w,nx=nx)
par(mfrow=c(1,3))
plotdetmodel(hr.name,hr.pars,ystart,xlim=c(0,w),what="px")
plot(seq(0,w,length=(nx-1)),deriv.x,type="l",xlab="x",ylab="Derivative")
abline(0,0,lty=2)
plotpimodel(pi.fun.name,logphi,w)
#dev.off()

arg.ip0_chnorm <- LT2D.fit(DataFrameInput = argdf,
                            hr = hr.name,
                            b = hr.pars,
                            ystart = ystart,
                            pi.x = pi.fun.name,
                            logphi = logphi,
                            w = w, 
                            hessian = TRUE,
                            control=list(trace=5))

plot(arg.h1.2_chnorm)
gof.LT2D(arg.h1.2_chnorm,plot=TRUE)
arg.h1.2_chnorm$ests
arg.h1.2_chnorm.deriv.x = dpdx(hr.name,arg.h1.2_chnorm$fit$par[1:2],ystart,w,nx=nx)
arg.h1.2_chnorm.deriv.x[1]<0
plot(seq(0,w,length=(nx-1)),arg.h1.2_chnorm.deriv.x,type="l")
abline(0,0,lty=2)

plot(arg.ip0_chnorm)
gof.LT2D(arg.ip0_chnorm,plot=TRUE)
arg.ip0_chnorm$ests
arg.ip0_chnorm.deriv.x = dpdx(hr.name,arg.ip0_chnorm$fit$par[1:2],ystart,w,nx=nx)
arg.ip0_chnorm.deriv.x[1]<0
plot(seq(0,w,length=(nx-1)),arg.ip0_chnorm.deriv.x,type="l")
abline(0,0,lty=2)

arg.h1.2_chnorm$fit$hessian; arg.ip0_chnorm$fit$hessian
arg.h1.2_chnorm$fit$AIC; arg.ip0_chnorm$fit$AIC



# ======
# Ibex
# ======

pi.fun.name <- 'pi.hnorm'; logphi <- log(c(700))

# Fit with hazard rate, constrained to have a shoulder
hr.name = "h1.2"; hr.pars=c(12, 0.25)
ystart = 2000
w=1000
par(mfrow=c(1,2))
plotdetmodel(hr.name,hr.pars,ystart,xlim=c(0,w),what="px")
plotpimodel(pi.fun.name,logphi,w)
#dev.off()

ib.h1.2_chnorm <- LT2D.fit(DataFrameInput = ibdf,
                            hr = hr.name,
                            b = hr.pars,
                            ystart = ystart,
                            pi.x = pi.fun.name,
                            logphi = logphi,
                            w = w, 
                            hessian = TRUE,
                            control=list(trace=5))

# Fit with model ip0
ystart = 2000
w=1000
hr.name = "ip0"; hr.pars=c(7, 0.25)
nx=101
deriv.x = dpdx(hr.name,hr.pars,ystart,w,nx=nx)
par(mfrow=c(1,3))
plotdetmodel(hr.name,hr.pars,ystart,xlim=c(0,w),what="px")
plot(seq(0,w,length=(nx-1)),deriv.x,type="l",xlab="x",ylab="Derivative")
abline(0,0,lty=2)
plotpimodel(pi.fun.name,logphi,w)
#dev.off()

ib.ip0_chnorm <- LT2D.fit(DataFrameInput = ibdf,
                           hr = hr.name,
                           b = hr.pars,
                           ystart = ystart,
                           pi.x = pi.fun.name,
                           logphi = logphi,
                           w = w, 
                           hessian = TRUE,
                           control=list(trace=5))

par(mfrow=c(2,2))
plot(ib.h1.2_chnorm)
gof.LT2D(ib.h1.2_chnorm,plot=TRUE)
ib.h1.2_chnorm$ests
ib.h1.2_chnorm.deriv.x = dpdx(h1.2,ib.h1.2_chnorm$fit$par[1:2],ystart,w,nx=nx)
ib.h1.2_chnorm.deriv.x[1]<0
plot(seq(0,w,length=(nx-1)),ib.h1.2_chnorm.deriv.x,type="l")
abline(0,0,lty=2)

par(mfrow=c(2,2))
plot(ib.ip0_chnorm)
gof.LT2D(ib.ip0_chnorm,plot=TRUE)
ib.ip0_chnorm$ests
ib.ip0_chnorm.deriv.x = dpdx(ip0,ib.ip0_chnorm$fit$par[1:2],ystart,w,nx=nx)
ib.ip0_chnorm.deriv.x[1]<0
plot(seq(0,w,length=(nx-1)),ib.ip0_chnorm.deriv.x,type="l")
abline(0,0,lty=2)

ib.h1.2_chnorm$fit$hessian; ib.ip0_chnorm$fit$hessian
ib.h1.2_chnorm$fit$AIC; ib.ip0_chnorm$fit$AIC



# ========
# Combined
# ========

# Fit with hazard rate, constrained to have a shoulder
pi.fun.name <- 'pi.hnorm'; logphi <- log(c(700))
hr.name = "h1.2"; hr.pars=c(12, 0.25)
ystart = 2000
w=1000
par(mfrow=c(1,2))
plotdetmodel(hr.name,hr.pars,ystart,xlim=c(0,w),what="px")
plotpimodel(pi.fun.name,logphi,w)
#dev.off()

comb.h1.2_chnorm <- LT2D.fit(DataFrameInput = combdf,
                           hr = "h1.2",
                           b = hr.pars,
                           ystart = ystart,
                           pi.x = pi.fun.name,
                           logphi = logphi,
                           w = w, 
                           hessian = TRUE,
                           control=list(trace=5))

# Fit with hazard rate, constrained to have a shoulder and UNIFORM x
pi.fun.name <- 'pi.const'; logphi <- NULL
hr.name = "h1.2"; hr.pars=c(12, 0.25)
ystart = 2000
w=1000
par(mfrow=c(1,2))
plotdetmodel(hr.name,hr.pars,ystart,xlim=c(0,w),what="px")
plotpimodel("pi.const",logphi,w)
#dev.off()

comb.h1.2_chnorm.unif <- LT2D.fit(DataFrameInput = combdf,
                             hr = hr.name,
                             b = hr.pars,
                             ystart = ystart,
                             pi.x =  "pi.const",
                             logphi = NULL,
                             w = w, 
                             hessian = TRUE,
                             control=list(trace=5))

# Fit with model ip0
ystart = 2000
w=1000
hr.name = "ip0"; hr.pars=c(7, 0.25)
nx=101
deriv.x = dpdx(hr.name,hr.pars,ystart,w,nx=nx)
par(mfrow=c(1,3))
plotdetmodel(hr.name,hr.pars,ystart,xlim=c(0,w),what="px")
plot(seq(0,w,length=(nx-1)),deriv.x,type="l",xlab="x",ylab="Derivative")
abline(0,0,lty=2)
plotpimodel(pi.fun.name,logphi,w)
#dev.off()

comb.ip0_chnorm <- LT2D.fit(DataFrameInput = combdf,
                          hr = "ip0",
                          b = hr.pars,
                          ystart = ystart,
                          pi.x = pi.fun.name,
                          logphi = logphi,
                          w = w, 
                          hessian = TRUE,
                          control=list(trace=5))


par(mfrow=c(2,2))
plot(comb.h1.2_chnorm)
comb.h1.2_chnorm.gof = gof.LT2D(comb.h1.2_chnorm,plot=TRUE)
comb.h1.2_chnorm$ests
comb.h1.2_chnorm.deriv.x = dpdx(h1.2,comb.h1.2_chnorm$fit$par[1:2],ystart,w,nx=nx)
comb.h1.2_chnorm.deriv.x[1]<0
plot(seq(0,w,length=(nx-1)),comb.h1.2_chnorm.deriv.x,type="l")
abline(0,0,lty=2)


par(mfrow=c(2,2))
plot(comb.h1.2_chnorm.unif)
comb.h1.2_chnorm.unif.gof = gof.LT2D(comb.h1.2_chnorm.unif,plot=TRUE)
comb.h1.2_chnorm.unif$ests
comb.h1.2_chnorm.unif.deriv.x = dpdx(h1.2,comb.h1.2_chnorm.unif$fit$par[1:2],ystart,w,nx=nx)
comb.h1.2_chnorm.unif.deriv.x[1]<0
plot(seq(0,w,length=(nx-1)),comb.h1.2_chnorm.unif.deriv.x,type="l")
abline(0,0,lty=2)


par(mfrow=c(2,2))
plot(comb.ip0_chnorm)
comb.ip0_chnorm.gof = gof.LT2D(comb.ip0_chnorm,plot=TRUE)
comb.ip0_chnorm$ests
comb.ip0_chnorm.deriv.x = dpdx(ip0,comb.ip0_chnorm$fit$par[1:2],ystart,w,nx=nx)
comb.ip0_chnorm.deriv.x[1]<0
plot(seq(0,w,length=(nx-1)),comb.ip0_chnorm.deriv.x,type="l")
abline(0,0,lty=2)

evaldf.2D = data.frame(dfdx.10 =c(comb.h1.2_chnorm.deriv.x[10], 
                                  comb.h1.2_chnorm.unif.deriv.x[10],
                                  comb.ip0_chnorm.deriv.x[10]),
                    Hessian = c(!is.null(comb.h1.2_chnorm$fit$hessian), 
                                !is.null(comb.h1.2_chnorm.unif$fit$hessian), 
                                !is.null(comb.ip0_chnorm$fit$hessian)),
                    AIC = c(comb.h1.2_chnorm$fit$AIC, 
                            comb.h1.2_chnorm.unif$fit$AIC, 
                            comb.ip0_chnorm$fit$AIC),
                    N = c(comb.h1.2_chnorm$ests[nrow(comb.h1.2_chnorm$ests),8],
                          comb.h1.2_chnorm.unif$ests[nrow(comb.h1.2_chnorm.unif$ests),8],
                          comb.ip0_chnorm$ests[nrow(comb.ip0_chnorm$ests),8]),
                    CvM.x = c(comb.h1.2_chnorm.gof$X[2],
                              comb.h1.2_chnorm.unif.gof$X[2],
                              comb.ip0_chnorm.gof$X[2]),
                    CvM.y = c(comb.h1.2_chnorm.gof$Y[2],
                              comb.h1.2_chnorm.unif.gof$Y[2],
                              comb.ip0_chnorm.gof$Y[2])
                    )
row.names(evaldf.2D) = c("h1.2","h1.2.unif","ip0")
evaldf.2D

# plot 95\% envelope for pi.hnorm:
xs = seq(0,w,length=100)
pivals = pi.hnorm(xs,comb.ip0_chnorm$fit$par[3],w=w)
logci = comb.ip0_chnorm$fit$par[3] + c(-1,1)*1.96*sqrt(comb.ip0_chnorm$fit$vcov[3,3])
lcl = pi.hnorm(xs,logci[1],w=w)
ucl = pi.hnorm(xs,logci[2],w=w)
plot(xs,pivals/max(pivals),type="l",ylim=c(0,1))
lines(xs,lcl/max(lcl),lty=2)
lines(xs,ucl/max(ucl),lty=2)

# Before moving on, just try best model (ip0) with w = 800:
ystart = 2000
w=800
hr.name = "ip0"; hr.pars=comb.ip0_chnorm$fit$par[1:2]
logphi = comb.ip0_chnorm$fit$par[3]
nx=101
deriv.x = dpdx(hr.name,hr.pars,ystart,w,nx=nx)
par(mfrow=c(1,3))
plotdetmodel(hr.name,hr.pars,ystart,xlim=c(0,w),what="px")
plot(seq(0,w,length=(nx-1)),deriv.x,type="l",xlab="x",ylab="Derivative")
abline(0,0,lty=2)
plotpimodel(pi.fun.name,logphi,w)
#dev.off()

comb.ip0_chnorm.800 <- LT2D.fit(DataFrameInput = combdf,
                            hr = "ip0",
                            b = hr.pars,
                            ystart = ystart,
                            pi.x = pi.fun.name,
                            logphi = logphi,
                            w = w, 
                            hessian = TRUE,
                            control=list(trace=5))

par(mfrow=c(2,2))
plot(comb.ip0_chnorm.800)
comb.ip0_chnorm.800.gof = gof.LT2D(comb.ip0_chnorm.800,plot=TRUE)
comb.ip0_chnorm.800$ests
comb.ip0_chnorm.800.deriv.x = dpdx(h1.2,comb.ip0_chnorm.800$fit$par[1:2],ystart,w,nx=nx)
comb.ip0_chnorm.800.deriv.x[1]<0
plot(seq(0,w,length=(nx-1)),comb.ip0_chnorm.800.deriv.x,type="l")
abline(0,0,lty=2)

saveRDS(list(comb.h1.2_chnorm=comb.h1.2_chnorm,
             comb.h1.2_chnorm.unif=comb.h1.2_chnorm.unif,
             comb.ip0_chnorm=comb.ip0_chnorm,
             comb.ip0_chnorm.800=comb.ip0_chnorm.800),
        file="./inst/comb_chnorm_fits.Rda")



#=========================================
# Double-observer
#=========================================
chapman <- function(Gobi){
  S1 <- sum(Gobi$Obs1)
  S2 <- sum(Gobi$Obs2)
  B <- sum(Gobi$Obs1[Gobi$Obs2 == 1])
  Ngroups <- (S1+1)*(S2+ 1)/(B+1) - 1
  varNhat <- (S1+1)*(S2+ 1)*(S1-B)*(S2-B)/((B+2)*(B+1)**2)
  d <- exp(1.96*sqrt(log(1 + varNhat/Ngroups**2)))
  ci <- c(Ngroups/d, Ngroups*d)
  return(c(Ngroups, ci))
}
combined <- chapman(Gobi)
ibex.cm <- chapman(Gobi[Gobi$Sightings == 'ibex',])
argali.cm <- chapman(Gobi[Gobi$Sightings == 'argali',])







#=========================================
# CDS
#=========================================
require(mrds)
require(Distance)
hist(CDSData$distance, main = "Perpendicular Distance", xlab = "Distance")
w.cds = 1000
hr_block <- ds(CDSdf, truncation = w.cds, formula = ~Region.Label, key = "hr")
hr_block$dht$clusters$N
#hr_ibex_block <- ds(ibexdf, truncation = w.cds, formula = ~Region.Label, key = "hr")
#hr_argali_block <- ds(argalidf, truncation = w.cds, formula = ~Region.Label, key = "hr")
#total.summary <- summary(hr_block)
#ibex.summary <- summary(hr_ibex_block)
#argali.summary <- summary(hr_argali_block)






#=========================================
# MRDS
#=========================================