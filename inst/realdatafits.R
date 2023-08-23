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

pi.fun.name <- 'pi.hnorm'; logphi <- log(c(700))

# Fit with hazard rate, constrained to have a shoulder
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
                          hr = hr.name,
                          b = hr.pars,
                          ystart = ystart,
                          pi.x = pi.fun.name,
                          logphi = logphi,
                          w = w, 
                          hessian = TRUE,
                          control=list(trace=5))

par(mfrow=c(2,2))
plot(comb.h1.2_chnorm)
gof.LT2D(comb.h1.2_chnorm,plot=TRUE)
comb.h1.2_chnorm$ests
comb.h1.2_chnorm.deriv.x = dpdx(h1.2,comb.h1.2_chnorm$fit$par[1:2],ystart,w,nx=nx)
comb.h1.2_chnorm.deriv.x[1]<0
plot(seq(0,w,length=(nx-1)),comb.h1.2_chnorm.deriv.x,type="l")
abline(0,0,lty=2)

par(mfrow=c(2,2))
plot(comb.ip0_chnorm)
gof.LT2D(comb.ip0_chnorm,plot=TRUE)
comb.ip0_chnorm$ests
comb.ip0_chnorm.deriv.x = dpdx(ip0,comb.ip0_chnorm$fit$par[1:2],ystart,w,nx=nx)
comb.ip0_chnorm.deriv.x[1]<0
plot(seq(0,w,length=(nx-1)),comb.ip0_chnorm.deriv.x,type="l")
abline(0,0,lty=2)

evaldf = data.frame(dfdx.10 =c(comb.h1.2_chnorm.deriv.x[10], comb.ip0_chnorm.deriv.x[10]),
                    Hessian = c(!is.null(comb.h1.2_chnorm$fit$hessian), !is.null(comb.ip0_chnorm$fit$hessian)),
                    AIC = c(comb.h1.2_chnorm$fit$AIC, comb.ip0_chnorm$fit$AIC))
row.names(evaldf) = c("h1.2","ip0")
evaldf

