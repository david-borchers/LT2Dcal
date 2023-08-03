td = read.csv("./inst/Gobi_gurvansaikhan_Detections.csv")
dat1 = td[td$Obs1==1,] # Use only detections by observer 1
names(dat1)
nd = dat1$Neardistance
r = dat1$Obs1.Distance
ang = (dat1$Obs1.AngleDetection-dat1$Obs1.AnglePath)
x = r*sin(ang*pi/180)
y = r*cos(ang*pi/180)
# get rid of data with NA angles
keep = !is.na(ang) 
ang = ang[keep]
x = x[keep]
y = y[keep]
r = r[keep]
# put all angles in (0, 180)
ang = abs(ang)
ang[ang>180] = 360 - ang[ang>180]
# convert to radians (because R works in radians)
rad = ang * pi/180
y = r*cos(rad)
x = r*sin(rad)

# Plot detection locations

# All origial data
pdf("./inst/xyplotAll.pdf")
plot(x,y,pch="+",xlab="Perpendicular distance",ylab="Forward distance")
points(x[y<0],y[y<0],pch="+",col="red")
abline(h=0,lty=2)
dev.off()
sum(y<0)/length(y) # proportion behind observer

# Extract only data with positive forward distance
w = 1500
ystart = 2000
keepdat = (y>0 & x<=w)
yobs = y[keepdat]
xobs = x[keepdat]

# All extracted data
pdf("./inst/xyplotForward.pdf")
plot(xobs,yobs,pch="+",xlab="Perpendicular distance",ylab="Forward distance")
dev.off()

pdf("./inst/histograms.pdf")
par(mfrow=c(2,2))
# perpendicular distance
hist(xobs,breaks=seq(0,1500,50),xlab="Perp distance intervals",main="Perp distance")
breaks=c(seq(0,100,length=11),seq(150,500,50),seq(600,1500,100))
hist(xobs,breaks=breaks,xlab="Perp distance intervals",main="Perp distance")
# forward distance
hist(yobs,breaks=seq(0,ystart,50),xlab="Forward distance intervals",main="Forward distance")
breaks=c(seq(0,100,length=11),seq(150,500,50),seq(600,2000,100))
hist(yobs,breaks=breaks,xlab="Forward distance intervals",main="Forward distance")
dev.off()


# Try fit a model
# fabricate some survey specification objects
L= sum(unique(dat1$Transect_Length))
A = 2*w*L
all.1s <- rep(1,length(xobs))
obj <- 1:length(xobs)
sim.df <- data.frame(x = xobs,
                     y = yobs,
                     stratum = all.1s,
                     transect = all.1s,
                     L = L,
                     area = A,
                     object = obj,
                     size = all.1s)

h.fun.name = "h1"
b = log(c(15,1.5))
h.fun = match.fun(h.fun.name) # make h.fun the function specified via a character variable
par(mfrow=c(1,3))
# plot detection function
nx = ny = 100
xs = seq(0,w,length=nx)
ys = seq(0,ystart,length=ny)
p.vals = p.approx(ys,xs,h.fun,b) # detection function values to plot
plot(xs,p.vals,type='l',ylim=range(0,p.vals),xlab='Perp. distance, x',ylab=expression(p(x)))
# and now the 2D detection function:
pmat = p.approx(ys,xs,h.fun,b,xy=TRUE)
persp(x=xs,y=ys,z=pmat,theta=120,phi=25, xlab="Perp dist (x)", ylab="Forward dist (y)",
      zlab="p(det by y)")
# and now the 2D pdf:
pdfmat = p.approx(ys,xs,h.fun,b,pdf=TRUE)
persp(x=xs,y=ys,z=pdfmat,theta=45,phi=25, xlab="Perp dist (x)", ylab="Forward dist (y)",
      zlab="f(y|x)")

# plot perp dist distribution 
pi.fun.name = "pi.norm"; phi = c(500,600); logphi = c(phi[1],log(phi[2])) # for pi.norm
pi.fun.name = "pi.hnorm"; phi = c(2000); logphi = log(phi[1]) # for pi.hnorm
pi.fun = match.fun(pi.fun.name) # make Dfun the function specified via a character variable
D.pdf = pi.fun(xs,logphi,w) # density pdf values to plot
plot(xs,D.pdf,type='l',ylim=range(0,D.pdf),xlab='Perp. distance, x',ylab=expression(pi(x)))


# fit an LT2D model
fit <- LT2D.fit(DataFrameInput = sim.df,
                hr = h.fun.name,
                # start values for b:
                b = b,
                ystart = ystart,
                pi.x = pi.fun.name,
                # start values for logphi:
                logphi = logphi,
                w = w,
                hessian = TRUE)
b.est = fit$fit$par[1:2]
p.vals = p.approx(ys,xs,h.fun,b.est) # detection function values to plot
plot(xs,p.vals,type='l',ylim=range(0,p.vals),xlab='Perp. distance, x',ylab=expression(p(x)))
# and now the 2D detection function:
pmat = p.approx(ys,xs,h.fun,b.est,xy=TRUE)
persp(x=xs,y=ys,z=pmat,theta=120,phi=25, xlab="Perp dist (x)", ylab="Forward dist (y)",
      zlab="p(det by y)")
# and now the 2D pdf:
pdfmat = p.approx(ys,xs,h.fun,b.est,pdf=TRUE)
persp(x=xs,y=ys,z=pdfmat,theta=45,phi=25, xlab="Perp dist (x)", ylab="Forward dist (y)",
      zlab="f(y|x)")

# plot perp dist distribution 
logphi.est = fit$fit$par[3:length(fit$fit$par)]
D.pdf = pi.fun(xs,logphi.est,w) # density pdf values to plot
plot(xs,D.pdf,type='l',ylim=range(0,D.pdf),xlab='Perp distance, x',ylab=expression(pi(x)))


fit$fit$AIC
par(mfrow=c(1,2))
gof.LT2D(fit, plot=TRUE)
par(mfrow=c(1,2))
plot(fit)

# Try fit to ad-hoc adjusted data
# ===============================
# Ad-hoc undoing of rounding to zero:
jitterdat = jitterzeros(xobs,yobs,xcut=10,ycut=10,anground.degree=5)
xjitter = jitterdat$x
yjitter = jitterdat$y
plot(xobs,xjitter)
plot(yobs,yjitter)
#pdf("./inst/xyplotJittered.pdf")
plot(xjitter,yjitter,pch="+",xlab="Perpendicular distance",ylab="Forward distance")
#dev.off()
#pdf("./inst/histograms.pdf")
par(mfrow=c(2,2))
# perpendicular distance
hist(xjitter,breaks=seq(0,1500,50),xlab="Perp distance intervals",main="Perp distance")
breaks=c(seq(0,100,length=11),seq(150,500,50),seq(600,1500,100))
hist(xjitter,breaks=breaks,xlab="Perp distance intervals",main="Perp distance")
# forward distance
hist(yjitter,breaks=seq(0,ystart,50),xlab="Forward distance intervals",main="Forward distance")
breaks=c(seq(0,100,length=11),seq(150,500,50),seq(600,2000,100))
hist(yjitter,breaks=breaks,xlab="Forward distance intervals",main="Forward distance")
#dev.off()

adj.df <- data.frame(x = xjitter,
                     y = yjitter,
                     stratum = all.1s,
                     transect = all.1s,
                     L = L,
                     area = A,
                     object = obj,
                     size = all.1s)

# half-normal
fitadj <- LT2D.fit(DataFrameInput = adj.df,
                hr = h.fun.name,
                # start values for b:
                b = b,
                ystart = ystart,
                pi.x = pi.fun.name,
                # start values for logphi:
                logphi = logphi,
                w = w,
                hessian = TRUE)
fitadj$ests
fit$ests

sigma.cv = exp(fitadj$fit$par[3] + c(-1,1)*1.96*sqrt(fitadj$fit$vcov[3,3]))
pi.est = match.fun(pi.fun.name)(xs,fitadj$fit$par[3],w)/match.fun(pi.fun.name)(0,fitadj$fit$par[3],w)
pi.lcl = match.fun(pi.fun.name)(xs,log(sigma.cv[1]),w)/match.fun(pi.fun.name)(0,log(sigma.cv[1]),w)
pi.ucl = match.fun(pi.fun.name)(xs,log(sigma.cv[2]),w)/match.fun(pi.fun.name)(0,log(sigma.cv[2]),w)
plot(xs,pi.est,type="l",ylim=c(0,1))
lines(xs,pi.lcl,lty=2)
lines(xs,pi.ucl,lty=2)

b.est = fitadj$fit$par[1:2]
p.vals = p.approx(ys,xs,h.fun,b.est) # detection function values to plot
plot(xs,p.vals,type='l',ylim=range(0,p.vals),xlab='Perp. distance, x',ylab=expression(p(x)))
# and now the 2D detection function:
pmat = p.approx(ys,xs,h.fun,b.est,xy=TRUE)
persp(x=xs,y=ys,z=pmat,theta=120,phi=25, xlab="Perp dist (x)", ylab="Forward dist (y)",
      zlab="p(det by y)")
# and now the 2D pdf:
pdfmat = p.approx(ys,xs,h.fun,b.est,pdf=TRUE)
persp(x=xs,y=ys,z=pdfmat,theta=45,phi=25, xlab="Perp dist (x)", ylab="Forward dist (y)",
      zlab="f(y|x)")

# plot perp dist distribution 
logphi.est = fitadj$fit$par[3:length(fitadj$fit$par)]
D.pdf = pi.fun(xs,logphi.est,w) # density pdf values to plot
plot(xs,D.pdf,type='l',ylim=range(0,D.pdf),xlab='Perp distance, x',ylab=expression(pi(x)))


fitadj$fit$AIC
par(mfrow=c(1,2))
gof.LT2D(fitadj, plot=TRUE)
par(mfrow=c(1,2))
plot(fitadj)


# Complimentary half-normal
pi.fun.name = "pi.chnorm"
logphi = c(0,log(2000))
fitadj.chn <- LT2D.fit(DataFrameInput = adj.df,
                   hr = h.fun.name,
                   # start values for b:
                   b = b,
                   ystart = ystart,
                   pi.x = pi.fun.name,
                   # start values for logphi:
                   logphi = logphi,
                   w = w,
                   hessian = TRUE)
fitadj.chn$ests
fitadj$ests
fit$ests

b.est = fitadj.chn$fit$par[1:2]
p.vals = p.approx(ys,xs,h.fun,b.est) # detection function values to plot
plot(xs,p.vals,type='l',ylim=range(0,p.vals),xlab='Perp. distance, x',ylab=expression(p(x)))
# and now the 2D detection function:
pmat = p.approx(ys,xs,h.fun,b.est,xy=TRUE)
persp(x=xs,y=ys,z=pmat,theta=120,phi=25, xlab="Perp dist (x)", ylab="Forward dist (y)",
      zlab="p(det by y)")
# and now the 2D pdf:
pdfmat = p.approx(ys,xs,h.fun,b.est,pdf=TRUE)
persp(x=xs,y=ys,z=pdfmat,theta=45,phi=25, xlab="Perp dist (x)", ylab="Forward dist (y)",
      zlab="f(y|x)")

# plot perp dist distribution 
logphi.est = fitadj.chn$fit$par[3:length(fitadj.chn$fit$par)]
D.pdf = match.fun(pi.fun.name)(xs,logphi.est,w) # density pdf values to plot
plot(xs,D.pdf,type='l',ylim=range(0,D.pdf),xlab='Perp distance, x',ylab=expression(pi(x)))


fitadj.chn$fit$AIC
par(mfrow=c(1,2))
gof.LT2D(fitadj.chn, plot=TRUE)
par(mfrow=c(1,2))
plot(fitadj.chn)
