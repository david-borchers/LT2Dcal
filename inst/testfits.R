td = read.csv("./inst/Gobi_gurvansaikhan_Detections.csv")
dat1 = td[td$Obs1==1,] # Use only detections by observer 1
names(dat1)
nd = dat1$Neardistance
r = dat1$Obs1.Distance
ang = (dat1$Obs1.AngleDetection-dat1$Obs1.AnglePath)
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
pdf("./inst/xyplot.pdf")
plot(x,y,pch="+",xlab="Perpendicular distance",ylab="Forward distance")
points(x[y<0],y[y<0],pch="+",col="red")
abline(h=0,lty=2)
dev.off()
sum(y<0)/length(y)

# Try fit a model
w = 1500
keepdat = (y>0 & x<=w)
yobs = y[keepdat]
xobs = x[keepdat]

L= 100
A = 1000
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

ystart = 2000
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
hist(xobs)

# plot perp dist distribution 
logphi.est = fit$fit$par[3:length(fit$fit$par)]
D.pdf = pi.fun(xs,logphi.est,w) # density pdf values to plot
plot(xs,D.pdf,type='l',ylim=range(0,D.pdf),xlab='Perp distance, x',ylab=expression(pi(x)))
hist(yobs)

fit$fit$AIC
par(mfrow=c(1,2))
gof.LT2D(fit, plot=TRUE)
par(mfrow=c(1,2))
plot(fit)

