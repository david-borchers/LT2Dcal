Gobi <- read.csv("./inst/report/Gobi.csv")
ystart = 1300; w = 1600
DistData <- Gobi[Gobi$Obs1 == 1,]
DistData$AngleDiff <- DistData$Obs1.AngleDetection-DistData$Obs1.AnglePath 
#convert to radians
DistData$AngleDiff <- DistData$AngleDiff * pi/180
#calculate distances
DistData$distance <- abs(DistData$Obs1.Distance*sin(DistData$AngleDiff)) 
DistData$forward <- DistData$Obs1.Distance*cos(DistData$AngleDiff)
#remove na's and negative y
DistData <- DistData[!is.na(DistData$forward)&DistData$forward>=0,]

#jitter zeros
xobs <- DistData$distance
yobs <- DistData$forward
#jitterdat = jitterzeros(xobs,yobs,xcut=10,ycut=10,anground.degree=5)
set.seed(1)
jitterdat = jitterzeros(xobs,yobs,xcut=30,ycut=30,anground.degree=10)
xjitter = jitterdat$x
yjitter = jitterdat$y
hist(xobs,breaks=c(seq(0,200,25),seq(300,2000,100)))
hist(xjitter,breaks=c(seq(0,200,25),seq(300,2000,100)))
hist(yobs,breaks=c(seq(0,200,25),seq(300,2000,100)))
hist(yjitter,breaks=c(seq(0,200,25),seq(300,2000,100)))

#add transects where we don't see anything
DistData$id <- 1:dim(DistData)[1]
missed <- Gobi[!Gobi$Site %in% DistData$Site,]
missed <- missed[!duplicated(missed$Site),]
missed$AngleDiff <- rep("", dim(missed)[1])
missed$Group.size <- rep(0, dim(missed)[1])
missed$distance <- rep(NA, dim(missed)[1])
missed$forward <- rep(NA, dim(missed)[1])
missed$id <- rep(NA, dim(missed)[1])
DistData <- rbind(DistData, missed)

#set length and area
L <- rep(0, dim(DistData)[1]) 
for(t in unique(DistData$Site)){ #probably a quicker way to do this
  L[DistData$Site == t][1] <- DistData$Transect_Length[DistData$Site == t][1]
}
A <- rep(0, dim(DistData)[1])
for(b in unique(Gobi$Block)){
  A[DistData$Block == b][1] <- 2*w*sum(unique(Gobi$Transect_Length[Gobi$Block == b]))
}
jitterdf <- data.frame(x = append(xjitter, missed$distance), y = append(yjitter,missed$forward),
                       stratum = DistData$Block,
                       transect = DistData$Site,
                       L = L,
                       area = A,
                       object = DistData$id,
                       size = DistData$Group.size)

#fit models
#combined.ip0.hnorm <- LT2D.fit(DataFrameInput = jitterdf,
#                               hr = 'ip0',
#                               b = c(4.9,0.035),
#                               ystart = 1300,      
#                               pi.x = 'pi.hnorm',   
#                               logphi = c(6.61348904, 4.84269089),
#                               w = 1600,
#                               hessian = TRUE,
#                               control=list(trace=5))

#combined.ep2.hnorm <- LT2D.fit(DataFrameInput = jitterdf,
#                               hr = 'ep2',
#                               b = c(qlogis(0.8), log(2), log(0.2), log(2)),
#                               ystart = 1300,      
#                               pi.x = 'pi.hnorm',   
#                               logphi = c(6.61348904, 4.84269089),
#                               w = 1500,
#                               hessian = TRUE,
#                               control=list(trace=5))
#combined.ip0.hnorm$fit$AIC; combined.ep2.hnorm$fit$AIC
#combined.ip0.hnorm$ests
#combined.ep2.hnorm$ests
#gof.LT2D(combined.ep2.hnorm, plot = T)
#gof.LT2D(combined.ip0.hnorm, plot = T)

#fit models
combined.ip0.hnorm2 <- LT2D.fit(DataFrameInput = jitterdf,
                               hr = 'ip0',
                               b = c(4.9,0.035),
                               ystart = 2000,      
                               pi.x = 'pi.hnorm',   
                               logphi = c(6.61348904, 4.84269089),
                               w = 1600,
                               hessian = TRUE,
                               control=list(trace=5))

combined.ep2.hnorm2 <- LT2D.fit(DataFrameInput = jitterdf,
                               hr = 'ep2',
                               b = c(qlogis(0.8), log(2), log(0.2), log(2)),
                               ystart = 2000,      
                               pi.x = 'pi.hnorm',   
                               logphi = c(6.61348904, 4.84269089),
                               w = 1500,
                               hessian = TRUE,
                               control=list(trace=5))
combined.ip0.hnorm2$fit$AIC; combined.ep2.hnorm2$fit$AIC
combined.ip0.hnorm2$ests
combined.ep2.hnorm2$ests
par(mfrow=c(1,2))
gof.LT2D(combined.ep2.hnorm2, plot = T)
gof.LT2D(combined.ip0.hnorm2, plot = T)
plot(combined.ep2.hnorm2)
plot(combined.ip0.hnorm2)
