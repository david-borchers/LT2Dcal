
# Set Up Functions --------------------------------------------------------
library(circular)
library(mrds)
library(Distance)
library(knitr)
library(LT2D)
# Gobi Analysis -----------------------------------------------------------
Gobi <- read.csv("Gobi.csv")
ystart = 1700; w = 1600
DistData <- Gobi[Gobi$Obs1 == 1,]
DistData$AngleDiff <- DistData$Obs1.AngleDetection-DistData$Obs1.AnglePath 
#convert to radians
DistData$AngleDiff <- DistData$AngleDiff * pi/180
#calculate distances
DistData$distance <- abs(DistData$Obs1.Distance*sin(DistData$AngleDiff)) 
DistData$forward <- DistData$Obs1.Distance*cos(DistData$AngleDiff)
#create copy for CDS
CDSData <- DistData 
#remove NAs and negative y
DistData <- DistData[!is.na(DistData$forward)&DistData$forward>=0,]

# CDS combined ------------------------------------------------------------
xobs <- CDSData$distance
yobs <- CDSData$forward
jitterdat = jitterzeros(xobs,yobs,xcut=10,ycut=10,anground.degree=5)
xjitter = jitterdat$x
missed <- Gobi[!Gobi$Site %in% CDSData$Site,]
missed <- missed[!duplicated(missed$Site),]
missed$AngleDiff <- rep(NA, dim(missed)[1])
missed$Group.size <- rep(0, dim(missed)[1])
missed$distance <- rep(NA, dim(missed)[1])
missed$forward <- rep(NA, dim(missed)[1])
missed$Sightings <- rep(NA, dim(missed)[1])
CDSData <- rbind(CDSData, missed)
L <- rep(0, dim(CDSData)[1]) 
#there's probably a way to do this without for loops
#but it would take longer for me to figure out than just running the loops
for(t in unique(CDSData$Site)){
  L[CDSData$Site == t] <- CDSData$Transect_Length[CDSData$Site == t][1]
}
A <- rep(0, dim(CDSData)[1])
for(b in unique(Gobi$Block)){
  A[CDSData$Block == b] <- 2*w*sum(unique(Gobi$Transect_Length[Gobi$Block == b]))
}

CDSdf <- data.frame(Region.Label = CDSData$Block,
                    Sample.Label = CDSData$Site,
                    distance = append(xjitter, missed$distance), 
                    size = CDSData$Group.size,
                    Effort = L,
                    Area = A,
                    species = CDSData$Sightings,
                    obs.name = CDSData$Choidogjamts)
# CDS Ibex only ---------------------------------------------------------------
ibexCDS <- CDSData[CDSData$Species == 'ibex',]
xobs <- ibexCDS$distance
yobs <- ibexCDS$forward
jitterdat = jitterzeros(xobs,yobs,xcut=10,ycut=10,anground.degree=5)
xjitter = jitterdat$x
missed <- Gobi[!Gobi$Site %in% ibexCDS$Site,]
missed <- missed[!duplicated(missed$Site),]
missed$AngleDiff <- rep(NA, dim(missed)[1])
missed$Group.size <- rep(0, dim(missed)[1])
missed$distance <- rep(NA, dim(missed)[1])
missed$forward <- rep(NA, dim(missed)[1])
missed$Sightings <- rep(NA, dim(missed)[1])
ibexCDS <- rbind(ibexCDS, missed)
L <- rep(0, dim(ibexCDS)[1]) 
for(t in unique(ibexCDS$Site)){
  L[ibexCDS$Site == t] <- ibexCDS$Transect_Length[ibexCDS$Site == t][1]
}
A <- rep(0, dim(ibexCDS)[1])
for(b in unique(Gobi$Block)){
  A[ibexCDS$Block == b] <- 2*w*sum(unique(Gobi$Transect_Length[Gobi$Block == b]))
}

ibexdf <- data.frame(Region.Label = ibexCDS$Block,
                     Sample.Label = ibexCDS$Site,
                     distance = append(xjitter, missed$distance), 
                     size = ibexCDS$Group.size,
                     Effort = L,
                     Area = A,
                     obs.name = ibexCDS$Choidogjamts)
#CDS Argali only -------------------------------------------------------------
argaliCDS <- CDSData[CDSData$Species == 'argali',]
xobs <- argaliCDS$distance
yobs <- argaliCDS$forward
jitterdat = jitterzeros(xobs,yobs,xcut=10,ycut=10,anground.degree=5)
xjitter = jitterdat$x
missed <- Gobi[!Gobi$Site %in% argaliCDS$Site,]
missed <- missed[!duplicated(missed$Site),]
missed$AngleDiff <- rep(NA, dim(missed)[1])
missed$Group.size <- rep(0, dim(missed)[1])
missed$distance <- rep(NA, dim(missed)[1])
missed$forward <- rep(NA, dim(missed)[1])
missed$Sightings <- rep(NA, dim(missed)[1])
argaliCDS <- rbind(argaliCDS, missed)
L <- rep(0, dim(argaliCDS)[1]) 
for(t in unique(argaliCDS$Site)){
  L[argaliCDS$Site == t] <- argaliCDS$Transect_Length[argaliCDS$Site == t][1]
}
A <- rep(0, dim(argaliCDS)[1])
for(b in unique(Gobi$Block)){
  A[argaliCDS$Block == b] <- 2*w*sum(unique(Gobi$Transect_Length[Gobi$Block == b]))
}

argalidf <- data.frame(Region.Label = argaliCDS$Block,
                       Sample.Label = argaliCDS$Site,
                       distance = append(xjitter, missed$distance), 
                       size = argaliCDS$Group.size,
                       Effort = L,
                       Area = A,
                       obs.name = argaliCDS$Choidogjamts)

# 2D combined -------------------------------------------------------------
xobs <- DistData$distance
yobs <- DistData$forward
jitterdat = jitterzeros(xobs,yobs,xcut=10,ycut=10,anground.degree=5)
xjitter = jitterdat$x
yjitter = jitterdat$y
#add missing
DistData$id <- 1:dim(DistData)[1]
missed <- Gobi[!Gobi$Site %in% DistData$Site,]
missed <- missed[!duplicated(missed$Site),]
missed$AngleDiff <- rep("", dim(missed)[1])
missed$Group.size <- rep(0, dim(missed)[1])
missed$distance <- rep(NA, dim(missed)[1])
missed$forward <- rep(NA, dim(missed)[1])
missed$id <- rep(NA, dim(missed)[1])
DistData <- rbind(DistData, missed)
L <- rep(0, dim(DistData)[1]) 
#there's probably a way to do this without for loops
#but it would take longer for me to figure out than just running the loops
for(t in unique(DistData$Site)){
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

# Ibex --------------------------------------------------------------------
ibex <- DistData[DistData$Obs1 == 1 & DistData$Sightings == "ibex",]
xobs <- ibex$distance
yobs <- ibex$forward
jitterdat = jitterzeros(xobs,yobs,xcut=10,ycut=10,anground.degree=5)
xjitter = jitterdat$x
yjitter = jitterdat$y
#adding missed
ibex$id <- 1:dim(ibex)[1]
missed <- Gobi[!Gobi$Site %in% ibex$Site,]
missed <- missed[!duplicated(missed$Site),]
missed$AngleDiff <- rep("", dim(missed)[1])
missed$Group.size <- rep(0, dim(missed)[1])
missed$distance <- rep(NA, dim(missed)[1])
missed$forward <- rep(NA, dim(missed)[1])
missed$id <- rep(NA, dim(missed)[1])
ibex <- rbind(ibex, missed)

Li <- rep(0, dim(ibex)[1]) 
for(t in unique(Gobi$Site)){
  Li[ibex$Site == t][1] <- Gobi$Transect_Length[Gobi$Site == t][1]
}
Ai <- rep(0, length(xjitter))
for(b in unique(Gobi$Block)){
  Ai[ibex$Block == b][1] <- 2*w*sum(unique(Gobi$Transect_Length[Gobi$Block == b]))
}
ibexjitter <- data.frame(x = append(xjitter,missed$distance), y = append(yjitter,missed$forward),
                         stratum = ibex$Block,
                         transect = ibex$Site,
                         L = Li,
                         area = Ai,
                         object = ibex$id,
                         size = ibex$Group.size)

# Argali ------------------------------------------------------------------
argali <- DistData[DistData$Obs1 == 1 & DistData$Sightings == "argali",]
xobs <- argali$distance
yobs <- argali$forward
jitterdat = jitterzeros(xobs,yobs,xcut=10,ycut=10,anground.degree=5)
xjitter = jitterdat$x
yjitter = jitterdat$y

#adding missed
argali$id <- 1:dim(argali)[1]
missed <- Gobi[!Gobi$Site %in% argali$Site,]
missed <- missed[!duplicated(missed$Site),]
missed$AngleDiff <- rep("", dim(missed)[1])
missed$Group.size <- rep(0, dim(missed)[1])
missed$distance <- rep(NA, dim(missed)[1])
missed$forward <- rep(NA, dim(missed)[1])
missed$id <- rep(NA, dim(missed)[1])
argali <- rbind(argali, missed)

La <- rep(0, dim(argali)[1]) 
for(t in unique(Gobi$Site)){
  La[argali$Site == t][1] <- argali$Transect_Length[argali$Site == t][1]
}
Aa <- rep(0, length(xjitter))
for(b in unique(Gobi$Block)){
  Aa[argali$Block == b][1] <- 2*w*sum(unique(Gobi$Transect_Length[Gobi$Block == b]))
}
argalijitter <- data.frame(x = append(xjitter,missed$distance), y = append(yjitter,missed$forward),
                         stratum = argali$Block,
                         transect = argali$Site,
                         L = La,
                         area = Aa,
                         object = argali$id,
                         size = argali$Group.size)

