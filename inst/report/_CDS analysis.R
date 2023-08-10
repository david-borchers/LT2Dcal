hn_basic <- ds(CDSdf, truncation = w, key = "hn")
hn_species <- ds(CDSdf, truncation = w, formula = ~species, key = "hn")
hn_size <- ds(CDSdf, truncation = w, formula = ~size, key = "hn")
hn_block <- ds(CDSdf, truncation = w, formula = ~Region.Label, key = "hn")
hn_obs <- ds(CDSdf, truncation = w, formula = ~obs.name, key = "hn")
summarize_ds_models(hn_basic, hn_species, hn_block, hn_size, hn_obs)#prefer basic

hr_basic <- ds(CDSdf, truncation = w, key = "hr")
hr_species <- ds(CDSdf, truncation = w, formula = ~species, key = "hr")
hr_size <- ds(CDSdf, truncation = w, formula = ~size, key = "hr")
hr_block <- ds(CDSdf, truncation = w, formula = ~Region.Label, key = "hr")
hr_obs <- ds(CDSdf, truncation = w, formula = ~obs.name, key = "hr")
summarize_ds_models(hr_basic, hr_species, hr_block, hr_size) #prefer by block

hr_blockspec <- ds(CDSdf, truncation = w, formula = ~Region.Label+species, key = "hr")
hr_int <- ds(CDSdf, truncation = w, formula = ~Region.Label*species, key = "hr")
summarize_ds_models(hn_basic, hr_block, hr_blockspec, hr_int) #hr block is best


# Ibex only ---------------------------------------------------------------
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

hn_ibex<- ds(ibexdf, truncation = w, key = "hn")
hn_ibex_size <- ds(ibexdf, truncation = w, formula = ~size, key = "hn")
hn_ibex_block <- ds(ibexdf, truncation = w, formula = ~Region.Label, key = "hn")
hn_ibex_obs <- ds(ibexdf, truncation = w, formula = ~obs.name, key = "hn")
summarize_ds_models(hn_ibex, hn_ibex_block, hn_ibex_size, hn_ibex_obs)#prefer basic

hr_ibex<- ds(ibexdf, truncation = w, key = "hr")
hr_ibex_size <- ds(ibexdf, truncation = w, formula = ~size, key = "hr")
hr_ibex_block <- ds(ibexdf, truncation = w, formula = ~Region.Label, key = "hr")
hr_ibex_obs <- ds(ibexdf, truncation = w, formula = ~obs.name, key = "hr")
summarize_ds_models(hr_ibex, hr_ibex_block, hr_ibex_size, hr_ibex_obs)#prefer basic
#prefer block

summarize_ds_models(hn_ibex, hr_ibex_block) #prefer block

# Argali only -------------------------------------------------------------
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

hn_argali<- ds(argalidf, truncation = w, key = "hn")
hn_argali_size <- ds(argalidf, truncation = w, formula = ~size, key = "hn")
hn_argali_block <- ds(argalidf, truncation = w, formula = ~Region.Label, key = "hn")
hn_argali_obs <- ds(argalidf, truncation = w, formula = ~obs.name, key = "hn")
summarize_ds_models(hn_argali, hn_argali_block, hn_argali_size, hn_argali_obs)#prefer observer

hr_argali<- ds(argalidf, truncation = w, key = "hr")
hr_argali_size <- ds(argalidf, truncation = w, formula = ~size, key = "hr")
hr_argali_block <- ds(argalidf, truncation = w, formula = ~Region.Label, key = "hr")
hr_argali_obs <- ds(argalidf, truncation = w, formula = ~obs.name, key = "hr")
summarize_ds_models(hr_argali, hr_argali_block, hr_argali_size, hr_argali_obs)#prefer basic
#prefer block

summarize_ds_models(hn_argali_obs, hn_argali, hr_argali_block)

hr_argali_block$ddf$criterion + hr_ibex_block$ddf$criterion
hr_block$ddf$criterion

plot(hr_ibex_block)
plot(hr_block)
gof_ds(hr_argali_block)
gof_ds(hr_ibex_block)
gof_ds(hr_block)
