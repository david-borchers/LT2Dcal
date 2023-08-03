
# Set Up Functions --------------------------------------------------------
library(circular)
library(dplyr)
library(mrds)
library(secr)
library(Distance)
movement <- function(N){ #random movement
  animals <- data.frame("id" = 1:N, "x" = runif(N,0,2), "y" = runif(N,0,2))
  angle <- rwrappedcauchy(N, mu = circular(0),rho = 0)
  distance <- abs(rnorm(N, 0.05, 0.05))
  animals$newx <- as.numeric(animals$x + distance*cos(angle))
  animals$newy <- as.numeric(animals$y + distance*sin(angle))
  return(animals)
}

reaction <- function(N, transect,theta){
  animals <- data.frame("id" = 1:N, "x" = runif(N,0,2), "y" = runif(N,0,2))
  animals$angle <- rep(0,N)
  rho <- invlogit(theta[1] - theta[2]*abs(animals$x-transect))
  dist <- exp(theta[3] - theta[4]*abs(animals$x-transect))
  for(i in animals$id[animals$x <= transect]){ #using for loops as rwrappedcauchy doesn't work with rho vectors 
    animals$angle[i] <- as.numeric(rwrappedcauchy(1, mu = circular(pi), rho = rho[i]))
  }
  for(j in animals$id[animals$x > transect]){
    animals$angle[j] <- as.numeric(rwrappedcauchy(1, mu = circular(0), rho = rho[j]))
  }
  animals$distancemoved <- abs(rnorm(N, mean = dist, sd = theta[5]))
  animals$newx <- animals$x + animals$distancemoved*cos(animals$angle)
  animals$newy <- animals$y + animals$distancemoved*sin(animals$angle)
  return(animals)
} 
#parameters for different avoidance levels
low <- c(logit(0.6),0.4,log(0.1),1, 0.05)
med <- c(logit(0.75),0.6,log(0.2),1.25, 0.05)
high <- c(logit(0.9),0.8,log(0.25),1.5, 0.05)

avoidance_demo <- function(Nanimals, transect){
  par(mfrow = c(2,2))
  low_avoid <- reaction(Nanimals, transect,low)
  med_avoid <- reaction(Nanimals, transect,med)
  high_avoid <- reaction(Nanimals, transect,high)
  plot(low_avoid$x, low_avoid$y, xlab = "x", ylab = "y", main = "Original Animal Positions")
  plot(low_avoid$newx, low_avoid$newy,xlab = "x", ylab = "y", main = "Low Avoidance Movement")
  abline(v = 1, col = "blue", lty = "dashed")
  plot(med_avoid$newx, med_avoid$newy,xlab = "x", ylab = "y", main = "Medium Avoidance Movement")
  abline(v = 1, col = "blue", lty = "dashed")
  plot(high_avoid$newx, high_avoid$newy,xlab = "x", ylab = "y", main = "High Avoidance Movement")
  abline(v = 1, col = "blue", lty = "dashed")
}
half_normal <- function(x,sigma){
  return(exp(-(x**2)/2*sigma**2))
}
hazard_rate <- function(x,sigma, beta){
  return(1-exp(-(x/sigma)**-beta))
}

# p(0) < 1 ----------------------------------------------------------------
recapture_p0 <- function(animals,transect = 1,
                                detectfn = "hr", hn_param = 5, hr_param = c(0.8,6),
                                p0 = 0.8){
  animals$distance<- abs(animals$x-transect)
  animals$newdistance <- abs(animals$newx-transect)
  if(detectfn =="hn"){
    detect_prob1 <- p0*half_normal(animals$distance,hn_param)
    detect_prob2 <- p0*half_normal(animals$newdistance,hn_param)
  }
  if(detectfn =="hr"){
    detect_prob1 <- p0*hazard_rate(animals$distance, hr_param[1], hr_param[2])
    detect_prob2 <- p0*hazard_rate(animals$newdistance, hr_param[1], hr_param[2])
  }
  animals$obs1 <- rbinom(length(animals$id),1,detect_prob1)
  animals$obs2 <- rbinom(length(animals$id),1,detect_prob2)
  return(animals)
}
get_p0_model <- function(animals){
  seen_animals <- animals[animals$obs1 == 1 |animals$obs2 == 1,]
  data <- data.frame("object" = rep(seen_animals$id, each = 2),
                     "observer" = rep(c(1,2), times = length(seen_animals$id)),
                     "detected" = rep(0, 2*length(seen_animals$id)),
                     "distance" = rep(0, 2*length(seen_animals$id)))
  
  data$detected[data$observer == 1] <- seen_animals$obs1
  data$detected[data$observer == 2] <- seen_animals$obs2
  
  data$distance[data$observer == 1] <- seen_animals$distance
  data$distance[data$observer == 2] <- seen_animals$newdistance
  data$distance[data$observer == 1][seen_animals$obs1 == 0 & seen_animals$obs2 == 1] <- seen_animals$newdistance[seen_animals$obs1 == 0 & seen_animals$obs2 == 1]
  data$distance[data$distance > 1] <- 1
  model <- fit_mrds_model(data)
  return(model)
}

get_p0_bias <- function(Nrep, Nanimals){
  estimates <- matrix(nrow = Nrep, ncol = 4, 
                      dimnames = list(1:Nrep, c("No","Low","Med","High")))
  for(i in 1:Nrep){
    for(j in 1:4){
      if(j == 1){
        sim <- animals <- data.frame("id" = 1:Nanimals, "x" = runif(Nanimals,0,2), "y" = runif(Nanimals,0,2))
        sim$newx <- sim$x
        sim$newy <- sim$y
        animals <- recapture_p0(sim,hr_param = hr_param,p0 = 0.8)
        model <- get_p0_model(animals)
        estimates[i,j] <- model$Nhat
      }else{
        if(j == 2){avoid = low}
        if(j == 3){avoid = med}
        if(j == 4){avoid = high}
        sim <- reaction(Nanimals,1,avoid)
        animals <- recapture_p0(sim,hr_param = hr_param,p0 = p0)
        model <- get_p0_model(animals)
        estimates[i,j] <- model$Nhat
      }
    }
  }
  bias <- (colMeans(estimates, na.rm = T) - Nanimals)/(Nanimals/100)
  return(bias)
}

# Imperfect Matching ------------------------------------------------------
recapture_imperfect <- function(animals,transect = 1,
                                detectfn = "hr", hn_param = 5, hr_param = c(0.8,6),
                                match_limits = c(0.1,0.5), match_param = c(0.04,4)){
  animals$distance<- abs(animals$x-transect)
  animals$newdistance <- abs(animals$newx-transect)
  if(detectfn =="hn"){
    detect_prob1 <- half_normal(animals$distance,hn_param)
    detect_prob2 <- half_normal(animals$newdistance,hn_param)
  }
  if(detectfn =="hr"){
    detect_prob1 <- hazard_rate(animals$distance, hr_param[1], hr_param[2])
    detect_prob2 <- hazard_rate(animals$newdistance, hr_param[1], hr_param[2])
  }
  animals$true.obs1 <- rbinom(length(animals$id),1,detect_prob1)
  animals$true.obs2 <- rbinom(length(animals$id),1,detect_prob2)
  #observers can't see animals that moved outside region
  animals$true.obs2[animals$newdistance > 1] <- 0
  
  animals$imperfect.obs1 <- animals$true.obs1
  animals$imperfect.obs2 <- rep(0,length(animals$id))
  #mismatched animals will have 2 new distances, 1 for true position of original animals
  #second for position of misidentified second animals
  animals$imperfect.distance <- rep(0,length(animals$id))
  animals$imperfect.x<- rep(0,length(animals$id))
  animals$imperfect.y<- rep(0,length(animals$id))
  
  distances <- edist(cbind(animals$x, animals$y), cbind(animals$newx[animals$true.obs2 == 1], animals$newy[animals$true.obs2 == 1]))
  firstobs <- animals$id[animals$true.obs1 == 1]
  
  for(i in 1:dim(animals[animals$true.obs2 == 1,])[1]){
    closest <- which.min(distances[firstobs,i])
    closest <- firstobs[closest]
    if(is.na(closest)){next}
    if(distances[closest,i] > match_limits[2]){ #animal is too far away, observers think two different animals
      animals <- add_row(animals, id = animals$id[animals$true.obs2 == 1][i] + 1000,
                         x = animals$x[animals$true.obs2 == 1][i], y = animals$y[animals$true.obs2 == 1][i],
                         newx = animals$newx[animals$true.obs2 == 1][i], newy = animals$newy[animals$true.obs2 == 1][i],
                         distance = animals$distance[animals$true.obs2 == 1][i], newdistance = animals$newdistance[animals$true.obs2 == 1][i],
                         true.obs1 = 0, true.obs2 = 0,
                         imperfect.obs1 = 0, imperfect.obs2 = 1)
    }else{
      if(distances[closest,i] <= match_limits[1]){
        animals$imperfect.obs2[animals$id == closest][1] <- 1
      }else{
      match_prob <- hazard_rate(distances[closest,i]-match_limits[1],match_param[1],match_param[2])
      animals$imperfect.obs2[animals$id == closest][1] <- rbinom(1,1, match_prob)
      }
      #if matched add imperfect distance
      if(animals$imperfect.obs2[animals$id == closest][1] == 1){
        animals$imperfect.distance[animals$id == closest][1] <- animals$newdistance[animals$true.obs2 == 1][i]
        animals$imperfect.x[animals$id == closest][1] <- animals$newx[animals$true.obs2 == 1][i]
        animals$imperfect.y[animals$id == closest][1] <- animals$newy[animals$true.obs2 == 1][i]
      }
      #if not matched add new observation
      if(animals$imperfect.obs2[animals$id == closest][1] == 0){
        animals <- add_row(animals, id = animals$id[animals$true.obs2 == 1][i] + 1000,
                           x = animals$x[animals$true.obs2 == 1][i], y = animals$y[animals$true.obs2 == 1][i],
                           newx = animals$newx[animals$true.obs2 == 1][i], newy = animals$newy[animals$true.obs2 == 1][i],
                           distance = animals$distance[animals$true.obs2 == 1][i], newdistance = animals$newdistance[animals$true.obs2 == 1][i],
                           true.obs1 = 0, true.obs2 = 0,
                           imperfect.obs1 = 0, imperfect.obs2 = 1)
      } 
      #remove closest from future checks
      firstobs <- firstobs[!firstobs == closest]
    }
  }
  animals <- arrange(animals, id)
  return(animals)
}
fit_mrds_model <- function(data){
  modelaics <- c()
  model1 <- ddf(method = "io", dsmodel =~cds(key ="hn"),
                mrmodel =~glm(link = "logit", formula = ~distance),
                data = data, meta.data = list(width = 1), control = list(refit = T, nrefit = 2, debug = T))
  modelaics[1] <- model1$criterion
  model2 <- ddf(method = "io", dsmodel =~cds(key ="hr"),
                mrmodel =~glm(link = "logit", formula = ~distance),
                data = data, meta.data = list(width = 1), control = list(refit = T, nrefit = 2, debug = T))
  modelaics[2] <- model2$criterion
  if(which.min(modelaics) == 1){
    return(model1)
  }else{
    return(model2)
  }
}


# MRDS --------------------------------------------------------------------
fit_mrds_model <- function(data){
  modelaics <- c()
  model1 <- ddf(method = "io", dsmodel =~cds(key ="hn"),
                mrmodel =~glm(link = "logit", formula = ~distance),
                data = data, meta.data = list(width = 1), control = list(refit = T, nrefit = 2, debug = T))
  modelaics[1] <- model1$criterion
  model2 <- ddf(method = "io", dsmodel =~cds(key ="hr"),
                mrmodel =~glm(link = "logit", formula = ~distance),
                data = data, meta.data = list(width = 1), control = list(refit = T, nrefit = 2, debug = T))
  modelaics[2] <- model2$criterion
  if(which.min(modelaics) == 1){
    return(model1)
  }else{
    return(model2)
  }
}

get_perfect_model <- function(animals){
  perf_animals <- animals[animals$true.obs1 == 1 |animals$true.obs2 == 1,]
  perf_data <- data.frame("object" = rep(perf_animals$id, each = 2),
                          "observer" = rep(c(1,2), times = length(perf_animals$id)),
                          "detected" = rep(0, 2*length(perf_animals$id)),
                          "distance" = rep(0, 2*length(perf_animals$id)))
  
  perf_data$detected[perf_data$observer == 1] <- perf_animals$true.obs1
  perf_data$detected[perf_data$observer == 2] <- perf_animals$true.obs2
  
  perf_data$distance[perf_data$observer == 1] <- perf_animals$distance
  perf_data$distance[perf_data$observer == 2] <- perf_animals$newdistance
  perf_data$distance[perf_data$observer == 1][perf_animals$true.obs1 == 0 & perf_animals$true.obs2 == 1] <- perf_animals$newdistance[perf_animals$true.obs1 == 0 & perf_animals$true.obs2 == 1]
  perf_data$distance[perf_data$distance > 1] <- 1
  model <- fit_mrds_model(perf_data)
  return(model)
}
get_imperfect_model <- function(animals){
  imperf_animals <- animals[animals$imperfect.obs1 == 1 |animals$imperfect.obs2 == 1,]
  imperf_data <- data.frame("object" = rep(imperf_animals$id, each = 2),
                            "observer" = rep(c(1,2), times = length(imperf_animals$id)),
                            "detected" = rep(0, 2*length(imperf_animals$id)),
                            "distance" = rep(0, 2*length(imperf_animals$id)))
  
  imperf_data$detected[imperf_data$observer == 1] <- imperf_animals$imperfect.obs1
  imperf_data$detected[imperf_data$observer == 2] <- imperf_animals$imperfect.obs2
  
  imperf_data$distance[imperf_data$observer == 1] <- imperf_animals$distance
  imperf_data$distance[imperf_data$observer == 2] <- imperf_animals$newdistance
  imperf_data$distance[imperf_data$observer == 1][imperf_animals$imperfect.obs1 == 0 & imperf_animals$imperfect.obs2 == 1] <- imperf_animals$newdistance[imperf_animals$imperfect.obs1 == 0 & imperf_animals$imperfect.obs2 == 1]
  imperf_data$distance[imperf_data$distance >1] <- 1
  model <- fit_mrds_model(imperf_data)
  return(model)
}
get_bias <- function(Nrep, Nanimals){
  estimates <- matrix(nrow = Nrep, ncol = 8, dimnames = list(1:Nrep, c("No Perfect","No Imperfect",
                                                                       "Low Perfect", "Low Imperfect",
                                                                       "Med Perfect", "Med Imperfect",
                                                                       "High Perfect", "High Imperfect")))
  for(i in 1:Nrep){
    for(j in 1:4){
      if(j == 1){
        sim <- movement(Nanimals)
        sim$newx <- sim$x
        sim$newy <- sim$y
        animals <- recapture_imperfect(sim,hr_param = hr_param,
                                       match_limits =  match_limits,match_param = match_param)
        perf_model <- get_perfect_model(animals)
        imperf_model <- get_imperfect_model(animals)
        estimates[i,2*j -1] <- perf_model$Nhat
        estimates[i,2*j] <- imperf_model$Nhat
      }else{
      if(j == 2){avoid = low}
      if(j == 3){avoid = med}
      if(j == 4){avoid = high}
      sim <- reaction(Nanimals,1,avoid)
      animals <- recapture_imperfect(sim, hr_param = hr_param,
                                     match_limits =  match_limits,match_param = match_param)
      perf_model <- get_perfect_model(animals)
      imperf_model <- get_imperfect_model(animals)
      estimates[i,2*j -1] <- perf_model$Nhat
      estimates[i,2*j] <- imperf_model$Nhat
      }
    }
  }
  bias <- (colMeans(estimates, na.rm = T) - Nanimals)/(Nanimals/100)
  return(bias)
}

# 2D Distance -------------------------------------------------------------
#default params are best model for real data
library(LT2D)
#Fits 2D model and no movment,random movement and random movement w/ imperfect matching
MCR_2D <- function(Nrep, Nanimals, 
                   b = c(4.945767, 0.714222), logphi = c(6.520254,4.754897),
                   w = 1500, ystart = 1300,
                   match_limits = c(200,500), match_params = c(250,3)){
  ests_2d <- c()
  chapman.nomvmnt <- c()
  chapman.mvmnt <- c()
  chapman.imperf <- c()
  for(j in 1:Nrep){
    #simulate animals
    simDat = simXY(Nanimals, 'pi.hnorm',
                   logphi, 'h1', 
                   b, w, 
                   ystart, discardNotSeen = FALSE)
    seen <- simDat$locs[simDat$locs$y >= 0,]
    all.1s <- rep(1,length(seen$x))
    obj <- 1:length(seen$x)
    sim.df <- data.frame(x = seen$x,
                         y = seen$y,
                         stratum = all.1s,
                         transect = all.1s,
                         L = 10,
                         area = 2*w*10,
                         object = obj,
                         size = all.1s)
    fit <- LT2D.fit(DataFrameInput = sim.df,
                    hr = 'h1',
                    b = b,
                    ystart = ystart,
                    pi.x = 'pi.hnorm',
                    logphi = logphi,
                    w = w,
                    hessian = TRUE)
    ests_2d[j] <- fit$ests[nrow(fit$ests),ncol(fit$ests)]
    #MCR
    ##format data
    mcr.df <- data.frame("id" =1:Nanimals, "x" = simDat$locs$x, "y" = runif(length(simDat$locs$x),0,10000),
                         "obs1" = rep(0, length(simDat$locs$x)),
                         "obs2" = rep(0, length(simDat$locs$x)),
                         "both" = rep(0, length(simDat$locs$x)))
    #y in simdat represents distance from observer, we need y co ordinate in space
    #assumes animals are uniform in y direction 
    mcr.df$obs1[simDat$locs$y >= 0] <- 1
    ##get detection probabilities for second observer
    xs <- mcr.df$x
    ys <- seq(0,10000, length.out = 200)
    detect_probs <- p.approx(ys, xs, h1, b, xy=TRUE)
    mcr.df$obs2 <- rbinom(Nanimals, 1, detect_probs)
    mcr.df$both[mcr.df$obs1 == 1 & mcr.df$obs2 == 1] <- 1
    chapman.nomvmnt[j] <- (sum(mcr.df$obs1)+1)*(sum(mcr.df$obs2)+ 1)/(sum(mcr.df$both)+1)
    
    #add movement and imperfect matching
    angle <- rwrappedcauchy(Nanimals, mu = circular(0),rho = 0)
    distance <- abs(rnorm(Nanimals, 100, 50))
    mcr.df$newx <- as.numeric(mcr.df$x + distance*cos(angle))
    mcr.df$newy <- as.numeric(mcr.df$y + distance*sin(angle))
    
    xs <- mcr.df$newx
    ys <- seq(0,10000, length.out = 200)
    detect_probs <- p.approx(ys, xs, h1, b, xy=TRUE)
    mcr.df$obs2.mvmnt <- rbinom(Nanimals, 1, detect_probs)
    mcr.df$both.mvmnt[mcr.df$obs1 == 1 & mcr.df$obs2.mvmnt == 1] <- 1
    chapman.mvmnt[j] <- (sum(mcr.df$obs1)+1)*(sum(mcr.df$obs2.mvmnt)+ 1)/(sum(mcr.df$both.mvmnt, na.rm = T)+1)
    
    #imperfect matching
    distances <- edist(cbind(mcr.df$x, mcr.df$y), cbind(mcr.df$newx[mcr.df$obs2.mvmnt == 1], mcr.df$newy[mcr.df$obs2.mvmnt == 1]))
    firstobs <- mcr.df$id[mcr.df$obs1 == 1]
    mcr.df$imperfect.obs1 <- mcr.df$obs1
    mcr.df$imperfect.obs2 <- rep(0,length(mcr.df$id))
    for(i in 1:dim(mcr.df[mcr.df$obs2.mvmnt == 1,])[1]){
      closest <- which.min(distances[firstobs,i])
      closest <- firstobs[closest]
      if(is.na(closest)|length(closest)!= 1){next}
      if(distances[closest,i] > match_limits[2]){ #animal is too far away, observers think two different mcr.df
        mcr.df <- add_row(mcr.df, id = mcr.df$id[mcr.df$obs2.mvmnt == 1][i] + 1000,
                           x = mcr.df$x[mcr.df$obs2.mvmnt == 1][i], y = mcr.df$y[mcr.df$obs2.mvmnt == 1][i],
                           newx = mcr.df$newx[mcr.df$obs2.mvmnt == 1][i], newy = mcr.df$newy[mcr.df$obs2.mvmnt == 1][i],
                           obs1 = 0, obs2.mvmnt = 0,
                           imperfect.obs1 = 0, imperfect.obs2 = 1)
      }else{
        if(distances[closest,i] <= match_limits[1]){
          mcr.df$imperfect.obs2[mcr.df$id == closest][1] <- 1
        }else{
          match_prob <- hazard_rate(distances[closest,i],match_param[1],match_param[2])
          mcr.df$imperfect.obs2[mcr.df$id == closest][1] <- rbinom(1,1, match_prob)
        }
        #if not matched add new observation
        if(mcr.df$imperfect.obs2[mcr.df$id == closest][1] == 0){
          mcr.df <- add_row(mcr.df, id = mcr.df$id[mcr.df$obs2.mvmnt == 1][i] + 1000,
                             x = mcr.df$x[mcr.df$obs2.mvmnt == 1][i], y = mcr.df$y[mcr.df$obs2.mvmnt == 1][i],
                             newx = mcr.df$newx[mcr.df$obs2.mvmnt == 1][i], newy = mcr.df$newy[mcr.df$obs2.mvmnt == 1][i],
                             obs1 = 0, obs2.mvmnt = 0,
                             imperfect.obs1 = 0, imperfect.obs2 = 1)
        } 
        #remove closest from future checks
        firstobs <- firstobs[!firstobs == closest]
      }
    } 
    mcr.df$imperf.both <- rep(0, dim(mcr.df)[1])
    mcr.df$imperf.both[mcr.df$imperfect.obs1 == 1 & mcr.df$imperfect.obs2 == 1] <- 1
    chapman.imperf[j] <- (sum(mcr.df$imperfect.obs1)+1)*(sum(mcr.df$imperfect.obs2, na.rm = T)+ 1)/(sum(mcr.df$imperf.both, na.rm = T)+1)
  }
  return(data.frame("2D" = ests_2d, "MCR.nomvmnt" = chapman.nomvmnt, 
                    "MCR.move" = chapman.mvmnt, "MCR.imperf" = chapman.imperf))
}


