
# Set Up Functions --------------------------------------------------------
library(circular)
library(dplyr)
library(mrds)
library(secr)
library(Distance)
library(knitr)
ystart = 1700; w = 1600

movement <- function(N){ #random movement
  animals <- data.frame("id" = 1:N, "x" = runif(N,-w,w), "y" = runif(N,0,10000))
  angle <- rwrappedcauchy(N, mu = circular(0),rho = 0)
  distance <- abs(rnorm(N, 300, 100))
  animals$newx <- as.numeric(animals$x + distance*cos(angle))
  animals$newy <- as.numeric(animals$y + distance*sin(angle))
  return(animals)
}

reaction <- function(N, transect,theta){
  animals <- data.frame("id" = 1:N, "x" = runif(N,-w,w), "y" = runif(N,0,10000))
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
low <- c(logit(0.6),0.0004,log(300),0.001, 100)
med <- c(logit(0.75),0.0006,log(400),0.00125, 100)
high <- c(logit(0.9),0.0008,log(500),0.0015, 100)

avoidance_demo <- function(Nanimals, transect = 0){
  par(mfrow = c(2,2))
  low_avoid <- reaction(Nanimals, transect,low)
  med_avoid <- reaction(Nanimals, transect,med)
  high_avoid <- reaction(Nanimals, transect,high)
  plot(low_avoid$x, low_avoid$y, xlab = "x", ylab = "y", main = "Original Animal Positions")
  plot(low_avoid$newx, low_avoid$newy,xlab = "x", ylab = "y", main = "Low Avoidance Movement")
  abline(v = transect, col = "blue", lty = "dashed")
  plot(med_avoid$newx, med_avoid$newy,xlab = "x", ylab = "y", main = "Medium Avoidance Movement")
  abline(v = transect, col = "blue", lty = "dashed")
  plot(high_avoid$newx, high_avoid$newy,xlab = "x", ylab = "y", main = "High Avoidance Movement")
  abline(v = transect, col = "blue", lty = "dashed")
}
half_normal <- function(x,sigma){
  return(exp(-(x**2)/2*sigma**2))
}
hazard_rate <- function(x,sigma, beta){
  return(1-exp(-(x/sigma)**-beta))
}

# p(0) < 1 ----------------------------------------------------------------
recapture_p0 <- function(animals,transect = 0,
                                detectfn = "hr", hn_param = 50, hr_param = c(374.7160,2.2416),
                                p0 = 0.8){
  animals$distance<- abs(animals$x-transect)
  animals$newdistance <- abs(animals$newx-transect)
  if(detectfn =="hn"){
    detect_prob1 <- p0*half_normal(animals$distance,hn_param)
    detect_prob2 <- p0*half_normal(animals$newdistance,hn_param)
  }else if(detectfn =="hr"){
    detect_prob1 <- p0*hazard_rate(animals$distance, hr_param[1], hr_param[2])
    detect_prob2 <- p0*hazard_rate(animals$newdistance, hr_param[1], hr_param[2])
  }
  animals$obs1 <- rbinom(length(animals$id),1,detect_prob1)
  animals$obs2 <- rbinom(length(animals$id),1,detect_prob2)
  animals$obs2[animals$newdistance > w] <- 0
  return(animals)
}
fit_mrds_model <- function(data){
  modelaics <- c()
  model1 <- ddf(method = "io", dsmodel =~cds(key ="hn"),
                mrmodel =~glm(link = "logit", formula = ~distance),
                data = data, control = list(refit = T, nrefit = 5, debug = T))
  modelaics[1] <- model1$criterion
  model2 <- ddf(method = "io", dsmodel =~cds(key ="hr"),
                mrmodel =~glm(link = "logit", formula = ~distance),
                data = data, control = list(refit = T, nrefit = 5, debug = T))
  modelaics[2] <- model2$criterion
  
  if(which.min(modelaics) == 1){
    return(model1)
  }else{
    return(model2)
  }
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
  data$distance[data$distance > w] <- w
  model <- fit_mrds_model(data)
  return(model)
}

get_p0_bias <- function(Nrep, Nanimals, detectfn = "hr",
                        hn_param = 5, hr_param = c(374.7160,2.2416),
                        p0 = 0.8){
  estimates <- matrix(nrow = Nrep, ncol = 4, 
                      dimnames = list(1:Nrep, c("No","Low","Med","High")))
  for(i in 1:Nrep){
    for(j in 1:4){
      if(j == 1){
        sim <- data.frame("id" = 1:Nanimals, "x" = runif(Nanimals,-w,w), "y" = runif(Nanimals,0,ystart))
        sim$newx <- sim$x
        sim$newy <- sim$y
        animals <- recapture_p0(sim,detectfn = detectfn, hr_param = hr_param,hn_param = hn_param, p0 = p0)
        try({
          model <- get_p0_model(animals)
          estimates[i,j] <- model$Nhat
        })
      }else{
        if(j == 2){avoid = low}
        if(j == 3){avoid = med}
        if(j == 4){avoid = high}
        sim <- reaction(Nanimals,0,avoid)
        animals <- recapture_p0(sim,detectfn = detectfn, hr_param = hr_param, hn_param = hn_param, p0 = p0)
        try({
          model <- get_p0_model(animals)
          estimates[i,j] <- model$Nhat
        })
      }
    }
  }
  bias <- 100*(colMeans(estimates, na.rm = T) - Nanimals)/(Nanimals)
  return(bias)
}

# Unmodelled Heterogeneity ------------------------------------------------
heterogeneity <- function(Nanimals, transect, hr_param = c(150,2.2416)){
  animals <- data.frame("id" = 1:Nanimals, "x" = runif(Nanimals,-w,w),
                        "y" = runif(Nanimals,0,10000), "size" = sample(1:5, Nanimals, replace = T))
  animals$distance <- abs(animals$x - transect)
  detectprob <- hazard_rate(animals$distance, hr_param[1] + 100*animals$size, hr_param[2])
  animals$obs1 <- rbinom(Nanimals, 1, detectprob)
  animals$obs2 <- rbinom(Nanimals, 1, detectprob)
  seen_animals <- animals[animals$obs1 == 1 |animals$obs2 == 1,]
  data <- data.frame("object" = rep(seen_animals$id, each = 2),
                     "observer" = rep(c(1,2), times = length(seen_animals$id)),
                     "detected" = rep(0, 2*length(seen_animals$id)),
                     "distance" = rep(0, 2*length(seen_animals$id)))
  
  data$detected[data$observer == 1] <- seen_animals$obs1
  data$detected[data$observer == 2] <- seen_animals$obs2
  
  data$distance[data$observer == 1] <- seen_animals$distance
  data$distance[data$observer == 2] <- seen_animals$distance
  data$distance[data$distance > w] <- w
  model <- fit_mrds_model(data)
  return(model)
}

# Imperfect Matching ------------------------------------------------------
recapture_imperfect <- function(animals,transect = 0,
                                detectfn = "hr", hn_param = 5, hr_param = c(374.7160,2.2416),
                                b = c(4.89553757, 0.03561984),
                                match_limits = c(0.3,1), match_param = c(0.04,4)){
  animals$distance<- abs(animals$x-transect)
  animals$newdistance <- abs(animals$newx-transect)
  if(detectfn =="hn"){
    detect_prob1 <- half_normal(animals$distance,hn_param)
    detect_prob2 <- half_normal(animals$newdistance,hn_param)
  }else if(detectfn =="hr"){
    detect_prob1 <- hazard_rate(animals$distance, hr_param[1], hr_param[2])
    detect_prob2 <- hazard_rate(animals$newdistance, hr_param[1], hr_param[2])
  }else{
    detect_prob1 <- p.approx(y = seq(0,1700, length.out = 200), x = animals$distance,
                                h.fun = detectfn, b = b, what = "px")
    detect_prob2 <- p.approx(y = seq(0,1700, length.out = 200), x = animals$newdistance,
                                h.fun = detectfn, b = b, what = "px")
  }
  animals$true.obs1 <- rbinom(length(animals$id),1,detect_prob1)
  animals$true.obs2 <- rbinom(length(animals$id),1,detect_prob2)
  #observers can't see animals that moved outside region
  animals$true.obs2[animals$newdistance > w] <- 0
  
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
    if(length(closest) != 1){next}
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

# MRDS --------------------------------------------------------------------

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
  perf_data$distance[perf_data$distance > w] <- w
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
  imperf_data$distance[imperf_data$distance >w] <- w
  model <- fit_mrds_model(imperf_data)
  return(model)
}
get_bias <- function(Nrep, Nanimals){
  estimates <- matrix(nrow = Nrep, ncol = 8, 
                      dimnames = list(1:Nrep, c("No Perfect","No Imperfect",
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
  return(estimates)
}

# 2D Distance -------------------------------------------------------------
library(LT2D)
#Fits 2D model and no movement,random movement and random movement w/ imperfect matching
MCR_2D <- function(Nrep, Nanimals, 
                   b = c(4.89553757, 0.03561984), logphi = c(6.61348904, 4.84269089),
                   w = w, ystart = ystart, L = 100,
                   match_limits = c(300,1000), match_param = c(250,3)){
  ests_2d <- c();ci_l <- rep(NA, Nrep);ci_u <- rep(NA, Nrep)
  chapman.nomvmnt <- c();nomvmnt_l <- c();nomvmnt_u <- c()
  chapman.mvmnt <- c();mvmnt_l <- c(); mvmnt_u <- c()
  chapman.imperf <- c(); imperf_l <- c(); imperf_u <- c()
  for(j in 1:Nrep){
    #simulate animals
    positions <- simpop2DLT(L,w = w, pi.x = pi.hnorm, logphi = logphi, 
                         En = 200, fixed.n = T) 
    simDat <- detect2DLT(positions$x, hr = ip0, b= b, ystart = ystart, ny = 1000)
    
    
    all.1s <- rep(1,length(simDat$x))
    obj <- 1:length(simDat$x)
    sim.df <- data.frame(x = simDat$x,
                         y = simDat$y,
                         stratum = all.1s,
                         transect = all.1s,
                         L = L,
                         area = 2*w*L,
                         object = obj,
                         size = all.1s)
    fit <- LT2D.fit(DataFrameInput = sim.df,
                    hr = 'ip0',
                    b = b,
                    ystart = ystart,
                    pi.x = 'pi.hnorm',
                    logphi = logphi,
                    w = w,
                    hessian = TRUE)
    ests_2d[j] <- fit$ests[nrow(fit$ests),ncol(fit$ests)]
    try({boot <- LT2D.bootstrap(fit)
    ci_l[j] <- boot$ci[1]
    ci_u[j] <- boot$ci[2]
    }, silent = TRUE)
    #MCR
    ##format data
    mcr.df <- data.frame("id" =1:Nanimals, "x" = positions$x, "y" = positions$y,
                         "obs1" = rep(0, Nanimals),
                         "obs2" = rep(0, Nanimals),
                         "both" = rep(0, Nanimals))
    
    #assumes animals are uniform in y direction 
    mcr.df$obs1[mcr.df$x %in% simDat$x] <- 1
    ##get detection probabilities for second observer
    xs <- mcr.df$x
    ys <- seq(0,ystart, length.out = 200)
    detect_probs <- p.approx(ys, xs, ip0, b, what = "px")
    mcr.df$obs2 <- rbinom(Nanimals, 1, detect_probs)
    mcr.df$both[mcr.df$obs1 == 1 & mcr.df$obs2 == 1] <- 1
    n1 = sum(mcr.df$obs1); n2 = sum(mcr.df$obs2); m2 = sum(mcr.df$both)
    chapman.nomvmnt[j] <- (n1+1)*(n2+ 1)/(m2+1) - 1
    varNhat <- (n1+1)*(n2+ 1)*(n1-m2)*(n2-m2)/((m2+2)*(m2+1)**2)
    d <- exp(1.96*sqrt(log(1 + varNhat/chapman.nomvmnt[j]**2)))
    nomvmnt_l[j] <- chapman.nomvmnt[j]/d
    nomvmnt_u[j] <- chapman.nomvmnt[j]*d
    
    #add movement
    angle <- rwrappedcauchy(Nanimals, mu = circular(0),rho = 0)
    distance <- abs(rnorm(Nanimals, 300, 100))
    mcr.df$newx <- as.numeric(mcr.df$x + distance*cos(angle))
    mcr.df$newy <- as.numeric(mcr.df$y + distance*sin(angle))
    
    #get observations after movement with perfect matching
    xs <- mcr.df$newx
    ys <- seq(0,10000, length.out = 200)
    detect_probs <- p.approx(ys, xs, ip0, b, what = "px")
    mcr.df$obs2.mvmnt <- rbinom(Nanimals, 1, detect_probs)
    mcr.df$both.mvmnt[mcr.df$obs1 == 1 & mcr.df$obs2.mvmnt == 1] <- 1
    
    n1 = sum(mcr.df$obs1); n2 = sum(mcr.df$obs2.mvmnt, na.rm = T); m2 = sum(mcr.df$both.mvmnt, na.rm = T)
    chapman.mvmnt[j] <- (n1+1)*(n2+ 1)/(m2+1)
    varNhat <- (n1+1)*(n2+ 1)*(n1-m2)*(n2-m2)/((m2+2)*(m2+1)**2)
    d <- exp(1.96*sqrt(log(1 + varNhat/chapman.mvmnt[j]**2)))
    mvmnt_l[j] <- chapman.mvmnt[j]/d
    mvmnt_u[j] <- chapman.mvmnt[j]*d
    
    #imperfect matching
    distances <- edist(cbind(mcr.df$x, mcr.df$y), cbind(mcr.df$newx[mcr.df$obs2.mvmnt == 1], mcr.df$newy[mcr.df$obs2.mvmnt == 1]))
    firstobs <- mcr.df$id[mcr.df$obs1 == 1]
    mcr.df$imperfect.obs1 <- mcr.df$obs1
    mcr.df$imperfect.obs2 <- rep(0,length(mcr.df$id))
    for(i in 1:dim(mcr.df[mcr.df$obs2.mvmnt == 1,])[1]){
      closest <- which.min(distances[firstobs,i])
      closest <- firstobs[closest]
      if(length(closest)!= 1){next}
      if(distances[closest,i] > match_limits[2]){ #moved too far away, observers think two different animals
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
    
    n1 = sum(mcr.df$imperfect.obs1, na.rm = T); n2 = sum(mcr.df$imperfect.obs2, na.rm = T); m2 = sum(mcr.df$imperf.both, na.rm = T)
    chapman.imperf[j] <- (n1+1)*(n2+ 1)/(m2+1)
    varNhat <- (n1+1)*(n2+ 1)*(n1-m2)*(n2-m2)/((m2+2)*(m2+1)**2)
    d <- exp(1.96*sqrt(log(1 + varNhat/chapman.imperf[j]**2)))
    imperf_l[j] <- chapman.imperf[j]/d
    imperf_u[j] <- chapman.imperf[j]*d
}
  return(data.frame("2D" = ests_2d, "2D.lower" = ci_l,
                    "2D.upper" = ci_u, "MCR.nomvmnt" = chapman.nomvmnt,
                    "nomvmnt.lower" = nomvmnt_l, "nomvmnt.upper" = nomvmnt_u,
                    "MCR.move" = chapman.mvmnt, "mvmnt.lower" = mvmnt_l,
                    "mvmnt.upper" = mvmnt_u, "MCR.imperf" = chapman.imperf,
                    "imperf.lower" = imperf_l, "imperf.upper" = imperf_u))
}


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

