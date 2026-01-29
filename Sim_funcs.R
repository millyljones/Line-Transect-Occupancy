# This is the file with all the simulation functions in it
# We will call this file at the beginning of every simulation
# in order to get the correct functions.

# For each simulation, we will then separately call a file with
# all the relevant simulation parameters for that simulation

library(nimble)

##################### simData #################################

simulateData <- function(SimParams){
  
  list2env(SimParams,environment()) 
  # Output from SimParams
  # num_sites = number of zones/site
  # num_trans = number of transects per zone
  # num_segments = number of segments per transect
  # psi = probability of site occupancy  
  # init = occupancy probability for first segment of transect 
  # lambda = rate of detection exponentially distributed  
  # R = length of segment
  # theta11 = P(occupied to occupied segments)
  # theta01 = P(unoccupied to occupied segments)  
  # max_det_site = just needed a default for max number of detections per segment
  #               to make code faster, currently set at 1000.
  
  # Returns:
  # site.occ = vector length num_trans, true site occupancy 
  # Occ = matrix(num_trans x num_segments), true segment occupancy
  # det_ = data for observation processes
  # All data are conditional on same Occ, site.occ. First generate
  # all detections in segment ~ exp(lambda1). Then discretise/simplify
  # as needed for the other processes. All values (other than the 
  # independent secondary processes) are based on this first data set.
  
  #---------------------------------------------------
  # Setting defaults and warnings
  
  lambda = lambda1 # set default
  
  # --------------------------------------------------
  
  # Site and segment occupancy processes
  site.occ = rep(NA, num_sites)
  Occ = matrix(NA, nrow = num_sites*num_trans, ncol = num_segments)
  for (i in 1:num_sites){
    site.occ[i] = rbinom(1, 1, psi)
    for (t in 1:num_trans){
      Occ[(i-1)*num_trans + t, 1] = site.occ[i] * rbinom(1, 1, init)
      for (j in 2:num_segments){
        Occ[(i-1)*num_trans + t,j] = site.occ[i] * rbinom(1, 1, 
                                        theta11*Occ[(i-1)*num_trans + t,j-1] + 
                                          theta01*(1-Occ[(i-1)*num_trans + t, j-1]))
      }
    }
  }
  
  
  # ---------------------------------------------------------------
  
  # generate exponential detections along each occupied segment
  site_det_hist = array(data = NA, dim = c(num_sites*num_trans, num_segments, max_det_site))
  n_detections_segment = c()
  for (i in 1:num_sites){ # looping through sites
    if (site.occ[i] == 1){ # if site is occupied
      for (t in 1:num_trans){ # loop through transects
        for (j in 1:num_segments){ # loop through segments
          if (Occ[(i-1)*num_trans + t,j] == 1){ # if segment is occupied
            current_num_det = 0 # current number of detections in this segment
            current_distance = 0
            while (current_distance < R){
              # generate a detection distance
              dist = rexp(1, lambda)
              current_distance = current_distance + dist
              if (current_distance < R){ # if inside segment
                current_num_det = current_num_det + 1
                site_det_hist[(i-1)*num_trans + t,j,current_num_det] = dist
              }
            } # end looping through detections
            n_detections_segment = c(n_detections_segment, current_num_det)
          } # end segment detections
        } # end looping through segments
      } # end transect
    } # end site
  } # end looping through sites
            
  # max number of detections per segment
  max_n = max(n_detections_segment)
  site_det_hist = site_det_hist[,,1:max_n] # truncate un-needed sections
  
  #------------------------------------------------------
  
  # Now go through and do various detection methods on all these results
  
  # Observation process DTD/1:-----------------------------------
  # taking only first detection distance
  det3 = site_det_hist[,,1] 
  
  # Observation process DND/1:-----------------------------------
  # 0,1 detection per segment
  det1 = 1*(!is.na(det3)) 
  
  # Observation process DTD/2:-----------------------------------
  # two continuous observations per segment
  
  # Using same detection rate both times
  det5 = array(NA, dim = c(num_sites*num_trans, num_segments, 2))
  det5[,,1] = det3
  for (i in 1:num_sites){
    for (t in 1:num_trans){
      for (j in 1:num_segments){
        if (Occ[(i-1)*num_trans + t,j] == 1){ # if we have an occupied segment
          det5[(i-1)*num_trans + t,j, 2] = rexp(1, lambda) # observe distance to second independent detection
          if (det5[(i-1)*num_trans + t,j,2] > R){
            det5[(i-1)*num_trans + t,j, 2]  = NA
          }
        }
      }
    }
  }
    
  # Observation process DND/2:-------------------------------------
  # two 0,1 detections per segment
  
  # If using same detection rate both times
  det2 = 1*(!is.na(det5)) 
  
  # Observation process C/1:------------------------------------------
  # number of detections per segment
  det4 = rowSums(1*(!is.na(site_det_hist)), dims = 2)
  
  # Return:---------------------------------------------------------
  return(list('Occ' = Occ,
              'site.occ' = site.occ,
              'det_alpha' = site_det_hist,
              'DTD1' = det3,
              'DND1' = det1,
              'DTD2' = det5,
              'DND2' = det2,
              'C1' = det4))
  
  
}

#################### dCT ######################################

# Have the density and random generating functions here

dexp_t <- nimbleFunction(
  # truncated exponential distribution
  # we already condition on the observation being before the truncation 
  # R.
  # Define the inputs:
  run = function(x = double(0), # scalar distance to first detection
                 lambda = double(0), # exponential rate
                 R = integer(0), # the truncation distance
                 occ = integer(0), # indicator (0,1) if site occupied
                 seg.occ = integer(0), # indicator (0,1) if segment is occupied
                 det.ind = integer(0), # indicator (0,1) if distance is a detection or tail
                 log = integer(0, default = 0)){
    
    returnType(double(0))
    
    if (seg.occ == 0 & det.ind == 1){
      stop('In nimbleFunction dexp_t : Simultaneously non occupied and detections')
    }
    
    if (x > R){
      stop('In nimbleFunction dexp_t : Observed distance is greater than segment length')
    }
    
    if (det.ind == 1){ # if this is a standard inter-detection distance
      
      # Prob detection occurs in less than R 
      prob_det = 1 - exp(-lambda*R)
      d_exp_t <- lambda*exp(-lambda*x) / prob_det
      
    } else { # if this is an observation with no detection 
      # i.e (distance to end of segment with no detections)
      
      d_exp_t <- 1  
      
    }
    
    if (log) {
      return(log(d_exp_t))
    }
    return(d_exp_t)
    
  }
)

rexp_t <- nimbleFunction(
  # random generating function. 
  run = function(n = integer(0), 
                 lambda = double(0),
                 R = integer(0),
                 occ = integer(0),
                 seg.occ = integer(0),
                 det.ind = integer(0)) {
  
    returnType(double(0))
  
  if(n != 1) print("rexp_t only allows n = 1; using n = 1.")
  
  if (occ == 0){ # if we have an un-occupied site
    x <- R
  } else { # if we do have an occupied site
    if (seg.occ == 0){ # if we have an unoccupied segment
      x <- R
    } else {# if we have an occupied segment
      if (det.ind == 1){ # if this is an observation
        x <- R + 1
        while (x > R){# need an observation that is less than R from truncated exponential
          x <- rexp(1, lambda)
        } 
      } else {
        x <- R
      }
    }
  }
  
  return(x)
  
})

##################### nimble models ###########################

Obs1d <- nimbleCode({
  # For DND1 - 1 observation discrete
  
  # priors
  theta11 ~ dunif(0,1)
  theta01 ~ dunif(0,1)
  p ~ dunif(0,1)
  psi ~ dunif(0,1)
  
  # different initial conditions
  if (stationary){
    init <- theta01/(theta01 + 1 - theta11)
  } else {
    init <- theta01
  }
  
  # latent segment occupancy process
  for (i in 1:num_sites){
    occ[i] ~ dbern(psi)
  }
  for (i in 1:(num_trans*num_sites)){
    seg.occ[i, 1] ~ dbern(init * occ[zone[i]])
    for (j in 2:num_segments){
      seg.occ[i,j] ~ dbern(occ[zone[i]] * (theta11*seg.occ[i,j-1] + 
                                             theta01*(1-seg.occ[i,j-1])))
    }
  }
  
  # Observation process
  for (i in 1:(num_trans*num_sites)){
    for (j in 1:num_segments){
      det[i,j] ~ dbern(seg.occ[i,j] * p)
    }
  }
  
})

Obs1c <- nimbleCode({
  # for DTD1
  
  # priors
  theta11 ~ dunif(0,1)
  theta01 ~ dunif(0,1)
  lambda ~ dgamma(scale = l.scale, shape = l.shape)
  p <- 1 - exp(-lambda*R) # prob of making an obs in segment
  psi ~ dunif(0,1)
  
  # different initial conditions
  if (stationary){
    init <- theta01/(theta01 + 1 - theta11)
  } else {
    init <- theta01
  }
  
  # latent segment occupancy process
  for (i in 1:num_sites){
    occ[i] ~ dbern(psi)
  }
  for (i in 1:(num_trans*num_sites)){
    seg.occ[i, 1] ~ dbern(init * occ[zone[i]])
    for (j in 2:num_segments){
      seg.occ[i,j] ~ dbern(occ[zone[i]] * (theta11*seg.occ[i,j-1] + 
                                             theta01*(1-seg.occ[i,j-1])))
    }
  }
  
  # Likelihood
  for (i in 1:(num_sites*num_trans)){
    for (j in 1:num_segments){
      det.ind[i,j] ~ dbern(seg.occ[i,j] * p) # we condition on probability of detection
      # on segment here, hence why dexp_t returns 1 if det.ind = 0
      det[i,j] ~ dexp_t(lambda = lambda, R = R, occ = occ[zone[i]],
                        seg.occ = seg.occ[i,j], det.ind = det.ind[i,j])
    }
  }
  
})

Obs2d <- nimbleCode({
  # for DND2
  
  # priors
  theta11 ~ dunif(0,1)
  theta01 ~ dunif(0,1)
  p ~ dunif(0,1)
  psi ~ dunif(0,1)
  
  # different initial conditions
  if (stationary){
    init <- theta01/(theta01 + 1 - theta11)
  } else {
    init <- theta01
  }
  
  # latent segment occupancy process
  for (i in 1:num_sites){
    occ[i] ~ dbern(psi)
  }
  for (i in 1:(num_trans*num_sites)){
    seg.occ[i, 1] ~ dbern(init * occ[zone[i]])
    for (j in 2:num_segments){
      seg.occ[i,j] ~ dbern(occ[zone[i]] * (theta11*seg.occ[i,j-1] + 
                                             theta01*(1-seg.occ[i,j-1])))
    }
  }
  
  # observation process
  for (i in 1:(num_sites*num_trans)){
    for (j in 1:num_segments){
      for (t in 1:2){
        det[i,j,t] ~ dbern(p * seg.occ[i,j])
      }
    }
  }
  
})

Obs2c <- nimbleCode({
  
  # priors
  theta11 ~ dunif(0,1)
  theta01 ~ dunif(0,1)
  lambda ~ dgamma(scale = l.scale, shape = l.shape) 
  p <- 1 - exp(-lambda*R) # prob of making an obs in segment
  psi ~ dunif(0,1)
  
  # different initial conditions
  if (stationary){
    init <- theta01/(theta01 + 1 - theta11)
  } else {
    init <- theta01
  }

  
  # latent segment occupancy process
  for (i in 1:num_sites){
    occ[i] ~ dbern(psi)
  }
  for (i in 1:(num_trans*num_sites)){
    seg.occ[i, 1] ~ dbern(init * occ[zone[i]])
    for (j in 2:num_segments){
      seg.occ[i,j] ~ dbern(occ[zone[i]] * (theta11*seg.occ[i,j-1] + 
                                             theta01*(1-seg.occ[i,j-1])))
    }
  }
  
  # Likelihood
  for (i in 1:(num_sites*num_trans)){
    for (j in 1:num_segments){
      for (t in 1:2){
        det.ind[i,j,t] ~ dbern(seg.occ[i,j] * p)
        det[i,j,t] ~ dexp_t(lambda = lambda, R = R, occ = occ[zone[i]],
                            seg.occ = seg.occ[i,j], det.ind = det.ind[i,j,t])
      }
    }
  }
  
})

Obs4d <- nimbleCode({
  # for C1
  
  # priors
  theta11 ~ dunif(0,1)
  theta01 ~ dunif(0,1)
  lambda ~ dgamma(shape = l.shape, scale = l.scale)
  psi ~ dunif(0,1)
  p <- 1 - exp(-lambda*R)
  
  # different initial conditions
  if (stationary){
    init <- theta01/(theta01 + 1 - theta11)
  } else {
    init <- theta01
  }

  # latent segment occupancy process
  for (i in 1:num_sites){
    occ[i] ~ dbern(psi)
  }
  for (i in 1:(num_trans*num_sites)){
    seg.occ[i, 1] ~ dbern(init * occ[zone[i]])
    for (j in 2:num_segments){
      seg.occ[i,j] ~ dbern(occ[zone[i]] * (theta11*seg.occ[i,j-1] + 
                                             theta01*(1-seg.occ[i,j-1])))
    }
  }
  
  # Observation process
  for (i in 1:(num_sites*num_trans)){
    for (j in 1:num_segments){
      det[i,j] ~ dpois(seg.occ[i,j] * lambda * R)
    }
  }
  
})

#################### Set model Parameters #####################

# set up the constants, data, and initial values for each model
# Each function will take only 'data', and then all other parameters
# will be taken from the global environment

set_model1d <- function(data, SimParams){
  
  list2env(SimParams,environment())
  
  zone <- rep(1:num_sites, each = num_trans)
  
  obs_occ_trans = 1*(rowSums(data$DND1) >= 1)
  obs_occ <- rep(NA, num_sites)
  for (i in 1:num_sites){
    obs_occ[i] <- 1*(sum(obs_occ_trans[which(zone == i)])>0)
  }
  obs_occ[obs_occ == 0] = NA # what sites do we know are occupied
  
  obs_seg.occ = data$DND1
  obs_seg.occ[obs_seg.occ == 0] = NA # what segments do we know are occupied
  
  Obs1dconstants <- list(num_sites = num_sites, num_segments = num_segments, R = R,
                         stationary = T, num_trans = num_trans, zone = zone)
  Obs1ddata = list(det = data$DND1, occ = obs_occ, seg.occ = obs_seg.occ)
  
  occ.init = rbinom(num_sites, 1, 0.5)
  occ.init[!is.na(obs_occ)] = NA
  
  seg.occ.init = data$Occ
  for (i in 1:(num_sites*num_trans)){
    if (!is.na(occ.init[zone[i]])){
      if (occ.init[zone[i]] == 0){
        seg.occ.init[i,] = 0 
      }
    }
  }
  seg.occ.init[!is.na(obs_seg.occ)] = NA
  
  Obs1dinits <- function(){
    list(theta11 = theta11, theta01 = theta01, 
         p = p, seg.occ = seg.occ.init, 
         psi = psi, occ = occ.init)
  }
  
  return(list(data = Obs1ddata,
              constants = Obs1dconstants,
              inits = Obs1dinits))
  
}

set_model1c <- function(data, SimParams){
  
  list2env(SimParams,environment())
  
  DTD1 = data$DTD1
  
  zone <- rep(1:num_sites, each = num_trans)
  
  # what sites can we see are occupied
  obs_occ_trans = apply(DTD1, 1, function(x) {1*any(!is.na(x))})
  obs_occ <- rep(NA, num_sites)
  for (i in 1:num_sites){
    obs_occ[i] <- 1*(sum(obs_occ_trans[which(zone == i)]) > 0)
  }
  obs_occ[obs_occ == 0] = NA
  
  # what segments can we see are occupied
  obs_seg.occ = 1*(!is.na(DTD1))
  obs_seg.occ[obs_seg.occ == 0] = NA
  
  det.ind = 1*(!is.na(data$DTD1))
  DTD1[is.na(DTD1)] = R # we let segments with no detections
  # simply have observation the segment length (the model knows this
  # is not a true observation because of det.ind, and will handle this
  # accordingly)

  Obs1cconstants <- list(num_sites = num_sites, num_segments = num_segments, R = R,
                         stationary = T, num_trans = num_trans,
                         zone = zone, l.shape = 1, l.scale = 1/R)
  Obs1cdata = list(det = DTD1, det.ind = det.ind, occ = obs_occ, seg.occ = obs_seg.occ)
  
  occ.init = rbinom(num_sites, 1, 0.5)
  occ.init[!is.na(obs_occ)] = NA
  
  seg.occ.init = data$Occ
  for (i in 1:(num_sites*num_trans)){
    if (!is.na(occ.init[zone[i]])){
      if (occ.init[zone[i]] == 0){
        seg.occ.init[i,] = 0 
      }
    }
  }
  seg.occ.init[!is.na(obs_seg.occ)] = NA
  
  Obs1cinits <- function(){
    list(theta11 = theta11, theta01 = theta01, 
         lambda = lambda, seg.occ = seg.occ.init, 
         psi = psi, occ = occ.init)
  }
  
  return(list(data = Obs1cdata,
              constants = Obs1cconstants,
              inits = Obs1cinits))
         
}

set_model2d <- function(data, SimParams){
  
  list2env(SimParams,environment())
  
  zone = rep(1:num_sites, each = num_trans)
  
  DND2 = data$DND2
  
  # what sites do we know are occupied
  obs_occ_trans = 1*(rowSums(DND2) >= 1)
  obs_occ <- rep(NA, num_sites)
  for (i in 1:num_sites){
    obs_occ[i] <- 1*(sum(obs_occ_trans[which(zone == i)]) > 0)
  }
  obs_occ[obs_occ == 0] = NA
  
  # what segments do we know are occupied
  obs_seg.occ = matrix(NA, nrow = num_sites*num_trans, ncol = num_segments)
  for (i in 1:(num_sites*num_trans)){
    for (j in 1:num_segments){
      if (sum(DND2[i,j,]) >= 1){
        obs_seg.occ[i,j] = 1
      }
    }
  }
  
  Obs2dconstants <- list(num_sites = num_sites, num_segments = num_segments, R = R,
                         num_trans = num_trans, zone = zone, stationary = T)
  Obs2ddata = list(det = DND2, occ = obs_occ, seg.occ = obs_seg.occ)
  
  occ.init = data$site.occ
  occ.init[!is.na(obs_occ)] = NA
  
  seg.occ.init = data$Occ
  for (i in 1:(num_sites*num_trans)){
    if (!is.na(occ.init[zone[i]])){
      if (occ.init[zone[i]] == 0){
        seg.occ.init[i,] = 0
      }
    }
  }
  seg.occ.init[!is.na(obs_seg.occ)] = NA
  
  Obs2dinits <- function(){
    list(theta11 = theta11, theta01 = theta01, 
         p = p, seg.occ = seg.occ.init, 
         psi = psi, occ = occ.init)
  }
  
  return(list(data = Obs2ddata,
              constants = Obs2dconstants,
              inits = Obs2dinits))
}
  
set_model2c <- function(data, SimParams){
  
  list2env(SimParams,environment())
  
  zone = rep(1:num_sites, each = num_trans)
  
  DTD2 = data$DTD2
  
  # what sites do we know are occupied
  obs_occ1_trans = apply(DTD2[,,1], 1, function(x) {1*any(!is.na(x))})
  obs_occ2_trans = apply(DTD2[,,2], 1, function(x) {1*any(!is.na(x))})
  obs_occ_trans <- obs_occ1_trans + obs_occ2_trans
  obs_occ <- rep(NA, num_sites)
  for (i in 1:num_sites){
    obs_occ[i] <- 1*(sum(obs_occ_trans[which(zone == i)]) > 0)
  }
  obs_occ[obs_occ == 0] = NA
  
  # what segments do we know are occupied
  obs_seg.occ = matrix(NA, nrow = num_sites*num_trans, ncol = num_segments)
  for (i in 1:(num_sites*num_trans)){
    for (j in 1:num_segments){
      if (any(!is.na(DTD2[i,j,]))){
        obs_seg.occ[i,j] = 1
      }
    }
  }
  obs_seg.occ[obs_seg.occ == 0] = NA
  
  # same logic as with Obs1c for the section below
  det.ind = 1*(!is.na(DTD2))
  DTD2[is.na(DTD2)] = R
  
  Obs2cconstants <- list(num_sites = num_sites, num_segments = num_segments, R = R,
                         stationary = T, num_trans = num_trans, zone = zone, l.shape = 1,
                         l.scale = 1/R)
  Obs2cdata = list(det = DTD2, det.ind = det.ind, occ = obs_occ, seg.occ = obs_seg.occ)
  
  occ.init = data$site.occ
  occ.init[!is.na(obs_occ)] = NA
  
  seg.occ.init = data$Occ
  for (i in 1:(num_sites*num_trans)){
    if (!is.na(occ.init[zone[i]])){
      if (occ.init[zone[i]] == 0){
        seg.occ.init[i,] = 0
      }
    }
  }
  seg.occ.init[!is.na(obs_seg.occ)] = NA
  
  Obs2cinits <- function(){
    list(theta11 = theta11, theta01 = theta01, lambda = lambda,
         seg.occ = seg.occ.init, psi = psi, occ = occ.init)
  }
  
  return(list(data = Obs2cdata,
              constants = Obs2cconstants,
              inits = Obs2cinits))
  
}
  
set_model4d <- function(data, SimParams){
  
  list2env(SimParams,environment())
  
  zone = rep(1:num_sites, each = num_trans)
  
  C1 = data$C1
  
  # what sites do we know are occupied
  obs_occ_trans = 1*(rowSums(C1) >= 1)
  obs_occ <- rep(NA, num_sites)
  for (i in 1:num_sites){
    obs_occ[i] <- 1*(sum(obs_occ_trans[which(zone == i)]) > 0)
  }
  obs_occ[obs_occ == 0] = NA
  
  # what segments do we know are occupied
  obs_seg.occ = C1
  obs_seg.occ = 1*(obs_seg.occ >= 1)
  obs_seg.occ[obs_seg.occ == 0] = NA
  
  Obs4dconstants <- list(num_sites = num_sites, num_segments = num_segments, R = R,
                         stationary = T, num_trans = num_trans, zone = zone, l.shape = 1, 
                         l.scale = 1/R)
  Obs4ddata = list(det = C1, occ = obs_occ, seg.occ = obs_seg.occ)
  
  occ.init = data$site.occ
  occ.init[!is.na(obs_occ)] = NA
  
  seg.occ.init = data$Occ
  for (i in 1:(num_sites*num_trans)){
    if (!is.na(occ.init[zone[i]])){
      if (occ.init[zone[i]] == 0){
        seg.occ.init[i,] = 0
      }
    }
  }
  seg.occ.init[!is.na(obs_seg.occ)] = NA
  
  Obs4dinits <- function(){
    list(theta11 = theta11, theta01 = theta01, 
         seg.occ = seg.occ.init, psi = psi, 
         occ = occ.init, lambda = lambda)
  }
  
  return(list(data = Obs4ddata,
              constants = Obs4dconstants,
              inits = Obs4dinits))
}


#################### Run Simulation Test #####################

# before running this, need to load in the correct Experiment 
# cloning data, also need to load the correct Parameters file

run_test = function(Sim_name = 'Simulation', 
                    nchains = 1, niter = 210000, nthin = 100,
                    nburn = 10000, N=40,
                    save.every = 10, add_AF_slice = TRUE){
  # Sim_name (str): the file name the simulation results will be saved under
  # nchains (int): number of chains to run
  # niters (int): number of iterations per chain
  # nburn (int): number of burn in iterations
  # nthin (int): number for thinning of chain
  # N (int): the number of data sets to simulate
  # save.every (int): save all chains after save.every simulations
  
  #######################################
  #######################################
  # Generate N data sets
  #######################################
  #######################################
  
  data_list <- list()
  for (i in 1:N){
    data <- simulateData(SimParams)
    data_list[[length(data_list)+1]] <- data
  }
  
  ########################################
  ########################################
  # Run analysis on each set
  ########################################
  ########################################
  
  # Initial data set for first in N
  data <- data_list[[1]]
  
  # Set up data and constants for each model
  model1d <- set_model1d(data, SimParams)
  model1c <- set_model1c(data, SimParams)
  model2d <- set_model2d(data, SimParams)
  model2c <- set_model2c(data, SimParams)
  model4d <- set_model4d(data, SimParams)

  DEconstants1d = model1d$constants
  DEconstants1c = model1c$constants
  DEconstants2d = model2d$constants
  DEconstants2c = model2c$constants
  DEconstants4d = model4d$constants

  DEdata1d = model1d$data
  DEdata1c = model1c$data
  DEdata2d = model2d$data
  DEdata2c = model2c$data
  DEdata4d = model4d$data

  DEinits1d = model1d$inits
  DEinits1c = model1c$inits
  DEinits2d = model2d$inits
  DEinits2c = model2c$inits
  DEinits4d = model4d$inits

  DEmodel1d = nimbleModel(Obs1d, 
                          constants = DEconstants1d, 
                          data = DEdata1d, 
                          inits = DEinits1d())
  
  DEmodel1c = nimbleModel(Obs1c, 
                          constants = DEconstants1c, 
                          data = DEdata1c, 
                          inits = DEinits1c())
  
  DEmodel2d = nimbleModel(Obs2d, 
                          constants = DEconstants2d, 
                          data = DEdata2d, 
                          inits = DEinits2d())
  
  DEmodel2c = nimbleModel(Obs2c, 
                          constants = DEconstants2c, 
                          data = DEdata2c, 
                          inits = DEinits2c())
  
  DEmodel4d = nimbleModel(Obs4d, 
                          constants = DEconstants4d, 
                          data = DEdata4d, 
                          inits = DEinits4d())
  
  for (i in 1:N){
    print(i)
    
    if (i != 1){ # re-set data and inits
      
      data <- data_list[[i]]
      
      DEdata1d = set_model1d(data, SimParams)$data
      DEmodel1d$resetData()
      DEmodel1d$setData(DEdata1d)
      
      DEdata1c = set_model1c(data, SimParams)$data
      DEmodel1c$resetData()
      DEmodel1c$setData(DEdata1c)
      
      DEdata2d = set_model2d(data, SimParams)$data
      DEmodel2d$resetData()
      DEmodel2d$setData(DEdata2d)
      
      DEdata2c = set_model2c(data, SimParams)$data
      DEmodel2c$resetData()
      DEmodel2c$setData(DEdata2c)
      
      DEdata4d = set_model4d(data, SimParams)$data
      DEmodel4d$resetData()
      DEmodel4d$setData(DEdata4d)
      
      DEinits1d = set_model1d(data, SimParams)$inits
      DEinits1c = set_model1c(data, SimParams)$inits
      DEinits2d = set_model2d(data, SimParams)$inits
      DEinits2c = set_model2c(data, SimParams)$inits
      DEinits4d = set_model4d(data, SimParams)$inits
      
      DEmodel1d$setInits(DEinits1d())
      DEmodel1c$setInits(DEinits1c())
      DEmodel2d$setInits(DEinits2d())
      DEmodel2c$setInits(DEinits2c())
      DEmodel4d$setInits(DEinits4d())
      
    }
    
    # Run the model
    mcmcConf1d = configureMCMC(DEmodel1d, monitors = c('p', 'theta01', 'theta11', 'psi'), print = FALSE, enableWAIC = TRUE)
    mcmcConf1c = configureMCMC(DEmodel1c, monitors = c('p', 'theta01', 'theta11', 'psi'), print = FALSE, enableWAIC = TRUE)
    mcmcConf2d = configureMCMC(DEmodel2d, monitors = c('p', 'theta01', 'theta11', 'psi'), print = FALSE, enableWAIC = TRUE)
    mcmcConf2c = configureMCMC(DEmodel2c, monitors = c('p', 'theta01', 'theta11', 'psi'), print = FALSE, enableWAIC = TRUE)
    mcmcConf4d = configureMCMC(DEmodel4d, monitors = c('p', 'theta01', 'theta11', 'psi'), print = FALSE, enableWAIC = TRUE)
    
    if (add_AF_slice){
      mcmcConf1d$removeSamplers(c('p', 'theta11', 'theta01'))
      mcmcConf1d$addSampler(target = c('p','theta11', 'theta01'), type = 'AF_slice')
    }
    
    DEmcmc1d <- buildMCMC(mcmcConf1d, print = FALSE) 
    DEmcmc1c <- buildMCMC(mcmcConf1c, print = FALSE) 
    DEmcmc2d <- buildMCMC(mcmcConf2d, print = FALSE) 
    DEmcmc2c <- buildMCMC(mcmcConf2c, print = FALSE) 
    DEmcmc4d <- buildMCMC(mcmcConf4d, print = FALSE) 
    
    cDEmodel1d <- compileNimble(DEmodel1d)
    cDEmodel1c <- compileNimble(DEmodel1c)
    cDEmodel2d <- compileNimble(DEmodel2d)
    cDEmodel2c <- compileNimble(DEmodel2c)
    cDEmodel4d <- compileNimble(DEmodel4d)
    
    cDEmcmc1d <- compileNimble(DEmcmc1d, project = cDEmodel1d)
    cDEmcmc1c <- compileNimble(DEmcmc1c, project = cDEmodel1c)
    cDEmcmc2d <- compileNimble(DEmcmc2d, project = cDEmodel2d)
    cDEmcmc2c <- compileNimble(DEmcmc2c, project = cDEmodel2c)
    cDEmcmc4d <- compileNimble(DEmcmc4d, project = cDEmodel4d)
    
    DEresults1d <- runMCMC(cDEmcmc1d, niter = niter, nburnin = nburn, 
                           thin = nthin, nchains = nchains, WAIC = TRUE)
    DEresults1c <- runMCMC(cDEmcmc1c, niter = niter, nburnin = nburn, 
                           thin = nthin, nchains = nchains, WAIC = TRUE)
    DEresults2d <- runMCMC(cDEmcmc2d, niter = niter, nburnin = nburn, 
                           thin = nthin, nchains = nchains, WAIC = TRUE)
    DEresults2c <- runMCMC(cDEmcmc2c, niter = niter, nburnin = nburn, 
                           thin = nthin, nchains = nchains, WAIC = TRUE)
    DEresults4d <- runMCMC(cDEmcmc4d, niter = niter, nburnin = nburn, 
                           thin = nthin, nchains = nchains, WAIC = TRUE)
    
    # store results
    if (i == 1){
      
      par_chain1d = list()
      par_chain1c = list()
      par_chain2d = list()
      par_chain2c = list()
      par_chain4d = list()
      
      par_names1d = colnames(DEresults1d$samples)
      par_names1c = colnames(DEresults1c$samples)
      par_names2d = colnames(DEresults2d$samples)
      par_names2c = colnames(DEresults2c$samples)
      par_names4d = colnames(DEresults4d$samples)
      
      nsaved.iters <- nrow(DEresults1d$samples)
      
      for (par in par_names1d){
        par_chain1d[[par]] = matrix(numeric(0), nrow = nsaved.iters, ncol = N)
      }
      for (par in par_names1c){
        par_chain1c[[par]] = matrix(numeric(0), nrow = nsaved.iters, ncol = N)
      }
      for (par in par_names2d){
        par_chain2d[[par]] = matrix(numeric(0), nrow = nsaved.iters, ncol = N)
      }
      for (par in par_names2c){
        par_chain2c[[par]] = matrix(numeric(0), nrow = nsaved.iters, ncol = N)
      }
      for (par in par_names4d){
        par_chain4d[[par]] = matrix(numeric(0), nrow = nsaved.iters, ncol = N)
      }
    }
    
    for (par in par_names1d){
      par_chain1d[[par]][,i] = DEresults1d$samples[,par]
    }
    for (par in par_names1c){
      par_chain1c[[par]][,i] = DEresults1c$samples[,par]
    }
    for (par in par_names2d){
      par_chain2d[[par]][,i] = DEresults2d$samples[,par]
    }
    for (par in par_names2c){
      par_chain2c[[par]][,i] = DEresults2c$samples[,par]
    }
    for (par in par_names4d){
      par_chain4d[[par]][,i] = DEresults4d$samples[,par]
    }
    
    if (i%%save.every == 0){
      filename <- paste0(Sim_name, '_', as.character(i%/%save.every), '.RData')
      save(par_chain1d, par_chain1c, par_chain2d, par_chain2c, par_chain4d,
           data_list, file = filename)
    } # iteratively save the data (just in case something goes wrong half way through)
    
  }
}