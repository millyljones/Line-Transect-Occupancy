# Compute metric for a given data set

Delta_metric <- function(J, data, sites){
  # J (integer): is the number of replicate segment level observations
  # data: a matrix (e.g. if J=1) or an array
    # If a matrix: nrow(data) = number of sites and ncol(data) = number of segments
    # If an array: dim(data)= number of sites, number of segments, number of replicate segment-level observations
    # data can be passed as an array when J=1 if dim(data)[3] = 1
    # Entries are 0,1 entries denoting detection/non-detection
  # sites: a vector of integers where i-th entry denotes in which site transect i is located
  
  if (length(dim(data))==2){
    data <- array(data, dim = c(dim(data)[1], dim(data)[2], 1))
  }
  
  # narrow data we are focusing on:
  data <- data[,,1:J, drop = F]
  num_sites <- max(sites)
  nsegs <- dim(data)[2]
  
  # get sites with at least one detection
  obs_occ_trans = 1*(rowSums(data, na.rm=T) >= 1)
  obs_occ <- rep(NA, num_sites)
  for (i in 1:num_sites){
    obs_occ[i] <- 1*(sum(obs_occ_trans[which(sites == i)]) > 0)
  }

  occ.sites <- which(obs_occ == 1)
  data.obs <- data[sites %in% occ.sites,,1:J, drop=F] # just get the occupied sites

  nsegs <- dim(data.obs)[2]

  # now compute the transition probs
  n10 <- 0 # number of 10's
  n11 <- 0 # number of 11's
  n00 <- 0 # number of 00's
  n01 <- 0 # number of 01's

  for (k in 1:dim(data.obs)[1]){ # how many times see pattern across sites
    for (j in 1:J){ # loop through the replicates
      n01 <- n01 + sum(paste0(as.character(data.obs[k,1:(nsegs-1),j]),
                                as.character(data.obs[k,2:nsegs,j])) == '01')
      n11 <- n11 + sum(paste0(as.character(data.obs[k,1:(nsegs-1),j]),
                                as.character(data.obs[k,2:nsegs,j])) == '11')
      n10 <- n10 + sum(paste0(as.character(data.obs[k,1:(nsegs-1),j]),
                                as.character(data.obs[k,2:nsegs,j])) == '10')
      n00 <- n00 + sum(paste0(as.character(data.obs[k,1:(nsegs-1),j]),
                                as.character(data.obs[k,2:nsegs,j])) == '00')
    }
  }
    
  p1.1 <- n11/(n11+n10) # prob of 1 given 1
  p1.0 <- n01/(n01+n00) # prob of 1 given 0
  diff <- p1.1 - p1.0 

  return(diff)
  
}

