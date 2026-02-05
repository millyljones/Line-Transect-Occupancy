# Compute metric for each data set

Delta_metric <- function(data, J.max){
  # J.max (integer): is the number of replicate segment level observations
  # data: If J.max = 1 (matrix) and an array otherwise.
    # If J.max = 1: nrow(data) = number of sites and ncol(data) = number of segments
    # If J.max > 1: dim(data)= number of sites, number of segments, number of replicate segment-level observations
    # Entries are 0,1 entries denoting detection/non-detection
  # Computes the Delta value for data

  # This needs generalising to data of form described above
  
  if (J.max == 1){
    
    num_zones <- length(data$site.occ)
    zones_to_trans <- rep(1:num_zones, each = nrow(data$Occ)/num_zones)
    
    obs_occ_trans = 1*(rowSums(data$DND1) >= 1)
    obs_occ <- rep(NA, num_zones)
    for (i in 1:num_zones){
      obs_occ[i] <- 1*(sum(obs_occ_trans[which(zones_to_trans == i)]) > 0)
    }
    
    occ.zones <- which(obs_occ == 1)
    DND1.obs <- data$DND1[zones_to_trans %in% occ.zones,] # just get the occupied zones
    dets.occ <- DND1.obs
    nsegs<- ncol(dets.occ)
    
    # now compute the transition probs
    n10 <- 0 # number of 10's
    n11 <- 0 # number of 11's
    n00 <- 0 # number of 00's
    n01 <- 0 # number of 01's
    
    for (k in 1:nrow(dets.occ)){ # how many times see pattern across sites
      n01 <- n01 + sum(paste0(as.character(dets.occ[k,1:(nsegs-1)]),
                              as.character(dets.occ[k,2:nsegs])) == '01')
      n11 <- n11 + sum(paste0(as.character(dets.occ[k,1:(nsegs-1)]),
                              as.character(dets.occ[k,2:nsegs])) == '11')
      n10 <- n10 + sum(paste0(as.character(dets.occ[k,1:(nsegs-1)]),
                              as.character(dets.occ[k,2:nsegs])) == '10')
      n00 <- n00 + sum(paste0(as.character(dets.occ[k,1:(nsegs-1)]),
                              as.character(dets.occ[k,2:nsegs])) == '00')
    }
    
    p1.1 <- n11/(n11+n10) # prob of 1 given 1
    p1.0 <- n01/(n01+n00) # prob of 1 given 0
    diff <- p1.1 - p1.0 
    
  } else if (J.max > 1) {
    
    num_zones <- length(data$site.occ)
    zones_to_trans <- rep(1:num_zones, each = nrow(data$Occ)/num_zones)
    
    obs_occ_trans = 1*(rowSums(data$DND2) >= 1)
    obs_occ <- rep(NA, num_zones)
    for (i in 1:num_zones){
      obs_occ[i] <- 1*(sum(obs_occ_trans[which(zones_to_trans == i)]) > 0)
    }
    
    occ.zones <- which(obs_occ == 1)
    DND2.obs <- data$DND2[zones_to_trans %in% occ.zones,,] # just get the occupied zones
    dets.occ <- DND2.obs
    nsegs<- ncol(dets.occ)
    
    # now compute the transition probs
    n10 <- 0 # number of 10's
    n11 <- 0 # number of 11's
    n00 <- 0 # number of 00's
    n01 <- 0 # number of 01's
    
    for (k in 1:nrow(dets.occ)){ # how many times see pattern across sites
      for (j in 1:J.max){ # loop through the replicates
        n01 <- n01 + sum(paste0(as.character(dets.occ[k,1:(nsegs-1),j]),
                                as.character(dets.occ[k,2:nsegs,j])) == '01')
        n11 <- n11 + sum(paste0(as.character(dets.occ[k,1:(nsegs-1),j]),
                                as.character(dets.occ[k,2:nsegs,j])) == '11')
        n10 <- n10 + sum(paste0(as.character(dets.occ[k,1:(nsegs-1),j]),
                                as.character(dets.occ[k,2:nsegs,j])) == '10')
        n00 <- n00 + sum(paste0(as.character(dets.occ[k,1:(nsegs-1),j]),
                                as.character(dets.occ[k,2:nsegs,j])) == '00')
      }
    }
    
    p1.1 <- n11/(n11+n10) # prob of 1 given 1
    p1.0 <- n01/(n01+n00) # prob of 1 given 0
    diff <- p1.1 - p1.0 
    
  }
  
  return(diff)
  
}

