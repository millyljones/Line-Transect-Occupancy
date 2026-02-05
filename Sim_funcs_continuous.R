run_occupancy = function(num_sites, num_trans, num_segments, psi, R, mu12, mu21, init){
  
  site.occ = rbinom(num_sites, 1, psi) # site occupancies
  
  switch_history = data.frame(matrix(NA, nrow = 0, ncol = 4))
  colnames(switch_history) = c('site', 'trans', 'locations', 'state')
  
  L = num_segments*R
  
  for (i in 1:num_sites){ # loop through sites
    for (m in 1:num_trans){
      if (site.occ[i] == 1){ # if site is occupied
      
      current_state = 1 + rbinom(1, 1, init) # state 2 is occupied state
      current_switch_rate = c(mu12, mu21)[current_state]
      current_location = 0
      
      # Run through switching process
      switch_locs = c(0)
      switch_states = c(current_state)
      
      while (current_location <= L){
        
        dist = rexp(1, current_switch_rate) # distance travelled in state
        current_location = current_location + dist
        current_state = c(1,2)[current_state != c(1:2)] # change state
        current_switch_rate = c(mu12, mu21)[current_state] # change switch rate
        
        if(current_location <= L){
          switch_locs = c(switch_locs, current_location)
          switch_states = c(switch_states, current_state)
        }
      }
      
      switch.i = data.frame('site' = i, 
                            'trans' = m,
                            'locations' = switch_locs, 
                            'state' = switch_states)
      switch_history = rbind(switch_history, switch.i)
      
    } else { # not occupied so no transect occupancy
      switch.i = data.frame('site' = i, 
                            'trans' = m,
                            'locations' = NA, 
                            'state' = NA)
      switch_history = rbind(switch_history, switch.i)
    }
    }# end loop through transects
  } # end loop through sites
  
  return(list('site.occ' = site.occ,
              'switch_history' = switch_history))
}

run_detections = function(num_sites, num_trans, num_segments, R, lambda, switch_history, site.occ){
  # places detections along the transects given the site occupancy and switching
  # processes along each transect.
  # This is the sort of data we would expect for 2MMPP analysis.
  # run_det4c then turns this into the form we would expect if we were doing method
  # 4c along segments.
  
  L = num_segments*R
  
  det_lambda = data.frame(matrix(NA, nrow = 0, ncol = 3))
  colnames(det_lambda) = c('site', 'trans', 'locations')
  
  for (i in 1:num_sites){ # looping through the new sites
    if(site.occ[i] == 1){ # if site is occupied
      for (m in 1:num_trans){
        switch.i = switch_history[switch_history$site == i & switch_history$trans == m, ]
        det_locs = c()
        n_sections = nrow(switch.i)
        for (n in 1:n_sections){ # looping through the sections
          if (switch.i$state[n] == 2){ # if in detection state
            # what are the bounds of this state?
            current_location = switch.i$locations[n] # start of section
            max_point = ifelse(is.na(switch.i$locations[n+1]),
                             L, switch.i$locations[n+1])
            # generate detections
            while (current_location <= max_point){
              
              dist = rexp(1, lambda)
              current_location = current_location + dist
              if (current_location <= max_point){
                det_locs = c(det_locs, current_location)
              }
            } # end detections for this occupied section
          } # end section with detections
        } # end detections for transect
        
        # if made no detections at transect (despite occupancy)
        if (length(det_locs) == 0){
          det_locs = NA
        }
        
        det.i = data.frame('site' = i, 'trans' = m, 'locations' = det_locs)
        det_lambda = rbind(det_lambda, det.i)
        
      } # end detections for all transects at site
        
    } else { # unoccupied sites
      for (m in 1:num_trans){
        det.i = data.frame('site' = i, 'trans' = m, 'locations' = NA)
        det_lambda = rbind(det_lambda, det.i)
      }
    }
    } # end looping through sites
    
    return(det_lambda)
}
        
# For run_det4c, just need det_lambda, cols are site, transect, and det_locs (i-th entry is distance
# between beginning on transect and detection)
# run_det4c assumes that the number of transects and segments are equal everywhere, but can easily be
# generalised.

run_det4c = function(num_sites, num_trans, num_segments, R, det_lambda){
  
  # turn det_lambda into det4c equivalent
  max_dets = max(table(det_lambda$trans)) # maximum number of detections per site
  det4c = array(data = NA, dim = c(num_sites*num_trans, num_segments, max_dets))
  for (i in 1:num_sites){
    for (m in 1:num_trans){
      # this sites detection locations
      det.i = det_lambda[det_lambda$site == i & det_lambda$trans == m,]
      for (j in 1:num_segments){
        # what are the bounds of this segment?
        seg_start = R*(j-1)
        seg_end = R*(j)
        # detections in this segment
        det.ij = det.i$locations[det.i$locations >= seg_start & det.i$locations <= seg_end]
        if (length(det.ij) > 0){ # if we have detections in this segment
          for (t in 1:length(det.ij)){ # loop through the detections
            if (t == 1){ # for first detection
              det4c[num_trans*(i-1) + m,j,t] = det.ij[t] - seg_start
            } else {
              det4c[num_trans*(i-1) + m,j,t] = det.ij[t] - det.ij[t-1]
            }
          }
        }
      } # end loop through segments
    } # end loop through trans
  } #end loop through sites
  
  max_dets = max(rowSums(!is.na(det4c), dim = 2, na.rm = T))
  det4c = det4c[,,1:max(1, max_dets), drop=F]
  
  return(det4c)
  
  # det4c has each [i,,] as a transect,
  # [,j,] as a segment,
  # [,,d] as an inter-detection distance (not including distance to end of segment)
  
}

simulate_data_continuous = function(num_sites, num_trans, num_segments, psi, init, 
                                    mu12, mu21, lambda, R){
  #num_sites
  #num_trans
  #num_segments - although not using segments in this continuous process
  #               I will use them with R to compute transect length
  #psi - site occupancy prob
  #init - prob of being in occupancy state
  #mu12 - transition rate state 1 to 2
  #mu21 - transition rate state 2 to 1
  #lambda - detection rate in occupancy state
  #lambda2 - detection rate of second process
  #R - length of segments
  
  # We assume that states switch according to a continuous time markov
  # process - and one state has no detections, whilst another does have
  # detections at exponential rate lambda.
  # This is the continuous analogue of the discrete segment occupancy methods
  # State 2 is the detecting state
  
  trans_occ = run_occupancy(num_sites, num_trans, num_segments, psi, R, mu12, mu21, init)
  switch_history = trans_occ$switch_history
  site.occ = trans_occ$site.occ
  
  det_lambda = run_detections(num_sites, num_trans, num_segments, R, lambda, switch_history, site.occ)
  
  # Turn this into the equivalent of the observation process
  # from the discrete methods
  # Obs det4c ----------------------------------------------------
  
  site_det_hist = run_det4c(num_sites, num_trans, num_segments, R, det_lambda)
  
  # Obs1c ---------------------------------------------------------
  # one continuous measurement
  det3 = site_det_hist[,,1] 
  
  # Obs1d----------------------------------------------------------
  # one discrete measurement
  det1 = 1*(!is.na(det3)) 
  
  # Obs2c----------------------------------------------------------
  # two continuous measurements (same rate)
  det5 = array(NA, dim = c(num_sites*num_trans, num_segments, 2))
  det5[,,1] = det3
  det2_lambda = run_detections(num_sites, num_trans, num_segments, R, lambda, switch_history, site.occ)
  site_det2_hist = run_det4c(num_sites, num_trans, num_segments, R, det2_lambda)
  det5[,,2] = site_det2_hist[,,1] 
  
  # Obs2d ---------------------------------------------------------
  # two discrete measurements (same rate)
  det2 = 1*(!is.na(det5))
  
  # Obs4d ----------------------------------------------------------
  # count per segment
  det4 = rowSums(1*(!is.na(site_det_hist)), dims = 2)
  
  return(list('site.occ' = site.occ,
              'DTD1' = det3,
              'DND1' = det1,
              'DTD2' = det5,
              'DND2' = det2,
              'det4c' = site_det_hist,
              'C1' = det4))
  
}
