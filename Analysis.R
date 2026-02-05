# Compute metric for each data set

Delta_metric <- function(data, J.max){
  # data is the data set
  # J.max is the number of replicates to consider in the computation
  # of the metric
  # computes the metric for this particular data set
  
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
    
  } else if (J.max == 2) {
    
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
      for (j in 1:2){ # loop through the replicates
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

p.diff<- function(p, theta01, theta11){
  theta11*p - (p*theta01*(1-theta11*p))/(1-theta11+theta01*(1-p))
}


#### Run through an analyse data sets

filenames <- paste0('P12/', c('P12_1_N1', 'P12_2_N1', 'P12_3_N1', 
                             'P12_4_N1', 'P12_5_N1'), '.RData')
df <- data.frame(matrix(nrow = 0, ncol = 5))
colnames(df) <- c('d1', 'd2', 'par', 'rate', 'method')
#models <- c('MP/DND/1', 'MP/DND/2', 'MP/DTD/1', 'MP/DTD/2', 'MP/C/1')
models <- c('MP/DND/2', 'MP/DTD/1', 'MP/DTD/2', 'MP/C/1')
D <- c(1,2,5,10,15,20)

for (i in 1:length(filenames)){
  
  load(filenames[i]) # load in data set  

  # get variances
  #DND1 <- as.data.frame(lapply(par_chain1d, function(x) {apply(x, 2, var)}))
  DTD1 <- as.data.frame(lapply(par_chain1c, function(x) {apply(x, 2, var)}))
  DND2 <- as.data.frame(lapply(par_chain2d, function(x) {apply(x, 2, var)}))
  DTD2 <- as.data.frame(lapply(par_chain2c, function(x) {apply(x, 2, var)}))
  C1 <- as.data.frame(lapply(par_chain4d, function(x) {apply(x, 2, var)}))
  
  # scale the variances
  #DND1 <- as.data.frame(matrix(unlist(apply(DND1, 1, function(x){ x/DND1[1,] })), ncol = length(par_chain1d), byrow = T))
  DTD1 <- as.data.frame(matrix(unlist(apply(DTD1, 1, function(x){ x/DTD1[1,] })), ncol = length(par_chain1c), byrow = T))
  DND2 <- as.data.frame(matrix(unlist(apply(DND2, 1, function(x){ x/DND2[1,] })), ncol = length(par_chain2d), byrow = T))
  DTD2 <- as.data.frame(matrix(unlist(apply(DTD2, 1, function(x){ x/DTD2[1,] })), ncol = length(par_chain2c), byrow = T))
  C1 <- as.data.frame(matrix(unlist(apply(C1, 1, function(x){ x/C1[1,] })), ncol = length(par_chain4d), byrow = T))
  
  # rate.dnd1 <- apply(DND1, 2, function(x){lm(log(x) ~ 0 + log(D))$coef[[1]]})
  # rate.dtd1 <- apply(DTD1, 2, function(x){lm(log(x) ~ 0 + log(D))$coef[[1]]})
  # rate.dnd2 <- apply(DND2, 2, function(x){lm(log(x) ~ 0 + log(D))$coef[[1]]})
  # rate.dtd2 <- apply(DTD2, 2, function(x){lm(log(x) ~ 0 + log(D))$coef[[1]]})
  # rate.c1 <- apply(C1, 2, function(x){lm(log(x) ~ 0 + log(D))$coef[[1]]})

  #mse.dnd1 <- apply(DND1, 2, function(x){mean((x-(1/D))^2)})
  mse.dtd1 <- apply(DTD1, 2, function(x){mean((x-(1/D))^2)})
  mse.dnd2 <- apply(DND2, 2, function(x){mean((x-(1/D))^2)})
  mse.dtd2 <- apply(DTD2, 2, function(x){mean((x-(1/D))^2)})
  mse.c1 <- apply(C1, 2, function(x){mean((x-(1/D))^2)})

  # last number in the filename is the data set from data_list
  set <- as.numeric(substr(filenames[i], nchar(filenames[i])-6, nchar(filenames[i])-6))
  
  df.i <- data.frame(matrix(nrow = 4*length(models), ncol = 5))
  colnames(df.i) <- c('d1', 'd2', 'par', 'rate', 'method')
  
  df.i$d1 <- Delta_metric(data_list[[set]], 1)
  df.i$d2 <- Delta_metric(data_list[[set]], 2)
  
  df.i$par <- rep(names(par_chain2d), length(models))
  # df.i$rate <- c(rate.dnd1, rate.dnd2, rate.dtd1,
  #               rate.dtd2, rate.c1)
  df.i$rate <- c(mse.dnd2, mse.dtd1, mse.dtd2, mse.c1) #c(mse.dnd1, mse.dnd2, mse.dtd1, mse.dtd2, mse.c1)
  df.i$method <- rep(models, each  = length(unique(df.i$par)))
  
  df <- rbind(df, df.i)

}

# for ignoring DND1 only
#df <- df[!(df$method %in% 'MP/DND/1'),]

save(df, file = 'df_P12.RData')


# Merge the results together to plot:

rm(list = ls())
filenames <- c(paste0('Non Parametric Results/', 'df_P', as.character(1:7), '.RData'),
               paste0('Non Parametric Results/', 'df_P', as.character(8), '*.RData'),
               paste0('Non Parametric Results/', 'df_P', as.character(9:11), '.RData'),
               paste0('Non Parametric Results/', 'df_P', as.character(12), '*.RData'))

df.plot <- data.frame(matrix(nrow = 0, ncol = 5))
colnames(df.plot) <- c('d1', 'd2', 'par', 'rate', 'method')

for (i in 1:length(filenames)){
  load(filenames[i])
  
  df.plot <- rbind(df.plot, df)
}

df.plot$method <- factor(df.plot$method, levels= c('MP/DND/1', 'MP/DND/2', 
                                                      'MP/DTD/1', 'MP/DTD/2', 'MP/C/1'))

save(df.plot, file = 'df_plot.RData')

library(ggplot2)

for (par in c('p', 'theta01', 'theta11')){
  df.par <- df.plot[df.plot$par == par,]
  p <- ggplot(df.par, aes(x = d1, y = rate, color = method)) + 
    geom_point() + #geom_hline(yintercept = -1) +  
    facet_wrap(~method, nrow = 1) + ggtitle(par)+
    theme_bw() + scale_y_continuous(trans='log10') +
    ylab('log(10) MSE') + 
    xlab('Metric') +
    scale_x_continuous(breaks=c(0.1,0.3, 0.5))
  ggsave(paste0(par, '.pdf'), width=6, height=3, units='in')
  print(p)
}
