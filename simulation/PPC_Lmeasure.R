# Compared SMMB with fitting a separate model for each user pair
rm(list=ls())
library(survival)
library(Matrix)

proj <- "simulation"
v <- 1
set.seed(1234)


if(!dir.exists("PPC_results")){
  dir.create("PPC_results")
}

# load observed data
load(paste0("ObsData/Obs_SMMB_v",v,".RData"))

# Set dimension information
## the number of objects
N <- max(Dat[,1], Dat[,2])
tol_rep <- nrow(Dat)
tol_pair <- nrow(num.pair)
p.cov <- ncol(Dat) - 6

# extract surivival observations and covariate matrix
nu <- Dat[,4]
X <- Dat[,5:(5+p.cov)]
Y <- Dat[,ncol(Dat)]

## Number of time knots
P <- 5
knots<-c(6,18,30,48,90)
diff.knots<-diff(c(0,knots))
knots0 <- c(0,knots)

# iteration number
num_iter <- 50000
burnin <- 25000
thin <- 1
num_replications <- (num_iter - burnin)/thin

# auxiliary functions
rCure <- function(num, lam, kn, the){
  
  # "the" variable can be a scalar, all survival time sampling with the same theta,
  # or "the" variable is a vector with length "num" 
  if(length(the) == 1){
    vec_the <- rep(the, num)
  }else{
    vec_the <- the
  }
  result <- rep(NA, num)
  kn0 <- c(0,kn)
  sur_ratio <- exp(-lam * diff(kn0))
  num_int <- length(lam)
  
  # calculate the survival probabilities of PCH distributions at each knot
  sur_p <- rep(1,num_int + 1)
  for(p in 2:(num_int +1)){
    sur_p[p:(num_int+1)] <- sur_p[p:(num_int+1)] * sur_ratio[p-1]
  }
  
  u <- runif(num)
  
  for(i in 1:num){
    # temp = 1 + log(u) / theta
    temp <- 1 + log(u[i]) / vec_the[i]
    
    p_i <- sum(temp < sur_p)
    if(p_i == num_int + 1){
      result[i] <- Inf
    }else{
      result[i] <- - log(temp / sur_p[p_i]) / lam[p_i] + kn0[p_i] 
    }
  }
  return(result)
}

rCure_truncated <- function(num, lam, kn, the, thres){
  
  # "the" variable can be a scalar, all survival time sampling with the same theta,
  # or "the" variable is a vector with length "num" 
  if(length(the) == 1){
    vec_the <- rep(the, num)
  }else{
    vec_the <- the
  }
  result <- rep(NA, num)
  kn0 <- c(0,kn)
  sur_ratio <- exp(-lam * diff(kn0))
  num_int <- length(lam)
  
  # calculate the survival probabilities of PCH distributions at each knot
  sur_p <- rep(1,num_int + 1)
  for(p in 2:(num_int +1)){
    sur_p[p:(num_int+1)] <- sur_p[p:(num_int+1)] * sur_ratio[p-1]
  }
  
  u <- runif(num, min = 0, max = 1 - thres)
  
  for(i in 1:num){
    # temp = 1 + log(u) / theta
    temp <- 1 + log(u[i]) / vec_the[i]
    
    p_i <- sum(temp < sur_p)
    if(p_i == num_int + 1){
      result[i] <- Inf
    }else{
      result[i] <- - log(temp / sur_p[p_i]) / lam[p_i] + kn0[p_i] 
    }
  }
  return(result)
}

# calculate the survival probability of PCH distribution at each knot
qCure <- function(t, lam, kn, the){
  #lam: the constant hazard in each interval, 
  #kn: knots
  #the: exp(x^T beta)
  
  # determine the interval
  
  P <- length(kn)
  p <- min(sum(t > kn) + 1, P)
  
  if(p == 1){
    sur_base <- exp(-lam[1] * t)
  }else{
    dif.kn <- diff(c(0,kn))
    sur_base <- exp(-sum(c(dif.kn[1:(p-1)],t-kn[p-1]) * lam[1:p]))
  }
  
  quan <- 1 - exp(-the * (1-sur_base))
  
  return(quan)
}

###################################
# Calculate L measure of the SMMB #
###################################
load(paste0("Inference/Inference_",proj,"_SMMB_K3_v",v,"_r1.RData"))

# divide all user pairs into 10 parts to run in parallel
large_value <- exp(10)

# recording
sum_logyrep_SMMB <- rep(0, tol_rep)
sum_logyrepsq_SMMB <- rep(0, tol_rep)
sum_logyobssq_SMMB <- rep(0, tol_rep)
mean_logyrep_SMMB <- rep(0, tol_rep)
mean_logyobssq_SMMB <- rep(0,tol_rep)

time_start_SMMB <- Sys.time()
index_rep <- 0

for(i in 1:tol_pair){
  
  nij <- num.pair[i,3]
  
  for(iter in 1:num_replications){
    
    # generate the replications of real dataset with the same censoring pattern as the simulated data
    # parameter and latent variable values at the current iteration
    lambda_cur <- lambda_record[[iter * thin + burnin]]
    beta_cur <- beta_record[[iter * thin + burnin]]
    role_pair_cur <- mp_record[[iter * thin + burnin]] + 1
    
    # Inverse sampling
    y_rep <- rep(NA, nij)
    mp_cur <- role_pair_cur[i]
    
    for(g in 1:nij){
      
      if(nu[index_rep + g] == 1){
        y_rep[g] <- rCure(1, lambda_cur[mp_cur,], knots, exp(Dat[index_rep + g, 5:(p.cov + 5)] %*% beta_cur[mp_cur, ]))
      }else{
        thres_temp <- qCure(Y[index_rep + g], lambda_cur[mp_cur,], knots, exp(Dat[index_rep + g, 5:(p.cov + 5)] %*% beta_cur[mp_cur, ]))
        y_rep[g] <- rCure_truncated(1, lambda_cur[mp_cur,], knots, exp(Dat[index_rep + g, 5:(p.cov + 5)] %*% beta_cur[mp_cur, ]), thres_temp)
      }
    }
    
    y_rep[which(y_rep > large_value)] <- large_value 
    sum_logyrep_SMMB[index_rep + 1:nij] <- sum_logyrep_SMMB[index_rep + 1:nij] + log(y_rep)
    sum_logyrepsq_SMMB[index_rep + 1:nij] <- sum_logyrepsq_SMMB[index_rep + 1:nij] + log(y_rep)^2
    
  }
  
  mean_logyrep_SMMB[index_rep + 1:nij] <- sum_logyrep_SMMB[index_rep + 1:nij]/num_replications
  
  
  for(iter in 1:num_replications){
    
    # generate the replications of real dataset with the same censoring pattern as the simulated data
    # parameter and latent variable values at the current iteration
    lambda_cur <- lambda_record[[iter + burnin]]
    beta_cur <- beta_record[[iter + burnin]]
    role_pair_cur <- mp_record[[iter * thin + burnin]] + 1
    
    # Inverse sampling
    y_cen <- rep(NA, nij)
    index_cen <- which(nu[index_rep + 1:nij] == 0)
    mp_cur <- role_pair_cur[i]
    
    for(g in index_cen){
      
      thres_temp <- qCure(Y[index_rep + g], lambda_cur[mp_cur,], knots, exp(Dat[index_rep + g, 5:(p.cov + 5)] %*% beta_cur[mp_cur, ]))
      y_cen[g] <- rCure_truncated(1, lambda_cur[mp_cur,], knots, exp(Dat[index_rep + g, 5:(p.cov + 5)] %*% beta_cur[mp_cur, ]), thres_temp)
      
    }
    
    y_cen[which(y_cen > large_value)] <- large_value 
    sum_logyobssq_SMMB[index_rep + index_cen] <- sum_logyobssq_SMMB[index_rep + index_cen] + (log(y_cen[index_cen]) - mean_logyrep_SMMB[index_rep + index_cen])^2
  }
  
  mean_logyobssq_SMMB[index_rep + index_cen] <- sum_logyobssq_SMMB[index_rep + index_cen]/num_replications
  
  
  for(g in which(nu[index_rep + 1:nij] == 1)){
    
    mean_logyobssq_SMMB[index_rep + g] <- mean_logyobssq_SMMB[index_rep + g] + (log(Y[index_rep + g]) - mean_logyrep_SMMB[index_rep + g])^2
    
  }
  
  index_rep <- index_rep + nij
  
  time_cur_SMMB <- Sys.time()
  message(paste0("It has taken ",difftime(time_cur_SMMB, time_start_SMMB, units = "mins")," minutes to finish the ",i,"-th user pair"))
}

###########################################
# Calculate L measure of the overall SCRM #
###########################################
load(paste0("Inference/Inference_",proj,"_SMMB_K1_v",v,"_r1.RData"))

# divide all user pairs into 10 parts to run in parallel
large_value <- exp(10)

# recording
sum_logyrep_overall <- rep(0, tol_rep)
sum_logyrepsq_overall <- rep(0, tol_rep)
sum_logyobssq_overall <- rep(0, tol_rep)
mean_logyrep_overall <- rep(0, tol_rep)
mean_logyobssq_overall <- rep(0,tol_rep)

time_start_overall <- Sys.time()

index_rep <- 0
for(i in 1:tol_pair){
  nij <- num.pair[i,3]
  
  for(iter in 1:num_replications){
    
    # generate the replications of real dataset with the same censoring pattern as the simulated data
    # parameter and latent variable values at the current iteration
    lambda_cur <- lambda_record[[iter * thin + burnin]]
    beta_cur <- beta_record[[iter * thin + burnin]]
    
    # Inverse sampling
    y_rep <- rep(NA, nij)
    for(g in 1:nij){
      
      if(nu[index_rep + g] == 1){
        y_rep[g] <- rCure(1, lambda_cur, knots, exp(Dat[index_rep + g, 5:(p.cov + 5)] %*% beta_cur))
      }else{
        thres_temp <- qCure(Y[index_rep + g], lambda_cur,knots, exp(Dat[index_rep + g, 5:(p.cov + 5)] %*% beta_cur))
        y_rep[g] <- rCure_truncated(1, lambda_cur, knots, exp(Dat[index_rep + g, 5:(p.cov + 5)] %*% beta_cur), thres_temp)
      }
    }
    
    y_rep[which(y_rep > large_value)] <- large_value 
    sum_logyrep_overall[index_rep + 1:nij] <- sum_logyrep_overall[index_rep + 1:nij] + log(y_rep)
    sum_logyrepsq_overall[index_rep + 1:nij] <- sum_logyrepsq_overall[index_rep + 1:nij] + log(y_rep)^2
    
  }
  
  mean_logyrep_overall[index_rep + 1:nij] <- sum_logyrep_overall[index_rep + 1:nij]/num_replications
  
  
  for(iter in 1:num_replications){
    
    # generate the replications of real dataset with the same censoring pattern as the simulated data
    # parameter and latent variable values at the current iteration
    lambda_cur <- lambda_record[[iter + burnin]]
    beta_cur <- beta_record[[iter + burnin]]
    
    # Inverse sampling
    y_cen <- rep(NA, nij)
    index_cen <- which(nu[index_rep + 1:nij] == 0)
    for(g in index_cen){
      
      thres_temp <- qCure(Y[index_rep + g], lambda_cur,knots, exp(Dat[index_rep + g, 5:(p.cov + 5)] %*% beta_cur))
      y_cen[g] <- rCure_truncated(1, lambda_cur, knots, exp(Dat[index_rep + g, 5:(p.cov + 5)] %*% beta_cur), thres_temp)
    }
    
    y_cen[which(y_cen > large_value)] <- large_value 
    sum_logyobssq_overall[index_rep + index_cen] <- sum_logyobssq_overall[index_rep + index_cen] + (log(y_cen[index_cen]) - mean_logyrep_overall[index_rep + index_cen])^2
  }
  
  mean_logyobssq_overall[index_rep + index_cen] <- sum_logyobssq_overall[index_rep + index_cen]/num_replications
  
  for(g in which(nu[index_rep + 1:nij] == 1)){
    
    mean_logyobssq_overall[index_rep + g] <- mean_logyobssq_overall[index_rep + g] + (log(Y[index_rep + g]) - mean_logyrep_overall[index_rep + g])^2
    
  }
  
  index_rep <- index_rep + nij
  
  time_cur_overall <- Sys.time()
  message(paste0("It has taken ",difftime(time_cur_overall, time_start_overall, units = "mins")," minutes to finish the ",i,"-th user pair"))
}

############################################
# Calculate L measure of the pairwise SCRM #
############################################
# recording
sum_logyrep_Pairwise <- rep(0, tol_rep)
sum_logyrepsq_Pairwise <- rep(0, tol_rep)
sum_logyobssq_Pairwise <- rep(0, tol_rep)
mean_logyrep_Pairwise <- rep(0, tol_rep)
mean_logyobssq_Pairwise <- rep(0,tol_rep)

time_start_Pairwise <- Sys.time()
index_rep <- 0

for(i in 1:tol_pair){
  
  load(paste0("Inference/Inference_simulation_pairwiseSCRM_pair",i,".RData"))
  nij <- num.pair[i,3]
  
  for(iter in 1:num_replications){
    
    # generate the replications of real dataset with the same censoring pattern as the simulated data
    # parameter and latent variable values at the current iteration
    lambda_cur <- lambda_record[iter,]
    beta_cur <- beta_record[iter,]
    
    # Inverse sampling
    y_rep <- rep(NA, nij)
    
    for(g in 1:nij){
      
      if(nu[index_rep + g] == 1){
        y_rep[g] <- rCure(1, lambda_cur, knots, exp(Dat[index_rep + g, 5:(p.cov + 5)] %*% beta_cur))
      }else{
        thres_temp <- qCure(Y[index_rep + g], lambda_cur,knots, exp(Dat[i, 5:(p.cov + 5)] %*% beta_cur))
        y_rep[g] <- rCure_truncated(1, lambda_cur, knots, exp(Dat[i, 5:(p.cov + 5)] %*% beta_cur), thres_temp)
      }
    }
    
    y_rep[which(y_rep > large_value)] <- large_value
    sum_logyrep_Pairwise[index_rep + 1:nij] <- sum_logyrep_Pairwise[index_rep + 1:nij] + log(y_rep)
    sum_logyrepsq_Pairwise[index_rep + 1:nij] <- sum_logyrepsq_Pairwise[index_rep + 1:nij] + log(y_rep)^2
  }
  
  mean_logyrep_Pairwise[index_rep + 1:nij] <- sum_logyrep_Pairwise[index_rep + 1:nij]/num_replications
  
  
  for(iter in 1:num_replications){
    
    # generate the replications of real dataset with the same censoring pattern as the simulated data
    # parameter and latent variable values at the current iteration
    lambda_cur <- lambda_record[iter,]
    beta_cur <- beta_record[iter,]
    
    # Inverse sampling
    y_cen <- rep(NA, nij)
    index_cen <- which(nu[index_rep + 1:nij] == 0)
    
    for(g in index_cen){
      
      thres_temp <- qCure(Y[index_rep + g], lambda_cur,knots, exp(Dat[i, 5:(p.cov + 5)] %*% beta_cur))
      y_cen[g] <- rCure_truncated(1, lambda_cur, knots, exp(Dat[i, 5:(p.cov + 5)] %*% beta_cur), thres_temp)
      
    }
    
    y_cen[which(y_cen > large_value)] <- large_value 
    sum_logyobssq_Pairwise[index_rep + index_cen] <- sum_logyobssq_Pairwise[index_rep + index_cen] + (log(y_cen[index_cen]) - mean_logyrep_Pairwise[index_rep + index_cen])^2
  }
  
  mean_logyobssq_Pairwise[index_rep + index_cen] <- sum_logyobssq_Pairwise[index_rep + index_cen]/num_replications
  
  
  for(g in which(nu[index_rep + 1:nij] == 1)){
    
    mean_logyobssq_Pairwise[index_rep + g] <- mean_logyobssq_Pairwise[index_rep + g] + (log(Y[index_rep + g]) - mean_logyrep_Pairwise[index_rep + g])^2
    
  }
  
  index_rep <- index_rep + nij
  
  time_cur_Pairwise <- Sys.time()
  message(paste0("It has taken ",difftime(time_cur_Pairwise, time_start_Pairwise, units = "mins")," minutes to finish the ",i,"-th user pair"))
}


#######################
# Summarize L measure #
#######################
Lmeasure_Collect <- matrix(NA, 3, 4)
rownames(Lmeasure_Collect) <- c("SMMB", "Overall SCRM", "Pairwise SCRM")
colnames(Lmeasure_Collect) <- c("Variance","Squared Bias","L measure delta = 1/2","L measure delta = 1/3")

#SMMB
Lmeasure_Collect[1,1] <- sum(sum_logyrepsq_SMMB/num_replications - mean_logyrep_SMMB^2)
Lmeasure_Collect[1,2] <- sum(mean_logyobssq_SMMB)
Lmeasure_Collect[1,3] <- Lmeasure_Collect[1,1] + Lmeasure_Collect[1,2]/2
Lmeasure_Collect[1,4] <- Lmeasure_Collect[1,1] + Lmeasure_Collect[1,2]/3

#OVerall
Lmeasure_Collect[2,1] <- sum(sum_logyrepsq_overall/num_replications - mean_logyrep_overall^2)
Lmeasure_Collect[2,2] <- sum(mean_logyobssq_overall)
Lmeasure_Collect[2,3] <- Lmeasure_Collect[2,1] + Lmeasure_Collect[2,2]/2
Lmeasure_Collect[2,4] <- Lmeasure_Collect[2,1] + Lmeasure_Collect[2,2]/3

#Pairwise
Lmeasure_Collect[3,1] <- sum(sum_logyrepsq_Pairwise/num_replications - mean_logyrep_Pairwise^2)
Lmeasure_Collect[3,2] <- sum(mean_logyobssq_Pairwise)
Lmeasure_Collect[3,3] <- Lmeasure_Collect[3,1] + Lmeasure_Collect[3,2]/2
Lmeasure_Collect[3,4] <- Lmeasure_Collect[3,1] + Lmeasure_Collect[3,2]/3

