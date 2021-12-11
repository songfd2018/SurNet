rm(list=ls())
library(survival)
library(Matrix)
set.seed(1234)

proj <-"simulation"
v <- 1
P <- 5

num_iter <- 50000
burnin <- 25000
thinning <- 1
pri_beta <- 5

# calculate the survival probability of PCH distribution at each knot
cal_sur_PCH <- function(lam,len_int){
  #lam the constant hazard in each interval, 
  #len_int the length of each interval
  num<-length(len_int)
  result <- rep(1,num)
  for(p in 1:num){
    result[p:num] <- result[p:num] * exp(-lam[p] * len_int[p])
  }
  return(result)
}

# calculate the survival probability of cure rate model at each knot
cal_sur_cure <- function(lam, len_int, the){
  sur_pch <- cal_sur_PCH(lam,len_int)
  sur_cure <- exp(-the * (1- sur_pch))
  return(sur_cure)
}

# sample from dirichlet distribution
rdiri <- function(par){
    vec_len <- length(par)
    res <- rgamma(vec_len, shape = par)
    res <- res / sum(res)
    return(res)
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
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
  
#####################################
# MCMC algorithm for pairwise SCRMs #
#####################################
# Auxillary variables
## Generate auxillary value p_ij indicates the interval Y_ij belongs to
Ind <- rep(0, tol_rep)
for(i in 1:tol_rep){
  Ind[i] <- sum(Y[i] > knots0[1:P])
}

## Construct two auxillary matrix for accumulating nu and Y - s[p-1] or s[p] - s[p-1]
A <- matrix(NA,tol_pair,P)
Mat_Y <- list()
index_rep <- 0
for(i in 1:tol_pair){
  nij <- num.pair[i,3]
  temp <- matrix(0, nij, P)
  for(j in 1:nij){
    pij <- Ind[index_rep + j]
    temp[j, 1:(pij - 1)] <- diff.knots[1:(pij - 1)]
    temp[j, pij] <- Y[index_rep + j] - knots0[pij]
  }
  Mat_Y[[i]] <- as(temp, "dgCMatrix")
  A[i,] <- apply(Mat_Y[[i]],2,sum)
  index_rep <- index_rep + nij
}

##B_{ijkp} = I(v_ijk = 1, p_ijk = p)
B <- matrix(NA,tol_pair,P)
Mat_nu <- list()
index_rep <- 0
for(i in 1:tol_pair){
  nij <- num.pair[i,3]
  temp <- matrix(0, nij, P)
  for(j in 1:nij){
    if(nu[index_rep + j] == 1){
      pij <- Ind[index_rep + j]
      temp[j, pij] <- 1
    }
  }
  Mat_nu[[i]] <- as(temp, "dgCMatrix")
  B[i,] <- apply(Mat_nu[[i]],2,sum)
  index_rep <- index_rep + nij
}

# Calculate the empirical estimator as prior distribution
## regression coefficients beta
beta.init <- rep(0, p.cov + 1) 
fit_all <- survfit(Surv(Y,nu) ~ 1)
sur_all <- cbind(fit_all$time, fit_all$surv)
cure_rate_rough <- sur_all[nrow(sur_all),2]
beta.init[1] <- log(-log(cure_rate_rough))

## Sample latent variable N_ij from Poisson(S_cure(y_ij|lambda)) + nu_ij
N.init <- rep(NA, tol_rep)
index_rep <- 0 # index (i,j,1) - 1 of Y matrix 
for(i in 1:tol_pair){
  nij <-  num.pair[i, 3]
  
  temp <- rep(NA, nij)
  for(k in 1:nij){
    temp[k] <- sum(sur_all[,1] <= Y[index_rep + k])
  }
  
  sur_prob <- sur_all[temp,2]
  sur_prob[which(is.na(sur_prob))] <- 0
  
  N.init[index_rep + 1:nij] <- rpois(nij, log(sur_prob/cure_rate_rough)) + nu[index_rep + 1:nij]
  index_rep <- index_rep + nij
}

## set the MLE of lambda of complete likelihood as the intial values
num_lam <- 0
den_lam <- 0
index_rep <- 0
for(i in 1:tol_pair){
  nij <-  num.pair[i, 3]
  num_lam <- num_lam + B[i,]
  den_lam <- den_lam + as.vector(N.init[index_rep + 1:nij] %*% Mat_Y[[i]])
  index_rep <- index_rep + nij
}
omega <- num_lam / den_lam

lambda.iter <- matrix(NA, tol_pair, P)
index_rep <- 0
for(i in 1:tol_pair){
  
  nij <-  num.pair[i, 3]
  num_lam <- B[i,]
  den_lam <- as.vector(N.init[index_rep + 1:nij] %*% Mat_Y[[i]])
  lambda.iter[i,] <- num_lam / den_lam
  lambda.iter[i,is.nan(lambda.iter[i,])] <- 0
  index_rep <- index_rep + nij
}

beta.iter<- matrix(beta.init, tol_pair, p.cov + 1, byrow =TRUE)

# add a penalty on beta
beta.prior <- list()
beta.prior$mean <- beta.init
beta.prior$var <- pri_beta
eta <- 1/2/beta.prior$var

# record the parameter values every 100 iterations
num_record <- (num_iter - burnin) / thinning
lambda_record <- matrix(NA, num_record, P)
beta_record <- matrix(NA, num_record, p.cov + 1)

########################
# Apply MCMC algorithm #
########################
if(!dir.exists("Inference")){
  dir.create("Inference")
}

start_time <- Sys.time()
index_rep <- 0
for(i in 1:tol_pair){
  
  nij <-  num.pair[i, 3]
  nu_temp <- nu[index_rep + 1:nij]
  Y_temp <- Y[index_rep + 1:nij]
  X_temp <- X[index_rep + 1:nij,]
  Ind_temp <- Ind[index_rep + 1:nij]
  lambda_temp <- lambda.iter[i,]
  beta_temp <- beta.iter[i,]
  
  re <- 1
  for(iter in 1:num_iter){
    # sample U_{ijg}
    U_temp <- nu_temp + rpois(nij, lambda = as.vector(exp(-Mat_Y[[i]] %*%lambda_temp + X_temp %*% beta_temp)))
    
    # sample lambda_m
    lam_num <- B[i,]
    lam_dem <- as.vector(U_temp %*% Mat_Y[[i]])
    lam_num[which(lam_dem>0)] <- lam_num[which(lam_dem>0)] + 1
    lambda_temp <- rgamma(P, lam_num,lam_dem)
    
    # sample beta_r
    for(r in 1:(p.cov + 1)){
      beta_mh <- rnorm(1, beta_temp[r], sd = 0.1)
      
      log_rho <- sum(X_temp[,r] * U_temp) * (beta_mh - beta_temp[r])
      log_rho <- log_rho - sum(exp(X_temp %*% beta_temp) * (exp(X_temp[,r] * (beta_mh - beta_temp[r])) - 1))
      if(r > 1){
        log_rho <- log_rho - eta * (beta_mh^2 - beta_temp[r]^2)
      }
      
      if(log(runif(1)) < log_rho){
        beta_temp[r] <- beta_mh
      }
    }
    
    # calculate likelihood
    betax <- X_temp %*% beta_temp
    
    # record the sampling
    if(iter > burnin){
      lambda_record[iter - burnin, ] <- lambda_temp
      beta_record[iter - burnin, ] <- beta_temp
    }
    
    
  }
  
  save(lambda_record,beta_record,file = paste0("Inference/Inference_",proj,"_pairwiseSCRM_pair",i,".RData"))
  
  index_rep <- index_rep + nij
  
  
  if(i %% 20 == 0){
    cur_time <- Sys.time()
    message(paste0("Finish ", i, " edges."))
    message(paste0("It has taken ",difftime(cur_time, start_time, units = "mins"), " mins."))
  }
}

end_time <- Sys.time()
message(paste0("The total running time is ",difftime(end_time, start_time, units = "mins"), " mins."))

save.image(paste0("Inference/Inference_",proj,"_pairwiseSCRM_v",v,".RData"))