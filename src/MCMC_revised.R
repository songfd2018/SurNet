rm(list=ls())
library(Matrix)
library(survival)
library(getopt)
library(RColorBrewer)
library(bayestestR)

spec = matrix(c(
  'project', 'p', 1, "character",
  'version'  , 'v', 1, "integer",
  'replication'   , 'r', 1, "integer",
  'num_membership'     , 'K', 1, "integer",
  'seed' , 's', 1, "integer",
  'num_knot', 'M', 2, "integer",
  'num_iteration', 'i', 2, "integer",
  'num_burnin', 'b', 2, "integer",
  'num_thin', 't', 2, "integer",
  'num_swap', 'w', 2, "integer",
  'num_core', 'c', 2, "integer",
  'prior_lambda', 'l', 2, "double",
  'prior_beta', 'm', 2, "double",
  'prior_xi', 'n', 2, "double",
  're_iter', 'o', 2, "integer",
  'DIC', 'd', 2, "integer",
  'help'   , 'h', 0, "logical"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

proj <- opt$project
v <- opt$version
r <- opt$replication
K <- opt$num_membership
set.seed(opt$seed)

# set the number of knots for real data analysis
if (is.null(opt$num_knot)) { 
  M <- 5
}else{
  M <- opt$num_knot
}

if (is.null(opt$num_iteration)) { 
  num_iter <- 1000
}else{
  num_iter <- opt$num_iteration
}

if (is.null(opt$num_thin)) { 
  num_thin <- 10
}else{
  num_thin <- opt$num_thin
}

if (is.null(opt$num_swap)){
  num_iter_swap <- 100
}else{
  num_iter_swap <- opt$num_swap
}

if (is.null(opt$num_burnin) ) { 
  burnin = num_iter/2    
}else{
  burnin = opt$num_burnin
}

if ( is.null(opt$num_core) ) { 
  nc = 1    
}else{
  nc = opt$num_core
}

if ( is.null(opt$prior_lambda) ) { 
  pri_lambda = 100    
}else{
  pri_lambda = opt$prior_lambda
}

if ( is.null(opt$prior_beta) ) { 
  pri_beta = 1    
}else{
  pri_beta = opt$prior_beta
}

if ( is.null(opt$prior_xi) ) { 
  pri_xi = 2    
}else{
  pri_xi = opt$prior_xi
}

if ( is.null(opt$re_iter) ) { 
  re_iter = 5    
}else{
  re_iter = opt$re_iter
}

if ( is.null(opt$DIC) ) { 
  compute_DIC <- FALSE    
}else{
  compute_DIC <- opt$DIC != 0
}

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

print(compute_DIC)

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

print("Load the observed data...")
if(proj == "sim"){
  
  ver.output<-paste0("K",K,"_",proj,"_v",v,"_r",r)
  
  # load observed data
  load(paste0("ObsData/Obs_SMMB_",proj,"_v",v,".RData"))
  
  # Set dimension information
  ## the number of objects
  N <- max(Dat[,1], Dat[,2])
  tol_rep <- nrow(Dat)
  tol_pair <- nrow(num.pair)
  p.cov <- ncol(Dat) - 5
  
  # extract surivival observations and covariate matrix
  nu <- Dat[,4]
  X <- Dat[,4+1:p.cov]
  Y <- Dat[,ncol(Dat)]
  
  ## Number of time knots
  knots<-c(6,24,60,90)
  knots0 <- c(0, knots)
  diff.knots<-diff(knots0)
  M <- length(knots)
  
}

if(proj == "simori"){
  
  ver.output<-paste0("K",K,"_",proj,"_v",v,"_r",r)
  
  # load observed data
  load(paste0("ObsData/Obs_SMMB_",proj,"_v",v,".RData"))
  
  # Set dimension information
  ## the number of objects
  N <- max(Dat[,1], Dat[,2])
  tol_rep <- nrow(Dat)
  tol_pair <- nrow(num.pair)
  p.cov <- ncol(Dat) - 5
  
  # extract surivival observations and covariate matrix
  nu <- Dat[,4]
  X <- Dat[,4+1:p.cov]
  Y <- Dat[,ncol(Dat)]
  
  ## Number of time knots
  knots<-c(6,18,30,48,90)
  knots0 <- c(0, knots)
  diff.knots<-diff(knots0)
  M <- length(knots)
  
}

if(proj == "enron"){
  
  ver.output<-paste0("K",K,"_",proj,"_v",v,"_r",r)
  
  if(v == 1){
    load("Data/Response_time_collection.RData")
    
  }
  
  if(v== 2){
    load("Data/Respons_confidential.RData")
    sur_collect <- sur_content
  }
  
  if(v==3){
    load("Data/Response_time.RData")
    sur_collect <- sur_time
  }

  tol_rep <- nrow(sur_collect)
  p.cov <- ncol(sur_collect) - 5
  usr_no_list <- as.numeric(names(table(sur_collect[,1:2])))
  N <-length(usr_no_list)
  
  # order the response time by senders and receivers
  Dat <- sur_collect
  Dat[,1] <- factor(sur_collect[,1], levels = usr_no_list)
  Dat[,2] <- factor(sur_collect[,2], levels = usr_no_list)
  ord_rep <- order(Dat[,1],Dat[,2])
  Dat <- Dat[ord_rep,]
  
  num.pair<-NULL# record the sender, receiver and num of reply 
  for (i in 2:tol_rep){
    if((Dat[i,1]!=Dat[i-1,1])||(Dat[i,2]!=Dat[i-1,2])){
      num.pair<-rbind(num.pair,Dat[i-1,1:3])
    }
  }
  num.pair<-rbind(num.pair,Dat[tol_rep,1:3])
  tol_pair<-nrow(num.pair)
  
  # extract surivival observations and covariate matrix
  nu <- Dat[,4]
  X <- Dat[,4 + 1:p.cov]
  Y <- Dat[,ncol(Dat)]
  
  #knots selection
  Y_res <- Y[which(nu==1)]
  # M<-10
  knots<-rep(NA,M)
  ordering<-order(Y_res)
  nj<-floor(length(Y_res)/M)
  for(p in 1:(M-1)){
    knots[p]<-Y_res[ordering[nj*p]]
  }
  knots[M]<- 21 * 24
  diff.knots<-diff(c(0,knots))
  knots0 <- c(0,knots)
  
}

print("Compute the auxiliary variables...")
## Generate auxiliary value p_ij indicates the interval Y_ij belongs to
Ind <- rep(0, tol_rep)
for(i in 1:tol_rep){
  Ind[i] <- sum(Y[i] > knots0[1:M])
}

## Construct two auxillary matrix for accumulating nu and Y - s[p-1] or s[p] - s[p-1]
A <- matrix(NA,tol_pair,M)
Mat_Y <- list()
index_rep <- 0
for(i in 1:tol_pair){
  nij <- num.pair[i,3]
  temp <- matrix(0, nij, M)
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
B <- matrix(NA,tol_pair,M)
Mat_nu <- list()
index_rep <- 0
for(i in 1:tol_pair){
  nij <- num.pair[i,3]
  temp <- matrix(0, nij, M)
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

## the third auxiliary variable tilde A_{ij} = sum_k=1^n_{ij} nu_{ijk} A_{ijk}
til_A <- matrix(NA,tol_pair,M)
index_rep <- 0
for(i in 1:tol_pair){
  nij <- num.pair[i,3]
  til_A[i, ] <- apply(Mat_Y[[i]] * nu[index_rep + 1:nij], 2, sum)
  index_rep <- index_rep + nij
}

## the forth auxiliary variable tilde B_{ij} = sum_k=1^n_{ij} nu_{ijk} X_{ijk}
til_B <- matrix(NA,tol_pair,p.cov)
index_rep <- 0
for(i in 1:tol_pair){
  nij <- num.pair[i,3]
  if(nij == 1){
    til_B[i, ] <- X[index_rep + 1:nij,] * nu[index_rep + 1:nij]
  }else{
    til_B[i, ] <- apply(X[index_rep + 1:nij,] * nu[index_rep + 1:nij], 2, sum)
  }
  index_rep <- index_rep + nij
}

# Set initial values
print("Set hyperparameters and the initial values of parameters...")
## the overall mixture proportion xi
xi.init <- rep(0.5,K)

## regression coefficients beta_lm
beta.init <- array(0,dim = c(K, K, p.cov)) 

## the prior for lambda by empirical estimator
fit_all <- survfit(Surv(Y,nu) ~ 1)
sur_all <- cbind(fit_all$time, fit_all$surv)

cure_rate_rough <- sur_all[nrow(sur_all),2]
beta.init[1] <- log(-log(cure_rate_rough))


#### Sample latent variable N_ij from Poisson(S_cure(y_ij|lambda)) + nu_ij
N_prior <- rep(NA, tol_rep)
index_rep <- 0 # index (i,j,1) - 1 of Y matrix 
for(i in 1:tol_pair){
  nij <-  num.pair[i, 3]
  
  temp <- rep(NA, nij)
  for(k in 1:nij){
    temp[k] <- sum(sur_all[,1] <= Y[index_rep + k])
  }
  
  sur_prob <- sur_all[temp,2]
  sur_prob[which(is.na(sur_prob))] <- 0
  
  N_prior[index_rep + 1:nij] <- rpois(nij, log(sur_prob/cure_rate_rough)) + nu[index_rep + 1:nij]
  index_rep <- index_rep + nij
}

### set the MLE of lambda of complete likelihood as the intial values
num_lam <- 0
den_lam <- 0
index_rep <- 0
for(i in 1:tol_pair){
  nij <-  num.pair[i, 3]
  num_lam <- num_lam + B[i,]
  den_lam <- den_lam + as.vector(N_prior[index_rep + 1:nij] %*% Mat_Y[[i]])
  index_rep <- index_rep + nij
}
lambda.prior <- num_lam / den_lam

## component-specific piecewise constant hazards lambda_lm
l_guess <- rep(NA,tol_pair)
m_guess <- rep(NA,tol_pair)

for(i in 1:tol_pair){
  l_guess[i] <- sample(1:K, size = 1, replace = TRUE, prob = xi.init/sum(xi.init))
  m_guess[i] <- sample(1:K, size = 1, replace = TRUE, prob = xi.init/sum(xi.init))
}

### sample N_ij | w_i = k from Poisson(S_cure(y_ij|lambda_k)) + nu_ij
#### fit a Kaplan Meier curve without covariates for each component (l,m)
Y_k <- list()
Y_k[[K * K]] <- NA
nu_k <- list()
nu_k[[K * K]] <- NA

index_rep <- 0
for(i in 1:tol_pair){
  
  l <- l_guess[i]
  m <- m_guess[i]
  nij <-  num.pair[i, 3]
  mem_index <- (l-1) * K + m
  
  Y_k[[mem_index]] <- c(Y_k[[mem_index]], Y[index_rep + 1:nij])
  nu_k[[mem_index]] <- c(nu_k[[mem_index]], nu[index_rep + 1:nij])
  index_rep <- index_rep + nij
  
}

Y_k[[K * K]] <- Y_k[[K * K]][-1]
nu_k[[K * K]] <- nu_k[[K * K]][-1]


sur_infor <- list()
for(k in 1:(K * K)){
  fit_k <- survfit(Surv(Y_k[[k]],nu_k[[k]]) ~ 1)
  sur_infor[[k]] <- cbind(fit_k$time, fit_k$surv)
}

#### Sample latent variable N_ij from Poisson(S_cure(y_ij|lambda_k)) + nu_ij
N_init <- rep(NA, tol_rep)
index_rep <- 0 # index (i,j,1) - 1 of Y matrix 
for(i in 1:tol_pair){
  nij <-  num.pair[i, 3]
  l <- l_guess[i]
  m <- m_guess[i]
  nij <-  num.pair[i, 3]
  mem_index <- (l-1) * K + m
  temp <- match(Y[index_rep + 1:nij],sur_infor[[mem_index]][,1])
  sur_prob <- sur_infor[[mem_index]][temp,2]
  sur_prob[which(is.na(sur_prob))] <- 0
  
  N_init[index_rep + 1:nij] <- rpois(nij, sur_prob) + nu[index_rep + 1:nij]
  index_rep <- index_rep + nij
}

### set the MLE of lambda of complete likelihood as the intial values
sum_nu <- array(0, dim = c(K, K, M))
sum_y <- array(0, dim = c(K, K, M))
index_rep <- 0
for(i in 1:tol_pair){
  l <- l_guess[i]
  m <- m_guess[i]
  nij <-  num.pair[i, 3]
  mem_index <- (l-1) * K + m
  
  sum_nu[l,m, ] <- sum_nu[l,m, ] + B[i,]
  sum_y[l,m, ] <- sum_y[l,m, ] + as.vector(N_init[index_rep + 1:nij] %*% Mat_Y[[i]])
  index_rep <- index_rep + nij
}
lambda.init <- sum_nu / sum_y

for(l in 1:K){
  for(m in 1:K){
    for(p in 1:M){
      lambda.init[l,m,p] <- max(lambda.init[l,m,p],exp(-100))
    }
  }
}

xi.iter <- t(xi.init)

lambda.iter <- matrix(NA, K^2, M)
for(l in 1:K){
  for(m in 1:K){
    lambda.iter[(l-1) * K + m,] <- lambda.init[l,m,]
  }
}

beta.iter <- matrix(NA, K^2, p.cov)
for(l in 1:K){
  for(m in 1:K){
    beta.iter[(l-1) * K + m,] <- beta.init[l,m,]
  }
}

mp.iter <- (l_guess - 1) * K + m_guess - 1
pi.iter <- matrix(1/K, K, N)
N.iter <- N_init

print("Start MCMC algorithm...")
# load the c++ dynamics library
dyn.load("/home/songfangda/SurNet/src/SMMB_revised.so")

########################
# Apply MCMC algorithm #
########################
num_record <- floor(num_iter/num_thin)
num_record_after_burnin <- floor((num_iter - burnin)/num_thin)
num_record_befor_burnin <- floor(burnin/num_thin)

log_like_record <- rep(0, num_record)
WAIC_record <- rep(NA, num_record_after_burnin)
ind_update_xi <- 1

if(compute_DIC){
  DIC_record <- matrix(NA, num_record_after_burnin, 2)
  RunTime <- rep(0, 9)
}else{
  RunTime <- rep(0, 7)
}

# iter_update_xi <- 0

# N_record <- list()
mp_record <- matrix(NA, num_record, tol_pair)
beta_record <- array(NA, dim = c(num_record, K * K, p.cov))
lambda_record <- array(NA, dim = c(num_record, K * K, M))
pi_record <- array(NA, dim = c(num_record, K, N))
xi_record <- matrix(NA, num_record, K)
time_record <- NULL

# print(compute_DIC)

start_time <- Sys.time()
ind_record <- 1
for(iter in 1:num_iter){
  
  MCMC_seed <- round(runif(1,1,10000))
  # print(paste0("Run on ",iter,"-th iteration."))
  if(compute_DIC){
    results <-.C("MCMC_iteration_DIC", 
                 # seed
                 seed = as.integer(MCMC_seed),
                 nc = as.integer(nc),
                 MCMC_iter = as.integer(iter),
                 thin_iter = as.integer(num_thin),
                 ind_swap = as.integer(iter %% num_iter_swap == 0),
                 # dim information
                 num_usr = as.integer(N),
                 num_mem = as.integer(K),
                 num_com = as.integer(tol_pair),
                 num_eml = as.integer(tol_rep),
                 num_int = as.integer(M),
                 num_cov = as.integer(p.cov),
                 # knots
                 knots = as.double(knots0),
                 # communication information
                 com_list <- as.integer(t(num.pair)),
                 sur_time <- as.double(Y),
                 cen_ind <- as.integer(nu),
                 cov <- as.double(t(X)),
                 interval <- as.integer(Ind), 
                 # parameters
                 xi = as.double(xi.iter), 
                 lambda = as.double(lambda.iter),
                 beta = as.double(beta.iter),
                 # prior / penatly
                 lambda.prior = as.double(lambda.prior/pri_lambda),
                 lambda_rate = as.double(pri_lambda),
                 beta.prior = as.double(c(beta.init[1], pri_beta)),
                 xi.prior = as.double(rep(pri_xi,2)),
                 update.xi = as.integer(ind_update_xi),
                 # auxiliary variables
                 sur_interval = as.double(t(til_A)),
                 fail_ind = as.integer(t(B)),
                 sur_cov = as.double(t(til_B)),
                 # output
                 mem_pair = as.integer(mp.iter),
                 pi = as.double(pi.iter),
                 intentions = as.integer(N.iter),
                 # loglikelihood
                 log_like = as.double(0),
                 ind_after_burnin = as.integer(iter > burnin),
                 resample_iter = as.integer(re_iter),
                 DIC = as.double(rep(0,2)),
                 WAIC = as.double(0),
                 Time = as.double(RunTime))
  }else{
    # print("Run MCMC but not compute DIC")
    results <-.C("MCMC_iteration", 
                 # seed
                 seed = as.integer(MCMC_seed),
                 nc = as.integer(nc),
                 MCMC_iter = as.integer(iter),
                 thin_iter = as.integer(num_thin),
                 ind_swap = as.integer(iter %% num_iter_swap == 0),
                 # dim information
                 num_usr = as.integer(N),
                 num_mem = as.integer(K),
                 num_com = as.integer(tol_pair),
                 num_eml = as.integer(tol_rep),
                 num_int = as.integer(M),
                 num_cov = as.integer(p.cov),
                 # knots
                 knots = as.double(knots0),
                 # communication information
                 com_list <- as.integer(t(num.pair)),
                 sur_time <- as.double(Y),
                 cen_ind <- as.integer(nu),
                 cov <- as.double(t(X)),
                 interval <- as.integer(Ind), 
                 # parameters
                 xi = as.double(xi.iter), 
                 lambda = as.double(lambda.iter),
                 beta = as.double(beta.iter),
                 # prior / penatly
                 lambda.prior = as.double(lambda.prior/pri_lambda),
                 lambda_rate = as.double(pri_lambda),
                 beta.prior = as.double(c(beta.init[1], pri_beta)),
                 xi.prior = as.double(rep(pri_xi,2)),
                 update.xi = as.integer(ind_update_xi),
                 # auxiliary variables
                 sur_interval = as.double(t(til_A)),
                 fail_ind = as.integer(t(B)),
                 sur_cov = as.double(t(til_B)),
                 # output
                 mem_pair = as.integer(mp.iter),
                 pi = as.double(pi.iter),
                 intentions = as.integer(N.iter),
                 # loglikelihood
                 log_like = as.double(0),
                 ind_after_burnin = as.integer(iter > burnin),
                 WAIC = as.double(0),
                 Time = as.double(RunTime))
    
  }
  
  mp.temp <- results$mem_pair
  pi.temp <- matrix(results$pi, K, N)
  N.temp <- results$intentions
  lambda.temp <- matrix(results$lambda, K^2, M)
  beta.temp <- matrix(results$beta, K^2, p.cov)
  xi.temp <- results$xi
  
  mp.iter <- mp.temp
  pi.iter <- pi.temp
  N.iter <- N.temp
  lambda.iter <- lambda.temp
  beta.iter <- beta.temp
  xi.iter <- xi.temp
  
  if(iter %% num_thin == 0){
    log_like_record[ind_record] <- results$log_like
    if(iter > burnin){
      if(compute_DIC){
        DIC_record[ind_record - num_record_befor_burnin, 1] <- results$DIC[1]
        DIC_record[ind_record - num_record_befor_burnin, 2] <- results$DIC[2]
      }
      WAIC_record[ind_record - num_record_befor_burnin] <- results$WAIC[1]
    }
    time_record <- rbind(time_record, results$Time)
    
    # record the posterior sampling
    mp_record[ind_record,] <- mp.iter
    beta_record[ind_record,,] <- beta.iter
    lambda_record[ind_record,,] <- lambda.iter
    pi_record[ind_record,,] <- pi.iter
    xi_record[ind_record,] <- xi.iter
    ind_record <- ind_record + 1
  }
  
  if(iter %% (100 * num_thin) == 0){
    cur_time <- Sys.time()
    print(paste0("Finish ", iter, " iterations."))
    print(paste0("The observed log-likelihood is ",log_like_record[ind_record - 1],"."))
    print(paste0("It has taken ",difftime(cur_time, start_time, units = "mins"), " mins."))
  }
  
}

if(compute_DIC){
  colnames(time_record) <- c("Allocate memory", "Update N", "Update beta",
                             "Update lambda", "Update pi", "Swapping step", 
                             "Compute likelihood",
                             "Compute expected deviance", 
                             "Compute fixed-point deviance")
}else{
  colnames(time_record) <- c("Allocate memory", "Update N", "Update beta",
                             "Update lambda", "Update pi", "Swapping step", 
                             "Compute likelihood")
}

if(!dir.exists("Inference")){
  dir.create("Inference")
}

# save.image(paste0("Inference/MCMC_sampling_",proj,"_SMMB_K",K,"_v",v,"_r",r,".RData"))

# posterior inference
print("Conducting posterior inference...")

iter_after_burnin <- num_record_befor_burnin + 1:num_record_after_burnin
if(K > 1){
  # posterior mean
  lambda.est <- apply(lambda_record[iter_after_burnin,,], c(2,3), mean)
  beta.est <- apply(beta_record[iter_after_burnin,,], c(2,3), mean)
  xi.est <- apply(xi_record[iter_after_burnin,], 2, mean)
  pi.est <- apply(pi_record[iter_after_burnin,,], c(2,3), mean)
  mp.est <- apply(mp_record[iter_after_burnin,], 2, getmode) + 1
  
  # posterior sd
  lambda.sd <- apply(lambda_record[iter_after_burnin,,], c(2,3), sd)
  beta.sd <- apply(beta_record[iter_after_burnin,,], c(2,3), sd)
  xi.sd <- apply(xi_record[iter_after_burnin,], 2, sd)
  pi.sd <- apply(pi_record[iter_after_burnin,,], c(2,3), sd)
  
  # posterior CI
  lambda.ci <- array(NA, dim = c(K * K, M, 2))
  beta.ci <- array(NA, dim = c(K * K, p.cov, 2))
  xi.ci <- matrix(NA, K, 2)
  pi.ci <- array(NA, dim = c(K, N, 2))
  
  for(l in 1:K){
    for(m in 1:K){
      ind_mem <- (l-1) * K + m
      for(p in 1:M){
        lambda.ci[ind_mem,p,] <- as.numeric(ci(lambda_record[iter_after_burnin,ind_mem,p], method = "ETI", ci = .95))[2:3]
      }
      for(p in 1:p.cov){
        beta.ci[ind_mem,p,] <- as.numeric(ci(beta_record[iter_after_burnin,ind_mem,p], method = "ETI", ci = .95))[2:3]
      }
    }
    xi.ci[l, ] <- as.numeric(ci(xi_record[iter_after_burnin,l], method = "ETI", ci = .95))[2:3]
    for(i in 1:N){
      pi.ci[l,i,] <- as.numeric(ci(pi_record[iter_after_burnin,l,i], method = "ETI", ci = .95))[2:3]
    }
  }
}else{
  # posterior mean
  lambda.est <- matrix(apply(lambda_record[iter_after_burnin,,], 2, mean), 1, M)
  beta.est <- matrix(apply(beta_record[iter_after_burnin,,], 2, mean),1, p.cov)
  xi.est <- mean(xi_record[iter_after_burnin,])
  pi.est <- matrix(apply(pi_record[iter_after_burnin,,], 2, mean),1,N)
  mp.est <- apply(mp_record[iter_after_burnin,], 2, getmode) + 1
  
  # posterior sd
  lambda.sd <- apply(lambda_record[iter_after_burnin,,], 2, sd)
  beta.sd <- apply(beta_record[iter_after_burnin,,], 2, sd)
  xi.sd <- sd(xi_record[iter_after_burnin,])
  pi.sd <- apply(pi_record[iter_after_burnin,,], 2, sd)
  
  # posterior CI
  lambda.ci <- array(NA, dim = c(K * K, M, 2))
  beta.ci <- array(NA, dim = c(K * K, p.cov, 2))
  xi.ci <- matrix(NA, K, 2)
  pi.ci <- array(NA, dim = c(K, N, 2))
  
  for(l in 1:K){
    for(m in 1:K){
      ind_mem <- (l-1) * K + m
      for(p in 1:M){
        lambda.ci[ind_mem,p,] <- as.numeric(ci(lambda_record[iter_after_burnin,ind_mem,p], method = "ETI", ci = .95))[2:3]
      }
      for(p in 1:p.cov){
        beta.ci[ind_mem,p,] <- as.numeric(ci(beta_record[iter_after_burnin,ind_mem,p], method = "ETI", ci = .95))[2:3]
      }
    }
    xi.ci[l, ] <- as.numeric(ci(xi_record[iter_after_burnin,l], method = "ETI", ci = .95))[2:3]
    for(i in 1:N){
      pi.ci[l,i,] <- as.numeric(ci(pi_record[iter_after_burnin,l,i], method = "ETI", ci = .95))[2:3]
    }
  }
}


# BIC
log_like <- 0
log_p <- array(0, dim = c(tol_pair, K, K))
ind_eml <- 0

for(i in 1:tol_pair){
  sender <- num.pair[i,1]
  receiver <- num.pair[i,2]
  n_eml <- num.pair[i,3]
  for(l in 1:K){
    for(m in 1:K){
      
      ind_mem <- (l - 1) * K + m
      log_p[i,l,m] <- log(pi.est[l,sender]) + log(pi.est[m,receiver])
      
      #### Corresponding to nu_ijk log(lambda_{lmp_{ijk}}) = (log(lambda_{lm}))^T B_{ij}
      log_p[i,l,m] <- log_p[i,l,m] + sum(B[i,] * log(lambda.est[ind_mem,]))
      
      #### Corresponding to nu_ij (- lambda_p(y_ij - s_{p-1}) - sum_q = 1 ^ p-1 lambda_q (s_q-s_{q-1}))
      ####                             = lambda_{lm}^T tilde A_{ij}
      log_p[i,l,m] <- log_p[i,l,m] - sum(til_A[i,] * lambda.est[ind_mem,])
      
      # Corresponding to nu_ij x_ijk^T beta_lm
      log_p[i,l,m] <- log_p[i,l,m] + sum(til_B[i,] * beta.est[ind_mem,])
      
      for(k in 1:n_eml){
        
        int_time <- Ind[ind_eml + k]
        
        ## exp(x_{ijk}^T beta_{lm}) 
        exp_xbeta <- exp(sum(X[ind_eml + k,] * beta.est[ind_mem,]))
        
        ## ( 1- S_pch(y_ij | lambda_lm) )
        if(int_time > 1){
          fail_prob <- exp(-sum(lambda.est[ind_mem, 1:(int_time-1)] * diff.knots[1:(int_time-1)]))
          fail_prob <- fail_prob * exp(-lambda.est[ind_mem, int_time] * (Y[ind_eml + k] - knots0[int_time]))
        }else{
          fail_prob <- exp(-lambda.est[ind_mem, int_time] * (Y[ind_eml + k] - knots0[int_time]))
        }
        
        fail_prob <- 1.0 - fail_prob
        
        log_p[i,l,m] <- log_p[i,l,m] - exp_xbeta * fail_prob
        
      }
    }
  }
  max_log <- max(log_p[i,,])
  log_like <- log_like + max_log + log(sum(exp(log_p[i,,] - max_log)))
  ind_eml <- ind_eml + n_eml
}

sum_xi <- sum(xi.est)
log_like <- log_like + N * lgamma(sum(xi.est))
log_like <- log_like - N * sum(lgamma(xi.est))
log_like <- log_like + sum((xi.est - 1) * apply(log(pi.est), 1, sum))

BIC <- - 2 * log_like + (K + K^2 * (M + p.cov)) * log(tol_rep)

# DIC 
if(compute_DIC){
  print("Calculating DIC...")
  pD <- 2 * (mean(DIC_record[,2]) - mean(DIC_record[,1]))
  DIC <- -2 * mean(DIC_record[,1]) + pD
  message(paste0("The DIC is ",DIC," with p_D = ", pD,"."))
}

# WAIC
max_WAIC <- max(WAIC_record)
WAIC <- -4 * (max_WAIC + log(mean(exp(WAIC_record - max_WAIC))))  + 2 * mean(WAIC_record)
end_time <- Sys.time()
message(paste0("The total running time is ",difftime(end_time, start_time, units = "mins"), " mins."))

print("Saving the whole workspace...")
save.image(paste0("Inference/Inference_",proj,"_SMMB_K",K,"_v",v,"_r",r,".RData"))