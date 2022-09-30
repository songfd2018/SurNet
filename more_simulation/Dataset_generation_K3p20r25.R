rm(list=ls())
set.seed(6332)
proj <- "sim"
data_rep <- 100
v <- "K3p20r25"

setwd(v)

library("igraph") # generate random graph
library("ggplot2")

# calculate lambda from the survival probability at each knot
cal.lambda <- function(sur_prob,len_int){
  
  #sur_prob the survival probablity at each knot, 
  #len_int the time length of intervals
  
  n<-length(sur_prob)
  sur_pre <- 1 #survival funciton at time 0
  
  lambda<-rep(0,n)
  for (i in 1:n){
    #S(s[p])/S(s[p-1]) = exp(-lambda[p] * (s[p] - s[p-1]))
    lambda[i] <- -log(sur_prob[i] / sur_pre) / len_int[i]
    sur_pre <- sur_prob[i]
  }
  return(lambda)
}

# generate samples from Dirichelet distribution
rDiri<-function(num,alpha){
  
  #num is # of sample, alpha is the parameter
  
  k<-length(alpha)
  result<-matrix(rep(0,k*num),num,k)
  for (j in 1:k){
    result[,j]<-rgamma(num,alpha[j])#sample from gamma(alpha[i],1)
  }
  result<-result/apply(result,1,sum)
  return(result)
}

# sampling survival time from cure rate model
# by inverse sampling
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
  
  # calculate the survival probabilities of PCH distributions at knots
  
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

# Create a folder to save the figures
if(!dir.exists(v)){
  dir.create(v)
}

# Number of components
K <- 3
xi.syn<-c(0.2,0.3,0.4)

#each edge has the same prob to be chosen
prob_edge <- 0.20 

# Number of time knots
knots<-c(6,24,60,90)
P <- length(knots)
diff.knots<-diff(c(0,knots))

# Setting the survival probabilities at the end of each time interval for PCH distribution
sur_interval <- array(NA, dim = c(K, K, P))
lambda <- array(NA, dim = c(K, K, P))

for(l in 1:K){
  for(m in 1:K){
    rand_sur <- runif(P)
    sur_interval[l, m, ] <- rand_sur[order(rand_sur, decreasing = TRUE)]
    if(sur_interval[l, m, P] < 1 - 10e-4){
      relative_failure <- (1 - sur_interval[l, m,]) / (1 - sur_interval[l, m, P]) * (1 - 10e-4)
      sur_interval[l, m, ] <- 1 - relative_failure
    }
    lambda[l, m, ] <- cal.lambda(sur_interval[l, m, ], diff.knots)
  }
}

# Setting cure rate parameter beta
p.cov <- 2
beta <- array(runif(K * K * p.cov, -1, 1), dim = c(K, K, p.cov))

# Calculate the survival probabilities of cure rate model
sur_knot <- array(NA, dim= c(K, K, P, 2))
for(l in 1:K){
  for(m in 1:K){
    sur_knot[l, m, ,1] <- exp(-exp(beta[l, m, 1]) * (1 - sur_interval[l, m, ]))
    sur_knot[l, m, ,2] <- exp(-exp(beta[l, m, 1] + beta[l, m, 2]) * (1 - sur_interval[l, m, ]))
  }
}

# Generate observations
N <- 150 # number of objects

#synthetic observed data
# num_community <- 4
# community_nodes <- c(25, 30, 40, 55)

graph_sparse <- erdos.renyi.game(n = N, p.or.m = prob_edge, directed = TRUE)
num.pair <- get.edgelist(graph_sparse)

#graph_dense <- list()
# node_cum <- 0
# for(c in 1:num_community){
#   graph_dense[[c]] <- erdos.renyi.game(n = community_nodes[c], p.or.m = prob_same, directed = TRUE)
#   num.pair <- rbind(num.pair, get.edgelist(graph_dense[[c]]) + node_cum)
#   node_cum <- node_cum + community_nodes[c]
# }

# exclude the repeat edges
num.pair <- num.pair[order(num.pair[,1],num.pair[,2]),]
num.pair <- as.data.frame(num.pair)
colnames(num.pair) <- c("i","j")
# num.pair <- as.matrix(unique(num.pair[c("i", "j")]))
tol_pair <- nrow(num.pair)

num.pair <- cbind(num.pair, rpois(tol_pair, 25) + 5)
colnames(num.pair) <- c("i", "j", "nij")

# covariate
tol_rep<- sum(num.pair[,3])
x<-matrix(1,tol_rep,p.cov)
x[,2]<-rnorm(tol_rep,0,0.5)

# draw the adjacent matrix of the synthetic network
# dat_numpair_sim <- data.frame(Sender = factor(num.pair[,1]),
#                               Receiver = factor(num.pair[,2]),
#                               Weights = num.pair[,3])
# 
# pdf(paste0("Images/Adjacent_matrix_grey_v",v,".pdf"),width = 6, height = 6)
# p_grey <- ggplot(dat_numpair_sim, aes(x = Sender, y = Receiver)) +
#   geom_tile() +
#   theme_minimal() +
#   theme(axis.text.x = element_blank(),
#         axis.text.y = element_blank(),
#         text = element_text(size=32),
#         axis.ticks = element_blank())
# p_grey
# dev.off()

################################
# Generate replicated datasets #
################################
seed_vec <- round(runif(data_rep) * 10000)

for(rep in 1:data_rep){
  
  set.seed(seed_vec[rep])
  dir_name <- paste0(v,"/",proj,"_v",rep)
  if(!dir.exists(dir_name)){
    dir.create(dir_name)
  }
  setwd(dir_name)
  
  # generate synthetic dataset
  pi.syn<-rDiri(N,xi.syn)
  
  l.syn<-rep(0,tol_pair)
  m.syn<-rep(0,tol_pair)
  
  Dat.l <- NULL
  Dat.m <- NULL
  for (i in 1:tol_pair){
    l.syn[i]<-which(rmultinom(1,1,pi.syn[num.pair[i,1],])==1)
    m.syn[i]<-which(rmultinom(1,1,pi.syn[num.pair[i,2],])==1)
    Dat.l<-c(Dat.l,rep(l.syn[i],num.pair[i,3]))
    Dat.m<-c(Dat.m,rep(m.syn[i],num.pair[i,3]))
  }
  w.syn <- cbind(l.syn,m.syn)
  
  # Contruct an observation matrix with p + 6 columns, where p denotes the number of covariates
  # the first column is the sender label
  # the second column is the receiver label
  # the third column is the index of response time for a given pair of users
  # the forth column is the censoring indicator
  # the fifth to (p+5)th columns include the covariate matrix
  # the last column denotes the survival time 
  Dat <- NULL
  
  for(i in 1:tol_pair){
    sender <- num.pair[i,1]
    receiver <- num.pair[i,2]
    nij <- num.pair[i,3]
    pair.reply <- cbind(rep(sender,nij), rep(receiver,nij), 1:nij)
    Dat <- rbind(Dat, pair.reply)
  }
  
  # nu
  Dat <- cbind(Dat,NA)

  #1 containing weekend and 0 no weekend
  #x[,3]<-rnorm(tol_rep,0,0.5)
  Dat<-cbind(Dat,x)
  
  # Inverse sampling for survival time
  Failure_Time <- NULL
  index_rep <- 0
  for(i in 1:tol_pair){
    l <- l.syn[i]
    m <- m.syn[i]
    nij <-  num.pair[i, 3]
    Failure_Time <- c(Failure_Time,rCure(nij, lambda[l, m, ], knots, #exp(beta[l,m,1])))
                                         exp(Dat[index_rep + 1:nij, 4 + 1:p.cov] %*% beta[l, m, ])))
    index_rep <- index_rep + nij
  }
  
  # Censoring uniform on (0,C) 
  # choosing different C for different censoring rates
  Cen_Time <- runif(tol_rep, 50, 150)
  
  Dat[,4] <- Failure_Time < Cen_Time
  
  Sur_obs <- rep(NA, tol_rep)
  for(i in 1:tol_rep){
    if(Dat[i,4]){
      Sur_obs[i] <- Failure_Time[i]
    }else{
      Sur_obs[i] <- Cen_Time[i]
    }
  }
  Dat <- cbind(Dat, Sur_obs)
  
  colnames(Dat) <- c("i","j","k","nu",paste0("x_",1:p.cov),"Y")
  
  # summary
  xi_output <- "("
  for(k in 1:K){
    if(k <K){
      xi_output <- paste0(xi_output,xi.syn[k],",")
    }else{
      xi_output <- paste0(xi_output,xi.syn[k],")")
    }
  }
  
  print(paste0("The ",rep,"-th random network G(",N,",",prob_edge,") is generated by igraph R package with ",
               tol_pair," directed edges, ", nrow(Dat), " events, ",sum(Dat[,4])," failure time and ",
               nrow(Dat) - sum(Dat[,4]), " censoring time and xi = ",xi_output,"."))
  
  # save the workspace
  if(!dir.exists("ObsData")){
    dir.create("ObsData")
  }
  save(Dat,num.pair,file = paste0("ObsData/Obs_SMMB_",proj,"_v",rep,".RData"))
  save.image(file = paste0("ObsData/workspace_SMMB_",proj,"_v",rep,".RData"))
}

