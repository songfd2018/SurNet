rm(list=ls())
set.seed(6332)
v <- 1

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
if(!dir.exists("Images")){
  dir.create("Images")
}


# Number of components
K <- 3

# Number of time knots
knots<-c(6,18,30,48,90)
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
p.cov <- 3
beta <- array(NA, dim = c(K, K, p.cov))

beta[1,1,1] <- 0.8
beta[1,2,1] <- 0.3
beta[1,3,1] <- -0.2
beta[2,1,1] <- 1.2
beta[2,2,1] <- 0.6
beta[2,3,1] <- -0.6
beta[3,1,1] <- 1.6
beta[3,2,1] <- 0
beta[3,3,1] <- -1

beta[1,1,2] <- 0.3 
beta[1,2,2] <- 0
beta[1,3,2] <- 0.6
beta[2,1,2] <- 0.9
beta[2,2,2] <- 0
beta[2,3,2] <- -0.3
beta[3,1,2] <- -0.9
beta[3,2,2] <- -0.6
beta[3,3,2] <- 0

beta[1,1,3] <- 1 
beta[1,2,3] <- 0.5
beta[1,3,3] <- 0
beta[2,1,3] <- -0.5
beta[2,2,3] <- -1
beta[2,3,3] <- 1.5
beta[3,1,3] <- -1.5
beta[3,2,3] <- 0
beta[3,3,3] <- -2

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
num_community <- 4
community_nodes <- c(25, 30, 40, 55)
prob_diff<- 0.08 #each edges has the same prob to be chosen
graph_sparse <- erdos.renyi.game(n = N, p.or.m = prob_diff, directed = TRUE)
prob_same <- 0.32
num.pair <- get.edgelist(graph_sparse)
graph_dense <- list()
node_cum <- 0
for(c in 1:num_community){
  graph_dense[[c]] <- erdos.renyi.game(n = community_nodes[c], p.or.m = prob_same, directed = TRUE)
  num.pair <- rbind(num.pair, get.edgelist(graph_dense[[c]]) + node_cum)
  node_cum <- node_cum + community_nodes[c]
}

# exclude the repeat edges
num.pair <- num.pair[order(num.pair[,1],num.pair[,2]),]
num.pair <- as.data.frame(num.pair)
colnames(num.pair) <- c("i","j")
num.pair <- as.matrix(unique(num.pair[c("i", "j")]))
tol_pair <- nrow(num.pair)

num.pair <- cbind(num.pair, rnbinom(tol_pair, size = 3, prob = 0.15) + 2)
colnames(num.pair) <- c("i", "j", "nij")

Dat <- NULL

# Contruct an observation matrix with p + 6 columns, where p denotes the number of covariates
# the first column is the sender label
# the second column is the receiver label
# the third column is the index of response time for a given pair of users
# the forth column is the censoring indicator
# the fifth to (p+5)th columns include the covariate matrix
# the last column denotes the survival time 

for(i in 1:tol_pair){
  sender <- num.pair[i,1]
  receiver <- num.pair[i,2]
  nij <- num.pair[i,3]
  pair.reply <- cbind(rep(sender,nij), rep(receiver,nij), 1:nij)
  Dat <- rbind(Dat, pair.reply)
}
tol_rep<- nrow(Dat)

# nu
Dat <- cbind(Dat,NA)

#covariate
x<-matrix(1,tol_rep,p.cov)
x[,2]<-rbinom(tol_rep,1,0.3)#1 containing weekend and 0 no weekend
x[,3]<-rnorm(tol_rep,0,0.5)
Dat<-cbind(Dat,x)

#synthetic parameter
xi.syn<-c(0.2,0.4,0.6)
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
Cen_Time <- runif(tol_rep, 0, 100)

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

# draw the adjacent matrix of the synthetic network
dat_numpair_sim <- data.frame(Sender = factor(num.pair[,1]),
                              Receiver = factor(num.pair[,2]),
                              Weights = num.pair[,3])

pdf(paste0("Images/Adjacent_matrix_grey_v",v,".pdf"),width = 6, height = 6)
p_grey <- ggplot(dat_numpair_sim, aes(x = Sender, y = Receiver)) +
  geom_tile() +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size=32),
        axis.ticks = element_blank())
p_grey
dev.off()

# summary
xi_output <- paste0("(",xi.syn[1],",",xi.syn[2],",",xi.syn[3],")")
print(paste0("A random network G(",N,",",prob_diff,") for different communities and G(",N,",",prob_same + prob_diff,
             ") for the same communities are generated by igraph R package with ",
             tol_pair," directed edges, ", nrow(Dat), " events, ",sum(Dat[,4])," failure time and ",
             nrow(Dat) - sum(Dat[,4]), " censoring time and xi = ",xi_output,"."))

# check the censoring rate
# 1 - mean(Dat[,4])
# check the non-cure censoring rate
# non_cure_censor <- Failure_Time > Cen_Time & is.finite(Failure_Time)
# mean(non_cure_censor)

# save the workspace
if(!dir.exists("ObsData")){
  dir.create("ObsData")
}
save(Dat,num.pair,file = paste0("ObsData/Obs_SMMB_v",v,".RData"))
save.image(file = paste0("ObsData/workspace_SMMB_v",v,".RData"))