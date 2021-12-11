rm(list=ls())
library(mclust)
library(xtable)

# load the estimated role pairs under each setting
K <- 3

wkspace_name <- paste0("Inference/Inference_simulation_SMMB_K",K,"_v1_r1.RData")
if(file.exists(wkspace_name)){
  load(wkspace_name)
}else{
  message(paste0("The inference of 1-th replication with K=",K," does not exist!"))
}

rolepair_record <- matrix(NA, tol_pair, 19)
rolepair_record[,1] <- unlist(mp.est)

for(r in 2:13){
  wkspace_name <- paste0("Inference/Inference_simulation_SMMB_K",K,"_v1_r",r,".RData")
  if(file.exists(wkspace_name)){
    load(wkspace_name)
    rolepair_record[,r] <- unlist(mp.est)
    
  }else{
    message(paste0("The inference of ",r,"-th replication with K=",K," does not exist!"))
  }
}

for(r in 1:6){
  wkspace_name <- paste0("Inference/Inference_simknots_SMMB_K",K,"_v1_r",r,".RData")
  if(file.exists(wkspace_name)){
    load(wkspace_name)
    rolepair_record[,r + 13] <- unlist(mp.est)
    
  }else{
    message(paste0("The inference of ",r,"-th replication with K=",K," does not exist!"))
  }
}

# Check the ARI for different hyperparameter values
# the hyperparameter kappa for lambda
kappa_sel <- c(30,40,60,70)
for(i in 1:4){
  print(paste0("The ARI between the estiamted role pairs with kappa = 50 and those with kappa = ",kappa_sel[i]," is ", adjustedRandIndex(rolepair_record[,1],rolepair_record[,i + 1])))
}


# the hyperparameter sigma^2 for beta
sigmasq_sel <- c(1,3,7,9)
for(i in 1:4){
  print(paste0("The ARI between the estiamted role pairs with sigma^2 = 5 and those with sigma^2 = ",sigmasq_sel[i]," is ", adjustedRandIndex(rolepair_record[,1],rolepair_record[,i + 5])))
}


# the hyperparameter a for xi
a_sel <- c(0.2,0.5,2,5)
for(i in 1:4){
  print(paste0("The ARI between the estiamted role pairs with a = 1 and those with a = ",sigmasq_sel[i]," is ", adjustedRandIndex(rolepair_record[,1],rolepair_record[,i + 9])))
}

# the number of knots J
J_sel <- c(3,4,5,6,7,8)
for(i in 1:6){
  print(paste0("The ARI between the estiamted role pairs with the true knots and those with the selected ",J_sel[i]," knots is ", adjustedRandIndex(rolepair_record[,1],rolepair_record[,i + 13])))
}
