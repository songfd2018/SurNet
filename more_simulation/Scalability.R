# Study the scalability of our proposed MCMC algorithm
load("RunningTime.Rdata")
if(!dir.exists("Images")){
  dir.create("Images")
}

# vary K from 2 to 5
K_sel <- 2:5 
prob <- 0.3
ave_res <- 25

nk <- length(K_sel)

jpeg(paste0("Images/cpu_time_by_role.jpg"),width = 800, height = 600)
par(mar=c(5.1,5.1,4.1,2.1))
plot(K_sel,time_K,xlab= "The number of roles K",ylab = "Time (hours)",type="n",cex.axis=2.5,cex.lab=3)
points(K_sel,time_K,type="b",pch=19,cex=3)
dev.off()

# vary connectivity probability from 0.2 to 0.5
K <- 3 
prob_sel <- (2:5)/10
ave_res <- 25

np <- length(prob_sel)

jpeg(paste0("Images/cpu_time_by_prob.jpg"),width = 800, height = 600)
par(mar=c(5.1,5.1,4.1,2.1))
plot(prob_sel,time_prob,xlab= "The connectivity probabilities",ylab = "Time (hours)",type="n",cex.axis=2.5,cex.lab=3)
points(prob_sel,time_prob,type="b",pch=19,cex=3)
dev.off()

# vary average number of response from 0.2 to 0.5
K <- 3 
prob <- 0.3
ave_res_sel <- seq(10,35,5)

nr <- length(ave_res_sel)

jpeg(paste0("Images/cpu_time_by_ave_res.jpg"),width = 800, height = 600)
par(mar=c(5.1,5.1,4.1,2.1))
plot(ave_res_sel,time_response,xlab= "The average number of response times",ylab = "Time (hours)",type="n",cex.axis=2.5,cex.lab=3)
points(ave_res_sel,time_response,type="b",pch=19,cex=3)
dev.off()