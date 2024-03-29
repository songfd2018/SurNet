rm(list=ls())
library(RColorBrewer)
library(ggplot2)
library(mclust)
library(xtable)
library(bayestestR)

proj <- "enron"


if(!dir.exists("Images")){
  dir.create("Images")
}

#################################
# Determine the number of knots #
#################################
proj <- "enron"
v <- 1
K.sel <- 1:5

if(!dir.exists("Images")){
  dir.create("Images")
}

pD.record <- rep(NA,length(K.sel))
DIC.record <- rep(NA,length(K.sel))

for(ind in 1:length(K.sel)){
  K <- K.sel[ind]
  # Analyze the obtained results
  wkspace_name <- paste0("Inference/Inference_",proj,"_SMMB_K",K,"_v",v,"_r1.RData")
  if(file.exists(wkspace_name)){
    load(wkspace_name)
    
    # DIC
    pD.record[ind] <- 2 * (mean(DIC_record[,2]) - mean(DIC_record[,1]))
    DIC.record[ind] <- -2 * mean(DIC_record[,1]) + pD.record[ind]
    # message(paste0("The DIC is ",DIC.record[ind],"."))
    
  }else{
    # message(paste0("The inference of 1-th replication on dataset \"",proj,"\" with K=",K," does not exist!"))
  }
}

# record the DIC
names(pD.record) <- paste0("K=",K.sel)

names(DIC.record) <- paste0("K=",K.sel)

print(paste0("The optimal number of roles is ",names(DIC.record[which.min(DIC.record)]),"."))

############################################
# Analyze without confidential information #
############################################
ver <- 1
K_opt <- 2
load(paste0("Inference/Inference_",proj,"_SMMB_K",K_opt,"_v",v,"_r1.RData"))

######
# Pi #
######
# Analyze the relationship bewteen membership probabilities and positions of all employee
emp.list <- read.csv("RawData/employeelist.csv")

num_reply_employee <- rep(0,N)
for(i in 1:tol_pair){
  sender <- num.pair[i,1]
  receiver <- num.pair[i,2]
  nij <- num.pair[i,3]
  
  num_reply_employee[sender] <- num_reply_employee[sender] + nij
  num_reply_employee[receiver] <- num_reply_employee[receiver] + nij
}

# reorder the position from high to low
pos_factor <- factor(emp.list$position, levels = c("CEO", "President","Vice President","Managing Director",
                                                   "Director","Manager",
                                                   "Trader","Employee","Lawyer","N/A"))

# load role probability by MMSB from files
pi_MMSB <- as.matrix(read.csv("MMSB/Enron_email_VI_pi_N148C2.csv", header = TRUE))
B_MMSB <- as.matrix(read.csv("MMSB/Enron_email_VI_B_N148C2.csv", header =FALSE))


dat_emp <- data.frame(Name = emp.list$folder,
                      Pi = pi.est[1,],
                      Pi_MMSB = pi_MMSB[,1],
                      Position = pos_factor,
                      NumEmail = num_reply_employee)

pos_interest <- c("CEO", "President","Managing Director",
                  "Vice President","Director","Manager",
                  "Lawyer","Trader","Employee")

pdf(paste0("Images/leadership_analysis_K",K_opt,"_v",ver,"_r1.pdf"),width = 12, height = 8)
p <- ggplot(data=dat_emp[dat_emp$Position %in% pos_interest,], mapping=aes(x=Pi, y=3 * log(1+NumEmail), col = Position, size = log(1+NumEmail))) +
  geom_point() + scale_color_manual(values = c(brewer.pal(8,"Spectral"),"black")) + 
  xlim(0,1) +
  xlab("The first membership proportion") + ylab(expression(log(1+"NumEmail"))) + theme_bw() +
  guides(col = guide_legend(override.aes = list(size = 5), order = 1),  size = guide_legend(order = 2)) + 
  theme(text = element_text(size=36),
        axis.text.x =element_text(#face = "bold", #color = "#993333", 
          size = 32,angle = 90), 
        axis.text.y = element_text(#face = "bold", #color = "#993333", 
          size = 32),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2))
p
dev.off()

pdf(paste0("Images/leadership_analysis_MMSB.pdf"),width = 12, height = 8)
p <- ggplot(data=dat_emp[dat_emp$Position %in% pos_interest,], mapping=aes(x=Pi_MMSB, y=3 * log(1+NumEmail), col = Position, size = log(1+NumEmail))) +
  geom_point() + scale_color_manual(values = c(brewer.pal(8,"Spectral"),"black")) + 
  xlim(0,1) +
  xlab("The first membership proportion") + ylab(expression(log(1+"NumEmail"))) + theme_bw() +
  guides(col = guide_legend(override.aes = list(size = 5), order = 1),  size = guide_legend(order = 2)) + 
  theme(text = element_text(size=36),
        axis.text.x =element_text(#face = "bold", #color = "#993333", 
          size = 32,angle = 90), 
        axis.text.y = element_text(#face = "bold", #color = "#993333", 
          size = 32),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2))
p
dev.off()

# return the membership of CEOs
CEO_emp <- which(emp.list$position=="CEO")
n_CEO <- length(CEO_emp)
CEO_membership <- matrix(NA,n_CEO,2 * K_opt + 2)
CEO_membership[,1] <- emp.list$eid[CEO_emp]
CEO_membership[,2] <- as.character(emp.list$folder)[CEO_emp]
CEO_membership[,2 + 1:K_opt] <- round(t(pi.iter[,CEO_emp]),digits = 3)
CEO_membership[,2 + K_opt + 1:K_opt] <- round(pi_MMSB[CEO_emp,],digits = 3)
colnames(CEO_membership) <- c("Employee ID","User Name",paste0("SMMB Role ",1:K_opt),paste0("MMSB Role ",1:K_opt))
rownames(CEO_membership) <- NULL

print("The role porportions of four CEO under SMMB and MMSB is given by")
print(CEO_membership)

##########
# Lambda #
##########
SurCure <- function(time, lam, kn, the){
  n_t <- length(time)
  kn0 <- c(0,kn)
  res <- rep(NA, n_t)
  for(ind_t in 1:n_t){
    int_t <- sum(time[ind_t] > kn0)
    if(int_t > length(kn)){
      cum_hazard <- lam[int_t - 1] * (time[ind_t] - kn0[int_t - 1]) 
      int_t <- int_t - 1
    }else{
      cum_hazard <- lam[int_t] * (time[ind_t] - kn0[int_t]) 
    }
    if(int_t > 1){
      dif_kn0 <- diff(kn0)
      cum_hazard <- cum_hazard + sum(lam[1:(int_t-1)] * dif_kn0[1:(int_t-1)])
    }
    res[ind_t] <- exp(- the * (1 - exp(-cum_hazard))) 
  }
  return(res)
}

# Calculate the survival probability at each knot
Sur_prob_mem <- matrix(NA, P, 2 + K_opt^2)
Sur_prob_mem[,1] <- 1:P
Sur_prob_mem[,2] <- round(knots,digits = 4)

for(l in 1:K_opt){
  for(m in 1:K_opt){
    if(proj == "simulation" | proj == "simknots"){
      ind_mem <- (switch[l]-1) * K_opt + switch[m]
    }else{
      ind_mem <- (l - 1) * K_opt + m
    }
    # Forward = FALSE, only 1 receiver
    Sur_prob_mem[,ind_mem + 2] <- cal_sur_cure(lambda.est[(l - 1) * K_opt + m, ],diff.knots, exp(beta.est[(l - 1) * K_opt + m, 1]))
  }
}

colnames(Sur_prob_mem) <- c("Knot #", "Knot (Hours)", paste0("(",rep(1:K_opt, each = K_opt),",",rep(1:K_opt, K_opt),")"))

# draw the survival curve
color_memb <- brewer.pal(K_opt * K_opt, "Set1")
color_name <- NULL
seniority<- c("junior","senior")
for(l in 1:K_opt){
  for(m in 1:K_opt){
    color_name <- c(color_name, paste0("(",l,",",m,"): Role ",m," (",seniority[m],") replies to the emails from role ",l," (",seniority[l],")")) 
  }
}

t_seq = seq(0.1, knots[P], by = 0.1)

sur_function_mem <- matrix(NA, length(t_seq), K_opt * K_opt)
for(l in 1:K_opt){
  for(m in 1:K_opt){
    if(proj == "simulation" | proj == "simknots"){
      ind_mem <- (switch[l]-1) * K_opt + switch[m]
    }else{
      ind_mem <- (l - 1) * K_opt + m
    }
    
    sur_function_mem[,ind_mem] <- SurCure(t_seq, lambda.est[(l-1) * K_opt + m,], knots, exp(beta.est[(l-1) * K_opt + m,1]))
  }
}
dat_sur <- data.frame(Time = rep(t_seq, K_opt * K_opt),
                      SurFunc = as.vector(sur_function_mem),
                      MemPair = rep(color_name, each = length(t_seq)))


pdf(paste0("Images/",proj,"_sur_prob_X0_est_K",K_opt,"_v",ver,".pdf"),width = 18, height = 8)
p <- ggplot(data=dat_sur, mapping=aes(x=Time, y=SurFunc, group = MemPair)) +
  geom_line(aes(linetype=MemPair,col = MemPair), size = 2) + scale_color_manual(values=c("#2b83ba", "#d7191c", "#abdda4", "#fdae61")) +
  scale_linetype_manual(values=c("solid", "longdash", "dotdash","dotted")) +
  xlab("Time") + ylab("Survival probability") + ylim(0,1) + theme_bw() +
  labs(col="Role Pair", linetype="Role Pair") +
  # guides(col = guide_legend(override.aes = list(size = 5))) + 
  theme(text = element_text(size=28),
        axis.text.x =element_text(#face = "bold", #color = "#993333", 
          size = 24,angle = 90), 
        axis.text.y = element_text(#face = "bold", #color = "#993333", 
          size = 24),
        legend.key.size = unit(6,"line"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2))
p
dev.off()

########
# Beta #
########
sum_betasq <- matrix(0, K_opt * K_opt, p.cov + 1)
for(iter in burnin + 1:num_after_iter){
  sum_betasq <- sum_betasq + beta_record[[iter]]^2
}

beta.mat <- matrix(NA, p.cov + 1, K_opt^2)
beta.sd <- matrix(NA, p.cov + 1, K_opt^2)
for(l in 1:K_opt){
  for(m in 1:K_opt){
    if(proj == "simulation" | proj == "simknots"){
      ind_mem <- (switch[l]-1) * K_opt + switch[m]
    }else{
      ind_mem <- (l - 1) * K_opt + m
    }
    beta.mat[,ind_mem] <- beta.est[(l-1)*K_opt+m,]
    beta.sd[,ind_mem] <- sqrt((sum_betasq[(l-1)*K_opt+m,] - num_after_iter * beta.est[(l-1)*K_opt+m, ]^2)/(num_after_iter - 1))
  }
}
colnames(beta.mat) <- paste0("(",rep(1:K_opt, each = K_opt),",",rep(1:K_opt, K_opt),")")
colnames(beta.sd) <- paste0("(",rep(1:K_opt, each = K_opt),",",rep(1:K_opt, K_opt),")")

# credible interval of beta
# convergence
beta.post <- array(NA, dim = c(K_opt * K_opt, p.cov + 1, num_iter))
for(iter in 1:num_iter){
  beta.post[,,iter] <- beta_record[[iter]]
}
iter_after_burin <- burnin + 1:num_after_iter
beta.interval.ETI <- array(NA, dim = c(K_opt * K_opt, p.cov + 1, 2))

for(l in 1:K_opt){
  for(m in 1:K_opt){
    if(proj == "simulation" | proj == "simknots"){
      ind_mem <- (switch[l]-1) * K_opt + switch[m]
    }else{
      ind_mem <- (l - 1) * K_opt + m
    }
    for(r in 1:(p.cov + 1)){
      cre_int_ETI <- ci(beta.post[(l-1)*K_opt+m,r,iter_after_burin], ci = 0.95, method = "ETI")
      beta.interval.ETI[ind_mem, r, 1] <- cre_int_ETI$CI_low
      beta.interval.ETI[ind_mem, r, 2] <- cre_int_ETI$CI_high
    }
  }
}

beta_out <- matrix(NA, p.cov * K_opt * K_opt, 4)
ind_row <- 1
for(r in 2:(p.cov+1)){
  for(l in 1:K_opt){
    for(m in 1:K_opt){
      ind_mem <- (l - 1) * K_opt + m
      beta_out[ind_row, 1] <- round(beta.mat[r, ind_mem],4)
      beta_out[ind_row, 2] <- round(beta.sd[r, ind_mem],4)
      beta_out[ind_row, 3] <- paste0("[",round(beta.interval.ETI[ind_mem, r, 1],4),", ",
                                     round(beta.interval.ETI[ind_mem, r, 2],4),"]")
      ind_row <- ind_row + 1 
    }
  }
}

rownames(beta_out) <- paste0("$\\beta^{",rep(rep(1:K_opt, each = K_opt),p.cov),
                             rep(1:K_opt, K_opt * (p.cov)),"}_",
                             rep(2:(p.cov+1), each = K_opt^2),"$")
colnames(beta_out) <- c("Posterior Mean", "Posterior SD","95% CI")
print(xtable(beta_out),sanitize.text.function=function(x){x})

#########################################
# Analyze with confidential information #
#########################################
ver <- 2
K_opt <- 2
load(paste0("Inference/Inference_",proj,"_SMMB_K",K_opt,"_v",ver,"_r1.RData"))

switch <- c(2,1)

######
# Pi #
######
# Analyze the relationship bewteen membership probabilities and positions of all employee
emp.list <- read.csv("RawData/employeelist.csv")

num_reply_employee <- rep(0,N)
for(i in 1:tol_pair){
  sender <- num.pair[i,1]
  receiver <- num.pair[i,2]
  nij <- num.pair[i,3]
  
  num_reply_employee[sender] <- num_reply_employee[sender] + nij
  num_reply_employee[receiver] <- num_reply_employee[receiver] + nij
}

# reorder the position from high to low
pos_factor <- factor(emp.list$position, levels = c("CEO", "President","Vice President","Managing Director",
                                                   "Director","Manager",
                                                   "Trader","Employee","Lawyer","N/A"))


dat_emp <- data.frame(Name = emp.list$folder,
                      Pi = pi.est[switch[1],],
                      Position = pos_factor,
                      NumEmail = num_reply_employee)

pos_interest <- c("CEO", "President","Managing Director",
                  "Vice President","Director","Manager",
                  "Lawyer","Trader","Employee")

pdf(paste0("Images/leadership_analysis_K",K_opt,"_v",ver,"_r1.pdf"),width = 12, height = 8)
# jpeg(paste0("Images/leadership_analysis_K",K_opt,"_v",v,"_r",r,".jpg"),width = 1440, height = 1080)
p <- ggplot(data=dat_emp[dat_emp$Position %in% pos_interest,], mapping=aes(x=Pi, y=3 * log(1+NumEmail), col = Position, size = log(1+NumEmail))) +
  geom_point() + scale_color_manual(values = c(brewer.pal(8,"Spectral"),"black")) + 
  xlim(0,1) + 
  xlab("The first membership proportion") + ylab(expression(log(1+"NumEmail"))) + theme_bw() +
  guides(col = guide_legend(override.aes = list(size = 5), order = 1),  size = guide_legend(order = 2)) + 
  theme(text = element_text(size=36),
        axis.text.x =element_text(#face = "bold", #color = "#993333", 
          size = 32,angle = 90), 
        axis.text.y = element_text(#face = "bold", #color = "#993333", 
          size = 32),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2))
p
dev.off()

# return the membership of CEOs
CEO_emp <- which(emp.list$position=="CEO")
n_CEO <- length(CEO_emp)
CEO_membership <- matrix(NA,n_CEO,K_opt + 2)
CEO_membership[,1] <- emp.list$eid[CEO_emp]
CEO_membership[,2] <- as.character(emp.list$folder)[CEO_emp]
CEO_membership[,2 + 1:K_opt] <- round(t(pi.est[switch[1:2],CEO_emp]),digits = 4)
# CEO_membership[,2 + K_opt + 1:K_opt] <- round(pi_MMSB[CEO_emp,],digits = 4)
colnames(CEO_membership) <- c("Employee ID","User Name",paste0("Role ",1:K_opt))
rownames(CEO_membership) <- NULL

print("The role porportions of four CEOs under SMMB and MMSB is given by")
print(CEO_membership)
write.table(CEO_membership[,1:4],file = "clipboard",sep = "\t",row.names = FALSE)


# Draw boxplot for each position
pdf(paste0("Images/boxplot_leadership_K",K_opt,"_v",v,"_r",r,".pdf"),width = 12, height = 12)
p <- ggplot(data=dat_emp[dat_emp$Position %in% pos_interest,], mapping=aes(x=Position, y=Pi, fill = Position)) +
  geom_boxplot() + scale_fill_manual(values = c(brewer.pal(8,"Spectral"),"black")) + 
  ylim(0,1) + 
  xlab("Position") + ylab("The first membership proportion") + theme_bw() +
  guides(fill = "none") + 
  theme(text = element_text(size=36),
        axis.text.x =element_text(#face = "bold", #color = "#993333", 
          size = 32,angle = 90), 
        axis.text.y = element_text(#face = "bold", #color = "#993333", 
          size = 32),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2))
p
dev.off()

##########
# Lambda #
##########
SurCure <- function(time, lam, kn, the){
  n_t <- length(time)
  kn0 <- c(0,kn)
  res <- rep(NA, n_t)
  for(ind_t in 1:n_t){
    int_t <- sum(time[ind_t] > kn0)
    if(int_t > length(kn)){
      cum_hazard <- lam[int_t - 1] * (time[ind_t] - kn0[int_t - 1]) 
      int_t <- int_t - 1
    }else{
      cum_hazard <- lam[int_t] * (time[ind_t] - kn0[int_t]) 
    }
    if(int_t > 1){
      dif_kn0 <- diff(kn0)
      cum_hazard <- cum_hazard + sum(lam[1:(int_t-1)] * dif_kn0[1:(int_t-1)])
    }
    res[ind_t] <- exp(- the * (1 - exp(-cum_hazard))) 
  }
  return(res)
}

# Calculate the survival probability at each knot
Sur_prob_mem <- matrix(NA, M, 2 + K_opt^2)
Sur_prob_mem[,1] <- 1:M
Sur_prob_mem[,2] <- round(knots,digits = 4)

for(l in 1:K_opt){
  for(m in 1:K_opt){
    if(proj == "simulation" | proj == "simknots"){
      ind_mem <- (switch[l]-1) * K_opt + switch[m]
    }else{
      ind_mem <- (l - 1) * K_opt + m
    }
    # Forward = FALSE, only 1 receiver
    Sur_prob_mem[,ind_mem + 2] <- cal_sur_cure(lambda.est[(l - 1) * K_opt + m, ],diff.knots, exp(beta.est[(l - 1) * K_opt + m, 1]))
  }
}

colnames(Sur_prob_mem) <- c("Knot #", "Knot (Hours)", paste0("(",rep(1:K_opt, each = K_opt),",",rep(1:K_opt, K_opt),")"))

# draw the survival curve
color_memb <- brewer.pal(K_opt * K_opt, "Set1")
color_name <- NULL
for(l in 1:K_opt){
  for(m in 1:K_opt){
    color_name <- c(color_name, paste0("(",l,",",m,")")) 
  }
}

t_seq = seq(0.1, knots[M], by = 0.1)

sur_function_mem <- matrix(NA, length(t_seq), K_opt * K_opt)
for(l in 1:K_opt){
  for(m in 1:K_opt){

    ind_mem <- (switch[l] - 1) * K_opt + switch[m]

    
    sur_function_mem[,ind_mem] <- SurCure(t_seq, lambda.est[(l-1) * K_opt + m,], knots, exp(beta.est[(l-1) * K_opt + m,1]))
  }
}
dat_sur <- data.frame(Time = rep(t_seq, K_opt * K_opt),
                      SurFunc = as.vector(sur_function_mem),
                      MemPair = rep(color_name, each = length(t_seq)))

legend_label <- paste0("(",rep(1:2,each = 2), ",",rep(1:2, 2),"): Role " , c(1,2,1,2),"(", c("junior","senior","junior","senior") ,
                       ") replies to the emails from role ", c(1,1,2,2),"(",c("junior","junior","senior","senior"),")")

pdf(paste0("Images/",proj,"_sur_prob_X0_est_K",K_opt,"_v",ver,".pdf"),width = 18, height = 8)
p <- ggplot(data=dat_sur, mapping=aes(x=Time, y=SurFunc, group = MemPair)) +
  geom_line(aes(linetype=MemPair,col = MemPair), size = 1.5) + 
  scale_colour_manual(values = brewer.pal(4, "Spectral")[c(4,1,3,2)], labels = legend_label) +
  labs(group = "Membership pair") +
  scale_linetype_manual(values=c("solid", "longdash", "dotdash","dotted"),labels = legend_label) +
  xlab("Time (hours)") + ylab("Survival probability") + ylim(0,1) + theme_bw() + 
  # guides(col = guide_legend(override.aes = list(size = 5))) + 
  theme(text = element_text(size=28),
        axis.text.x =element_text(#face = "bold", #color = "#993333", 
          size = 24,angle = 90), 
        axis.text.y = element_text(#face = "bold", #color = "#993333", 
          size = 24),
        legend.key.size = unit(6,"line"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2))
p
dev.off()

########
# Beta #
########
beta_out <- matrix(NA, (p.cov - 1) * K_opt * K_opt, 3)
ind_row <- 1
for(r in 2:p.cov){
  for(l in 1:K_opt){
    for(m in 1:K_opt){
      ind_mem <- (switch[l]-1) * K_opt + switch[m]
      beta_out[ind_row, 1] <- round(beta.est[ind_mem, r],4)
      beta_out[ind_row, 2] <- round(beta.sd[ind_mem, r],4)
      beta_out[ind_row, 3] <- paste0("[",round(beta.ci[ind_mem, r, 1],4),", ",
                                     round(beta.ci[ind_mem, r, 2],4),"]")
      ind_row <- ind_row + 1 
    }
  }
}

rownames(beta_out) <- paste0("$\\beta^{",rep(rep(1:K_opt, each = K_opt),p.cov-1),
                             rep(1:K_opt, K_opt * (p.cov - 1)),"}_",
                             rep(2:p.cov, each = K_opt^2),"$")
colnames(beta_out) <- c("Posterior Mean", "Posterior SD","95% CI")
print(xtable(beta_out),sanitize.text.function=function(x){x})

#########################################
# Analyze with confidential information #
#########################################
ver <- 3
K_opt <- 2
load(paste0("Inference/Inference_",proj,"_SMMB_K",K_opt,"_v",ver,"_r1.RData"))

switch <- c(1,2)

######
# Pi #
######
# Analyze the relationship bewteen membership probabilities and positions of all employee
emp.list <- read.csv("RawData/employeelist.csv")

num_reply_employee <- rep(0,N)
for(i in 1:tol_pair){
  sender <- num.pair[i,1]
  receiver <- num.pair[i,2]
  nij <- num.pair[i,3]
  
  num_reply_employee[sender] <- num_reply_employee[sender] + nij
  num_reply_employee[receiver] <- num_reply_employee[receiver] + nij
}

# reorder the position from high to low
pos_factor <- factor(emp.list$position, levels = c("CEO", "President","Vice President","Managing Director",
                                                   "Director","Manager",
                                                   "Trader","Employee","Lawyer","N/A"))


dat_emp <- data.frame(Name = emp.list$folder,
                      Pi = pi.est[switch[1],],
                      Position = pos_factor,
                      NumEmail = num_reply_employee)

pos_interest <- c("CEO", "President","Managing Director",
                  "Vice President","Director","Manager",
                  "Lawyer","Trader","Employee")

pdf(paste0("Images/leadership_analysis_K",K_opt,"_v",ver,".pdf"),width = 12, height = 8)
# jpeg(paste0("Images/leadership_analysis_K",K_opt,"_v",v,"_r",r,".jpg"),width = 1440, height = 1080)
p <- ggplot(data=dat_emp[dat_emp$Position %in% pos_interest,], mapping=aes(x=Pi, y=3 * log(1+NumEmail), col = Position, size = log(1+NumEmail))) +
  geom_point() + scale_color_manual(values = c(brewer.pal(8,"Spectral"),"black")) + 
  xlim(0,1) + 
  xlab("The first membership proportion") + ylab(expression(log(1+"NumEmail"))) + theme_bw() +
  guides(col = guide_legend(override.aes = list(size = 5), order = 1),  size = guide_legend(order = 2)) + 
  theme(text = element_text(size=36),
        axis.text.x =element_text(#face = "bold", #color = "#993333", 
          size = 32,angle = 90), 
        axis.text.y = element_text(#face = "bold", #color = "#993333", 
          size = 32),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2))
p
dev.off()

# return the membership of CEOs
CEO_emp <- which(emp.list$position=="CEO")
n_CEO <- length(CEO_emp)
CEO_membership <- matrix(NA,n_CEO,K_opt + 2)
CEO_membership[,1] <- emp.list$eid[CEO_emp]
CEO_membership[,2] <- as.character(emp.list$folder)[CEO_emp]
CEO_membership[,2 + 1:K_opt] <- round(t(pi.est[switch[1:2],CEO_emp]),digits = 4)
# CEO_membership[,2 + K_opt + 1:K_opt] <- round(pi_MMSB[CEO_emp,],digits = 4)
colnames(CEO_membership) <- c("Employee ID","User Name",paste0("Role ",1:K_opt))
rownames(CEO_membership) <- NULL

print("The role porportions of four CEOs under SMMB and MMSB is given by")
print(CEO_membership)
write.table(CEO_membership[,1:4],file = "clipboard",sep = "\t",row.names = FALSE)


# Draw boxplot for each position
pdf(paste0("Images/boxplot_leadership_K",K_opt,"_v",ver,".pdf"),width = 12, height = 12)
p <- ggplot(data=dat_emp[dat_emp$Position %in% pos_interest,], mapping=aes(x=Position, y=Pi, fill = Position)) +
  geom_boxplot() + scale_fill_manual(values = c(brewer.pal(8,"Spectral"),"black")) + 
  ylim(0,1) + 
  xlab("Position") + ylab("The first membership proportion") + theme_bw() +
  guides(fill = "none") + 
  theme(text = element_text(size=36),
        axis.text.x =element_text(#face = "bold", #color = "#993333", 
          size = 32,angle = 90), 
        axis.text.y = element_text(#face = "bold", #color = "#993333", 
          size = 32),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2))
p
dev.off()

##########
# Lambda #
##########
SurCure <- function(time, lam, kn, the){
  n_t <- length(time)
  kn0 <- c(0,kn)
  res <- rep(NA, n_t)
  for(ind_t in 1:n_t){
    int_t <- sum(time[ind_t] > kn0)
    if(int_t > length(kn)){
      cum_hazard <- lam[int_t - 1] * (time[ind_t] - kn0[int_t - 1]) 
      int_t <- int_t - 1
    }else{
      cum_hazard <- lam[int_t] * (time[ind_t] - kn0[int_t]) 
    }
    if(int_t > 1){
      dif_kn0 <- diff(kn0)
      cum_hazard <- cum_hazard + sum(lam[1:(int_t-1)] * dif_kn0[1:(int_t-1)])
    }
    res[ind_t] <- exp(- the * (1 - exp(-cum_hazard))) 
  }
  return(res)
}

# Calculate the survival probability at each knot
Sur_prob_mem <- matrix(NA, M, 2 + K_opt^2)
Sur_prob_mem[,1] <- 1:M
Sur_prob_mem[,2] <- round(knots,digits = 4)

for(l in 1:K_opt){
  for(m in 1:K_opt){

    ind_mem <- (switch[l]-1) * K_opt + switch[m]

    # Forward = FALSE, only 1 receiver
    Sur_prob_mem[,ind_mem + 2] <- cal_sur_cure(lambda.est[(l - 1) * K_opt + m, ],diff.knots, exp(beta.est[(l - 1) * K_opt + m, 1]))
  }
}

colnames(Sur_prob_mem) <- c("Knot #", "Knot (Hours)", paste0("(",rep(1:K_opt, each = K_opt),",",rep(1:K_opt, K_opt),")"))

# draw the survival curve
color_memb <- brewer.pal(K_opt * K_opt, "Set1")
color_name <- NULL
for(l in 1:K_opt){
  for(m in 1:K_opt){
    color_name <- c(color_name, paste0("(",l,",",m,")")) 
  }
}

t_seq = seq(0.1, knots[M], by = 0.1)

sur_function_mem <- matrix(NA, length(t_seq), K_opt * K_opt)
for(l in 1:K_opt){
  for(m in 1:K_opt){
    
    ind_mem <- (switch[l] - 1) * K_opt + switch[m]
    
    
    sur_function_mem[,ind_mem] <- SurCure(t_seq, lambda.est[(l-1) * K_opt + m,], knots, exp(beta.est[(l-1) * K_opt + m,1]))
  }
}
dat_sur <- data.frame(Time = rep(t_seq, K_opt * K_opt),
                      SurFunc = as.vector(sur_function_mem),
                      MemPair = rep(color_name, each = length(t_seq)))

legend_label <- paste0("(",rep(1:2,each = 2), ",",rep(1:2, 2),"): Role " , c(1,2,1,2),"(", c("junior","senior","junior","senior") ,
                       ") replies to the emails from role ", c(1,1,2,2),"(",c("junior","junior","senior","senior"),")")

pdf(paste0("Images/",proj,"_sur_prob_X0_est_K",K_opt,"_v",ver,".pdf"),width = 18, height = 8)
p <- ggplot(data=dat_sur, mapping=aes(x=Time, y=SurFunc, group = MemPair)) +
  geom_line(aes(linetype=MemPair,col = MemPair), size = 1.5) + 
  scale_colour_manual(values = brewer.pal(4, "Spectral")[c(4,1,3,2)], labels = legend_label) +
  labs(group = "Membership pair") +
  scale_linetype_manual(values=c("solid", "longdash", "dotdash","dotted"),labels = legend_label) +
  xlab("Time (hours)") + ylab("Survival probability") + ylim(0,1) + theme_bw() + 
  # guides(col = guide_legend(override.aes = list(size = 5))) + 
  theme(text = element_text(size=28),
        axis.text.x =element_text(#face = "bold", #color = "#993333", 
          size = 24,angle = 90), 
        axis.text.y = element_text(#face = "bold", #color = "#993333", 
          size = 24),
        legend.key.size = unit(6,"line"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2))
p
dev.off()

########
# Beta #
########
beta_out <- matrix(NA, (p.cov - 1) * K_opt * K_opt, 3)
ind_row <- 1
for(r in 2:p.cov){
  for(l in 1:K_opt){
    for(m in 1:K_opt){
      ind_mem <- (switch[l]-1) * K_opt + switch[m]
      beta_out[ind_row, 1] <- round(beta.est[ind_mem, r],4)
      beta_out[ind_row, 2] <- round(beta.sd[ind_mem, r],4)
      beta_out[ind_row, 3] <- paste0("[",round(beta.ci[ind_mem, r, 1],4),", ",
                                     round(beta.ci[ind_mem, r, 2],4),"]")
      ind_row <- ind_row + 1 
    }
  }
}

rownames(beta_out) <- paste0("$\\beta^{",rep(rep(1:K_opt, each = K_opt),p.cov-1),
                             rep(1:K_opt, K_opt * (p.cov - 1)),"}_",
                             rep(2:p.cov, each = K_opt^2),"$")
colnames(beta_out) <- c("Posterior Mean", "Posterior SD","95% CI")
print(xtable(beta_out),sanitize.text.function=function(x){x})
