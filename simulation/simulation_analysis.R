rm(list=ls())
library(RColorBrewer)
library(ggplot2)
library(mclust)
library(bayestestR)
proj <- "simulation"
v <- 1
K.sel <- 1:5

if(!dir.exists("Images")){
  dir.create("Images")
}

# load the true values of parameters
if(proj == "simulation" | proj == "simknots"){
  load(paste0("ObsData/workspace_SMMB_v",v,".RData"))
  mp.syn <- (l.syn - 1) * K + m.syn
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
      message(paste0("The inference of 1-th replication on dataset \"",proj,"\" with K=",K," does not exist!"))
    }
}

# record the DIC
names(pD.record) <- paste0("K=",K.sel)

names(DIC.record) <- paste0("K=",K.sel)

# draw the chart plot of DIC
dat_dic <- data.frame(K = 1:5,
                      DIC = DIC.record)

pdf(paste0("Images/",proj,"_DIC_v",v,".pdf"),width = 8, height = 6)
ggplot(dat_dic, aes(x=K)) + 
  geom_line(aes(y=DIC.record),size = 1.5) + 
  ylab("DIC") +
  theme_bw() +
  theme(axis.text = element_text(size = 24),  # rotate x axis text
        axis.title = element_text(size = 32))  # turn off minor grid
dev.off()

############################
# Analyze the optimal case #
############################
K_opt <- 3
load(paste0("Inference/Inference_",proj,"_SMMB_K",K_opt,"_v",v,"_r1.RData"))

if(proj == "simulation" | proj == "simknots"){
  # ad-hoc switch roles
  switch <- c(2,1,3)
}

######
# Pi #
######
if(proj == "simulation" | proj == "simknots"){
  l.est <- (mp.est-1) %/% 3 + 1
  m.est <- (mp.est-1) %% 3 + 1
  mp.syn <- (l.syn - 1) * K_opt + m.syn
  
  mp.switch <- (switch[l.est] - 1) * K_opt + switch[m.est]
  
  tab_latex <- table(mp.syn,mp.switch)
  rownames(tab_latex) <- paste0("(",rep(1:K_opt, each = K_opt),",",rep(1:K_opt, K_opt),")")
  colnames(tab_latex) <- paste0("(",rep(1:K_opt, each = K_opt),",",rep(1:K_opt, K_opt),")")
  
  print(tab_latex)
  
}

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

if(proj == "simulation" | proj == "simknots"){
  
  # Number of time knots
  knots.syn<-c(6,18,30,48,90)
  diff.knots.syn <-diff(c(0,knots.syn))
  P.syn <- 5
  
  Sur_prob_syn <- matrix(NA, P, 2 + K_opt^2)
  Sur_prob_syn[,1] <- 1:P
  Sur_prob_syn[,2] <- round(knots,digits = 4)
  
  
  
  for(l in 1:K_opt){
    for(m in 1:K_opt){
      ind_mem <- (l-1) * K_opt + m
      # Forward = FALSE, only 1 receiver
      Sur_prob_syn[,ind_mem + 2] <- SurCure(knots,lambda[l,m,],knots.syn, exp(beta[l, m, 1]))
    }
  }
  colnames(Sur_prob_syn) <- c("Knot #", "Knot (Hours)", paste0("(",rep(1:K_opt, each = K_opt),",",rep(1:K_opt, K_opt),")"))
}

# draw the survival curve
color_memb <- brewer.pal(K_opt * K_opt, "Set1")
color_name <- NULL
for(l in 1:K_opt){
  for(m in 1:K_opt){
    color_name <- c(color_name, paste0("(",l,",",m,")")) 
  }
}

# Draw the heatmap to compare the survival probability at each knot
library("grid")
library("gridExtra")
library(gtable)
# library("pdp")
index.l<-c(1,1,1,2,2,2,3,3,3)
index.m<-c(1,2,3,1,2,3,1,2,3)

p_syn <- list()
p_est <- list()
for(p in 1:P){
  data.syn<-data.frame(l=index.l,
                       m=index.m,
                       prob=Sur_prob_syn[p, 2 + 1:(K_opt*K_opt)])
  p_syn[[p]] <- ggplot(data.syn,aes(x=m, y=l, fill=prob)) + geom_tile(color="white", size=0.1) +
    ylab('Sender')+xlab("Receiver") +
    scale_fill_gradientn(colours = rev(heat.colors(10)),limits = c(0,1))+
    theme_bw() +
    #scale_x_discrete(limits=c(1,2,3))+
    geom_text(aes(label=format(round(prob,3),nsmall = 3,scientific = FALSE)), angle=45)+
    theme(panel.grid =element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          axis.text=element_text(size=18),
          axis.title=element_text(size=24))
  
  data.est<-data.frame(l=index.l,
                       m=index.m,
                       prob=Sur_prob_mem[p, 2 + 1:(K_opt*K_opt)])
  p_est[[p]] <- ggplot(data.est, aes(x=m, y=l, fill=prob))+
    ylab('Sender')+xlab("Receiver")+
    geom_tile(color="white", size=0.1)+
    scale_fill_gradientn(colours = rev(heat.colors(10)),limits = c(0,1))+
    guides(fill=FALSE)+
    theme_bw() +
    #scale_x_discrete(limits=c(1,2,3))+
    geom_text(aes(label=format(round(prob,3),nsmall = 3,scientific = FALSE)), angle=45)+
    theme(panel.grid =element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          axis.text=element_text(size=18),
          axis.title=element_text(size=24))
}

legend <- gtable_filter(ggplot_gtable(ggplot_build(
  p_syn[[1]] + labs(fill = "Prob") +
    theme(legend.position="right",
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 24)))), "guide-box")

pdf(paste0("Images/",proj,"_baseline_sur_prob_comparison.pdf"),width = 14.4, height = 6.3)
label_size <- 28
p_label <- list()
p_label[[1]] <- arrangeGrob(p_syn[[1]] + theme(legend.position="none"),
                            top = textGrob("a", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=label_size, fontfamily = "serif", fontface = "bold")))

p_label[[2]] <- arrangeGrob(p_syn[[2]] + theme(legend.position="none"),
                            top = textGrob("b", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=label_size, fontfamily = "serif", fontface = "bold")))

p_label[[3]] <- arrangeGrob(p_syn[[3]] + theme(legend.position="none"),
                            top = textGrob("c", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=label_size, fontfamily = "serif", fontface = "bold")))

p_label[[4]] <- arrangeGrob(p_syn[[4]] + theme(legend.position="none"),
                            top = textGrob("d", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=label_size, fontfamily = "serif", fontface = "bold")))

p_label[[5]] <- arrangeGrob(p_syn[[5]] + theme(legend.position="none"),
                            top = textGrob("e", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=label_size, fontfamily = "serif", fontface = "bold")))

p_label[[6]] <- arrangeGrob(p_est[[1]] + theme(legend.position="none"),
                            top = textGrob("f", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=label_size, fontfamily = "serif", fontface = "bold")))

p_label[[7]] <- arrangeGrob(p_est[[2]] + theme(legend.position="none"),
                            top = textGrob("g", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=label_size, fontfamily = "serif", fontface = "bold")))

p_label[[8]] <- arrangeGrob(p_est[[3]] + theme(legend.position="none"),
                            top = textGrob("h", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=label_size, fontfamily = "serif", fontface = "bold")))

p_label[[9]] <- arrangeGrob(p_est[[4]] + theme(legend.position="none"),
                            top = textGrob("i", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=label_size, fontfamily = "serif", fontface = "bold")))

p_label[[10]] <- arrangeGrob(p_est[[5]] + theme(legend.position="none"),
                             top = textGrob("j", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
                                            gp=gpar(col="black", fontsize=label_size, fontfamily = "serif", fontface = "bold")))

grid.arrange(arrangeGrob(p_label[[1]],p_label[[2]],p_label[[3]],p_label[[4]],p_label[[5]],
                         p_label[[6]],p_label[[7]],p_label[[8]],p_label[[9]],p_label[[10]],
                         nrow=2),legend,widths=c(12.5,1.5),nrow = 1)
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


if(proj == "simulation" | proj == "simknots"){
  beta.true <- matrix(NA, p.cov + 1, K_opt^2)
  for(l in 1:K_opt){
    for(m in 1:K_opt){
      ind_mem <- (l-1) * K_opt + m
      beta.true[,ind_mem] <- beta[l,m,]
    }
  }
  colnames(beta.true) <- paste0("(",rep(1:K_opt, each = K_opt),",",rep(1:K_opt, K_opt),")")
}

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

beta_latex <- matrix(NA, (p.cov + 1) * K_opt * K_opt, 5)
ind_row <- 1
for(r in 1:(p.cov+1)){
  for(l in 1:K_opt){
    for(m in 1:K_opt){
      ind_mem <- (l - 1) * K_opt + m
      beta_latex[ind_row, 1] <- beta.true[r, ind_mem]
      beta_latex[ind_row, 2] <- beta.mat[r, ind_mem]
      beta_latex[ind_row, 3] <- beta.sd[r, ind_mem]
      beta_latex[ind_row, 4] <- beta.interval.ETI[ind_mem, r, 1]
      beta_latex[ind_row, 5] <- beta.interval.ETI[ind_mem, r, 2]
      ind_row <- ind_row + 1 
    }
  }
}

rownames(beta_latex) <- rep(paste0("(",rep(1:K_opt, each = K_opt),",",rep(1:K_opt, K_opt),")"),p.cov + 1)
colnames(beta_latex) <- c("True Value","Posterior Mean", "Posterior SD","Lower 95% CI", "Upper 95% CI")
print(round(beta_latex,digits = 3))

save.image(paste0("Inference/Analysis_",proj,"_K",K_opt,"_v",v,"_r1.RData"))