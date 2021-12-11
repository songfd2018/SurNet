rm(list=ls())
library(survival)
load("Reponse_time_collection.RData")

#########################
# summarize the dataset #
#########################
# sort reply emails by the user order and include pairs without failure time
tol_rep <- nrow(sur_collect)
p.cov <- ncol(sur_collect) - 6
usr_no_list <- as.numeric(names(table(sur_collect[,1:2])))
N <-length(usr_no_list)

# order the response time by senders and receivers
Dat <- sur_collect
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

num.pair <- cbind(num.pair, 0)
index_rep <- 0
for(i in 1:tol_pair){
  nij <- num.pair[i,3]
  num.pair[i,4] <- sum(Dat[index_rep + 1:nij, 4])
  index_rep <- index_rep + nij
}

# Status:
print(paste0("The number of response time is ", sum(Dat[,4]),"."))
print(paste0("The number of censoring time is ", nrow(Dat) - sum(Dat[,4]),"."))


# covariate 1: Weekend effect
print(paste0("The number of emails sent on weekends is ", sum(Dat[,6]),"."))
print(paste0("The number of emails sent on weekdays is ", nrow(Dat) - sum(Dat[,6]),"."))

# covariate 2: Forwarding emails
print(paste0("The number of forwarding emails is ", sum(Dat[,7]),"."))
print(paste0("The number of non-forwarding emails is ", nrow(Dat) - sum(Dat[,7]),"."))

# covariate 3: Log-scale recevier number
print(paste0("The mean of log-scale reivecer numbers is ", mean(Dat[,8]),"."))
print(paste0("The standard deviation of log-scale reivecer numbers is ", sd(Dat[,8]),"."))

# covariate 4: Log-scale word count
print(paste0("The mean of log-scale word counts is ", mean(Dat[,9]),"."))
print(paste0("The standard deviation of log-scale word counts is ", sd(Dat[,9]),"."))

#####################################
# draw heatmap after ordering users #
#####################################
library(igraph)
employees<-read.table("RawData/employeelist.txt",stringsAsFactors = FALSE, header = T)
emp.list <- employees[usr_no_list,]

# order positions for coloring
emp.list$status_lab <- factor(emp.list$status, levels = c("CEO", "President","Trader","Vice President",
                                                          "Managing Director","Director","Manager",
                                                          "In House Lawyer","Employee","N/A"))

# match the eid after exclude a user
emp.list$eid <- 1:N

num.pair[,1:2] <- factor(num.pair[,1:2], levels = usr_no_list)
g <- graph_from_data_frame(d=num.pair, vertices=emp.list, directed=TRUE)

mem <- cluster_walktrap(g)

# order by cluster size
re_order_vertice <- NULL
mem_index <- membership(mem)
cluster_size <- table(mem_index)
cluster_order <- names(sort(cluster_size))

for(cl in 1:length(cluster_size)){
  re_order_vertice <- c(re_order_vertice, which(mem_index == cluster_order[cl]))
}

# draw the adjacent matrix
dat_numpair <- data.frame(Sender = factor(num.pair[,1],levels = re_order_vertice),
                          Receiver = factor(num.pair[,2], levels = re_order_vertice),
                          Weights = num.pair[,3])


if(!dir.exists("Images")){
  dir.create("Images")
}
pdf(paste0("Images/Adjacent_matrix_enron_grey_v1.pdf"),width = 6, height = 6)
p_grey <- ggplot(dat_numpair, aes(x = Sender, y = Receiver)) +
  geom_tile() +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size=32),
        axis.ticks = element_blank())
p_grey
dev.off()
