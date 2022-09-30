rm(list=ls())
library(ggplot2)

load("RecoveryResults.RData")

if(!dir.exists("Images")){
  dir.create("Images")
}

################################
# vary role number from 2 to 5 #
################################
pdf(paste0("Images/Recovery_performance_by_role.pdf"),width = 8, height = 8)
par_name <- c(expression(lambda),expression(beta),expression(xi))
p <- ggplot(data=Recovery_role, mapping=aes(x=Roles, y=CP, fill = Parameter)) +
  geom_boxplot() + 
  # scale_fill_brewer(palette = "Spectral") + 
  ylim(0,1) + 
  xlab("Number of Roles") + ylab("Coverage probabilities") + theme_bw() +
  # scale_fill_discrete(name = eval(parse(text= par_name))) +
  scale_fill_discrete(labels = par_name) +
  geom_hline(yintercept=0.95, linetype="dashed")+
  theme(text = element_text(size=32),
        axis.text.x =element_text(#face = "bold", #color = "#993333", 
          size = 28,angle = 90), 
        axis.text.y = element_text(#face = "bold", #color = "#993333", 
          size = 28),
        legend.key.height = unit(2, "line"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2),
        legend.position="none")
p
dev.off()

#################################################
# vary connectivity probability from 0.2 to 0.5 #
#################################################
pdf(paste0("Images/Recovery_performance_by_prob.pdf"),width = 8, height = 8)
par_name <- c(expression(lambda),expression(beta),expression(xi))
p <- ggplot(data=Recovery_prob, mapping=aes(x=Prob, y=CP, fill = Parameter)) +
  geom_boxplot() + 
  # scale_fill_brewer(palette = "Spectral") + 
  ylim(0,1) + 
  xlab("Connectivity probability") + ylab("Coverage probabilities") + theme_bw() +
  # scale_fill_discrete(name = eval(parse(text= par_name))) +
  scale_fill_discrete(labels = par_name) +
  geom_hline(yintercept=0.95, linetype="dashed")+
  theme(text = element_text(size=32),
        axis.text.x =element_text(#face = "bold", #color = "#993333", 
          size = 28,angle = 90), 
        axis.text.y = element_text(#face = "bold", #color = "#993333", 
          size = 28),
        legend.key.height = unit(2, "line"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2),
        legend.position="none")
p
dev.off()

#################################################
# vary average number of response from 10 to 35 #
#################################################
pdf(paste0("Images/Recovery_performance_by_resp.pdf"),width = 8, height = 8)
par_name <- c(expression(lambda),expression(beta),expression(xi))
p <- ggplot(data=Recovery_resp, mapping=aes(x=Response, y=CP, fill = Parameter)) +
  geom_boxplot() + 
  # scale_fill_brewer(palette = "Spectral") + 
  ylim(0,1) + 
  xlab("Average response number") + ylab("Coverage probabilities") + theme_bw() +
  # scale_fill_discrete(name = eval(parse(text= par_name))) +
  scale_fill_discrete(labels = par_name) +
  geom_hline(yintercept=0.95, linetype="dashed")+
  theme(text = element_text(size=32),
        axis.text.x =element_text(#face = "bold", #color = "#993333", 
          size = 28,angle = 90), 
        axis.text.y = element_text(#face = "bold", #color = "#993333", 
          size = 28),
        legend.key.height = unit(2, "line"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2),legend.position="none")
p
dev.off()

################################################
# Recovery performance of the original setting #
################################################
pdf(paste0("Images/Recovery_performance_Originial_setting.pdf"),width = 12, height = 6)
par_name <- c(expression(lambda),expression(beta),expression(xi))
p <- ggplot(data=Original_df, mapping=aes(x = Parameter, y=CP, fill = Parameter)) +
  geom_boxplot() + 
  # scale_fill_brewer(palette = "Spectral") + 
  ylim(0,1) + 
  xlab("Parameter") + ylab("Coverage probabilities") + theme_bw() +
  scale_x_discrete(labels= par_name) +
  guides(fill="none") +
  geom_hline(yintercept=0.95, linetype="dashed")+
  theme(text = element_text(size=32),
        axis.text.x =element_text(#face = "bold", #color = "#993333", 
          size = 28), 
        axis.text.y = element_text(#face = "bold", #color = "#993333", 
          size = 28),
        legend.key.height = unit(2, "line"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black",size = 2))
p
dev.off()
