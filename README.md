# Genome-wide-in-vivo-identification-of-virulence-factors-for-streptococcal-infective-endocarditis
#In the paper entitled "Genome-wide in vivo identification of virulence factors for streptococcal infective endocarditis", we have developed codes for counting of mutant aboundance and visulizing data in R and unix. In total we have generated 5 displays including figure 2F, 2G, 3A, S2, S4A and S4B.
#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Figure 2F is "abundance ratio" volcano plot for controls"
jpeg("volcano plot for controls.jpg", width=24000, height=18000, res=1800)
adj.p.c<-read.csv("adj.P.controls.csv")
ggplot(adj.p.c, aes(x = log(mean.PR), y = -log10(P.value.by.one.samplet.t.test))) +
  xlim(-16, 5) +  
  ylim(0, 25)  +  
  geom_point(aes(color = Mutants), alpha = 0.8, size = 1.5,shape = 19) +
  scale_color_manual(values = c("red", "grey","green","purple","black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = log(1), linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(
    x = "Log2 Ratio",
    y = "-Log10 P-Value",
    color = "Mutants"
  ) +
  theme(
    text = element_text(face = "bold", color = "black",size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(face = "bold", color = "black", size = 20),  # Bold, black x-axis text
    axis.text.y = element_text(face = "bold", color = "black", size = 20),  # Bold, black y-axis text
    axis.title.x = element_text(face = "bold", color = "black", size = 25), # Bold, black x-axis title
    axis.title.y = element_text(face = "bold", color = "black", size = 25)  # Bold, black y-axis title
    
  ) 
dev.off()

#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Figure 2G is "abundance ratio" of different mutants"
#virulence.figure3
virulence.figure3<-read.csv("virulence.figure3.1.csv",header = T)
jpeg("virulence.figure3.jpg", width=36000, height=18000, res=1800)

ggplot(virulence.figure3, aes(x =mutant, y = log(virulence), fill = mutant)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  labs(
    x = "Mutants",
    y = "Log2 Ratio")+
  theme(legend.position = "",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,face = "bold", color = "black", size = 15),text = element_text(face = "bold", color = "black",size = 20),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(face = "bold", color = "black", size = 15),  # Bold, black y-axis text
        axis.title.x = element_text(face = "bold", color = "black", size = 20), # Bold, black x-axis title
        axis.title.y = element_text(face = "bold", color = "black", size = 20)  # Bold, black y-axis title
  ) 
dev.off()

#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Figure 3A is "volcano plot for all mutants"
jpeg("volcano plot for all mutants.jpg", width=24000, height=18000, res=1800)
p.adjusted.all$significant <- ifelse(p.adjusted.all$P.value.by.one.samplet.t.test < 0.05 & log(p.adjusted.all$mean.PR) < 0, "yes", "no")
ggplot(p.adjusted.all, aes(x = log(mean.PR), y = -log10(P.value.by.one.samplet.t.test))) +
  geom_point(aes(color = significant),  alpha = 0.8, size = 1.5,shape = 19) +
  xlim(-16, 5) +  
  ylim(0, 25)  +  
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = log(1), linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Ratio",
    y = "-Log10 P-Value",
    color = "Significant"
  ) +
  theme(
    text = element_text(face = "bold", color = "black",size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(face = "bold", color = "black", size = 20),  # Bold, black x-axis text
    axis.text.y = element_text(face = "bold", color = "black", size = 20),  # Bold, black y-axis text
    axis.title.x = element_text(face = "bold", color = "black", size = 25), # Bold, black x-axis title
    axis.title.y = element_text(face = "bold", color = "black", size = 25)  # Bold, black y-axis title
    
  ) 
dev.off()

#----------------------------------------------------------------------------------------------------------------------------------------------------------
#Figure S2 is the pair-wise comparison of technical replicates
#technical/biological reproducibility in vegetation for final1.
final1<-read.csv("final1-one inpout and five outputs.csv")
pdf("f1.serum.BHI2.pdf", width = 8, height = 6, pointsize = 12)
ggpairs(final1[,2:19],lower = list( continuous = wrap("points", size = .1)),upper = list(continuous = wrap("cor", size=2)))+
  theme(legend.position = "",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,face = "bold", color = "black", size = 5),text = element_text(face = "bold", color = "black",size = 5),plot.title = element_text(hjust = 0.5))
        

dev.off()

#----------------------------------------------------------------------------------------------------------------------------------------------------------
#Figure S4A is tested.once mutants
#volcano plot for tested.once 
jpeg("volcano plot for tested.once mutants.jpg", width=24000, height=18000, res=1800)
tested.once<-read.csv("tested.once.csv")
tested.once$significant <- ifelse(tested.once$P.value.by.one.samplet.t.test < 0.05& log(tested.once$mean.PR) < 0, "yes", "no")
ggplot(tested.once, aes(x = log(mean.PR), y = -log10(P.value.by.one.samplet.t.test))) +
  xlim(-16, 5) +  
  ylim(0, 25)  +  
  geom_point(aes(color = significant),  alpha = 0.8, size = 1.5,shape = 19) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = log(1), linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(
    x = "Log2 Ratio",
    y = "-Log10 P-Value",
    color = "Significant")+
  theme(
    text = element_text(face = "bold", color = "black",size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(face = "bold", color = "black", size = 20),  # Bold, black x-axis text
    axis.text.y = element_text(face = "bold", color = "black", size = 20),  # Bold, black y-axis text
    axis.title.x = element_text(face = "bold", color = "black", size = 25), # Bold, black x-axis title
    axis.title.y = element_text(face = "bold", color = "black", size = 25)  # Bold, black y-axis title
    
  ) 
dev.off()


#----------------------------------------------------------------------------------------------------------------------------------------------------------
#Figure S4B is "tested.twice.or.more mutants"
#volcano plot for tested.twice.or.more mutants
jpeg("volcano plot for tested.twice.or.more mutants.jpg", width=24000, height=18000, res=1800)
tested.twice.or.more<-read.csv("tested.twice.or.more.csv")
tested.twice.or.more$significant <- ifelse(tested.twice.or.more$P.value.by.one.samplet.t.test < 0.05& log(tested.twice.or.more$mean.PR) < 0, "yes", "no")
ggplot(tested.twice.or.more, aes(x = log(mean.PR), y = -log10(P.value.by.one.samplet.t.test))) +
  geom_point(aes(color = significant),  alpha = 0.8, size = 1.5,shape = 19) +
  xlim(-16, 5) +  
  ylim(0, 25)  +  
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = log(1), linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(
    x = "Log2 Ratio",
    y = "-Log10 P-Value",
    color = "Significant")+
  theme(
    text = element_text(face = "bold", color = "black",size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(face = "bold", color = "black", size = 20),  # Bold, black x-axis text
    axis.text.y = element_text(face = "bold", color = "black", size = 20),  # Bold, black y-axis text
    axis.title.x = element_text(face = "bold", color = "black", size = 25), # Bold, black x-axis title
    axis.title.y = element_text(face = "bold", color = "black", size = 25)  # Bold, black y-axis title
    
  ) 
dev.off()



