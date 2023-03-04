library(ggplot2)
bp <- read.csv("IL6_Phosphotyrosine_MediumVsHeavyDEX_BP",sep="\t",header = TRUE)
cc <- read.csv("IL6_Phosphotyrosine_MediumVsHeavyDEX_CC",sep="\t",header = TRUE)
mf <- read.csv("IL6_Phosphotyrosine_MediumVsHeavyDEX_MF",sep="\t",header = TRUE)

colnames(bp)
colnames(cc)
colnames(mf)
data1 <- bp[,c(1,2,3,5)]
data2 <- cc[,c(1,2,3,5)]
data3 <- mf[,c(1,2,3,5)]
data1
data2
data3

data <- rbind(data1,data2,data3)
head(data)

library(viridis)
library(ggrepel)
#setEPS()
#postscript("EnrichedGO.eps",width=15,height=10)
#ggplot(data, aes(-log(PValue),log(Count), label=Term, color=-log(PValue)))+geom_point(aes(size = log(Count)))+
#  facet_wrap(~Category, scales = "free")
#  scale_color_viridis(discrete = FALSE, option = "D",direction=1)+
#  theme(axis.text.x = element_text(angle = 90, hjust = 5))+
#  geom_text_repel(size=3)+
#  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
#  xlab("Gene Ontology term")
#dev.off() 
head(data)
ggplot(data, aes(Term,-log(PValue),label=Term))+
  geom_point(aes(color = Category, size = log(Count),alpha=0.5)) +
    scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"),labels=c("BP","CC","MF")) +
    scale_size(range = c(0.5, 12))+
geom_text_repel(data=subset(data, -log(PValue) > 32),aes(Term,-log(PValue),label=Term),size=4)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
  
#png("EnrichedGO.png",width=600,height=500)
#ggplot(data, aes(Term,-log(PValue), label=Name, color=-log(PValue)))+geom_point(aes(size = log(Count)))+
#  facet_wrap(~Category, scales = "free")+
#  scale_color_viridis(discrete = FALSE, option = "D",direction=1)+
#  theme(axis.text.x = element_text(angle = 90, hjust = 5))+
#  geom_text_repel(size=3)+
#  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
#  xlab("Gene Ontology term")
#dev.off()  

