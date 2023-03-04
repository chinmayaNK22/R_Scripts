library(reshape)
library(ggplot2)
library(DEGreport)
library(DESeq2)
setwd('D:\\Phagosomes_NTNU/PD_Results/Reanalysis/Proteins_of_Interest_Identified/')

data <- read.csv('LIPASE.txt', sep = '\t', header = T)
colnames(data)
normalization_data <- cbind(data[,1:2], data[,10:16])
normalization_final <- melt(normalization_data)
colnames(normalization_final) <- c('Gene', 'Protein', 'Time', 'Normalized_Abundance')

ggplot(normalization_final) +
  geom_point(aes(x = Gene, y = log2(Normalized_Abundance), color = Time), position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("Normalized Abundance") +
  #ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))

