library(ggplot2)

setwd("E:\\NTMs_Raw_Data/MSMS_Cluster/")

df2 <- read.csv("Species_specific_clusters_scans_new.txt", sep = '\t', header = T)

colnames(df2)

p <- ggplot(df2, aes(x=P.value, y=Scans, group=Species)) +
  geom_line(aes(color=Species, linetype=Species),size = 2)+
  geom_point(aes(color=Species))+
  labs(title="Species-specific scans",x="Species", y = "Scans")+
  theme_classic()
p + expand_limits(y=c(0, 3000000))


###
p <-ggplot(df2, aes(x=Species, y=Scans, group=P.value)) +
  geom_line(aes(color=P.value, linetype=P.value), size = 1)+
  geom_point(aes(color=P.value))+
  labs(title="Species-specific scans",x="Species", y = "Scans")+
  theme_classic()
p + expand_limits(y=c(0, 2500000))

