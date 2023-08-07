setwd("D:\\OSCC_DIA/DIA_Proteomics_April_2021/Post_Processing_Analysis/")
data <- read.csv("OSCC_DIA_Protein_Abundance_Quantile_Normalized_Min_2_abundance_ANOVA_Tukeys_Fold_Change_092321.txt", sep = '\t', header = T, row.names = 1)
class(data)
colnames(data)
data.all <- as.matrix(data)
data.mat <- data[29:34]
colnames(data.mat)
## Volcano plot
par(mfrow=c(1,2))
with(data.mat, plot(log2(Fold_Change..Smokers_Vs_Control.), -log10(p.value..Smokers.Control.), pch=20, main="Volcano plot Smokers vs Control", ylab="-log10(pvalue)"))
with(subset(data.mat, -log10(p.value..Smokers.Control.) > 2 & log2(Fold_Change..Smokers_Vs_Control.)>1), points(log2(Fold_Change..Smokers_Vs_Control.), -log10(p.value..Smokers.Control.), pch=25, col="green"))
with(subset(data.mat, -log10(p.value..Smokers.Control.) > 2 & log2(Fold_Change..Smokers_Vs_Control.)>1), points(log2(Fold_Change..Smokers_Vs_Control.), -log10(p.value..Smokers.Control.), pch=20, col="green"))
with(subset(data.mat, -log10(p.value..Smokers.Control.) > 2 & log2(Fold_Change..Smokers_Vs_Control.)<(-1)), points(log2(Fold_Change..Smokers_Vs_Control.), -log10(p.value..Smokers.Control.), pch=25, col="red"))
with(subset(data.mat, -log10(p.value..Smokers.Control.) > 2 & log2(Fold_Change..Smokers_Vs_Control.)<(-1)), points(log2(Fold_Change..Smokers_Vs_Control.), -log10(p.value..Smokers.Control.), pch=20, col="red"))
abline(h = 2, col = "blue", lty = 2, lwd = 1)
abline(v = c(-1,1), col = "blue", lty = 2, lwd = 1)

with(data.mat, plot(log2(Fold_Change..Chewers_Vs_Control.), -log10(p.value..Chewers.Control.), pch=20, main="Volcano plot Chewers vs Control", ylab="-log10(pvalue)"))
with(subset(data.mat, -log10(p.value..Chewers.Control.) > 2 & log2(Fold_Change..Chewers_Vs_Control.)>1), points(log2(Fold_Change..Chewers_Vs_Control.), -log10(p.value..Chewers.Control.), pch=25, col="green"))
with(subset(data.mat, -log10(p.value..Chewers.Control.) > 2 & log2(Fold_Change..Chewers_Vs_Control.)>1), points(log2(Fold_Change..Chewers_Vs_Control.), -log10(p.value..Chewers.Control.), pch=20, col="green"))
with(subset(data.mat, -log10(p.value..Chewers.Control.) > 2 & log2(Fold_Change..Chewers_Vs_Control.)<(-1)), points(log2(Fold_Change..Chewers_Vs_Control.), -log10(p.value..Chewers.Control.), pch=25, col="red"))
with(subset(data.mat, -log10(p.value..Chewers.Control.) > 2 & log2(Fold_Change..Chewers_Vs_Control.)<(-1)), points(log2(Fold_Change..Chewers_Vs_Control.), -log10(p.value..Chewers.Control.), pch=20, col="red"))
abline(h = 2, col = "blue", lty = 2, lwd = 1)
abline(v = c(-1,1), col = "blue", lty = 2, lwd = 1)

## Find the differentially expressed smokers vs control
MH.up <- rownames(data.mat[which(data.mat$Fold_Change..Smokers_Vs_Control. > 2.0  & data.mat$p.value..Smokers.Control. < 0.01),])
MH.down <- rownames(data.mat[which(data.mat$Fold_Change..Smokers_Vs_Control. < 0.50  & data.mat$p.value..Smokers.Control. < 0.01),])
MH.dex <- c(MH.up,MH.down)
MH.t1 <- cbind(data.all[match(row.names(data.mat),data.all),],data.mat)
MH.t2 <- MH.t1[match(MH.dex,data.all),]
write.table(MH.t2,file="OSCC_DIA_Smokers_vs_Control_DEX_092321.txt",sep="\t")

## Find the differentially expressed peptides Light vs Heavy
LH.up <- rownames(toVolcano.LH[which(toVolcano.LH$fc1 > 1.5  & toVolcano.LH$re1 < 0.05),])
LH.down <- rownames(toVolcano.LH[which(toVolcano.LH$fc1 < 0.67  & toVolcano.LH$re1 < 0.05),])
LH.dex <- c(LH.up,LH.down)
LH.t <- cbind(a.extracted[match(row.names(combat_edata), rownames(a.extracted)), ], combat_edata)
LH.t1 <- cbind(LH.t[match(row.names(toVolcano.LH),LH.t$uniq.id),],toVolcano.LH)
LH.t2 <- LH.t1[match(LH.dex,LH.t1$uniq.id),]
write.table(LH.t2,file="LH_dex_all_062619.txt",sep="\t")


#####volcano plot by ggplot ####
library(ggrepel)
# plot adding up all layers we have seen so far

# add a column of NAs
data.mat$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
data.mat$diffexpressed[log2(data.mat$Fold_Change..Smokers_Vs_Control.) > 1 & data.mat$p.value..Smokers.Control. < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
data.mat$diffexpressed[log2(data.mat$Fold_Change..Smokers_Vs_Control.) < -1 & data.mat$p.value..Smokers.Control. < 0.05] <- "DOWN"
data.mat$id <- row.names(data.mat)
ggplot(data=data.mat, aes(x=log2(Fold_Change..Smokers_Vs_Control.), y=-log10(p.value..Smokers.Control.), col=diffexpressed)) +
  geom_point(size = 3) + 
  theme_minimal() +
  geom_text_repel(data=subset(data.mat, data.mat$diffexpressed == "UP"|data.mat$diffexpressed =="DOWN"), aes(label=id),hjust=0, vjust=0, size = 5, col="black") +
  scale_color_manual(values=c("blue", "black", "red"))+
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

#### Volcano plot with different symbols ####
ggplot(data=hie3, aes(x=log2(`Fold Change (HIE3/Healthy)`), y=-log10(`P-value`), label = Gene)) +
  geom_point(data=subset(hie3, hie3$diffexpressed == "UP"), shape=24, fill="green", size = 3, stroke = 0) +
  geom_point(data=subset(hie3, hie3$diffexpressed == "DOWN"), shape=25, fill="red", size = 3, stroke = 0) +
  geom_point(data=subset(hie3, hie3$diffexpressed == "NO"), color="black", fill="#B8B7B5", size = 2, shape = 21, stroke = 0.1) +
  theme_minimal() +
  geom_text_repel(data=subset(hie3, hie3$diffexpressed == "UP"), aes(label=Gene),hjust=0, vjust=0, size = 5, col="black") +
  geom_text_repel(data=subset(hie3, hie3$diffexpressed =="DOWN"), aes(label=Gene),hjust=0, vjust=0, size = 5, col="black") +
  #scale_color_manual(values=c("red", "#B8B7B5", "green"))+
  geom_vline(xintercept=c(-0.58, 0.58), col="blue", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="blue", linetype = "dashed")
