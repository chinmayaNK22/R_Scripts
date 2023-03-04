library(tidyverse)
library(preprocessCore)
setwd("E:\\MTB_Multi-PTM\\Multi-PTM_PXD009549_MTB")
peptides <- read.csv("PXD009549_Multi-PTM_052220_PeptideGroups.txt", sep = '\t', header = T)
colnames(peptides)

peptides <- mutate(peptides, ID = paste(Annotated.Sequence,Modifications, X..Protein.Groups, Theo..MH...Da.))
row.names(peptides) <- peptides$ID

peptides.df <- peptides[,3:30]

#Replacing NA's with zero
peptides.df[is.na(peptides.df)] <- 0
a <- peptides.df
a$new <- rowSums(a[,18:28])

head(a)
#Extracting rows with all zero values
a.extracted <- a[which(a$new > 0),]
data.m <- data.matrix(a.extracted[,18:28])

qs.data <- normalize.quantiles(data.m)
colnames(qs.data) <- colnames(data.m)
colnames(qs.data)
rownames(qs.data)<- rownames(data.m)
combat_edata <- qs.data
back <- qs.data

par(mfrow=c(1,2))
boxplot(qs.data, ylim=c(0,1000000), las=2,col= rep(c("maroon","blue","yellow", "grey","green", "purple", "red"),each=3))
boxplot(data.m, ylim=c(0,1000000), las=2,col= rep(c("maroon","blue","yellow", "grey","green", "purple", "red"),each=3))

CS.t <- cbind(a.extracted[match(rownames(a.extracted),row.names(back)), ], back)
dim(back)
dim(a.extracted)

#Anova analysis
sampleInfo = read.table("E:\\MTB_Multi-PTM/Multi-PTM_PXD009549_MTB/ColumnNames.txt", sep="\t",stringsAsFactors=T,header=T)
TIME_ORDER2 <- c("Culture_Filtrate", "Cell_wall", "Cytosol", "Membrane")
cov.time2 = sampleInfo[match(colnames(back), sampleInfo$ColumnNames), ]$Subcell
cov.time2 = factor(as.character(cov.time2),levels= TIME_ORDER2)
model.fits = lm(t(back) ~  cov.time2)
res = summary(aov(model.fits))
head(res)

#Extracting respone and p values from "res" and giving responeid as rownames
names(res) = sapply(names(res),function(id){unlist(strsplit(id,"Response "))[2]})
pValues = t(as.data.frame(lapply(res, function(thisRes){ thisRes$"Pr(>F)"[1] })))
rownames(pValues) = names(res)
colnames(pValues) = c('p-Values')

#Adjusting p-values by Bonferroni method
pValues_adj = pValues
method="bonferroni"
pValues_adj[,1] = p.adjust(pValues[,1], method)
head(pValues_adj)

pValThreshold = 0.05

#Extracting significantly expressed peptides
sig <- which(pValues_adj[,1] < pValThreshold)
Si <- row.names(as.data.frame(sig))
toPlot <- data.m[match(Si, rownames(data.m)), ]

CS.t2 <- cbind(CS.t[match(row.names(pValues_adj), rownames(CS.t)), ], pValues_adj)
CS.t3 <- cbind(CS.t[match(row.names(new_hc), rownames(CS.t)), ], new_hc)
write.table(CS.t2, file="PXD009549_Peptides_ANOVA_Bonferroni_adjusted_pvalues_052720.txt",sep="\t")
write.table(CS.t3, file="PXD009549_Peptides_ANOVA_Bonferroni_adjusted_pvalues_Clusters_052720.txt",sep="\t")

normalization1<-function(x){
  dimm=dim(x)
  for(i in 1:dimm[1]){
    x[i,]=(x[i,]-mean(x[i,]))/sd(x[i,]) 
  }
  return(x)
}

toPlot.norm <- normalization1(toPlot)

#plotting
library(pheatmap)
H = pheatmap(toPlot.norm, cluster_cols = TRUE)
#rownames(toPlot.norm[H$tree_row[["order"]],])
hc_clust <- cutree(H$tree_row, k=4)
#rownames(toPlot.norm[H$tree_row[["order"]],])
hc_clust <- cutree(H$tree_row, k=4)
new_hc <- c()
for (i in 1:length(hc_clust)){
  new_hc[i] <- paste("Cluster ",hc_clust[i])
  names(new_hc)[i] <- names(hc_clust[i])
}
new_hc <- data.frame(new_hc)
colnames(new_hc) <- "Clusters"
newCols <- colorRampPalette(grDevices::rainbow(length(unique(new_hc$Clusters))))
mycolors <- newCols(length(unique(new_hc$Clusters)))
names(mycolors) <- unique(new_hc$Clusters)
mycolors <- list(Clusters = mycolors)
pheatmap(toPlot.norm[H$tree_row$order,],cluster_cols = FALSE,cluster_rows = FALSE,cellwidth = 4,cutree_rows = 6,show_rownames = F, annotation_row = new_hc,annotation_colors = mycolors)
pheatmap(toPlot.norm[H$tree_row$order,],cluster_cols = FALSE,cluster_rows = FALSE,cellwidth = 4,cutree_rows = 6,show_rownames = F, annotation_row = new_hc,annotation_colors = mycolors )
