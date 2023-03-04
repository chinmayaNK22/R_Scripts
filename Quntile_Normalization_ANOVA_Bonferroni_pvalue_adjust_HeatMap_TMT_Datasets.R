setwd("Desktop/Yenepoya/Sumrati_PhD/Manuscript3/Proteomics_Data_Analysis/Serum_MQ_Result/")

data2use = read.csv("proteinGroups_formatted.txt", sep="\t", header=T, row.names = 1)

colnames(data2use)
dim(data2use)

#Replacing NA's with zero
data2use[is.na(data2use)] <- 0
a <- data2use
a$new <- rowSums(a[,54:77])

head(a)
#Extracting rows with all zero values
a.extracted <- a[which(a$new > 0),]

#Add a unique number to each row as id and prepare a matrix
#uniq.id <-c(1:length(a.extracted$new ))
#a.extracted <- cbind(a.extracted,uniq.id)
#row.names(a.extracted) <- t(a.extracted[ncol(a.extracted)])
#View(a.extracted)
#a.mat <- a.extracted[,5:28]
#rownames(a.mat) <- rownames(a.extracted)

a.mat <- a.extracted[,54:77] 
#converting into matrix
data.m <- as.matrix(a.mat)

#Normalization  
library(preprocessCore)
qs.data <- normalize.quantiles(data.m)
colnames(qs.data) <- colnames(data.m)
colnames(qs.data)
rownames(qs.data)<- rownames(data.m)
combat_edata <- qs.data
back <- data.m

#Batch correction
#library(sva)
#mat_val_b <- read.table("batchR1.txt", header = T, sep="\t",row.names = 1)
#modcombat <- model.matrix(~1, data=mat_val_b)
#combat_edata <- ComBat(dat=qs.data, batch = mat_val_b$batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
#combat_edata <- qs.data

#plotting
par(mfrow=c(1,2))
boxplot(qs.data, ylim=c(0,15000), las=2,col= rep(c("maroon","blue","yellow", "grey","green", "purple", "red"),each=6))
boxplot(data.m, ylim=c(0,15000), las=2,col= rep(c("maroon","blue","yellow", "grey","green", "purple", "red"),each=6))

CS.t <- cbind(a.extracted[match(rownames(a.extracted),row.names(back)), ], back)
write.table(CS.t, file="BA_Serum_Proteomics_MQ_Protein_Abundance_Quantile_Normalized_070322.txt",sep="\t")
dim(back)
dim(a.extracted)

#Anova analysis
sampleInfo = read.table("Sample_info.txt", sep="\t",stringsAsFactors=T,header=T)
TIME_ORDER2 <- c("Control", "HIE1", "HIE2", "HIE3")
cov.time2 = sampleInfo[match(colnames(back), sampleInfo$Sample), ]$Condition
cov.time2 = factor(as.character(cov.time2),levels= TIME_ORDER2)
model.fits = lm(t(back) ~  cov.time2)
res = summary(aov(model.fits))
head(res)

#Extracting respone and p values from "res" and giving responeid as rownames
names(res) = sapply(names(res),function(id){unlist(strsplit(id,"Response "))[2]})
pValues = t(as.data.frame(lapply(res, function(thisRes){ thisRes$"Pr(>F)"[1] })))
rownames(pValues) = names(res)
colnames(pValues) = c('P-value')

#Adjusting p-values by Bonferroni method
pValues_adj = pValues
method="bonferroni"
pValues_adj[,1] = p.adjust(pValues[,1], method)
head(pValues_adj)

#### ANOVA - Tukey's Post hoc analysis ####
back_T <- t(back)

if(rownames(back)[1] == colnames(back_T)[1]) {
  dim(back_T)[2]
  l <- list()
  for(i in 1:dim(back_T)[2]){
    model.fits1 = lm(back_T[,i] ~  cov.time2)
    anova_pvalue1 <- aov(model.fits1)
    tus <-TukeyHSD(anova_pvalue1)
    row <- colnames(back_T)[i]
    val <- tus$cov.time2[,4]
    a <- list(val)
    names(a)[1] <- c(row)
    l <- append(l,a)
  }
} else {
  print("Wromg data matrix is used")
}

ANOVA_tukeys <-do.call(rbind, lapply(l, function(x) x[match(names(l[[1]]), names(x))]))

####Extracting significantly expressed peptides####
tukeys.m <- data.matrix(ANOVA_tukeys)
dim(tukeys.m)
pvaluethreshold = 0.05

Si <- c()
for(j in 1:dim(tukeys.m)[2]){
  ab <- which(tukeys.m[,j] < pvaluethreshold)
  Si <- c(Si,row.names(as.data.frame(tukeys.m[ab,])))
}
Si <- unique(Si)
toplot <- back[match(Si, rownames(back)), ]

CS.t2 <- cbind(CS.t[match(row.names(ANOVA_tukeys), rownames(CS.t)), ], ANOVA_tukeys)
CS.t3 <- cbind(CS.t2[match(row.names(pValues_adj), rownames(CS.t2)),],pValues_adj)
write.table(CS.t3, file="BA_Serum_Proteomics_MQ_Protein_Abundance_Quantile_Normalized_ANOVA_Tukeys_Post-hoc_pvalues_071722.txt",sep="\t")

####plotting HeatMap####
library(pheatmap)

normalization1<-function(x){
  dimm=dim(x)
  for(i in 1:dimm[1]){
    x[i,]=(x[i,]-mean(x[i,]))/sd(x[i,]) 
  }
  return(x)
}

toPlot.norm <- normalization1(toplot)

#plotting

H = pheatmap(toPlot.norm, cluster_cols = TRUE)
#rownames(toPlot.norm[H$tree_row[["order"]],])
hc_clust <- cutree(H$tree_row, k=8)
#rownames(toPlot.norm[H$tree_row[["order"]],])
hc_clust <- cutree(H$tree_row, k=8)
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
pheatmap(toPlot.norm[H$tree_row$order,], cluster_cols = F, cluster_rows = T, cellwidth = 10, cutree_rows = 8, show_rownames = F, annotation_row = new_hc,annotation_colors = mycolors)

cols = colorRampPalette(c("darkblue", "white", "darkred"))(100)
pheatmap(toPlot.norm[H$tree_row$order,],color=cols,cluster_cols = F,cluster_rows = T,cellwidth = 10,cutree_rows = 8,show_rownames = F, annotation_row = new_hc,annotation_colors = mycolors)

Significant_ids <- cbind(CS.t3[match(rownames(new_hc),row.names(CS.t3)), ], new_hc)
write.table(Significant_ids, file="BA_Serum_Proteomics_MQ_Proteins_ANOVA_Tukeys_Cluster_071722.txt",sep="\t",row.names = T)
