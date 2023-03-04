library('Hmisc')
library('corrplot')
library('vegan')
library('igraph')

setwd('D:\\Sumrati_PhD/Manuscript3/Proteomics_results')
proteins <- read.csv("BA_Serum_and_Urine_Proteomics_MQ_Proteins_Common_DEX_corr_input_010323_final.txt", sep = '\t', header = 1, row.names = 1)

colnames(proteins)
dim(proteins)

t.mat <- as.matrix(t(proteins[,1:6]))
normal.mat <- as.matrix(proteins[,1:6])

normalization1<-function(x){
  dimm=dim(x)
  for(i in 1:dimm[1]){
    x[i,]=(x[i,]-mean(x[i,]))/sd(x[i,]) 
  }
  return(x)
}

toPlot.norm <- normalization1(normal.mat)

## Correlation of values from matrix by pearson or spearman
pcorr.m <- rcorr(toPlot.norm, type=c("spearman"))
#scorr.m <- rcorr(pep.m, pro.m, type=c("spearman"))
pcorrP <- pcorr.m$P
pcorrR <- pcorr.m$r
#M <- cor(pcorr)
cols = colorRampPalette(c("red", "white", "darkgreen"))(100)
pearson.plot <- corrplot(pcorrR, method = "color", order = "hclust", hclust.method = "ward.D2", col=cols, 
                         tl.cex = 0.5, tl.col = "black", sig.level = 0.01)

pearson.plot <- corrplot(pcorrR, method = "color", order = "hclust", col=cols,
                         tl.cex = 0.5, tl.col = "black", sig.level = 0.01)

cols = colorRampPalette(c("red", "white", "darkgreen"))(100)
pearson.plot <- corrplot(pcorrR, method = "color", order = "hclust", col=cols, tl.cex = 0.5, tl.col = "black", sig.level = 0.01)

#new_hc <- c()
#for (i in 1:length(pearson)){
#  new_hc[i] <- paste("Cluster ",pearson[i])
#  names(new_hc)[i] <- names(pearson[i])
#}
#new_hc <- data.frame(new_hc)
#colnames(new_hc) <- "Clusters"
write.table(pcorrR, file="BA_Serum_and_Urine_Pearson_Correlated_significant_proteins_022523.txt",sep="\t")

# Multiple testing correction using Benjamini-Hochberg standard false discovery rate correction
#matrix.cor.p <- p.adjust(pcorrP, method = "BH")

# Consider positive cooccurence at given coefficient (cor.cutoff) and p-value cutoffs
matrix.cor1 <- pcorrR
matrix.cor1.p <- pcorrP
matrix.cor1[which(matrix.cor1 <= 0.75)] <- 0
matrix.cor1[which(matrix.cor1.p > 0.01)] <- 0
# delete those rows and columns with sum = 0
matrix.cor1 <- matrix.cor1[which(rowSums(matrix.cor1) != 1), ]
matrix.cor1 <- matrix.cor1[, which(colSums(matrix.cor1) != 0)]

# Consider netagive cooccurence at given coefficient (-cor.cutoff) and p-value cutoffs
matrix.cor2 <- pcorrR
matrix.cor2.p <- pcorrP
matrix.cor2[which(matrix.cor2 > (-0.75))] <- 0
matrix.cor2[which(matrix.cor2.p > 0.01)] <- 0
# delete those rows and columns with sum = 0
matrix.cor2 <- matrix.cor2[which(rowSums(matrix.cor2) != 0), ]
matrix.cor2 <- matrix.cor2[, which(colSums(matrix.cor2) != 0)]

# Consider both positive and netagive cooccurence at given coefficient (cor.cutoff) and p-value cutoffs
matrix.cor3 <- pcorrR
matrix.cor3.p <- pcorrP
matrix.cor3[which(matrix.cor3 >= (-0.75) & matrix.cor3 <= 0.75)] <- 0
matrix.cor3[which(matrix.cor3.p > 0.01)] <- 0

# delete those rows and columns with sum = 0
matrix.cor3 <- matrix.cor3[which(rowSums(matrix.cor3) != 1), ]
matrix.cor3 <- matrix.cor3[, which(colSums(matrix.cor3) != 0)]

# get pairs r
# This is to remove redundancy as upper correlation matrix == lower
ma1 <- matrix.cor1
ma2 <- matrix.cor2
ma3 <- matrix.cor3
ma1[upper.tri(matrix.cor1, diag = TRUE)] <- NA
pair.r1 <- reshape2::melt(ma1, na.rm = TRUE, value.name = "cor")
ma2[upper.tri(ma2, diag = TRUE)] <- NA
pair.r2 <- reshape2::melt(ma2, na.rm = TRUE, value.name = "cor")
ma3[upper.tri(ma3, diag = TRUE)] <- NA
pair.r3 <- reshape2::melt(ma3, na.rm = TRUE, value.name = "cor")
pair.r1<-pair.r1[which(pair.r1[,3]!=0),]
pair.r2<-pair.r2[which(pair.r2[,3]!=0),]
pair.r3<-pair.r3[which(pair.r3[,3]!=0),]
write.csv(pair.r1, file = "BA_Serum_and_Urine_Sig_Proteins_Pos_otu_022523.csv",quote = F, row.names = F)
write.csv(pair.r2, file = "BA_Serum_and_Urine_Sig_Proteins_Neg_otu_022523.csv",quote = F,row.names = F)
write.csv(pair.r3, file = "BA_Serum_and_Urine_Sig_Proteins_PosNeg_otu_022523.csv",quote = F,row.names = F)

# generating graph using igraph
g1 <- graph.adjacency(matrix.cor1, weight = T, mode = "undirected")
g1 <- simplify(g1)
V(g1)$label <- V(g1)$name
V(g1)$degree <- degree(g1)

g2 <- graph.adjacency(matrix.cor2, weight = T, mode = "undirected")
g2 <- simplify(g2)
V(g2)$label <- V(g2)$name
V(g2)$degree <- degree(g2)

g3 <- graph.adjacency(matrix.cor3, weight = T, mode = "undirected")
g3 <- simplify(g3)
V(g3)$label <- V(g3)$name
V(g3)$degree <- degree(g3)

# append the output into results
result <- list()
result$matrix.cor <- pcorrR
result$matrix.cor.p <- pcorrP

result$matrix.cor1 <- matrix.cor1
result$graph1 <- g1

result$matrix.cor2 <- matrix.cor2
result$graph2 <- g2

result$matrix.cor3 <- matrix.cor3
result$graph3 <- g3
return(result)

#install.packages("gtools")
library(gtools)
f <- foldchange(Abu[,1],Abu[,14])
cbind(Abu[,1],Abu[,14],a,b,f)

# Co-occurrence-network-analysis
## OTU filtering, network generation, topological analysis and export OTU table
library(igraph)
library(Hmisc)

setwd("director/path")

Abu <- read.table("oturelative.txt", header = T)
Abu <- read.table("otuabu.txt", header = T)
Abu <- as.matrix(Abu)

### Filtering OTUs
table <- Abu
table[table > 0] <- 1
table.generalist <- Abu[which(rowSums(table) >= 12), ]
Abu <- table.generalist

## Creating gml files of network (to be visulized in Gephi or Cytoscape)

## cutoffs for correlation coefficient and P-value
pattern <- coRnetwork(Abu, 0.6, 0.01)

#write.graph(pattern$graph1, "Pos0.6-rela.gml", format = "gml") # network file for positive association
#write.graph(pattern$graph2, "Neg0.6-rela.gml", format = "gml") # network file for negative association
#write.graph(pattern$graph3, "PosNeg0.6-rela.gml", format = "gml") # network file for all association

write.graph(result$graph1, "SARS_CoV-2_Sig_Mod_Peptides_Proteins_Pos_otu_052120.gml", format = "gml") # network file for positive association
write.graph(result$graph2, "SARS_CoV-2_Sig_Mod_Peptides_Proteins_Neg_otu_052120.gml", format = "gml") # network file for negative association
write.graph(result$graph3, "SARS_CoV-2_Sig_Mod_Peptides_Proteins_PosNeg_otu_052120.gml", format = "gml") # network file for all association
