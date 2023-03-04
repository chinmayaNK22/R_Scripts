setwd("~/Desktop/Yenepoya/Sumrati_PhD/Manuscript3/Proteomics_Data_Analysis/")

data <- read.csv("BA_Urine_DEX_Proteins_abundance.txt", sep = "\t", row.names = 1, header = T)

colnames(data)

toplot <- data.matrix(data[1:24])

###plotting HeatMap####
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
pheatmap(toPlot.norm, cluster_cols = F, cluster_rows = T, cellwidth = 10, cutree_rows = 8, show_rownames = F, annotation_row = new_hc, annotation_colors = mycolors)

cols = colorRampPalette(c("maroon", "white", "darkgreen"))(100)
pheatmap(toPlot.norm[H$tree_row$order,], color=cols, cluster_cols = F,cluster_rows = T,cellwidth = 8,show_rownames = F, annotation_row = new_hc, annotation_colors = mycolors)
pheatmap(toPlot.norm[H$tree_row$order,], cluster_cols = F,cluster_rows = T,cellwidth = 8,show_rownames = F, annotation_row = new_hc, annotation_colors = mycolors)

Significant_ids <- cbind(data[match(rownames(new_hc),row.names(data)), ], new_hc)
write.table(Significant_ids, file="BA_Urine_Proteomics_MQ_Proteins_ANOVA_Tukeys_Cluster_022823.txt",sep="\t",row.names = T)

####ComplexHeatmap####
library(ComplexHeatmap)
library(circlize)
library(dendextend)

dend = as.dendrogram(hclust(dist(toPlot.norm)))
dend = color_branches(dend, k = 8)

f1 = colorRamp2(seq(min(toPlot.norm), max(toPlot.norm), length = 3), c("maroon", "#EEEEEE", "darkgreen"), space = "RGB")
Heatmap(toPlot.norm, name = "DEX", col = f1, column_names_rot = 80, cluster_rows = dend, cluster_columns = FALSE, 
        row_split  = 8)
