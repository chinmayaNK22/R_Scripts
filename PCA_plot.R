setwd()

data <- read.csv("input.txt", sep = '\t', 
                 row.names = 1, header=TRUE)

colnames(data)

data_t <- sapply(data.frame(t(data)[(138:161),]), as.numeric)

rownames(data_t) <- row.names(data.frame(t(data)[(138:161),]))

data.t <- data.frame(data_t)

data.pca <- prcomp(data.t, center = TRUE,scale. = TRUE)

#A file consisting of two columns, first one for the condition name in the input.txt and the second one is the actual condition name to be displayed in the figure
SampleInfo <- read.table("SampleInfo.txt", sep = "\t", row.names = 1, header = T)

data.t <- cbind(data.t[match(rownames(SampleInfo), row.names(data.t)),],SampleInfo)

summary(data.pca)

biplot(data.pca, scale = 0)

library("factoextra")
fviz_pca_ind(data.pca, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = data.t$Condition, 
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Conditions") +
  ggtitle("2D PCA-plot") +
  theme(plot.title = element_text(hjust = 0.5))


library(ggfortify)
autoplot(data.pca)
autoplot(data.pca, data = data.t, colour = "Condition",  frame = TRUE, frame.type = 't',
         frame.colour = 'Condition', size = 5) + scale_size_continuous(range = c(1,3))

