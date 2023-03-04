# Library
library(fmsb)

# Create data: note in High school for several students
set.seed(99)
RaRv <- read.csv("E:\\MTB_Multi-PTM/Manuscript/H37Rv_Ra_Modified_Sites.txt", sep = '\t', header = T, row.names = 1)
#data <- as.data.frame(matrix( sample( 0:20 , 15 , replace=F) , ncol=5))
colnames(RaRv) <- c("Acetyla_K" , "Butyryl_K" , "Crotonyl_K" , "Citrullin_R" , "Glutaryl_K", "Glycosyl_N", "Glycosyl_P", "Hex_N", "Hex_S", "Malonyl_K",
                    "Methyl_K", "Methyl_R", "Phospho_S", "Phospho_T", "Phospho_Y")
rownames(RaRv) <- row.names(RaRv)

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
RaRv <- rbind(rep(3000,5) , rep(0,5) , RaRv)

# Color vector
colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )

# plot with default options:
radarchart(RaRv[,1:15], axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
            #custom the grid
            cglcol="grey", cglty=2, axislabcol="grey", caxislabels=seq(0,3000,1), cglwd=1,
            #custom labels
            vlcex=1 
)

# Add a legend
legend(x=1, y=1, legend = rownames(RaRv[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)

