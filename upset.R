library(UpSetR)
library(ggplotify)

setwd("E:\\PXD017710/Summarization_PTMs_062821/")
myData1 <- read.csv("ForUpset3.txt",sep="\t", header=TRUE)

head(myData1)

upset(myData1, keep.order = F, sets = c("Citrullination","Methylation","Phosphorylation","Succinylation","Crotonylation","Acetylation","Ubiquitination","Hex","Glutarylation","Glycosylation"), 
      sets.bar.color = c("blue"),
      att.color = c("red"),matrix.color="red", order.by = "freq", point.size = 6,line.size=2,text.scale=2)
