library(clusterProfiler)
library(org.Hs.eg.db)
dat.medium <- read.csv("IL6_MediumVsHeavyDEX.txt",header=FALSE)
dat.light <- read.csv("IL6_LightVsHeavyDEX.txt",header=FALSE)

medium <- rep("Medium",each=dim(dat.medium)[1])
light <- rep("Light",each=dim(dat.light)[1])

m.df <- cbind(dat.medium,medium)
colnames(m.df) <- c("name","Type")
l.df <- cbind(dat.light,light)
colnames(l.df) <- c("name","Type")

my.df <- rbind(m.df,l.df)
my.df
df.gene.entrez <- bitr(my.df$name, fromType = "SYMBOL",
                       toType = c("ENSEMBL",  "ENTREZID"),
                       OrgDb = org.Hs.eg.db)
df.gene.entrez
mydf.gene <- merge(my.df, df.gene.entrez, by.x="name", by.y="SYMBOL")
gene.cluster <- compareCluster(ENTREZID~Type, data=mydf.gene, fun='enrichGO', OrgDb='org.Hs.eg.db')
dotplot(gene.cluster)
