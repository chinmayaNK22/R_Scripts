####
setwd("D:\\Sumrati_PhD/Manuscript3/Proteomics_results/GO_Pathway_Results/")

CC <- read.csv("GO_Cellular_Component_Serum_pval_0.01_filtered.txt", sep = '\t', header = T)
MF <- read.csv("GO_Molecular_Function_Serum_pval_0.01_filtered.txt", sep = '\t', header = T)
BP <- read.csv("GO_Biological_Process_Serum_pval_0.01_sumrati_filtered.txt", sep = '\t', header = T)

#### Bubble plot ####
library(viridis)
library(ggrepel)
library(ggplot2)

CC <- cbind(CC,'CC')
MF <- cbind(MF,'MF')
BP <- cbind(BP,'BP')

colnames(CC)[10] <- 'Category'
colnames(BP)[10] <- 'Category'
colnames(MF)[10] <- 'Category'

data <- rbind(CC,MF,BP)

colnames(CC)

ggplot(data, aes(Term,-log(P.value),label=Term))+
  geom_point(aes(color = Category, size = log(Number_of_genes), alpha=0.5)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "violetred1", "darkslateblue", "blueviolet", "chartreuse", "navyblue", "hotpink4", "salmon", "mediumaquamarine"),
                     labels=c("Biological Process","Cellular Component","Molecular Function")) +
  scale_size(range = c(0.5, 12)) + facet_wrap(~Category,  ncol = 5, scales = "free") + 
  geom_text_repel(data=subset(data, -log(P.value) > 10),aes(Term,-log(P.value),label=Term), size=3)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())

#### Dotplot ####

pathway <- read.csv("Pathway_Serum_fdr_1_Sumrati_filtered.txt", sep = '\t', header = T)

colnames(pathway)

dotplot <- ggplot(pathway, aes(x=Entities.pValue, y=Pathway.name)) +
  geom_dotplot(binaxis='y', stackdir='center', stackratio=1, dotsize=0.2) +
  theme(axis.text.x = element_text(angle = 45,  hjust = 1), plot.title = element_text(size=10))

plot(dotplot)

ggplot(pathway, aes(x=Entities.ratio, y=reorder(Pathway.name, Entities.ratio))) +
  geom_point(aes(size = X.Entities.found, color = Entities.FDR)) +
  theme_bw(base_size = 15) +
  scale_colour_gradient(limits=c(0, 0.05), low="red", high = "blue") +
  ylab(NULL) +
  ggtitle("Reactome pathway enrichment")

library(ReactomePA)
enrichMap(pathway, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
dotplot(pathway, showCategory=15)
