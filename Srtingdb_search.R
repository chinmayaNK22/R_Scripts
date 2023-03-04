library(STRINGdb)
STRINGdb$methods()

#### Specify taxonomy ID from STRINGdb for species ####
string_db <- STRINGdb$new( version="11", species=83332, score_threshold=400, input_directory="E:\\MTB_Multi-PTM/")

#### StringDB_Annotation ####
#annotations <- string_db$get_annotations(proteins_vector)

#### Fetch the list of proteins or genes for STRINGdb search ####
setwd("E:\\MTB_Multi-PTM/Reanalysis_Results/PTM_Summary/")
proteins <- read.csv("STRINGdb/For_StringDB_Search.txt", sep = "\t", header = T)
colnames(proteins)

for(i in 1:dim(proteins)[2]){
  j <- paste0(proteins[1,i],"090420.txt")
  if(file.exists(j))
  {
  }
  else
  {
    file.create(j)
  }
  enrichment <- string_db$get_enrichment(proteins[,i])
  write.table(enrichment, file = j, sep = "\t", append = TRUE, quote = FALSE,col.names = FALSE, row.names = FALSE)
}

#Mod <- as.data.frame(Acetyl)
#proteins_vector <- as.vector(Mod)
#proteins_vector <- unlist(proteins_vector)
#head(proteins_vector)

Crotonyl <- proteins$Crotonylation

#### StringDB_enrichment ####
Crotonyl <- string_db$get_enrichment(Crotonyl)
Crotonyl$Modification <- rep("Crotonylation",each=dim(Crotonyl)[1])

#Crotonyl_CC <- Crotonyl[Crotonyl$category == "Component", ]


#write.table(enrichment,file = "STRINGdb_enrichment.txt", sep = "\t")

#### plot PPI network ####
string_db$plot_network(proteins_vector)

####
setwd("STRINGdb/")
Acetylation <- read.csv("Acetylation090420.txt", sep = '\t', header = T)
Butyryl <- read.csv("Butyrylation090420.txt", sep = '\t', header = T)
Citrullin <- read.csv("Citrullination090420.txt", sep = '\t', header = T)
Crotonyl <- read.csv("Crotonylation090420.txt", sep = '\t', header = T)
Glutaryl <- read.csv("Glutarylation090420.txt", sep = '\t', header = T)
Glycosyl <- read.csv("Glycosylation090420.txt", sep = '\t', header = T)
Hex <- read.csv("Hex090420.txt", sep = '\t', header = T)
Malonyl <- read.csv("Malonylation090420.txt", sep = '\t', header = T)
Methyl <- read.csv("Methylation090420.txt", sep = '\t', header = T)
Phospho <- read.csv("Phoshorylation090420.txt", sep = '\t', header = T)

GO_term <- "Component"
Acetyl_GO <- Acetylation[Acetylation$Category == GO_term, ]
Butyryl_GO <- Butyryl[Butyryl$Category == GO_term, ]
Citrullin_GO <- Citrullin[Citrullin$Category == GO_term, ]
Crotonyl_GO <- Crotonyl[Crotonyl$Category == GO_term, ]
Glutaryl_GO <- Glutaryl[Glutaryl$Category == GO_term, ]
Glycosyl_GO <- Glycosyl[Glycosyl$Category == GO_term, ]
Hex_GO <- Hex[Hex$Category == GO_term, ]
Methyl_GO <- Methyl[Methyl$Category == GO_term, ]
Malonyl_GO <- Malonyl[Malonyl$Category == GO_term, ]
Phospho_GO <- Phospho[Phospho$Category == GO_term, ]

#### Bubble plot ####
library(viridis)
library(ggrepel)

#String_enrich <- rbind(enrichment[1:41,],enrichment[98:155,], enrichment[42:53,])
#String_enrich$term<-as.factor(String_enrich$term)
#str(String_enrich$term)
#str(String_enrich)

data <- rbind(Acetyl_GO,Butyryl_GO,Citrullin_GO,Crotonyl_GO,Glutaryl_GO,Glycosyl_GO,Hex_GO,Malonyl_GO,Methyl_GO,Phospho_GO)

ggplot(data, aes(Term,-log(p.value),label=Term))+
  geom_point(aes(color = Modification, size = log(Number_of_genes), alpha=0.5)) +
    scale_color_manual(values = c("#00AFBB", "#E7B800", "darkslateblue", "blueviolet", "chartreuse", "violetred1", "navyblue", "hotpink4", "salmon", "mediumaquamarine"),
    labels=c("Acetylation","Butyrylation","Citrullination","Crotonylation","Glutarylation","Glycosylation","Hex","Malonylation","Methylation","Phosphorylation")) +
    scale_size(range = c(0.5, 12))+
    facet_wrap(~Modification,  ncol = 5, scales = "free") +
geom_text_repel(data=subset(data, -log(p.value) > 10),aes(Term,-log(p.value),label=Term), size=3)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())

summary <- string_db$get_summary(proteins_vector)
all <- string_db$load()
plot(all)


#### Dotplot ####
CC <- enrichment[1:15,]
dotplot <- ggplot(String_enrich, aes(x=category, y=term, fill=category)) +
  geom_dotplot(binaxis='y', stackdir='center', stackratio=1, dotsize=0.2) + coord_flip()+
  theme(axis.text.x = element_text(angle = 45,  hjust = 1), plot.title = element_text(size=10))

plot(dotplot)
