setwd("E:/NTMs_Raw_Data/")

data <- read.csv("NTM_Db_Search_Parameters_011223_proteins.txt", sep = "\t", header = T)

head(data)

# Plotting Charts and adding a secondary axis
library(ggplot2)

data$grp <- paste(data$Strain,data$ID,data$Search.Engine)

head(data)

ggp <- ggplot(data, aes(grp, Protein.Sequences, fill=Species, group = 1)) + 
  geom_bar(aes(x=grp, y=Proteins,fill=Species), position=position_dodge(width=0.5),stat="identity",colour="#006000")+
  geom_line(stat="identity", color="red", size=2) +
  labs(title= "Proteome Coverage of NTM species", x="NTM strains", y="Proteins") +
  scale_y_continuous(sec.axis=sec_axis(~.*1, name="Protein sequences"))+
  theme(axis.text.x = element_text(angle = 90, hjust=1))
ggp

ggp + geom_line(aes(x=grp, y=Protein.Sequences), stat="identity",color="red",size=2)

class(as.character(data$Adjusted.P.value))
class(as.character(data$Term))

library(ggpubr)
bp <- ggbarplot(data, x = "grp", y = "Proteins", fill = "Search.Engine", 
                color = "white", palette = "jco", sort.val = "desc",
                sort.by.groups = TRUE, x.text.angle = 90)

bp

ggbarplot(data, x = "Term", y = "Genes.Count",
           fill = "GO",                                
           palette = "jco", 
           sorting = "descending",                       
           add = "segments",                             # Add segments from y = 0 to dots
           rotate = TRUE,                                # Rotate vertically
           group = "GO",                                # Order by groups
           dot.size = 6,                                 # Large dot size
           label = round(data$Adjusted.P.value,8),                        # Add values as dot labels
           font.label = list(color = "black", size = 9, vjust = 1, hjust = 1),               # Adjust label parameters
           ggtheme = theme_pubr())


