library(ggplot2)
library(scales)

df <- data.frame(
  group = c("Male", "Female", "Child"),
  value = c(25, 25, 50)
)
head(df)

setwd('F:/NTNU_M_avium_Proteomics/')

data <- read.csv('Mavium_DDA_DIA_proteins_pie_chart_input.txt', sep = '\t', header = T)

# Barplot
bp<- ggplot(data, aes(x="", y=Proteins, fill=Group))+
  geom_bar(width = 1, stat = "identity")
bp
pie <- bp + coord_polar("y", start=0)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )


pie + scale_fill_brewer("Blues") + blank_theme +
  theme(axis.text.x=element_blank())+
  geom_text(aes(y = Proteins/3 + c(0, cumsum(Proteins)[-length(Proteins)]), 
                label = percent(Proteins/100)), size=5)
