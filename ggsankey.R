# install.packages("remotes")
# remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(dplyr)
library(ggplot2)

setwd("Desktop/Yenepoya/NTNU_Mavium_Proteomics/Sankey_Plot/")

data <- read.csv("Variants_VFDB_Output.csv", sep = ",", header = T)

class(data)
colnames(data)
df <- data %>%
  make_long(Virulence.factors.class, Virulence.factors, Related.genes, 
            Protein_Accession, SNP.Single.Letter.)

ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, 
               fill = factor(node), label = node)) + 
  geom_sankey(flow.alpha = .9) +
  geom_sankey_label(size = 3, color = "white") +
  scale_fill_viridis_d(option = "C", alpha = .8)+
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
  ggtitle("VFDB")

ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, 
               fill = factor(node), label = node)) + geom_alluvial(flow.alpha = .6) +
  geom_alluvial_text(size = 3, color = "white") +
  scale_fill_viridis_d() +
  theme_alluvial(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("VFDB")
