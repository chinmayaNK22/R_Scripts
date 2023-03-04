library(ggplot2)
library(tidyverse)
library(ggdendro)
library(phylogram)
library(dendextend)
library(ape)
setwd('F:/NTNU_M_avium_Proteomics/Genome/Mavium_complete_genomes/')

ggd <-  read.dendrogram('Mavium_genome_mashtree_new_060421.dnd')

ggd1 <- as.ggdend(ggd)

### Plot phylogeny ###
ggplot(ggd1,labels = F) + 
  scale_y_reverse(expand = c(0.2, 0)) +
  coord_polar(theta="x")

ggd %>% set("branches_k_color", k = 4) %>% plot(main = "Default colors")


dend_data <- dendro_data(ggd, type = "triangle")
names(dend_data)
head(dend_data$segments)
head(dend_data$labels)


p <- ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 3)+
  ylim(-3, 15)
plot(p)

write.dendrogram(ggd, file = "Mavium_genome_tree_newick_new_060421.txt", append = FALSE, edges = TRUE)
ape::write.tree(my_tree, file='Mavium_genome_tree_newick.txt')


