library('pathfindR')
library('tidyr')

PPI <- read.csv('SARS-CoV-2_PPI_BioGrid.txt', sep = '\t', header = F)
colnames(PPI)

PPI.df <- data.frame(PPI)
Log2FC.df <- as.data.frame(log2(PPI.df[,2]))
PPI.final.df <- cbind(PPI.df[,1], Log2FC.df, PPI.df[,3])
colnames(PPI.final.df)<- c('Gene_symbol', 'log2FC', 'pvalue')

PPI.output.df <- run_pathfindR(PPI.final.df, gene_sets = "Reactome", pin_name_path = "Biogrid", list_active_snw_genes = TRUE)
