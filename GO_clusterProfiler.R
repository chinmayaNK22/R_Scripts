library('clusterProfiler')
library("tidyverse")
library("org.Hs.eg.db")
proteins <- read.csv("SARS-CoV-2_Proteins_Q-Normalized_FC_pvalue_gene_symbol_042320.txt", sep = '\t', header = 1)
colnames(proteins)

##Signifcant_peptides from 24h
sig_24h_up <- proteins %>% filter(proteins$p.value_24h < 0.05 & proteins$FC_24h > 1.5)
sig_24h_down <- proteins %>% filter(proteins$p.value_24h < 0.05 & proteins$FC_24h < 0.66)

#accession_24h_up <- separate(sig_24h_up, Master.Protein.Accessions, into = c("A1", "A2"), sep = ";")
accession <- data.matrix(sig_24h_up$Gene_symbol)
#write.table(accession, file = "peptide_accession.txt", sep = "\t")

a <- compareCluster(accession, fun = "enrichGO", OrgDb='org.Hs.eg.db')

eg = bitr(accession, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = TRUE)

b <- enrichGO(accession, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "BH", 
              universe, 
              qvalueCutoff = 0.2, minGSSize = 5, maxGSSize = 500, readable = FALSE, pool = FALSE)
cnetplot(b, foldChange= "geneList")
p2 <- cnetplot(b, categorySize="pvalue", foldChange="geneList")
p3 <- cnetplot(b, foldChange="geneList", circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
#goplot(b)
plot(p3)
