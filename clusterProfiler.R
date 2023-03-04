library('clusterProfiler')
library("org.Hs.eg.db")
library("DOSE")
library('msigdbr')
library(tidyverse)
proteins <- read.csv("OSCC_DIA_Protein_Abundance_Quantile_Normalized_ANOVA_Tukeys_Fold_Change_Cluster_gene_symbol_052521.txt", sep = '\t', header = 1)
colnames(proteins)

class(proteins)
## Accession
FC_Upreg_Cutoff = 1.5
FC_Downreg_Cutoff = 0.666
p.value_cutoff = 0.01
Smokers_up <- proteins %>% filter(p.value_Smokers.Control < p.value_cutoff & Fold_Change..Smokers_Vs_Control. >= FC_Upreg_Cutoff)
Smokers_Down <- proteins %>% filter(p.value_Smokers.Control < p.value_cutoff & Fold_Change..Smokers_Vs_Control. >= FC_Downreg_Cutoff)
Chewers_up <- proteins %>% filter(p.value_Chewers.Control < p.value_cutoff & Fold_Change..Chewers_Vs_Control. >= FC_Upreg_Cutoff)
Chewers_Down <- proteins %>% filter(p.value_Chewers.Control < p.value_cutoff & Fold_Change..Chewers_Vs_Control. >= FC_Downreg_Cutoff)
Smokers_up.df <- as.data.frame(Smokers_up$Gene_Symbol)
Smokers_Down.df <- as.data.frame(Smokers_Down$Gene_Symbol)

Smokers_up.GeneID = bitr(Smokers_up$Gene_Symbol, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db", drop = TRUE)
Smokers_Down.GeneID = bitr(Smokers_Down$Gene_Symbol, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db", drop = TRUE)
Chewers_up.GeneID = bitr(Chewers_up$Gene_Symbol, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db", drop = TRUE)
Chewers_Down.GeneID = bitr(Chewers_Down$Gene_Symbol, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db", drop = TRUE)

b <- enrichGO(GeneID$SYMBOL, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL", pvalueCutoff = 0.01, pAdjustMethod = "BH", 
              universe, qvalueCutoff = 0.01, minGSSize = 5, maxGSSize = 500, readable = FALSE, pool = FALSE)

p1 <- emapplot(b)
plot(p1)


## MSiGdb enrichment analysis
m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame

m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% dplyr::select(gs_name, entrez_gene)
head(m_t2g)

MSiGdb <- read.gmt('msigdb.v7.1.entrez.gmt')
C6 <- read.gmt('c6.all.v7.1.entrez.gmt')
C7 <- read.gmt('c7.all.v7.1.entrez.gmt')

Immune_MSig_gmt <- enricher(GeneID$ENTREZID, pvalueCutoff = 0.01,
                 pAdjustMethod = "BH",
                 universe,
                 minGSSize = 5,
                 maxGSSize = 500,
                 qvalueCutoff = 0.01,
                 TERM2GENE = C7,
                 TERM2NAME = NA)
head(Immune_MSig_gmt)

library(enrichplot)
barplot(Immune_MSig_gmt, showCategory=20)
dotplot(Immune_MSig_gmt, showCategory=30)

## DisGeNET analysis
smokers_up_dgn <- enrichDGN(Smokers_up.GeneID$ENTREZID, pvalueCutoff = 0.01, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.01,)
Smokers_Down_dgn <- enrichDGN(Smokers_Down.GeneID$ENTREZID, pvalueCutoff = 0.01, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.01,)
Chewers_up_dgn <- enrichDGN(Chewers_up.GeneID$ENTREZID, pvalueCutoff = 0.01, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.01,)
Chewers_Down_dgn <- enrichDGN(Chewers_Down.GeneID$ENTREZID, pvalueCutoff = 0.01, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.01,)

head(dgn)

#### Barplot ####
library(cowplot)
p1 <- barplot(Chewers_up_dgn, showCategory=20) + ggtitle("Upregulated Proteins")
p2 <- barplot(Chewers_Down_dgn, showCategory=20) + ggtitle("Downregulated Proteins")
plot_grid(p1, p2, ncol=2)

#### Dotplot ####
p1 <- dotplot(Chewers_up_dgn, showCategory=20) + ggtitle("Upregulated Proteins")
p2 <- dotplot(Chewers_Down_dgn, showCategory=20) + ggtitle("Downregulated Proteins")
plot_grid(p1, p2, ncol=2)


write.table(Chewers_up_dgn, file = "OSCC_DIA_Chewers_Up_Proteins_Disease_Gene_Set_Enrichment_052521.txt", sep = '\t')
write.table(Chewers_Down_dgn, file = "OSCC_DIA_Chewers_Down_Proteins_Disease_Gene_Set_Enrichment_052521.txt", sep = '\t')

