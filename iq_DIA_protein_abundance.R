#install.packages("iq")
library("iq")

#### Set the working directory ####
setwd("E:/PGIMER_Chandighar/")

data <- read.csv("PGIMER_GastricLavage_TB_Infection_Proteomics_DIA_SpecLib_search_Report_300123.tsv_iq_input.txt", sep = '\t', header = T)

head(data)

#### Preprocess the data ####
normalized_list <- iq::preprocess(data, median_normalization = FALSE, log2_intensity_cutoff = 0, 
                                  pdf_out = "PGIMER_GL_TB_Infection_Proteomics_qc-plots_300123.pdf", pdf_width = 12, pdf_height = 8)

#### Creat protein list ####
protein_list <- create_protein_list(normalized_list)
head(protein_list)

#### Create protein table ####
protein_table <- create_protein_table(protein_list, method = "maxLFQ")
head(protein_table)

#### Write output file #####
write.table(protein_table, file="PGIMER_GL_TB_Infection_Proteomics_Protein_Abundance_values_300123.txt", sep = '\t')
