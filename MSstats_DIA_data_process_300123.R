setwd("E:/PGIMER_Chandighar/")

library(MSstats)

input <- read.csv(file="PGIMER_GastricLavage_TB_Infection_Proteomics_DIA_SpecLib_search_MSStats_Input_300123_formatted.tsv", sep = '\t')
raw <- SkylinetoMSstatsFormat(input)

#MSstatsConvert::MSstatsLogsSettings(FALSE)
#ready_for_process = MSstatsPrepareForDataProcess(raw, 2, NULL)
#normalized_data = MSstatsNormalize(ready_for_process, "EQUALIZEMEDIANS")

QuantData <-dataProcess(raw, use_log_file = TRUE, featureSubset = "highQuality",
                        remove_uninformative_feature_outlier = TRUE, min_feature_count = 2, 
                        remove50missing = FALSE)

Protein_norm_abundance <- as.data.frame(QuantData$ProteinLevelData)
Feature_norm_abundance <- as.data.frame(QuantData$FeatureLevelData)

write.table(Protein_norm_abundance, file ="PGIMER_GastricLavage_TB_Infection_Proteomics_DIA_Normalized_Protein_Abundance_300123.txt", sep = '\t', row.names = F, 
            col.names = T)

### feature plots ####
# Profile plot
dataProcessPlots(data=QuantData,type="ProfilePlot")
# Quality control plot 
dataProcessPlots(data=QuantData,type="QCPlot")	
# Quantification plot for conditions
dataProcessPlots(data=QuantData,type="ConditionPlot")

####Group comparison####
head(QuantData$FeatureLevelData)
levels(QuantData$ProteinLevelData$GROUP)
Control_TB <- matrix(c(1,-1),nrow=1)

comparison <- Control_TB
row.names(comparison)<-"TB vs Control"
groups = levels(QuantData$ProteinLevelData$GROUP)
colnames(comparison) <- groups[order(as.numeric(groups))]

# Tests for differentially abundant proteins with models:
# label-based SRM experiment with expanded scope of biological replication.
testResultMultiComparisons <- groupComparison(contrast.matrix=comparison, data=QuantData,
                                           use_log_file = FALSE)
# table for result
testResultMultiComparisons$ComparisonResult

# Tests for differentially abundant proteins with models:
# label-based SRM experiment with expanded scope of biological replication.
testResultOneComparison <- groupComparison(contrast.matrix=comparison, data=QuantData,
                                           use_log_file = FALSE)
####Group comparson plots####
# Volcano plot with FDR cutoff = 0.05 and no FC cutoff
groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult, sig = 0.05, type="VolcanoPlot",
                     FCcutoff = 2, logBase.pvalue=2, ylimUp=F, ProteinName=TRUE, address="PGIMER_GL_TB_Proteomics")

# Heatmap with FDR cutoff = 0.05
groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult, type="Heatmap",
                     FCcutoff=2, logBase.pvalue=2, address="PGIMER_GL_TB_Proteomics")

# Comparison plot
groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult, type="ComparisonPlot",
                     address="PGIMER_GL_TB_Proteomics")

# normal quantile-quantile plots
modelBasedQCPlots(data=testResultMultiComparisons, type="QQPlots", address="NMO_Serum_DIA")
# residual plots
modelBasedQCPlots(data=testResultMultiComparisons, type="ResidualPlots", address="NMO_Serum_DIA")

output <- testResultMultiComparisons$ComparisonResult

write.table(output,file ="NCCS_Mouse_Secretomics_DIA_MSstats_group_comparison_output_010423.txt", sep = '\t', row.names = F, 
            col.names = T)

#### Volcano Plot ####
data <- read.csv("Indo-Dutch_Plasma_DIA_MSstats_group_comparison_output_121222.txt", sep = "\t", header = T)

input <- as.data.frame(data)
with(input, plot(log2FC, -log10(pvalue), pch=20, main="Volcano plot HIV vs Control", ylab="-log10(pvalue)"))
with(subset(input, -log10(pvalue) > 1.3 & log2FC>1), points(log2FC, -log10(pvalue), pch=25, col="green"))
with(subset(input, -log10(pvalue) > 1.3 & log2FC>1), points(log2FC, -log10(pvalue), pch=20, col="green"))
with(subset(input, -log10(pvalue) > 1.3 & log2FC<(-1)), points(log2FC, -log10(pvalue), pch=25, col="red"))
with(subset(input, -log10(pvalue) > 1.3 & log2FC<(-1)), points(log2FC, -log10(pvalue), pch=20, col="red"))
abline(h = 1.3, col = "blue", lty = 2, lwd = 1)
abline(v = c(-1,1), col = "blue", lty = 2, lwd = 1)

library(calibrate)
with(subset(input, -log10(pvalue) > 1 & log2FC>1.0), textxy(log2FC, -log10(pvalue), labs=Gene, cex=0.8, pos=4))
with(subset(input, -log10(pvalue) > 1 & log2FC<(-1.0)), textxy(log2FC, -log10(pvalue), labs=Gene, cex=0.8, pos=2))


####MSstats group comparison####
group_comparison_input = MSstatsPrepareForGroupComparison(QuantData)
levels(QuantData$ProteinLevelData$GROUP)
AFL_healthy <- matrix(c(-1,0,0,1,0,0,0),nrow=1)
ALC_healthy <- matrix(c(0,-1,0,1,0,0,0),nrow=1)
ALF_healthy <- matrix(c(0,0,-1,1,0,0,0),nrow=1)
NAFL_healthy <- matrix(c(0,0,0,1,-1,0,0),nrow=1)
NALC_healthy <- matrix(c(0,0,0,1,0,-1,0),nrow=1)
NALF_healthy <- matrix(c(0,0,0,1,0,0,-1),nrow=1)

comparison<-rbind(AFL_healthy,ALC_healthy, ALF_healthy, NAFL_healthy, NALC_healthy, NALF_healthy)
row.names(comparison)<-c("AFL_healthy","ALC_healthy","ALF_healthy","NAFL_healthy","NALC_healthy","NALF_healthy")
groups = levels(QuantData$ProteinLevelData$GROUP)
colnames(comparison) <- groups[order(as.numeric(groups))]

samples_info = getSamplesInfo(QuantData)
repeated = checkRepeatedDesign(QuantData)

group_comparison = MSstatsGroupComparison(group_comparison_input, comparison, FALSE, repeated, samples_info)

group_comparison_final = MSstatsGroupComparisonOutput(group_comparison, QuantData)

output <- group_comparison_final[["ComparisonResult"]]

write.table(output,file ="Liver_Cirrhosis_serum_DIA_MSstats_group_comparison_output_111222.txt", sep = '\t', row.names = F, 
            col.names = T)

####Group comparson plots####
# Volcano plot with FDR cutoff = 0.05 and no FC cutoff
groupComparisonPlots(data=group_comparison_final$ComparisonResult, type="VolcanoPlot",
                     FCcutoff=2, logBase.pvalue=2, ylimUp=100, ProteinName=TRUE, address="LC_Serum_DIA")

# Heatmap with FDR cutoff = 0.05
groupComparisonPlots(data=group_comparison_final$ComparisonResult, type="Heatmap",
                     FCcutoff=2, logBase.pvalue=2, address="LC_Serum_DIA")

# Comparison plot
groupComparisonPlots(data=group_comparison_final$ComparisonResult, type="ComparisonPlot",
                     address="Ex1_")
