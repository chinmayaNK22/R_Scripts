setwd("D:\\Sumrati_PhD/ROC/")

library(ggplot2)
library(dplyr)

data <- read.csv("ELISA_All_Proteins_HIE3_for_ROC.txt", sep = "\t", header = T)

Up_data <- read.csv("ELISA_UpRegulated_Proteins_AGT_FABP1_HIE3_for_ROC.txt", sep = "\t", header = T)
Down_data <- read.csv("ELISA_DownRegulated_Proteins_APP_FN1_HIE3_for_ROC.txt", sep = "\t", header = T)


####pROC####
library(pROC)
## Type 'citation("pROC")' for a citation.
## 
## Attaching package: 'pROC'
## The following objects are masked from 'package:stats':
## 
##     cov, smooth, var
pROC_obj <- roc(data$Label,data$Concentration,
                 smoothed = TRUE,
                 # arguments for ci
                 ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                 # arguments for plot
                 plot=TRUE, auc.polygon=TRUE, max.auc.polygon=F, grid=F,
                 print.auc=TRUE, show.thres=TRUE)

par(mfrow= c(1,2))
pROC_Up <- roc(Up_data$Label,Up_data$Concentration,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=F, grid=F,
                print.auc=TRUE, show.thres=TRUE, main = "Up Regulated")

pROC_Down <- roc(Down_data$Label,Down_data$Concentration,
                 smoothed = TRUE,
                 # arguments for ci
                 ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                 # arguments for plot
                 plot=TRUE, auc.polygon=TRUE, max.auc.polygon=F, grid=F,
                 print.auc=TRUE, show.thres=TRUE, main = "Down Regulated")

