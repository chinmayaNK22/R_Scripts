setwd("D:\\Sumrati_PhD/ROC/")

library(ggplot2)
library(dplyr)

data <- read.csv("AGT_Healthy_HIE1_24.tsv", sep = "\t", header = T)

####ROCR####
library(ROCR)
## Loading required package: gplots
## 
## Attaching package: 'gplots'
## The following object is masked from 'package:stats':
## 
##     lowess
# plot a ROC curve for a single prediction run
# and color the curve according to cutoff.
data(ROCR.simple)
df <- data.frame(ROCR.simple)
pred <- prediction(data$Concentration, data$Label)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE)

####ROCit####
library(ROCit)
## Warning: package 'ROCit' was built under R version 3.5.2
ROCit_obj <- rocit(score=data$Concentration,class=data$Label)
plot(ROCit_obj)

####precrec####
library(precrec)
## 
## Attaching package: 'precrec'
## The following object is masked from 'package:pROC':
## 
##     auc
precrec_obj <- evalmod(scores = data$Concentration, labels = data$Label)
autoplot(precrec_obj)

####PRROC####
library(PRROC)

PRROC_obj <- roc.curve(scores.class0 = data$Concentration, weights.class0=data$Label,
                       curve=TRUE)
plot(PRROC_obj)

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


sens.ci <- ci.se(pROC_obj)
plot(sens.ci, type="shape", col="lightblue")
## Warning in plot.ci.se(sens.ci, type = "shape", col = "lightblue"): Low
## definition shape.
plot(sens.ci, type="bars")

####Multiple ROC in single plot####
HIE1_24 = roc(data$Label, data$HIE1_24_Conc, smoothed = TRUE)
HIE2_24 = roc(data$Label, data$HIE2_24_Conc)
HIE3_24 = roc(data$Label, data$HIE3_24_Conc)
HIE1_72 = roc(data$Label, data$HIE1_72_Conc)
HIE2_72 = roc(data$Label, data$HIE2_72_Conc)
HIE3_72 = roc(data$Label, data$HIE3_72_Conc)

str(HIE1_24)

plot(HIE2_24, col = 1, lty = 2, main = "AGT")

roc_plot <- plot(HIE1_24, col = 1, lty = 2, main = "AGT")
roc_plot <- plot(HIE2_24, col = 2, lty = 3, add = TRUE)
roc_plot <- plot(HIE3_24, col = 3, lty = 1, add = TRUE)
roc_plot <- plot(HIE1_72, col = 4, lty = 4, add = TRUE)
roc_plot <- plot(HIE2_72, col = 5, lty = 5, add = TRUE)
roc_plot <- plot(HIE3_72, col = 6, lty = 6, add = TRUE)

lines(smooth(roc_plot), method = "density")

####ModelGood####
library(ModelGood)

# generate som data
set.seed(40)
N=40
Y=rbinom(N,1,.5)
X1=rnorm(N)
X1[Y==1]=rnorm(sum(Y==1),mean=rbinom(sum(Y==1),1,.5))
X2=rnorm(N)
X2[Y==0]=rnorm(sum(Y==0),mean=rbinom(sum(Y==0),1,.5))
dat <- data.frame(Y=Y,X1=X1,X2=X2)

# fit two logistic regression models
lm1 <- glm(data$Label~data$HIE_1_24_Conc,data=data,family="binomial")
lm2 <- glm(data$Label~data$HIE_2_24_Conc+data$HIE_1_24_Conc,data=data,family="binomial")
plot(Roc(list(lm1,lm2),data=data, legend=T))

colnames(data)

####plotROC####
library(plotROC)

longtest <- melt_roc(data, "Label", c("HIE_1_24_Conc", "HIE_2_24_Conc", "HIE_3_24_Conc", 
                                      "HIE_1_72_Conc", "HIE_2_72_Conc", "HIE_3_72_Conc"))
head(longtest)
ggplot(longtest, aes(d = D, m = M, color = name)) + geom_roc(n.cuts = 0) + style_roc()

test <- melt_roc(data, "Label", c("HIE_3_72_Conc"))
gplot(test, aes(d = D, m = M, color = name)) + geom_roc() + style_roc()

####ROCR####
library(ROCR)
data(ROCR.hiv)
x   <- prediction(ROCR.hiv$hiv.nn$predictions, ROCR.hiv$hiv.nn$labels)
ROC <- performance(x, "tpr", "fpr")
plot(ROC, col = as.list(1:10))