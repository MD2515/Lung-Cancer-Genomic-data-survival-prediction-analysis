library(dplyr)
#if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("survcomp")
library(survcomp)
library(survival)
#importFrom(R6,R6Class)
#install.packages(R6)
library(R6)
#BiocManager::install("RegParallel", force = TRUE)
library(RegParallel)
#importFrom(caret,createFolds)
#importFrom(caret,train)
library(caret)
#importFrom(data.table,rbindlist)
#BiocManager::install('data.table')
library(data.table)
#importFrom(textshape,column_to_rownames)
#install.packages('textshape')
library(textshape)
#install.packages('ggplot2')
library(ggplot2)
#install.packages('tictoc')
library(tictoc)
#install.packages("randomForestSRC")
library(randomForestSRC)
#install.packages("rms")
library(rms)
#install.packages("glmnet")
library(glmnet)
#install.packages("ranger")
library(ranger)
#install.packages("gbm")
library(gbm)

#AllinOne <- CBRCONF_AllinOne(GeneSub, Num = 10, Percent = 0.5, Cores = 1, blocksize = 200, n.folds = 10, KM.Plot = FALSE)

knitr::opts_chunk$set(echo = TRUE)
#install.packages("C:/Users/mdesai/Desktop/699/FLAIRS 2020/R Packages/CBRCONF - New/CBRCONF_1.0.2.tar.gz", repos = NULL, type = "source")
library(CBRCONF)



##################################################################################################################################################
######################################################LUAD Methylation##########################################################
luadmethylation <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/methylation/amlaf.csv")
#Name rownames cseid
row.names(luadmethylation) <- luadmethylation[,1]
#Delete caseid column for analysis
luadmethylation <- luadmethylation[, -1]

#CBRCONF Analysis
CBRCONF <- create_CBRCONF(luadmethylation)
CBRCONF <- cox_module(CBRCONF, Num = 10, Percent = 0.5, Cores = 1)
CBRCONF <- risk_groups_module(CBRCONF, KM.Plot = FALSE)
CBRCONF <- prototype_module(CBRCONF)
CBRCONF <- prototype_distance_module(CBRCONF)

#write significant data by using CBRCONF
amsfd <- CBRCONF[["sig_feature_data"]]
write.csv(amsfd, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/methylation/amsfd.csv")

#write low risk feature dataset
amlrf <- CBRCONF[["low_risk_features"]]
write.csv(amlrf, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/methylation/amlrf.csv")

#write medium risk feature
ammrf <- CBRCONF[["medium_risk_features"]]
write.csv(ammrf, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/methylation/ammrf.csv")

#write high risk feature
amhrf <- CBRCONF[["high_risk_features"]]
write.csv(amhrf, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/methylation/amhrf.csv")


CBRCONF <- CBR_module(CBRCONF, n.folds = 10, KM.Plot = FALSE)

##################################################################################################################################
#########################################################LUAD Gene Expression###################################################################
luadgene <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/gene expression/aglaf.csv")
#Name rownames cseid
row.names(luadgene) <- luadgene[,1]
#Delete caseid column for analysis
luadgene <- luadgene[, -1]

#CBRCONF Analysis
CBRCONF <- create_CBRCONF(luadgene)
CBRCONF <- cox_module(CBRCONF, Num = 10, Percent = 0.5, Cores = 1)
CBRCONF <- risk_groups_module(CBRCONF, KM.Plot = FALSE)
CBRCONF <- prototype_module(CBRCONF)
CBRCONF <- prototype_distance_module(CBRCONF)

#write significant data by using CBRCONF
agsfd <- CBRCONF[["sig_feature_data"]]
write.csv(agsfd, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/gene expression/agsfd.csv")

#write low risk feature dataset
aglrf <- CBRCONF[["low_risk_features"]]
write.csv(aglrf, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/gene expression/aglrf.csv")

#write medium risk feature
agmrf <- CBRCONF[["medium_risk_features"]]
write.csv(agmrf, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/gene expression/agmrf.csv")

#write high risk feature
aghrf <- CBRCONF[["high_risk_features"]]
write.csv(aghrf, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/gene expression/aghrf.csv")


CBRCONF <- CBR_module(CBRCONF, n.folds = 10, KM.Plot = FALSE)

##################################################################################################################################
##################################################################################################################################

###################################################LUSC Methylation###############################################################
luscmethylation <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/Methylation/smlaf.csv")
#Name rownames cseid
row.names(luscmethylation) <- luscmethylation[,1]
#Delete caseid column for analysis
luscmethylation <- luscmethylation[, -1]

#CBRCONF Analysis
CBRCONF <- create_CBRCONF(luscmethylation)
CBRCONF <- cox_module(CBRCONF, Num = 10, Percent = 0.5, Cores = 1)
CBRCONF <- risk_groups_module(CBRCONF, KM.Plot = FALSE)
CBRCONF <- prototype_module(CBRCONF)
CBRCONF <- prototype_distance_module(CBRCONF)

#write significant data by using CBRCONF
smsfd <- CBRCONF[["sig_feature_data"]]
write.csv(smsfd, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/Methylation/smsfd.csv")

#write low risk feature dataset
smlrf <- CBRCONF[["low_risk_features"]]
write.csv(smlrf, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/Methylation/smlrf.csv")

#write medium risk feature
smmrf <- CBRCONF[["medium_risk_features"]]
write.csv(smmrf, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/Methylation/smmrf.csv")

#write high risk feature
smhrf <- CBRCONF[["high_risk_features"]]
write.csv(smhrf, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/Methylation/smhrf.csv")


CBRCONF <- CBR_module(CBRCONF, n.folds = 10, KM.Plot = FALSE)


##########################################################################################################################################
######################################################LUSC Gene Expression################################################################
luscgene <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/Gene Expression/sglaf.csv")
#Name rownames cseid
row.names(luscgene) <- luscgene[,1]
#Delete caseid column for analysis
luscgene <- luscgene[, -1]

#CBRCONF Analysis
CBRCONF <- create_CBRCONF(luscgene)
CBRCONF <- cox_module(CBRCONF, Num = 10, Percent = 0.5, Cores = 1)
CBRCONF <- risk_groups_module(CBRCONF, KM.Plot = FALSE)
CBRCONF <- prototype_module(CBRCONF)
CBRCONF <- prototype_distance_module(CBRCONF)

#write significant data by using CBRCONF
sgsfd <- CBRCONF[["sig_feature_data"]]
write.csv(sgsfd, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/Gene Expression/sgsfd.csv")

#write low risk feature dataset
sglrf <- CBRCONF[["low_risk_features"]]
write.csv(sglrf, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/Gene Expression/sglrf.csv")

#write medium risk feature
sgmrf <- CBRCONF[["medium_risk_features"]]
write.csv(sgmrf, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/Gene Expression/sgmrf.csv")

#write high risk feature
sghrf <- CBRCONF[["high_risk_features"]]
write.csv(sghrf, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/Gene Expression/sghrf.csv")


CBRCONF <- CBR_module(CBRCONF, n.folds = 10, KM.Plot = FALSE)










##########################################################################################################################################
##################################################################################################################################################
##################################################################################################################################
##########################################################Significant feature Data########################################################################


#Significant feature data LUAD Methylation
luadsigmethylation <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/methylation/amsfd.csv")
#Name rownames cseid
row.names(luadsigmethylation) <- luadsigmethylation[,1]
#Delete caseid column for analysis
luadsigmethylation <- luadsigmethylation[, -1]

#CBRCONF Analysis
CBRCONF <- create_CBRCONF(luadsigmethylation)
CBRCONF <- cox_module(CBRCONF, Num = 10, Percent = 0.5, Cores = 1)
CBRCONF <- risk_groups_module(CBRCONF, KM.Plot = FALSE)
CBRCONF <- prototype_module(CBRCONF)
CBRCONF <- prototype_distance_module(CBRCONF)

CBRCONF <- CBR_module(CBRCONF, n.folds = 10, KM.Plot = FALSE)

#Significant feature data LUAD gene expression
luadsiggene <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/gene expression/agsfd.csv")
#Name rownames cseid
row.names(luadsiggene) <- luadsiggene[,1]
#Delete caseid column for analysis
luadsiggene <- luadsiggene[, -1]

#CBRCONF Analysis
CBRCONF <- create_CBRCONF(luadsiggene)
CBRCONF <- cox_module(CBRCONF, Num = 10, Percent = 0.5, Cores = 1)
CBRCONF <- risk_groups_module(CBRCONF, KM.Plot = FALSE)
CBRCONF <- prototype_module(CBRCONF)
CBRCONF <- prototype_distance_module(CBRCONF)

CBRCONF <- CBR_module(CBRCONF, n.folds = 10, KM.Plot = FALSE)


#Significant feature data LUSC methylation
luscsigmethylation <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/Methylation/smsfd.csv")
#Name rownames cseid
row.names(luscsigmethylation) <- luscsigmethylation[,1]
#Delete caseid column for analysis
luscsigmethylation <- luscsigmethylation[, -1]

#CBRCONF Analysis
CBRCONF <- create_CBRCONF(luscsigmethylation)
CBRCONF <- cox_module(CBRCONF, Num = 10, Percent = 0.5, Cores = 1)
CBRCONF <- risk_groups_module(CBRCONF, KM.Plot = FALSE)
CBRCONF <- prototype_module(CBRCONF)
CBRCONF <- prototype_distance_module(CBRCONF)

CBRCONF <- CBR_module(CBRCONF, n.folds = 10, KM.Plot = FALSE)



#Significant feature data LUSC gene expression 
luscsiggene <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/Gene Expression/sgsfd.csv")
#Name rownames cseid
row.names(luscsiggene) <- luscsiggene[,1]
#Delete caseid column for analysis
luscsiggene <- luscsiggene[, -1]

#CBRCONF Analysis
CBRCONF <- create_CBRCONF(luscsiggene)
CBRCONF <- cox_module(CBRCONF, Num = 10, Percent = 0.5, Cores = 1)
CBRCONF <- risk_groups_module(CBRCONF, KM.Plot = FALSE)
CBRCONF <- prototype_module(CBRCONF)
CBRCONF <- prototype_distance_module(CBRCONF)

CBRCONF <- CBR_module(CBRCONF, n.folds = 10, KM.Plot = FALSE)

