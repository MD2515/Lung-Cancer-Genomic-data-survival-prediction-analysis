#install.packages("survival")
#install.packages("survminer")
#install.packages("DMwR2")

#update packages
#install.packages("stringi", dependencies = TRUE, INSTALL_opts = c('--no-lock'))
#install.packages("dendextend", dependencies = TRUE, INSTALL_opts = c('--no-lock'))
#install.packages("ps", dependencies = TRUE, INSTALL_opts = c('--no-lock'))


#load the bio-manager installer
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()


#install TCGA package
#BiocManager::install("RTCGA")

#install clinical TCGA package
#BiocManager::install("RTCGA.clinical")


#install mRNA TCGA package
#BiocManager::install("RTCGA.mRNA")


#load libraries
library(survival)
library(survminer)
library(RTCGA)
library(RTCGA.clinical)
library(ggplot2)
library(tidyverse)
library(DMwR2)
library(zoo)
library(mice)

#display informaation on TCGA packaage
infoTCGA()

#display informaation on the clinical daata of BRCA
browseVignettes("RTCGA.clinical")


#####################################################################################################################################################################################
#####################################Methylation#################################################################################################################################
#####################################################################################################################################################################################

#Clinical data
#Either Download clinical data from here(http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LUAD/20160128/gdac.broadinstitute.org_LUAD.Merge_Clinical.Level_1.2016012800.0.0.tar.gz)
#use LUAD-TP.merged_data.txt data file.
# Or follow steps below to retrieve it by R
clinical <- survivalTCGA(LUAD.clinical,
                         extract.cols="admin.disease_code")
clinical <- clinical[, -4]
colnames(clinical)[2] <- 'caseid'
colnames(clinical)[1] <- 'time'
colnames(clinical)[3] <- 'status'

clinical$time <- as.numeric(clinical$time)/30
clinical[,'time']=round(clinical[,'time'],2)
head(clinical)



#Methylation data
#Download methylation data from here (http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LUAD/20160128/gdac.broadinstitute.org_LUAD.Methylation_Preprocess.Level_3.2016012800.0.0.tar.gz)
#use LUAD-TP.meth.by_min_clin_corr.data.txt data file.
methylation <- read.table("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/methylation/luadmeth.txt", sep = "\t")
methylation2 <- rbind(colnames(methylation), methylation)
methylation2 = t(methylation2)
methylation2 <- methylation2[, -1]
methylation2 <- methylation2[,-2]
colnames(methylation2) <- methylation2[1, ]
methylation2 <- methylation2[-1, ]

# Check the class of methylation2
class(methylation2)

# If methylation2 is not a data frame, convert it into one
methylation2 <- as.data.frame(methylation2)

colnames(methylation2)[1] <- 'caseid'

#Remove the normal sample data
samples_normal <- grepl("\\-11$", methylation2$caseid)

methylation2 <- methylation2[!samples_normal, ]
# Remove .01 from sample ids for id match
methylation2$caseid <- sub('\\-01$', '', methylation2$caseid)

#Replace 0 with NA
methylation2[methylation2 == 0] <- NA


#Replace Missing values with KNN
#methylation2 <- knnImputation(methylation2)

#remove columns of genes which has 0 or na values
methylation2 <- methylation2[, colSums(is.na(methylation2)) ==0]

#impute missing data using MICE function
#methylation2 <- complete(mice(as.data.frame(methylation2)))

############Merge data sets################################
methylation2 <- merge(clinical, methylation2, by = 'caseid', all = FALSE)

#set row names as case ids
row.names(methylation2) <- methylation2[,1]
#remove case id column
methylation2 <- methylation2[, -1]
# Delete all the rows that has time less than a month
methylation2 <- methylation2[methylation2$time >=1.0, ]


#Save data frame
write.csv(methylation2, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/Methylation/amlaf.csv", row.names = TRUE)


#####################################################################################################################################################################################
#####################################Gene Expression#################################################################################################################################
#####################################################################################################################################################################################

#Clinical data
#Either Download clinical data from here(http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LUAD/20160128/gdac.broadinstitute.org_LUAD.Merge_Clinical.Level_1.2016012800.0.0.tar.gz)
#use LUAD-TP.merged_data.txt data file.
# Or follow steps below to retrieve it by R
clinical <- survivalTCGA(LUAD.clinical,
                         extract.cols="admin.disease_code")
clinical <- clinical[, -4]
colnames(clinical)[2] <- 'caseid'
colnames(clinical)[1] <- 'time'
colnames(clinical)[3] <- 'status'

clinical$time <- as.numeric(clinical$time)/30
clinical[,'time']=round(clinical[,'time'],2)

head(clinical)



#Gene Expression data
#Download Gene Expression data from here (http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LUAD/20160128/gdac.broadinstitute.org_LUAD.mRNAseq_Preprocess.Level_3.2016012800.0.0.tar.gz)
#use LUAD-TP.uncv2.mRNAseq_RSEM_normalized_log2.txt file.
gene2 <- read.table("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/gene expression/luadgene.txt", sep = "\t")
gene2 <- rbind(colnames(gene2), gene2)
gene3 = t(gene2)

gene3 <- gene3[,-1]
colnames(gene3) <- gene3[1, ]
colnames(gene3)[1] <- 'caseid'
gene3 <- gene3[-1, ]

# Check the class of methylation2
class(gene3)

# If methylation2 is not a data frame, convert it into one
gene3 <- as.data.frame(gene3)

#Remove normal sample data
samples_normal <- grepl("\\-11$", gene3$caseid)

gene3 <- gene3[!samples_normal, ]


# Now, try to access the column using the $ operator
gene3$caseid <- sub('\\-01$', '', gene3$caseid)

#Replace 0 with NA
gene3[gene3 == 0] <- NA


#Replace Missing values with KNN
#gene3 <- knnImputation(gene3)

#remove columns of genes which has 0 or na values
gene3 <- gene3[, colSums(is.na(gene3)) ==0]

#impute missing data using MICE function
#gene3 <- complete(mice(as.data.frame(gene3)))

############Merge data sets################################
gene3 <- merge(clinical, gene3, by = 'caseid', all = FALSE)

#set row names as case ids
row.names(gene3) <- gene3[,1]
#remove case id column
gene3 <- gene3[, -1]
# Delete all the rows that has time less than a month
gene3 <- gene3[gene3$time >=1.0, ]


#Save data frame
write.csv(gene3, "D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/gene expression/aglaf.csv", row.names = TRUE)

