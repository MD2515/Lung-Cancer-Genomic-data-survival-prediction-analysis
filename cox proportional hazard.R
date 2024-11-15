library(survival)
library(survminer)
library(tidyr)
#install.packages("glmnet")
#install.packages("mice")
library(glmnet)
library(mice)


##################################################LUAD Methylation#####################################################

luadmethylation <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis documents/Data/lung/LUAD/methylation/amlaf.csv")

#Name rownames cseid
row.names(luadmethylation) <- luadmethylation[,1]
#Delete caseid column for analysis
luadmethylation <- luadmethylation[, -1]

set.seed(42)
# creating Surv object
surv_obj <- with(luadmethylation, Surv(time, status))

# only gene data bhy removing time and status
gene_data <- luadmethylation[, -(1:2)]  

# Impute missing values
#imputed_data <- complete(mice(as.data.frame(gene_data)))

# Cox model with LASSO
lasso_model <- glmnet(as.matrix(gene_data), surv_obj, family = "cox")

# lamda using cross validation
cv_fit <- cv.glmnet(as.matrix(gene_data), surv_obj, family = "cox")
best_lambda <- cv_fit$lambda.min

# coefficients using lamda
lasso_coefs <- coef(lasso_model, s = best_lambda)

# coefficients zero
selected_features <- rownames(lasso_coefs)[which(lasso_coefs != 0)]

# data subset with features
selected_data <- gene_data[, selected_features]

# cox model with feature
cox_model_selected <- coxph(surv_obj ~ ., data = selected_data)

# summary
summary(cox_model_selected)


# summary in object
cox_summary <- summary(cox_model_selected)

# variable that has p value <0.05
significant_variables <- cox_summary$coefficients[cox_summary$coefficients[, "Pr(>|z|)"] < 0.05, ]

# significant variable
print(significant_variables)



######################################################LUAD Gene Expression##################################################
luadgene <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis documents/Data/lung/LUAD/gene expression/aglaf.csv")

#Name rownames cseid
row.names(luadgene) <- luadgene[,1]
#Delete caseid column for analysis
luadgene <- luadgene[, -1]

set.seed(42)
# creating Surv object
surv_obj <- with(luadgene, Surv(time, status))

# only gene data bhy removing time and status
gene_data <- luadgene[, -(1:2)]  

# Impute missing values
#imputed_data <- complete(mice(as.data.frame(gene_data)))

# Cox model with LASSO
lasso_model <- glmnet(as.matrix(gene_data), surv_obj, family = "cox")

# lamda using cross validation
cv_fit <- cv.glmnet(as.matrix(gene_data), surv_obj, family = "cox")
best_lambda <- cv_fit$lambda.min

# coefficients using lamda
lasso_coefs <- coef(lasso_model, s = best_lambda)

# coefficients zero
selected_features <- rownames(lasso_coefs)[which(lasso_coefs != 0)]

# data subset with features
selected_data <- gene_data[, selected_features]

# cox model with feature
cox_model_selected <- coxph(surv_obj ~ ., data = selected_data)

# summary
summary(cox_model_selected)


# summary in object
cox_summary <- summary(cox_model_selected)

# variable that has p value <0.05
significant_variables <- cox_summary$coefficients[cox_summary$coefficients[, "Pr(>|z|)"] < 0.05, ]

# significant variable
print(significant_variables)



##############################################LUSC Methylation##############################################################
luscmethylation <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis documents/Data/lung/LUSC/methylation/smlaf.csv")

#Name rownames cseid
row.names(luscmethylation) <- luscmethylation[,1]
#Delete caseid column for analysis
luscmethylation <- luscmethylation[, -1]

set.seed(42)
# creating Surv object
surv_obj <- with(luscmethylation, Surv(time, status))

# only gene data bhy removing time and status
gene_data <- luscmethylation[, -(1:2)]  

# Impute missing values
#imputed_data <- complete(mice(as.data.frame(gene_data)))

# Cox model with LASSO
lasso_model <- glmnet(as.matrix(gene_data), surv_obj, family = "cox")

# lamda using cross validation
cv_fit <- cv.glmnet(as.matrix(gene_data), surv_obj, family = "cox")
best_lambda <- cv_fit$lambda.min

# coefficients using lamda
lasso_coefs <- coef(lasso_model, s = best_lambda)

# coefficients zero
selected_features <- rownames(lasso_coefs)[which(lasso_coefs != 0)]

# data subset with features
selected_data <- gene_data[, selected_features]

# cox model with feature
cox_model_selected <- coxph(surv_obj ~ ., data = selected_data)

# summary
summary(cox_model_selected)


# summary in object
cox_summary <- summary(cox_model_selected)

# variable that has p value <0.05
significant_variables <- cox_summary$coefficients[cox_summary$coefficients[, "Pr(>|z|)"] < 0.05, ]

# significant variable
print(significant_variables)



#######################################################LUSC Gene Expression################################################
luscgene <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis documents/Data/lung/LUSC/gene expression/sglaf.csv")

#Name rownames cseid
row.names(luscgene) <- luscgene[,1]
#Delete caseid column for analysis
luscgene <- luscgene[, -1]

set.seed(42)
# creating Surv object
surv_obj <- with(luscgene, Surv(time, status))

# only gene data bhy removing time and status
gene_data <- luscgene[, -(1:2)]  

# Impute missing values
#imputed_data <- complete(mice(as.data.frame(gene_data)))

# Cox model with LASSO
lasso_model <- glmnet(as.matrix(gene_data), surv_obj, family = "cox")

# lamda using cross validation
cv_fit <- cv.glmnet(as.matrix(gene_data), surv_obj, family = "cox")
best_lambda <- cv_fit$lambda.min

# coefficients using lamda
lasso_coefs <- coef(lasso_model, s = best_lambda)

# coefficients zero
selected_features <- rownames(lasso_coefs)[which(lasso_coefs != 0)]

# data subset with features
selected_data <- gene_data[, selected_features]

# cox model with feature
cox_model_selected <- coxph(surv_obj ~ ., data = selected_data)

# summary
summary(cox_model_selected)


# summary in object
cox_summary <- summary(cox_model_selected)

# variable that has p value <0.05
significant_variables <- cox_summary$coefficients[cox_summary$coefficients[, "Pr(>|z|)"] < 0.05, ]

# significant variable
print(significant_variables)






##################################################################################################################################
###############################################################################################################################
###############################################Significant feature data#######################################################

#Significant feature data LUAD methylation
luadsigmeth <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/methylation/amsfd.csv")

#Name rownames cseid
row.names(luadsigmeth) <- luadsigmeth[,1]
#Delete caseid column for analysis
luadsigmeth <- luadsigmeth[, -1]

set.seed(42)
# creating Surv object
surv_obj <- with(luadsigmeth, Surv(time, status))

# only gene data bhy removing time and status
gene_data <- luadsigmeth[, -(1:2)]  

# Impute missing values
#imputed_data <- complete(mice(as.data.frame(gene_data)))

# Cox model with LASSO
lasso_model <- glmnet(as.matrix(gene_data), surv_obj, family = "cox")

# lamda using cross validation
cv_fit <- cv.glmnet(as.matrix(gene_data), surv_obj, family = "cox")
best_lambda <- cv_fit$lambda.min

# coefficients using lamda
lasso_coefs <- coef(lasso_model, s = best_lambda)

# coefficients zero
selected_features <- rownames(lasso_coefs)[which(lasso_coefs != 0)]

# data subset with features
selected_data <- gene_data[, selected_features]

# cox model with feature
cox_model_selected <- coxph(surv_obj ~ ., data = selected_data)

# summary
summary(cox_model_selected)


# summary in object
cox_summary <- summary(cox_model_selected)

#variable that has p value <0.05
significant_variables <- cox_summary$coefficients[cox_summary$coefficients[, "Pr(>|z|)"] < 0.05, ]

# significant variable
print(significant_variables)




#Significant feature data LUAD gene expression
luadsiggene <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/gene expression/agsfd.csv")

#Name rownames cseid
row.names(luadsiggene) <- luadsiggene[,1]
#Delete caseid column for analysis
luadsiggene <- luadsiggene[, -1]

set.seed(42)
# creating Surv object
surv_obj <- with(luadsiggene, Surv(time, status))

# only gene data bhy removing time and status
gene_data <- luadsiggene[, -(1:2)]  

# Impute missing values
#imputed_data <- complete(mice(as.data.frame(gene_data)))

# Cox model with LASSO
lasso_model <- glmnet(as.matrix(gene_data), surv_obj, family = "cox")

# lamda using cross validation
cv_fit <- cv.glmnet(as.matrix(gene_data), surv_obj, family = "cox")
best_lambda <- cv_fit$lambda.min

# coefficients using lamda
lasso_coefs <- coef(lasso_model, s = best_lambda)

# coefficients zero
selected_features <- rownames(lasso_coefs)[which(lasso_coefs != 0)]

# data subset with features
selected_data <- gene_data[, selected_features]

# cox model with feature
cox_model_selected <- coxph(surv_obj ~ ., data = selected_data)

# summary
summary(cox_model_selected)


# summary in object
cox_summary <- summary(cox_model_selected)

# variable that has p value <0.05
significant_variables <- cox_summary$coefficients[cox_summary$coefficients[, "Pr(>|z|)"] < 0.05, ]

# significant variable
print(significant_variables)






#Significant feature data LUSC methylation
luscsigmeth <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/methylation/smsfd.csv")

#Name rownames cseid
row.names(luscsigmeth) <- luscsigmeth[,1]
#Delete caseid column for analysis
luscsigmeth <- luscsigmeth[, -1]

set.seed(42)
# creating Surv object
surv_obj <- with(luscsigmeth, Surv(time, status))

# only gene data bhy removing time and status
gene_data <- luscsigmeth[, -(1:2)]  

# Impute missing values
#imputed_data <- complete(mice(as.data.frame(gene_data)))

# Cox model with LASSO
lasso_model <- glmnet(as.matrix(gene_data), surv_obj, family = "cox")

# lamda using cross validation
cv_fit <- cv.glmnet(as.matrix(gene_data), surv_obj, family = "cox")
best_lambda <- cv_fit$lambda.min

# coefficients using lamda
lasso_coefs <- coef(lasso_model, s = best_lambda)

# coefficients zero
selected_features <- rownames(lasso_coefs)[which(lasso_coefs != 0)]

# data subset with features
selected_data <- gene_data[, selected_features]

# cox model with feature
cox_model_selected <- coxph(surv_obj ~ ., data = selected_data)

# summary
summary(cox_model_selected)


# summary in object
cox_summary <- summary(cox_model_selected)

# variable that has p value <0.05
significant_variables <- cox_summary$coefficients[cox_summary$coefficients[, "Pr(>|z|)"] < 0.05, ]

# significant variable
print(significant_variables)



#Significant feature data LUSC gene expression
luscsiggene <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/gene expression/sgsfd.csv")

#Name rownames cseid
row.names(luscsiggene) <- luscsiggene[,1]
#Delete caseid column for analysis
luscsiggene <- luscsiggene[, -1]

set.seed(42)
# creating Surv object
surv_obj <- with(luscsiggene, Surv(time, status))

# only gene data bhy removing time and status
gene_data <- luscsiggene[, -(1:2)]  

# Impute missing values
#imputed_data <- complete(mice(as.data.frame(gene_data)))

# Cox model with LASSO
lasso_model <- glmnet(as.matrix(gene_data), surv_obj, family = "cox")

# lamda using cross validation
cv_fit <- cv.glmnet(as.matrix(gene_data), surv_obj, family = "cox")
best_lambda <- cv_fit$lambda.min

# coefficients using lamda
lasso_coefs <- coef(lasso_model, s = best_lambda)

# coefficients zero
selected_features <- rownames(lasso_coefs)[which(lasso_coefs != 0)]

# data subset with features
selected_data <- gene_data[, selected_features]

# cox model with feature
cox_model_selected <- coxph(surv_obj ~ ., data = selected_data)

# summary
summary(cox_model_selected)


# summary in object
cox_summary <- summary(cox_model_selected)

# variable that has p value <0.05
significant_variables <- cox_summary$coefficients[cox_summary$coefficients[, "Pr(>|z|)"] < 0.05, ]

# significant variable
print(significant_variables)
