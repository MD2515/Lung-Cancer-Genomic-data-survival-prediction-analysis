library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
#install.packages("ggfortify")
library(ggfortify)


#read file
methylation <- read.csv("D:/BHI/Spring 2024/BHI 699/lung/sfd.csv", header = TRUE)
#Name rownames cseid
row.names(methylation) <- methylation[,1]
#Delete caseid column for analysis
methylation <- methylation[, -1]

os <- with(methylation, Surv(methylation$time, methylation$status))
os_fit <- survfit(Surv(methylation$time, methylation$status) ~1, data = methylation)
autoplot(os_fit)

#survival fit for specific gene
os_gene_fit <- survfit(Surv(methylation$time, methylation$status) ~ methylation$gene, data = methylation)
autoplot(os_gene_fit)

#Cox proportional hazard 
cox <- coxph(Surv(methylation$time, methylation$status) ~ ., data = methylation)
summary(cox)

# put cox summary in object
cox_summary <- summary(cox)

# Filter important variables based on significance level (e.g., p-value < 0.05)
significant_variables <- cox_summary$coefficients[cox_summary$coefficients[, "Pr(>|z|)"] < 0.05, ]

# Print the summary of important variables
print(significant_variables)
