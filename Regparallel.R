library(RegParallel)
dataset <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis documents/data/lung/luad/methylation/amlaf.csv")
row.names(dataset) <- dataset[,1]
#Delete caseid column for analysis
dataset <- dataset[, -1]
#Running a multivariate cox regression to test each feature for significance
res <- RegParallel(
  data = dataset,
  formula = 'Surv(time, status) ~ [*]',
   FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'efron',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(dataset)[3:ncol(dataset)],
  blocksize = ncol(dataset) / 4,
  cores = 2,
  nestedParallel = FALSE,
  conflevel = 95)
  
res <- res[!is.na(res$P),]
res <- res[order(res$LogRank, decreasing = FALSE),]
summary(res)
# Getting only the features with a log rank less than 0.01
final <- subset(res, LogRank < 0.01)
print(final$P)  
  