library(e1071)

###############################################LUAD#########################################################################

#LUAD Methylation
luadmeth <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/methylation/amsfd.csv")
row.names(luadmeth) <- luadmeth[,1]
luadmeth <- luadmeth[,-1]


# Split the dataset into 70% and 30% in to training and testing sets
set.seed(42)
train_indices <- sample(1:nrow(luadmeth), 0.7 * nrow(luadmeth))  
train_data <- luadmeth[train_indices, ]
test_data <- luadmeth[-train_indices, ]

# Prepare training and testing data
X_train <- train_data[, !(names(train_data) %in% c("time", "status"))]
y_train <- as.factor(train_data$status)        

X_test <- test_data[, !(names(test_data) %in% c("time", "status"))]
y_test <- as.factor(test_data$status)          

# training SVM classifier
svm_model <- svm(x = X_train, y = y_train)

# Summary
summary(svm_model)

# perfoem on test data
predictions <- predict(svm_model, newdata = X_test)

# acuuracy
accuracy <- mean(predictions == y_test)
cat("Accuracy:", accuracy, "\n")


#LUAD Gene Expression
luadgene <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/gene expression/agsfd.csv")
row.names(luadgene) <- luadgene[,1]
luadgene <- luadgene[,-1]


# Split the dataset into 70% and 30% in to training and testing sets
set.seed(42)
train_indices <- sample(1:nrow(luadgene), 0.7 * nrow(luadgene))  
train_data <- luadgene[train_indices, ]
test_data <- luadgene[-train_indices, ]

# Prepare training and testing data
X_train <- train_data[, !(names(train_data) %in% c("time", "status"))]
y_train <- as.factor(train_data$status)        

X_test <- test_data[, !(names(test_data) %in% c("time", "status"))]
y_test <- as.factor(test_data$status)          

# training SVM classifier
svm_model <- svm(x = X_train, y = y_train)

# Summary
summary(svm_model)

# perfoem on test data
predictions <- predict(svm_model, newdata = X_test)

# acuuracy
accuracy <- mean(predictions == y_test)
cat("Accuracy:", accuracy, "\n")




################################################################LUSC########################################################

#LUSC Methylation
luscmeth <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/methylation/smsfd.csv")
row.names(luscmeth) <- luscmeth[,1]
luscmeth <- luscmeth[,-1]


# Split the dataset into 70% and 30% in to training and testing sets
set.seed(42)
train_indices <- sample(1:nrow(luscmeth), 0.7 * nrow(luscmeth))  
train_data <- luscmeth[train_indices, ]
test_data <- luscmeth[-train_indices, ]

# Prepare training and testing data
X_train <- train_data[, !(names(train_data) %in% c("time", "status"))]
y_train <- as.factor(train_data$status)        

X_test <- test_data[, !(names(test_data) %in% c("time", "status"))]
y_test <- as.factor(test_data$status)          

# training SVM classifier
svm_model <- svm(x = X_train, y = y_train)

# Summary
summary(svm_model)

# perfoem on test data
predictions <- predict(svm_model, newdata = X_test)

# acuuracy
accuracy <- mean(predictions == y_test)
cat("Accuracy:", accuracy, "\n")


#LUSC Gene Expression
luscgene <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/gene expression/sgsfd.csv")
row.names(luscgene) <- luscgene[,1]
luscgene <- luscgene[,-1]


# Split the dataset into 70% and 30% in to training and testing sets
set.seed(42)
train_indices <- sample(1:nrow(luscgene), 0.7 * nrow(luscgene))  
train_data <- luscgene[train_indices, ]
test_data <- luscgene[-train_indices, ]

# Prepare training and testing data
X_train <- train_data[, !(names(train_data) %in% c("time", "status"))]
y_train <- as.factor(train_data$status)        

X_test <- test_data[, !(names(test_data) %in% c("time", "status"))]
y_test <- as.factor(test_data$status)          

# training SVM classifier
svm_model <- svm(x = X_train, y = y_train)

# Summary
summary(svm_model)

# perfoem on test data
predictions <- predict(svm_model, newdata = X_test)

# acuuracy
accuracy <- mean(predictions == y_test)
cat("Accuracy:", accuracy, "\n")
