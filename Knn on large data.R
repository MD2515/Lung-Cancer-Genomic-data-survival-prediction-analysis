

library(class)
#######################################################LUAD#################################################################

#LUAD methylation

luadmeth <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/methylation/amsfd.csv")
row.names(luadmeth) <- luadmeth[,1]
luadmeth <- luadmeth[,-1]

# Split the dataset into training and testing sets
set.seed(42) 
train_indices <- sample(1:nrow(luadmeth), 0.7 * nrow(luadmeth))  # 70% for training
train_data <- luadmeth[train_indices, ]
test_data <- luadmeth[-train_indices, ]

# train dataset without time and status
X_train <- train_data[, !(names(train_data) %in% c("time", "status"))]
y_train <- train_data$status  

# test data
X_test <- test_data[, !(names(test_data) %in% c("time", "status"))]
y_test <- test_data$status    # Actual labels (status) for testing

# KNN
k <- 10  
knn_model <- knn(train = X_train, test = X_test, cl = y_train, k = k)

# Evaluate the accuracy of the predictions
accuracy <- mean(knn_model == y_test)
cat("Accuracy:", accuracy, "\n")


#LUAD Gene Expression

luadgene <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUAD/gene expression/agsfd.csv")
row.names(luadgene) <- luadgene[,1]
luadgene <- luadgene[,-1]

# Split the dataset into training and testing sets
set.seed(42)  
train_indices <- sample(1:nrow(luadgene), 0.7 * nrow(luadgene))  # 70% for training
train_data <- luadgene[train_indices, ]
test_data <- luadgene[-train_indices, ]

# train dataset without time and status
X_train <- train_data[, !(names(train_data) %in% c("time", "status"))]
y_train <- train_data$status  

# test data
X_test <- test_data[, !(names(test_data) %in% c("time", "status"))]
y_test <- test_data$status    # Actual labels (status) for testing

# KNN
k <- 10  
knn_model <- knn(train = X_train, test = X_test, cl = y_train, k = k)

# Evaluate the accuracy of the predictions
accuracy <- mean(knn_model == y_test)
cat("Accuracy:", accuracy, "\n")



####################################################LUSC#####################################################################


#LUSC Methylation

luscmeth <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/methylation/smsfd.csv")
row.names(luscmeth) <- luscmeth[,1]
luscmeth <- luscmeth[,-1]

# Split the dataset into training and testing sets
set.seed(42)  
train_indices <- sample(1:nrow(luscmeth), 0.7 * nrow(luscmeth))  # 70% for training
train_data <- luscmeth[train_indices, ]
test_data <- luscmeth[-train_indices, ]

# train dataset without time and status
X_train <- train_data[, !(names(train_data) %in% c("time", "status"))]
y_train <- train_data$status  

# test data
X_test <- test_data[, !(names(test_data) %in% c("time", "status"))]
y_test <- test_data$status    # Actual labels (status) for testing

# KNN
k <- 10  
knn_model <- knn(train = X_train, test = X_test, cl = y_train, k = k)

# Evaluate the accuracy of the predictions
accuracy <- mean(knn_model == y_test)
cat("Accuracy:", accuracy, "\n")


#LUSC Gene Expression

luscgene <- read.csv("D:/BHI/Spring 2024/BHI 699/Thesis Documents/Data/lung/LUSC/gene expression/sgsfd.csv")
row.names(luscgene) <- luscgene[,1]
luscgene <- luscgene[,-1]

# Split the dataset into training and testing sets
set.seed(42)  
train_indices <- sample(1:nrow(luscgene), 0.7 * nrow(luscgene))  # 70% for training
train_data <- luscgene[train_indices, ]
test_data <- luscgene[-train_indices, ]

# train dataset without time and status
X_train <- train_data[, !(names(train_data) %in% c("time", "status"))]
y_train <- train_data$status

# test data
X_test <- test_data[, !(names(test_data) %in% c("time", "status"))]
y_test <- test_data$status    # Actual labels (status) for testing

# KNN
k <- 10  
knn_model <- knn(train = X_train, test = X_test, cl = y_train, k = k)

# Evaluate the accuracy of the predictions
accuracy <- mean(knn_model == y_test)
cat("Accuracy:", accuracy, "\n")
