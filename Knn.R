# Load the necessary package for KNN
library(class)  # for KNN classifier

methylation_data <- read.csv("D:/BHI/Spring 2024/BHI 699/lung/LUAD/methylation/sfd.csv")
row.names(methylation_data) <- methylation_data[,1]
methylation_data <- methylation_data[,-1]

# Split the dataset into training and testing sets
set.seed(42)  # Set seed for reproducibility
train_indices <- sample(1:nrow(methylation_data), 0.7 * nrow(methylation_data))  # 70% for training
train_data <- methylation_data[train_indices, ]
test_data <- methylation_data[-train_indices, ]

# Extract features (methylation values) for training (excluding "time" and "status" columns)
X_train <- train_data[, !(names(train_data) %in% c("time", "status"))]
y_train <- train_data$status  # Labels (status) for training

# Extract features for testing (excluding "time" and "status" columns)
X_test <- test_data[, !(names(test_data) %in% c("time", "status"))]
y_test <- test_data$status    # Actual labels (status) for testing

# Train a KNN classifier
k <- 5  # Specify the number of neighbors (k)
knn_model <- knn(train = X_train, test = X_test, cl = y_train, k = k)

# Evaluate the accuracy of the predictions
accuracy <- mean(knn_model == y_test)
cat("Accuracy:", accuracy, "\n")
