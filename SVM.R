library(e1071)
methylation_data <- read.csv("D:/BHI/Spring 2024/BHI 699/lung/LUAD/methylation/sfd.csv")
row.names(methylation_data) <- methylation_data[,1]
methylation_data <- methylation_data[,-1]
# Split the dataset into training and testing sets
set.seed(42)  # Set seed for reproducibility
train_indices <- sample(1:nrow(methylation_data), 0.7 * nrow(methylation_data))  # 70% for training
train_data <- methylation_data[train_indices, ]
test_data <- methylation_data[-train_indices, ]

# Prepare training and testing data
X_train <- train_data[, !(names(train_data) %in% c("time", "status"))]
y_train <- as.factor(train_data$status)        # Labels (status) for training

X_test <- test_data[, !(names(test_data) %in% c("time", "status"))]
y_test <- as.factor(test_data$status)          # Actual labels (status) for testing

# Train an SVM classifier
svm_model <- svm(x = X_train, y = y_train)

# Summary of the trained SVM model
summary(svm_model)

# Make predictions on the test set
predictions <- predict(svm_model, newdata = X_test)

# Evaluate the accuracy of the predictions
accuracy <- mean(predictions == y_test)
cat("Accuracy:", accuracy, "\n")
