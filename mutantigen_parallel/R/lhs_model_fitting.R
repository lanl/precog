library(caret)
library(gbm)
set.seed(022924)

#### Fitting Model ####

# Compiled dataset of all simulations that we've run with their parameters and current set of outcomes
super_outcomes <- read.csv("../mutantigen_experiment/params_outcomes.csv")


# 1. Criteria1 = Did it finish?
# 2. Criteria2 = Did it finish & (turnover by the end OR median # of circulating variants)


# split up data into 80% training and 20% testing
super_outcomes$criteria2 <- factor(super_outcomes$criteria2)
levels(super_outcomes$criteria2) <- make.names(levels(super_outcomes$criteria2))  # Ensures valid names
print(levels(super_outcomes$criteria2))  # Check new levels
super_outcomes$criteria2 <- relevel(super_outcomes$criteria2, ref = "X0")


training_ids <- sample(x = seq(from = 1, to = nrow(super_outcomes)), size = .80 * nrow(super_outcomes), replace = FALSE)
colnames(super_outcomes)

# Remove other outcome columns
train_df <- super_outcomes[training_ids, c(7:16, 18)] 
test_df <- super_outcomes[-training_ids, c(7:16, 18)]

# Test different metrics for training
sensitivity_metric <- function(data, lev = NULL, model = NULL) {
  cm <- confusionMatrix(data$pred, data$obs, positive = "X1")
  #specificity <- cm$byClass["Specificity"]
  #c(Specificity = specificity)
  sensitivity <- cm$byClass["Sensitivity"]
  c(Sensitivity = sensitivity)
}


F1_metric <- function(data, lev = NULL, model = NULL) {
  if (!all(c("pred", "obs") %in% colnames(data))) {
    stop("Missing 'pred' or 'obs' column in the data")
  }
  
  precision <- posPredValue(data$pred, data$obs, positive = "X1", na.rm = TRUE)
  recall <- sensitivity(data$pred, data$obs, positive = "X1", na.rm = TRUE)
  
  # Ensure F1 is not NA
  if (is.na(precision) | is.na(recall) | (precision + recall) == 0) {
    F1 <- 0
  } else {
    F1 <- 2 * (precision * recall) / (precision + recall)
  }
  
  return(c(F1 = F1))
}

F1_weighted_metric <- function(data, lev = NULL, model = NULL, beta = 2) {
  if (!all(c("pred", "obs") %in% colnames(data))) {
    stop("Missing 'pred' or 'obs' column in the data")
  }
  
  precision <- posPredValue(data$pred, data$obs, positive = "X1", na.rm = TRUE)
  recall <- sensitivity(data$pred, data$obs, positive = "X1", na.rm = TRUE)
  
  # Ensure F1 is not NA
  if (is.na(precision) | is.na(recall) | (beta^2 * precision + recall) == 0) {
    F_beta <- 0
  } else {
    F_beta <- (1 + beta^2) * (precision * recall) / ((beta^2 * precision) + recall)
  }
  
  return(c(F1_beta = F_beta))
}


# Test different metrics for minimizing
# Went through and tried sampling = "down" & sampling = "smote" 

#### Test Sensitivity #### 
train_control <- trainControl(
  method = "cv",      # Cross-validation
  number = 5,         # Number of folds
  classProbs = TRUE,
  verboseIter = TRUE,  # Show progress
  summaryFunction = sensitivity_metric,
  sampling = "down" # Still have about a 1:2 ratio
)

trainGrid <- expand.grid(n.trees = c(15, 20, 25),
                         interaction.depth = c(2, 3, 5), 
                         shrinkage = c(0.1),
                         n.minobsinnode = 20)

gbm_model_sensitivity <- train(
  criteria2 ~ .,         # Formula (target ~ predictors)
  data = train_df,          # Data frame
  method = "gbm",      # Logistic Regression
  trControl = train_control,
  tuneGrid = trainGrid,
  #  weights = ifelse(train_df$category == "1", 2, 1),
  tuneLength = 5,
)

# Evaluate testing data model performance
gbm_test_predictions <- predict(gbm_model_sensitivity, test_df)
gbm_test_cm_sensitivity <- confusionMatrix(gbm_test_predictions, test_df$criteria2, positive = "X1")

##### Testing F1 ####
train_control <- trainControl(
  method = "cv",      # Cross-validation
  number = 5,         # Number of folds
  classProbs = TRUE,
  verboseIter = TRUE,  # Show progress
  summaryFunction = F1_metric,
  sampling = "down"
)

trainGrid <- expand.grid(n.trees = c(20, 35, 50, 100, 500),
                         interaction.depth = c(3, 5, 10, 20), 
                         shrinkage = 0.1,
                         n.minobsinnode = 20)
gbm_model_f1 <- train(
  criteria2 ~ .,         # Formula (target ~ predictors)
  data = train_df,          # Data frame
  method = "gbm",      # Logistic Regression
  trControl = train_control,
  tuneGrid = trainGrid,
  #  weights = ifelse(train_df$category == "1", 2, 1),
  tuneLength = 5,
  metric = "F1",
)

gbm_test_predictions <- predict(gbm_model_f1, test_df)
gbm_test_cm_f1 <- confusionMatrix(gbm_test_predictions, test_df$criteria2, positive = "X1")


#### Testing F2 ####
train_control <- trainControl(
  method = "cv",      # Cross-validation
  number = 5,         # Number of folds
  classProbs = TRUE,
  verboseIter = TRUE,  # Show progress
  summaryFunction = function(data, lev, model) F1_weighted_metric(data, lev, model, beta = 2),  # More weight on recall
  sampling = "down"
)

trainGrid <- expand.grid(n.trees = c(20, 35, 50, 100, 500),
                         interaction.depth = c(3, 5, 10, 20), 
                         shrinkage = 0.1,
                         n.minobsinnode = 20)

gbm_model_f1_weighted <- train(
  criteria2 ~ .,         # Formula (target ~ predictors)
  data = train_df,          # Data frame
  method = "gbm",      # Logistic Regression
  trControl = train_control,
  tuneGrid = trainGrid,
  #  weights = ifelse(train_df$category == "1", 2, 1),
  tuneLength = 5,
  metric = "F1_beta"
)

# Evaluate testing data model performance
gbm_test_predictions <- predict(gbm_model_f1_weighted, test_df)
gbm_test_cm_weighted <- confusionMatrix(gbm_test_predictions, test_df$criteria, positive = "X1")

## Looking at my three different metrics
down_results <- data.frame(sensitivity = gbm_test_cm_sensitivity$byClass,
                      f1 = gbm_test_cm_f1$byClass,
                      fbeta_2 = gbm_test_cm_weighted$byClass)


# Based on these results, it seems like down_f1 is the best
# F1 scores are really low - lots of false positives (low precision)

down_results
# smote_results
saveRDS(gbm_model_f1, "../mutantigen_experiment/gbm_fbeta_full.RDS")


# Visualization of predictions
gbm_test_predictions <- predict(gbm_model_f1, test_df)
test_df$predictions <- gbm_test_predictions
test_df$result <- test_df$predictions == test_df$criteria2

library(GGally)
ggpairs(test_df,                # Data frame
        columns = c(1:10),        # Columns
        aes(color = as.factor(result),  # Color by group (cat. variable)
            alpha = 0.5))

# Looks like we are getting lambda predictions wrong -- in which direction?
actual_split <- test_df %>% ggplot(aes(x = lambda, fill = criteria2)) + geom_density(alpha = .5)

predicted_split <- test_df %>% ggplot(aes(x = lambda, fill = predictions)) + geom_density(alpha = .5)

accuracy_split <- test_df %>% ggplot(aes(x = lambda, fill = result)) + geom_density(alpha = .5)


# when lambda is low, there is a tendency to think it will finish and have interesting dynamics, but it may not always
plot_grid(actual_split, predicted_split, accuracy_split, ncol = 1)


#### Example of applying the model to a new lhs ####
library(lhs)   # For running the parameter lhs
library(yaml)  # For working with yaml files 
library(readr) # For reading data


#### Step 1: Determine the parameter settings 

# Define the number of samples
n_samples <- 1000

# Define the parameter ranges
param_ranges <- list(
  initialNs = c(10^4, 10^8),
  demeAmplitudes = c(0, 0.25),
  lambdaAntigenic = c(8.57143e-05, 0.002571429),
  meanAntigenicSize = c(1.2e-3, 1.2e-1),
  lambda = c(0.0095, 4.08),
  mutCost = c(8e-4, 8e-2),
  beta = c(0.1428, 2.25),
  nu = c(1/4, 1/14),
  epsilon_mut = c(0.5, 1.5),
  initialI_prop = c(0.0001, .001)
)

#### Step 2:  Generate the Latin Hypercube Sample 
lhs_samples <- randomLHS(n_samples, length(param_ranges))

# Scale the LHS samples to the defined parameter ranges

# For (initialNs) - uniform across the log scale 
range       <- param_ranges[[1]]
log_samples <- (log10(range[2]) - log10(range[1])) * lhs_samples[,1] + log10(range[1])
log_samples <- 10^log_samples

# For rest of the samples
scaled_samples <- t(apply(lhs_samples, 1, function(x) {
  sapply(2:length(param_ranges), function(i) {
    range <- param_ranges[[i]]
    range[1] + x[i] * (range[2] - range[1])
  })
}))

scaled_samples <- cbind(log_samples, scaled_samples)


# Convert the scaled samples into a data frame for readability
colnames(scaled_samples) <- names(param_ranges)
lhs_df <- as.data.frame(scaled_samples)


#### Step 3: Parameters that have to be adjusted based on others 
lhs_df$thresholdAntigenicSize <- lhs_df$meanAntigenicSize
lhs_df$initialI <- lhs_df$initialNs * lhs_df$initialI_prop
lhs_df$externalMigration <- lhs_df$initialNs* 0.00001
lhs_df$epsilon <- (.16/40000000) * lhs_df$initialNs * lhs_df$epsilon_mut

#### Step 4: Apply predictions
lhs_end_predictions <- predict(gbm_model_f1, newdata = lhs_df)
lhs_df_sub <- lhs_df[which(lhs_end_predictions == "X1"), ]

# Remove extra column 
lhs_df_sub <- lhs_df[,-which(colnames(lhs_df_sub) == "initialI_prop")]





