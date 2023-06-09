# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-06-07
# Updated... 2023-06-08
#
# Description:
# Train a Ridge Regression and Lasso model from the TE-embeddings generated from
# from the Nucleotide transformer model
# Best demo for ridge and lasso regression
# https://www.pluralsight.com/guides/linear-lasso-and-ridge-regression-with-r
# ------------------------------------------------------------------------------

# Set directory to files
setwd('/workdir/mbb262/te/te_nucleotideTransformer_results/')

# Load packages
library(glmnet)
library(dplyr)
library(data.table)
# library(caret)

# Load in all chunks of results, rbind them
te_results <- list.files('/workdir/mbb262/te/te_nucleotideTransformer_results/',
                         pattern = "*chunk.txt", full.names = T)
te_results_list <- lapply(te_results, data.table::fread,  nThread= 60)
results_df <- data.table::rbindlist(te_results_list)
dim(results_df)

# Load in expression data, drop the sequence
exp_df <- read.delim(file = '../te_sequence_with_walley_expression.txt') %>% 
  select(-"teSeq")
dim(exp_df)

# Merge with expression data
exp_resuls_df <- cbind(exp_df, results_df)
dim(exp_resuls_df)


# Split into training and testing ----------------------------------------------

# Set seed for consistency
set.seed(100) 

# Do a 70/30 training/testing split and sample rows
index = sample(1:nrow(exp_resuls_df), 0.7*nrow(exp_resuls_df)) 
train = exp_resuls_df[index,] # Create the training data 
test = exp_resuls_df[-index,] # Create the test data

# Check if they're subsampled
dim(train)
dim(test)

# Turn into matrix
x_train <- train[,-c(1,2)] %>% as.matrix()
y_train <- train$exp %>% as.matrix()

x_test <- test[,-c(1,2)] %>% as.matrix()
y_test <- test$exp %>% as.matrix()


# ridge regression -------------------------------------------------------------

# Setting the range of lambda values
lambdas <- 10^seq(2, -3, by = -.1)

# Using glmnet function to build the ridge regression in r
ridge_reg = glmnet(x_train, y_train, nlambda = 25, alpha = 0, family = 'gaussian', lambda = lambdas)

# Checking the model
summary(ridge_reg)

# Using cross validation glmnet to find the optimal value for lambda
cv_ridge <- cv.glmnet(x_train, y_train, alpha = 0, lambda = lambdas)
optimal_lambda <- cv_ridge$lambda.min
optimal_lambda

# Compute R^2 from true and predicted values
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  
  # Model performance metrics
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
}

# Prediction on train data --------------
predictions_train_rr <- predict(ridge_reg, s = optimal_lambda, newx = x_train)
eval_results(y_train, predictions_train_rr, train)
plot_train_rr <- data.frame(y_train, predictions_train_rr)
write.table(plot_train_rr, "../predicted_observed_train_ridge_regression.csv",
            quote = FALSE, row.names = FALSE)

# Plot predicted vs observed
plot_value <- paste0("R^2 = ", round(eval_results(y_train, predictions_train_rr, train)[2], 3))
rr_train_plot <- ggplot(plot_train_rr, aes(x= s1, y = y_train))+
  geom_point() +
  geom_label(label = plot_value, x = 1.5, y = 13) +
  xlab("Predicted") +
  ylab("Observed") +
  ggtitle("Train: RR on B73 Mature Leaf ~ All TE Sequence within 10 Kb")
ggsave("/workdir/mbb262/te/rr_train_plot.png", rr_train_plot, 
       width = 4, height = 5, units = "in")


# Prediction on test data ----------------
predictions_test_rr <- predict(ridge_reg, s = optimal_lambda, newx = x_test)
eval_results(y_test, predictions_test_rr, test)
plot_test_rr <- data.frame(y_test, predictions_test_rr)
write.table(plot_test_rr, "../predicted_observed_test_ridge_regression.csv",
            quote = FALSE, row.names = FALSE)

# Plot test data 
plot_value <- paste0("R^2 = ", round(eval_results(y_test, predictions_test_rr, test)[2], 3))
rr_test_plot <- ggplot(plot_test_rr, aes(x= s1, y = y_test))+
  geom_point() +
  geom_label(label = plot_value, x = 1.5, y = 13) +
  xlab("Predicted") +
  ylab("Observed") +
  ggtitle("Test: RR on B73 Mature Leaf ~ All TE Sequence within 10 Kb")
ggsave("/workdir/mbb262/te/rr_test_plot.png", rr_test_plot)


# Run lasso --------------------------------------------------------------------

# Set range for lambdas
lambdas <- 10^seq(2, -3, by = -.1)

# Setting alpha = 1 implements lasso regression
lasso_reg <- cv.glmnet(x_train, y_train, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)

# Find best lambda 
lambda_best <- lasso_reg$lambda.min
plot(lambda_best)
lambda_best

# Train the lasso model
lasso_model <- glmnet(x_train, y_train, alpha = 1, lambda = lambda_best, standardize = TRUE)

# Predict on training data ------------
predictions_train_lasso <- predict(lasso_model, s = lambda_best, newx = x_train)
eval_results(y_train, predictions_train_lasso, train)
plot_train_lasso <- data.frame(y_train, predictions_train_lasso)
write.table(plot_train_lasso, "../predicted_observed_train_lasso.csv",
            quote = FALSE, row.names = FALSE)

# Plot predicted vs observed training data
plot_value <- paste0("R^2 = ", round(eval_results(y_train, predictions_train_lasso, train)[2], 3))
lasso_train_plot <- ggplot(plot_train_lasso, aes(x= s1, y = y_train))+
  geom_point() +
  geom_label(label = plot_value, x = 1.5, y = 13) +
  xlab("Predicted") +
  ylab("Observed") +
  ggtitle("Train: Lasso on B73 Mature Leaf ~ All TE Sequence within 10 Kb")
ggsave("/workdir/mbb262/te/lasso_train_plot.png", lasso_train_plot)


# Predict on testing data --------------
predictions_test_lasso <- predict(lasso_model, s = lambda_best, newx = x_test)
eval_results(y_test, predictions_test_lasso, test)
plot_test_lasso <- data.frame(y_test, predictions_test_lasso)
write.table(plot_test_lasso, "../predicted_observed_test_lasso.csv",
            quote = FALSE, row.names = FALSE)

# Plot predicted vs observed testing data
plot_value <- paste0("R^2 = ", round(eval_results(y_test, predictions_test_lasso, test)[2], 3))
lasso_test_plot <- ggplot(plot_test_lasso, aes(x= s1, y = y_test))+
  geom_point() +
  geom_label(label = plot_value, x = 1.5, y = 13) +
  xlab("Predicted") +
  ylab("Observed") +
  ggtitle("Test: Lasso on B73 Mature Leaf ~ All TE Sequence within 10 Kb")
ggsave("/workdir/mbb262/te/lasso_test_plot.png", lasso_test_plot)


