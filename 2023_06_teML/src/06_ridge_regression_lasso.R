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

# Set directory to b73 gff file
setwd('/Users/mbb262-admin/Library/CloudStorage/Box-Box/git_projects/hackathons/2023_06_teML/data/te_nucleotideTransformer_results/')

# Load packages
library(glmnet)
library(dplyr)
library(data.table)
library(caret)

# Load in all chunks of results, rbind them
te_results <- list.files('/Users/mbb262-admin/Library/CloudStorage/Box-Box/git_projects/hackathons/2023_06_teML/data/te_nucleotideTransformer_results',
                         pattern = "*.txt", full.names = T)
te_results_list <- lapply(te_results, data.table::fread)
results_df <- data.table::rbindlist(te_results_list)
dim(results_df)

# Load in expression data, drop the sequence
exp_df <- read.delim(file = 'te_sequence_with_walley_expression.txt') %>% 
  select(-"teSeq")
dim(exp_df)


# Merge with expression data
exp_resuls_df <- cbind(exp_df, results_df)

# Hack a thing -----------------------------------------------------------------
test1 <- data.table::fread('/Users/mbb262-admin/Library/CloudStorage/Box-Box/git_projects/hackathons/2023_06_teML/data/te_nucleotideTransformer_results/01_te_sequence_with_walley_expression.txt_result_output_chunk.txt')
test_exp_df <- exp_df[1:610,]
sub_exp_resuls_df <- cbind(test_exp_df, test1)

# Split into training and testing
set.seed(100) 
index = sample(1:nrow(test1), 0.7*nrow(test1)) 
train = test1[index,] # Create the training data 
test = test1[-index,] # Create the test data
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
ridge_reg = glmnet(x, y_train, nlambda = 25, alpha = 0, family = 'gaussian', lambda = lambdas)

# Checking the model
summary(ridge_reg)

# Using cross validation glmnet to find the optimal value for lambda
cv_ridge <- cv.glmnet(x_train, y_train, alpha = 0, lambda = lambdas)
optimal_lambda <- cv_ridge$lambda.min
optimal_lambda

# rebild model with best lambda
best_fit <- ridge_cv$glmnet.fit
head(best_fit)

# Rebuilding the model with optimal lambda value
best_ridge <- glmnet(x, y, alpha = 0, lambda = best_lambda)

# here x is the test dataset
pred <- predict(best_ridge, s = best_lambda, newx = x)

# R squared formula
actual <- test$Price
preds <- test$PreditedPrice
rss <- sum((preds - actual) ^ 2)
tss <- sum((actual - mean(actual)) ^ 2)
rsq <- 1 - rss/tss
rsq


# Run lasso --------------------------------------------------------------------
fit <- glmnet(x, y)
plot(fit)
print(fit)
cvfit <- cv.glmnet(x, y)
plot(cvfit)
predict(cvfit, newx = x[1:5,], s = "lambda.min")


