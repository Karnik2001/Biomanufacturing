# Load necessary libraries
library(dplyr)
library(tidyverse)
library(glmnet)

# Load the dataset
pro <- read.csv("C:\\Users\\visha\\Downloads\\100_Batches_IndPenSim_Statistics.csv")
View(pro)

# Linear Regression Analysis
# Scatter plot between harvested during batch and harvested end of batch
pro1 <- scatter.smooth(pro$Penicllin_harvested_during_batch.kg.,
                       pro$Penicllin_harvested_end_of_batch..kg., 
                       main = " Start vs End Batch")

# Linear regression model between harvested during batch and harvested end of batch
pro1m <- lm(pro$Penicllin_harvested_during_batch.kg.~ pro$Penicllin_harvested_end_of_batch..kg.)
pro2m <- lm(pro$Penicllin_yield_total..kg.~ pro$Penicllin_harvested_during_batch.kg.)
pro3m <- lm(pro$Penicllin_yield_total..kg.~ pro$Penicllin_harvested_end_of_batch..kg.)

# Perform correlation tests between variables
cor.test(pro$Penicllin_harvested_during_batch.kg.,pro$Penicllin_harvested_end_of_batch..kg.)
cor.test(pro$Penicllin_yield_total..kg.,pro$Penicllin_harvested_during_batch.kg.)
cor.test(pro$Penicllin_yield_total..kg.,pro$Penicllin_harvested_end_of_batch..kg.)

# Display summary statistics for each linear model
summary(pro1m)
summary(pro2m)
summary(pro3m)

# Scatter plot of harvested during batch vs harvested end of batch
par(mar = c(4,4,4,4))  # Set margins for the plot
scatter.smooth(pro$Penicllin_harvested_during_batch.kg,pro$Penicllin_harvested_end_of_batch..kg.,
               main = " Harvested During Batch vs Harvested End Of Batch",  pch = 19, col = 'deepskyblue')

# Scatter plot of total yield vs harvested end of batch
par(mar = c(4,4,4,4))  
scatter.smooth(pro$Penicllin_yield_total..kg.,pro$Penicllin_harvested_end_of_batch..kg.,
               main = " Total Yield vs Harvested End Of Batch", pch = 19, col = 'cornflowerblue')

# Scatter plot of total yield vs harvested during the batch
par(mar = c(4,4,4,4))  
scatter.smooth(pro$Penicllin_yield_total..kg.,pro$Penicllin_harvested_during_batch.kg.,
               main = " Total Yield vs Harvested During The Batch", pch = 19, col = 'cyan2')

# Residual analysis for each linear model
pres <- resid(pro1m)
plot(fitted(pro1m), pres)
abline(0,0)
qqnorm(pres)
qqline(pres)

pres2 <- resid(pro2m)
plot(fitted(pro2m), pres2)
abline(0,0)
qqnorm(pres2)
qqline(pres2)

pres3 <- resid(pro3m)
plot(fitted(pro3m), pres3)
abline(0,0)
qqnorm(pres3)
qqline(pres3)

# Logistic Regression

# Set seed for reproducibility
set.seed(100)

# Split data into training (80%) and testing (20%) sets
bat <- sample(c(TRUE,FALSE),nrow(pro), replace = TRUE, prob = c(0.8,0.2))
btr <- pro[bat,]  # Training data
bte <- pro[!bat, ]  # Testing data

# Logistic regression model for predicting fault based on harvested quantities and yield
bm <- glm(Fault.ref.0.NoFault.1.Fault. ~ Penicllin_harvested_during_batch.kg. + Penicllin_harvested_end_of_batch..kg.
          + Penicllin_yield_total..kg.,family = "binomial", data = btr)

# Model summary and evaluation metrics
summary(bm)
pscl::pR2(bm)["McFadden"]  # McFadden's pseudo R-squared
caret::varImp(bm)  # Variable importance
car::vif(bm)  # Variance inflation factor to check multicollinearity

# Ensure testing data is in a data frame format
bte <- as.data.frame(bte)

# Predict on the testing set
bmp <- predict(bm, newdata = bte, type = "response")

# Convert probabilities to binary predictions (0 or 1)
bnp <- ifelse(bmp > 0.5, 1, 0)

# Calculate accuracy of the model
bacc <- mean(bnp == bte$Fault.ref.0.NoFault.1.Fault.)

# Display the accuracy score
print(paste("Accuracy Score:", bacc))

# View the testing set
View(bte)

# Principal Component Analysis (PCA)

# Perform PCA on the data (scale data before PCA)
prop <- prcomp(pro, scale. = TRUE)
prop$rotation <- -1*prop$rotation  # Invert the rotation matrix
prop$rotation

prop$x <- -1*prop$x  # Invert the principal component scores
head(prop$x)

# Plot the biplot of the principal components
biplot(prop, scale = 0)

# Calculate and display the proportion of variance explained by each component
prop$sdev^2/sum(prop$sdev^2)

provar = prop$sdev^2/sum(prop$sdev^2)  # Variance explained
provar

# Create a scree plot to visualize the explained variance
qplot(c(1:5), provar) + geom_line() +
  xlab("Principal Component") +
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0,1)

# Lasso Regression

# Define the dependent variable and independent variables
prol <- pro$Penicllin_yield_total..kg.
prol1 <- data.matrix(pro[,c("Penicllin_harvested_during_batch.kg.","Penicllin_harvested_end_of_batch..kg.")])

# Perform cross-validation for Lasso regression
prolc <- cv.glmnet(prol1, prol, alpha = 1)

# Extract the optimal lambda value
prol2 <- prolc$lambda.min

# Plot the cross-validation results
plot(prolc)

# Fit the Lasso regression model with the optimal lambda
prolb <- glmnet(prol1,prol, alpha = 1, lambda = prol2)
coef(prolb)

# Predict using the fitted Lasso model
prolm = matrix(c(9.00e+04, 9.80e-01), nrow = 1, ncol = 2)
predict(prolb, s = prol2, newx = prolm)

# Calculate predicted values for the entire dataset
prolp <- predict(prolb, s = prol2, newx = prol1)

# Compute R-squared for the model
pros <- sum((prol - mean(prol))^2)
pross <- sum((prolp - prol)^2)
propr <- 1 - pross/pros
propr

# Ridge Regression

# Define the dependent variable and independent variables for Ridge regression
proly <- pro$Penicllin_yield_total..kg.
prolx <- data.matrix(pro[,c("Penicllin_harvested_during_batch.kg.","Penicllin_harvested_end_of_batch..kg.")])

# Fit the Ridge regression model
prolcx <- glmnet(prolx, proly, alpha = 0)

# Perform cross-validation for Ridge regression
prolcx1 <- cv.glmnet(prolx, proly, alpha = 0)

# Extract the optimal lambda value
prolb1 <- prolcx1$lambda.min

# Plot the cross-validation results
plot(prolcx1)

# Fit the Ridge regression model with the optimal lambda
propb <- glmnet(prolx, proly, alpha = 0, lambda = prolb1)
coef(propb)

# Plot the Ridge regression path
plot(prolcx, xvar = "lambda")

# Predict using the fitted Ridge model
prolpp <- predict(prolcx, s = prolb1, newx = prolx)

# Compute R-squared for the Ridge model
props <- sum((proly - mean(proly))^2)
propss <- sum((prolpp - proly)^2)
propr1 <- 1 - propss/props
propr1
