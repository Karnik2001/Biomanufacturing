# Loading required libraries
library(dplyr)
library(tidyverse)
library(readxl)
library(dbscan)
library(fastICA)
library(umap)
library(Rtsne)
library(MASS)
library(pls)
library(neuralnet)
library(glmnet)
library(fpc)

# Set the theme for plots to be black and white
theme_set(theme_bw(16))

# Load the data from an Excel file
drug <- read_excel("C:\\Users\\visha\\Downloads\\Bio-Drug_Production_Process_Simulation_DATA.xlsx")
View(drug)  # View the data in a tabular format

# Correlation Tests: Perform simple linear regression to check relationships between different variables
bio <- lm(drug$`Initial Biomass X_0` ~ drug$`Protein X_F`)  # Regression between Initial Biomass and Protein X_F
summary(bio)  # Display summary statistics for the regression model

bio1 <- lm(drug$`Initial Biomass X_0` ~ drug$`Impurity I_F`)  # Regression between Initial Biomass and Impurity I_F
summary(bio1)

bio2 <- lm(drug$`Initial Biomass X_0` ~ drug$Purity)  # Regression between Initial Biomass and Purity
summary(bio2)

bio3 <- lm(drug$`Protein X_F` ~ drug$Purity)  # Regression between Protein X_F and Purity
summary(bio3)

bio4 <- lm(drug$`Protein X_F` ~ drug$`Impurity I_F`)  # Regression between Protein X_F and Impurity I_F
summary(bio4)

# Correlation Tests: Checking the correlation between various variables
cor.test(drug$`Initial Biomass X_0`, drug$`Protein X_F`)  # Correlation between Initial Biomass and Protein X_F
cor.test(drug$`Initial Biomass X_0`, drug$`Impurity I_F`)  # Correlation between Initial Biomass and Impurity I_F
cor.test(drug$`Initial Biomass X_0`, drug$Purity)  # Correlation between Initial Biomass and Purity
cor.test(drug$`Protein X_F`, drug$Purity)  # Correlation between Protein X_F and Purity
cor.test(drug$`Protein X_F`, drug$`Impurity I_F`)  # Correlation between Protein X_F and Impurity I_F

# Scatter plots to visualize the relationships between variables
scatter.smooth(drug$`Initial Biomass X_0`, drug$`Protein X_F`,
               main = 'Biomass vs. Protein Production', pch = 19, col = 'blue')

scatter.smooth(drug$`Initial Biomass X_0`, drug$`Impurity I_F`,
               main = 'Biomass vs. Impurity', pch = 19, col = 'red')

scatter.smooth(drug$`Initial Biomass X_0`, drug$Purity, 
               main = 'Biomass vs. Purity', pch = 19, col = 'darkgreen')

scatter.smooth(drug$`Protein X_F`, drug$Purity, 
               main = 'Protein Production vs. Purity', pch = 19, col = 'purple')

scatter.smooth(drug$`Protein X_F`, drug$`Impurity I_F`, 
               main = 'Protein Production vs. Impurity', pch = 19, col = 'orange')

# Residual plots to check the fit of the linear models
br <- resid(bio)  # Residuals for the first linear model
br1 <- resid(bio1)  # Residuals for the second linear model
br2 <- resid(bio2)  # Residuals for the third linear model
br3 <- resid(bio3)  # Residuals for the fourth linear model
br4 <- resid(bio4)  # Residuals for the fifth linear model

# QQ-plots to check normality of residuals
qqnorm(br)  # QQ-plot for the first model
qqline(br)  # Add the reference line for normal distribution

qqnorm(br1)  # QQ-plot for the second model
qqline(br1)

qqnorm(br2)  # QQ-plot for the third model
qqline(br2)

qqnorm(br3)  # QQ-plot for the fourth model
qqline(br3)

qqnorm(br4)  # QQ-plot for the fifth model
qqline(br4)

# Linear Regression using Neural Networks

# MAX-MIN NORMALIZATION function to scale the data
drugn <- function(d) {
  return ((d - min(d)) / (max(d) - min(d) + 1e-6))  # Normalizing to a [0, 1] range
}

# Apply normalization to the relevant columns
drugdf <- as.data.frame(lapply(drug[2:5], drugn))  # Normalize columns 2 to 5

# Splitting the data into training (dmtr) and testing (dmte) sets
dmtr <- drugdf[1:15, ]  # Training set (first 15 rows)
dmte <- drugdf[16:20, ]  # Testing set (next 5 rows)

# Training the neural network model
dn <- neuralnet(dmtr$Protein.X_F ~., 
                data=dmtr, hidden=c(5,5), linear.output=TRUE, threshold=0.05)  # 5 neurons in 2 layers
dn$result.matrix  # Display the results from the neural network training

# Plot the neural network
par(plt = c(5, 5, 5, 5))  # Set plot region size
plot(dn)  # Neural network plot

# Extracting features for testing
dnt <- subset(dmte, select = c("Initial.Biomass.X_0", "Purity", "Impurity.I_F"))

# Predict using the trained neural network model
dn.results <- compute(dn, dnt)
dp <- dn.results$net.result  # Predicted values from the model

# Denormalizing the predictions to original scale
dp1 <- dp * (max(dmtr$Protein.X_F) - min(dmtr$Protein.X_F) + 1e-6) + min(dmtr$Protein.X_F)

# Create a dataframe to compare actual vs predicted values
dc <- data.frame(predicted = dp1, actual = dmte$Protein.X_F)

# Calculate deviation between actual and predicted values
dd <- (dc$actual - dc$predicted) / dc$actual

# Add deviation to the comparison dataframe
dc1 <- cbind(dc, dd)

# Calculate accuracy as 1 - average absolute deviation
dac <- 1 - abs(mean(deviation, na.rm = TRUE))

# Print the accuracy of the model
print(dac)

# Density-Based Clustering (DBSCAN)

bm <- as.matrix(drug[, -1])  # Convert data to matrix, excluding the first column (target variable)
kNNdistplot(bm, k=2)  # Plot k-nearest neighbor distances
abline(h=13, col="blue")  # Add a horizontal line at threshold distance

set.seed(200)  # Set seed for reproducibility
bd = dbscan(bm, 13, 2)  # Perform DBSCAN clustering (eps=13, minPts=2)
bd  # Display the clustering result

# Plot clusters
hullplot(bm, bd$cluster)
table(drug$`antigen A`, bd$cluster)  # Cross-tabulate clusters with antigen A values

# Evaluate clustering performance
bcs = cluster.stats(dist(drug[1:4]), bd$cluster)
bcs[c("within.cluster.ss", "avg.silwidth")]  # Display within-cluster sum of squares and average silhouette width

# Principal Component Regression (PCR)

set.seed(100)  # Set seed for reproducibility

# Perform Principal Component Regression with Cross-Validation
dr <- pcr(drug$`Initial Biomass X_0` ~ ., data = drug, scale = TRUE, validation = "CV")

# Display summary and validation plots for PCR
summary(dr)
validationplot(dr)
validationplot(dr, val.type="MSEP")  # Mean Squared Error of Prediction
validationplot(dr, val.type="R2")  # R-squared

# Prepare training and testing data for PCR
dtrx <- drug[1:15, c("Initial Biomass X_0", "Protein X_F", "Impurity I_F", "Purity")]  # Predictor variables
dtry <- drug[16:20, c("Protein X_F")]  # Response variable for testing
dtry1 <- drug[16:20, c("Initial Biomass X_0", "Impurity I_F", "Purity")]  # Testing predictors

# Train PCR model
dr1 <- pcr(dtrx$`Protein X_F` ~ ., data = dtrx, scale = TRUE, validation = "CV")

# Make predictions with the PCR model
dpcr <- predict(dr1, dtry1, ncomp = 2)  # Use 2 principal components for prediction

# Calculate RMSE (Root Mean Squared Error)
dtry <- as.numeric(dtry$`Protein X_F`)
sqrt(mean((dpcr - dtry)^2))  # RMSE

# Lasso Regression

dri <- drug$`Protein X_F`  # Response variable for Lasso

# Prepare predictor variables
drii <- data.matrix(drug[, c('Initial Biomass X_0', 'Impurity I_F', 'Purity')])

# Perform Cross-Validation for Lasso Regression
dric <- cv.glmnet(drii, dri, alpha = 1)

# Get the best lambda (regularization parameter)
dricb <- dric$lambda.min
plot(dric)

# Train the Lasso Regression model with the best lambda
dricbm <- glmnet(drii, dri, alpha = 1, lambda = dricb)
coef(dricbm)  # Display the coefficients

# Make predictions with the trained Lasso model
dricby <- predict(dricbm, s = dricb, newx = drii)

# Calculate R-squared for model performance
drugs <- sum((dri - mean(dri))^2)
druge <- sum((dricby - dri)^2)
drugr <- 1 - druge/drugs  # R-squared value

print(drugr)

# Ridge Regression

driy <- drug$`Protein X_F`  # Response variable for Ridge

# Prepare predictor variables for Ridge
driix <- data.matrix(drug[, c('Initial Biomass X_0', 'Impurity I_F', 'Purity')])

# Perform Ridge Regression
dricc <- glmnet(driix, driy, alpha = 0)  # alpha = 0 for Ridge

# Perform Cross-Validation for Ridge Regression
driccc <- cv.glmnet(driix, driy, alpha = 0)

# Get the best lambda for Ridge
dricbb <- dric$lambda.min
plot(driccc)

# Train the Ridge Regression model with the best lambda
dricbmm <- glmnet(driix, driy, alpha = 0, lambda = dricbb)
coef(dricbmm)  # Display the coefficients

# Make predictions with the trained Ridge model
dricbyy <- predict(dricc, s = dricbb, newx = driix)

# Plot the Ridge model results
plot(driccc, xvar = "lambda")

# Calculate R-squared for Ridge model performance
drugs1 <- sum((driy - mean(driy))^2)
druge1 <- sum((dricbyy - driy)^2)
drugr1 <- 1 - druge1/drugs1  # R-squared value for Ridge

print(drugr1)
