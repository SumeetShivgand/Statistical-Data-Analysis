
library(tidyverse)
library(car)
library(caret)
require("tidyverse")
require("MASS")



################ Question 1  ##############

# Read the Data
amino_acid <- read.csv("D:\\CIT\\Semester 2\\Statistical Data Analysis\\Assignment 1\\amino_acidF.csv",header=TRUE)
View(amino_acid)

################# Question A ##############

plot(amino ~ folate, data = amino_acid)
cor(amino_acid$amino, amino_acid$folate)

summary(amino_acid)

# Histogram
hist(amino_acid$folate)

# Scatterplot
scatter.smooth(amino_acid$amino, amino_acid$folate, main = " Folate ~ Amino")

# Check for Outliers
#Boxplot
par(mfrow=c(1, 2))
boxplot(amino_acid$amino, main = "Amino", sub=paste("Outlier rows: ", boxplot.stats(amino_acid$amino)$out))
boxplot(amino_acid$folate, main = "Folate", sub=paste("Outlier rows: ", boxplot.stats(amino_acid$folate)$out))


############## Question B ###############

#fit a linear model to the 'folate' variable  and examine the residuals
model1<- lm(folate ~ amino, data= amino_acid)
summary(model1)
plot(folate ~ amino, data = amino_acid)
abline(model1)
# The fitted model was: 
#y = 3.5581 - 0.3718x
# Interpretation of β1 coefficient
#β1 = - 0.3718


########## Question C ########################
#?confint
confint(model1,level = 0.95,"amino")

################### Question D ##############
# Test the hypothesis:
# H0: β1 = 0 (There is no linear relationship between amino and folate)
# HA: β1 ≠ 0 (There is linear relationship between amino and folate)
anova(model1)

################ Question E #################
## create a sequence of values for the explanatory variable (x)
Amino_grid = seq(min(amino_acid$amino), max(amino_acid$amino), by = 0.01)


## create a set of 95% prediction intervals at each point in Amino_grid
dist_pi_band = predict(model1, newdata = data.frame(amino = Amino_grid), interval = "prediction", level = 0.95) 

## create a set of 95% prediction intervals at each point in Amino_grid
dist_ci_band = predict(model1, newdata = data.frame(amino = Amino_grid), interval = "confidence", level = 0.95) 

## plot scatter graph of Amino against Folate
plot(amino_acid$folate ~ amino_acid$amino, xlab = "Amino ((μ mol/L) )",  ylab = "Folate ((μg/L) )",  pch  = 20,
     cex  = 0.75, ylim = c(min(dist_pi_band), max(dist_pi_band)))

# Add regression Line
abline(model1)

## plot 95% prediction and confidence bands
lines(Amino_grid, dist_ci_band[,"lwr"], col = "blue", lwd = 1, lty = 2)
lines(Amino_grid, dist_ci_band[,"upr"], col = "blue", lwd = 1, lty = 2)
lines(Amino_grid, dist_pi_band[,"lwr"], col = "red", lwd = 1, lty = 3)
lines(Amino_grid, dist_pi_band[,"upr"], col = "red", lwd = 1, lty = 3)

############## Question F ####################

rstud <- rstudent(model1)

plot(rstud~fitted(model1), xlab = "Fitted values", ylab =" studentised residuals")
abline(0,0)


############ Question G #######################
# plot leverage and find any abservation with high leverage.
# see the plot no 3. 
plot(model1)

############ Question H ################

# Refer report for answer.



######################################
########### Question 2 ###############

# Read the Data
divusaF <- read.csv("D:\\CIT\\Semester 2\\Statistical Data Analysis\\Assignment 1\\divusaF.csv",header=TRUE)
View(divusaF)
# Year is a discrete variable and we will not be using it as a predictor in this assignment.
divusaF[,1]<-NULL

########### Question A #################

########### Examine the data ###########
## Pairs plot
##function to put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}


## function to put correlations on the upper panels,
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- (cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 2
  text(0.5, 0.5, txt, cex = cex.cor)
}
## pairs plot 
windows(10,10)
pairs(divusaF, diag.panel = panel.hist, upper.panel = panel.cor)
## pairs plot ##
# Note that the unemployed variable and military variables are skewed
# There is a positive correlation between divorce and femlab (r =0.83 )
# There is a positive correlation between unemployed and birth (r =0.67)
# There is a positive correlation between marriage and birth (r =0.56)


windows(10,10)
par(mfrow=c(2,3))
boxplot(divusaF$divorce, xlab="divorce")
boxplot(divusaF$unemployed, xlab="unemployed")
boxplot(divusaF$femlab, xlab="femlab")
boxplot(divusaF$marriage, xlab="marriage")
boxplot(divusaF$birth, xlab="birth")
boxplot(divusaF$military, xlab="military")


## To examine simple linear models for each predictor ##
m_unemployed <- lm(divorce ~ unemployed, data = divusaF)
m_femlab <- lm(divorce ~ femlab, data = divusaF)
m_marriage <- lm(divorce ~ marriage, data = divusaF)
m_birth <- lm(divorce ~ birth, data = divusaF)
m_military <- lm(divorce ~ military, data = divusaF)

## plot each of the regression lines
# Note the relationship between divorce and each of the predictors
windows(10,7)
par(mfrow = c(2,3)) 

plot(divorce ~ unemployed, data = divusaF)
abline(m_unemployed)
plot(divorce ~ femlab, data = divusaF)
abline(m_femlab)
plot(divorce ~ marriage, data = divusaF)
abline(m_marriage)
plot(divorce ~ birth, data = divusaF)
abline(m_birth)
plot(divorce ~ military, data = divusaF)
abline(m_military)


#################### Question B ########################
# Fit the model
# y = β_0 + β_1 unemployed + β_2 femlab + β_3 marriage + β_4 birth + β_3 military + e

############## (i)#############
m1<-lm(divorce ~ unemployed + femlab + marriage + birth + military, data = divusaF)
summary(m1)

# Interpret the coefficient for 'femlab'
# y = 33.793826 + 0.009136 unemployed + 0.284877 femlab - 0.281241 marriage - 0.101615 birth - 0.027808 military + e
# coefficient for 'femlab'is 0.284877 

########### (ii) ############
####Calculate the variance inflation factors for this model and discuss their implications for collinearity in the model.
# The VIFs lie between 1 and 2.1 indicating that collinearity is not 
# having a large impact on the coefficient estimates for this model.
vif(m1)

########## (iii) ########
### Fitting alternative model to check collinearity
m2 <-lm(divorce ~ femlab + marriage + birth, data = divusaF)
summary(m2)
vif(m2)
# For alternative model,the VIFs lie between 1 and 2 indicating that collinearity is not 
# having a large impact on the coefficient estimates for this model.

######### (iv) ########
# Create a partial regression plot showing the relationship between 
# birth and divorce adjusted for unemployed, femlab, marriage and military. 

# first create the model birth ~ unemployed + femlab + marriage + military
m_pr_birth <-lm(birth ~ unemployed + femlab + marriage + military, data = divusaF)

# second create the model divorce ~ unemployed + femlab + marriage + military
m_pr_divorce <-lm(divorce ~ unemployed + femlab + marriage + military, data = divusaF)

# plot residuals
windows(5,5)
plot(m_pr_divorce$res ~ m_pr_birth$res, xlab = "residuals birth ~ unemployed + femlab + marriage + military", 
     ylab = "residuals divorce ~ unemployed + femlab + marriage + military")

# fit a regression model to the residuals
m_res <-lm(m_pr_divorce$res ~ m_pr_birth$res)
abline(m_res)
summary(m_res)

################# (v) #############
# H0: β_1=β_2=β_3=β_4=0
# HA: at least one of the β_i≠0

summary(m1)


################### (vi) ############

### 7.	Assess the fit of the model using diagnostic plots and try to improve the 
### fit of the model.
windows(10,10)
par(mfrow=c(2,2))
plot(m1)

# try a log transformation 
divusaF$log_military<-log(divusaF$military)
windows(10,10)
pairs(divusaF, diag.panel = panel.hist, upper.panel = panel.cor)

# refit model
m3<-lm(divorce ~ unemployed + femlab + marriage + birth + log_military, data = divusaF)
summary(m3)
windows(10,10)
par(mfrow=c(2,2))
plot(m2)

# The pattern in the plot of residuals vs fitted appears to be gone however
# we note some observations with high leverage
# Try transforming the unemployed variable

divusaF$log_unemployed<-log(divusaF$unemployed)
windows(10,10)
pairs(divusaF, diag.panel = panel.hist, upper.panel = panel.cor)

# refit model
m4<-lm(divorce ~ log_unemployed + femlab + marriage + birth + log_military, data = divusaF)
summary(m4)
windows(10,10)
par(mfrow=c(2,2))
plot(m4)


########## Question C #############

full_model <-lm(divorce ~ unemployed + femlab + marriage + birth + military, data = divusaF)
summary(full_model)

reduced_model <-lm(divorce ~ femlab + marriage + birth + military, data = divusaF)
summary(reduced_model)

anova(reduced_model, full_model)


########## Question D ##############
# Cross Validation
train.control <- trainControl(method = "repeatedcv", number = 10, repeats = 50)

full_model_cv <-train(divorce ~ unemployed + femlab + marriage + birth + military, data = divusaF,
                      method = "lm", trControl = train.control)
print(full_model_cv)

reduced_model_cv <-train(divorce ~ femlab + marriage + birth + military, data = divusaF,
                         method = "lm", trControl = train.control)
print(reduced_model_cv)


############ Question E ##############

# Refer report for answer.



######################################
########### Question 3 ##############

########### Question A #############

# Refer report for answer.

########### Question B ############

# Stepwise regression function
step(full_model,direction = "backward")


########### Question C #############
set.seed(850) # so can reproduce results
k = 20 # 20 explanatory variables
k_ln = 5 # 5 explanatory variables are linearly related to Y
B_mag = 0.5 # the magnitude of the explanatory variables that are linearly related to Y
n = 1000 # number of observations in data set

#create a vector containing the population beta coefficients
B <- c(rep(B_mag,k_ln),rep(0,k-k_ln))
X <- matrix(rnorm(n*k), nrow=n) # create the explanatory variables
Y <- X%*%B + matrix(rnorm(n),nrow=n) # create the response variable

# combine the response and explanatory variables in a single dataframe
DF <- cbind(Y,X)
DF<-as.data.frame(DF)

# assign column names
colnames(DF)<-c("Y","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11",
                "X12","X13","X14","X15","X16","X17","X18","X19","X20")

# Fit linear regression model using all the predictors
full_model_sim <- lm(Y ~ .,data = DF)

# the results of the process are stored in step.model
step.model<-stepAIC(full_model_sim, direction = "backward")


# Let’s examine the model selected using stepAIC()
summary(step.model)

















