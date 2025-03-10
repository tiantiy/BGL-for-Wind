#### This R script provide following codes used in the BGL for Wind data paper:
## 1. Load wind data
## 2. Descriptive statistics of wind data
## 3. Define PDF and CDF functions for BGL distribution
## 4. Calculate Maximum Likelihood Estimation (MLE)
## 5. Draw Histogram of wind data and BGL density fit
## 6. Calculate Kolmogorov-Smirnov (KS) statistic
## 7. Calculate Anderson-Darling (AD) statistic
## 8. Calculate estimated 95th and 99th percentiles and Bias


rm(list=ls())

# set up local folder path
setwd("/restricted/projectnb/adiposity/tiantiy/GNN/engGNN/")
getwd()

###################### 1. Load Wind Data

library(data.table)

M2_wind_spead_2010_2020 <- fread("M2_wind_spead_2010_2020.csv")
dim(M2_wind_spead_2010_2020)

## long-term (Year 2010-2020) wind speed data at heights 10, 20, 50, 80 m
M2_wind_spead_2010_2020_height_10 <- M2_wind_spead_2010_2020[[5]]
length(M2_wind_spead_2010_2020_height_10)
M2_wind_spead_2010_2020_height_20 <- M2_wind_spead_2010_2020[[6]]
length(M2_wind_spead_2010_2020_height_20)
M2_wind_spead_2010_2020_height_50 <- M2_wind_spead_2010_2020[[7]]
length(M2_wind_spead_2010_2020_height_50)
M2_wind_spead_2010_2020_height_80 <- M2_wind_spead_2010_2020[[8]]
length(M2_wind_spead_2010_2020_height_80)

## annual wind speed data for Year 2010, 2015, 2020 at height 80
M2_wind_spead_2010_height_80 <- M2_wind_spead_2010_2020_height_80[1:8760]
length(M2_wind_spead_2010_height_80)
M2_wind_spead_2015_height_80 <- M2_wind_spead_2010_2020_height_80[43825:52584]
length(M2_wind_spead_2015_height_80)
M2_wind_spead_2020_height_80 <- M2_wind_spead_2010_2020_height_80[87649:96432]
length(M2_wind_spead_2020_height_80)


#### Now we use the annual year wind speed data at height 80 for a test run 
xvec <- M2_wind_spead_2010_height_80
(n <- length(xvec))  


###################### 2. Descriptive statistics of wind data

library(moments) # for skewness and kurtosis

## Define a function to calculate and return statistics 
calculate_stats <- function(data) {
  minimum <- min(data)
  maximum <- max(data)
  median <- median(data)
  mean <- mean(data)
  variance <- var(data)
  skewness <- skewness(data)
  kurtosis <- kurtosis(data) # excess kurtosis
  percentile95 <- quantile(data, 0.95)
  percentile99 <- quantile(data, 0.99)
  return(c(minimum, maximum, median, mean, variance, skewness, kurtosis, percentile95, percentile99))
}


DS.results1 <-calculate_stats(xvec)
DS.results1 <- t(DS.results1)
colnames(DS.results1) <- c("Minimum", "Maximum", "Median", "Mean", "Variance", "Skewness", "Kurtosis", "95th", "99th") # Set column names
print(DS.results1)


###################### 3. Define PDF and CDF functions for BGL distribution 


## Define PDF
f1=function(x,a,b,alpha,lambda){
  y = alpha*lambda^2*(1+x)*exp(-lambda*x)*(1-(1+lambda+lambda*x)*exp(-lambda*x)/(1+lambda))^(a*alpha-1)*(1-(1-(1+lambda+lambda*x)*exp(-lambda*x)/(1+lambda))^(alpha))^(b-1)/((1+lambda)*beta(a,b))
  return(y)
}

## Define CDF
F1=function(x,a,b,alpha,lambda){
  y=pbeta((1-(1+lambda+lambda*x)*exp(-lambda*x)/(1+lambda))^alpha,a,b) # pbeta calculates CDF of the beta distribution.
  return(y)
}


###################### 4. Calculate Maximum Likelihood Estimation (MLE)

library(bbmle)

## BGL 
fn1 <- function(alpha,lambda,a,b) {
  -sum( log(alpha*lambda^2/((1+lambda)*beta(a,b)))+log(1+xvec)-lambda*xvec + (a*alpha-1)*log(1-(1+lambda+lambda*xvec)*exp(-lambda*xvec)/(1+lambda)) + (b-1)*log(1-(1-(1+lambda+lambda*xvec)*exp(-lambda*xvec)/(1+lambda))^(alpha)) )
}

## MLE results
mle.results1<-mle2(fn1,start=list(alpha=30.83,lambda=0.2606,a=0.0335,b=0.3))
summary(mle.results1)  
(mle.results1_estimates <- coef(mle.results1))
(p.results1 <- length(coef(mle.results1)))  # Number of parameters
(logL.results1 <- logLik(mle.results1))

## list results for -2log(L), AIC, AICc, and BIC
(AIC.results1 <- c(
  "-2logL" = -2*logL.results1,
  AIC = 2*p.results1 -2*logL.results1,
  AICc = 2*p.results1 -2*logL.results1 + (2*p.results1*(p.results1+1))/(n-p.results1-1),
  BIC = p.results1*log(n)-2*logL.results1
))


###################### 5. Draw Histogram of wind data and BGL density fit

year <- "2010"
height <- 80

## Draw histogram
hist(xvec, breaks=80, 
     xlim=c(0,max(xvec)), ylim=c(0, max(density(xvec)$y) * 1.1), 
     xlab = "Wind Speed (m/s)", ylab="PDF", 
     main= paste("BGL fit for Year", year, "at Height", height, "m"), # Use paste()
     las=1, freq=0)

## BGL density fit
(alpha.results1 = mle.results1_estimates["alpha"])    
(lambda.results1 = mle.results1_estimates["lambda"])   
(a.results1 = mle.results1_estimates["a"])   
(b.results1 = mle.results1_estimates["b"])   
xfit=seq(0,30,length=300)
yfit=f1(xfit,a.results1,b.results1,alpha.results1,lambda.results1)
lines(xfit,yfit,col="red",lwd=1)



###################### 6. Calculate Kolmogorov-Smirnov (KS) statistic

#install.packages("nortest")
library(nortest)

## BGL
(KS.results1 <- ks.test(xvec, F1, a.results1, b.results1, alpha.results1, lambda.results1))
KS.results1$statistic


###################### 7. Calculate Anderson-Darling (AD) statistic

## Define a function for AD statistic
ad.test <- function(data, distr, ...) {
  n <- length(data)
  sorted_data <- sort(data) # Sort the data
  F_x <- distr(sorted_data, ...) # Compute the CDF of the specified distribution
  
  # Calculate the Anderson-Darling test statistic
  statistic <- -n - (1 / n) * sum((2 * (1:n) - 1) * (log(F_x + 10^{-12}) + log(1 + 10^{-12} - rev(F_x)))) 
  
  return(statistic = statistic) 
}

## BGL
(ad_v2.results1 <- ad.test(xvec, F1, a.results1, b.results1, alpha.results1, lambda.results1))


###################### 8. Calculate estimated 95th and 99th percentiles and Bias

## Define the find_quantile for BGL
find_quantile1 <- function(cdf_func, p, a, b, alpha, lambda, lower = 0, upper = 100) {
  uniroot(function(x) cdf_func(x, a, b, alpha, lambda) - p, lower = lower, upper = upper, tol = 1e-8)$root
}

## BGL
(p95.results1 <- round(find_quantile1(F1, 0.95, a.results1, b.results1, alpha.results1, lambda.results1), 4))
(p99.results1 <- round(find_quantile1(F1, 0.99, a.results1, b.results1, alpha.results1, lambda.results1), 4))
(bias95.results1 <- p95.results1-DS.results1[8])
(bias99.results1 <- p99.results1-DS.results1[9])

(percentiles.results1 <- c(DS.results1[8], p95.results1, bias95.results1, DS.results1[9], p99.results1, bias99.results1))
percentiles.results1 <- t(percentiles.results1)
colnames(percentiles.results1) <- c("95th obs", "95th est", "Bias", "99th obs", "99th est", "Bias")
print(percentiles.results1)
