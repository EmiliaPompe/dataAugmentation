rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("probit_algorithms.R")

#Generate test data
set.seed(13)
beta = matrix(c(-2, 7, 3), ncol=1) # a column vector of regression coefficients
nr = 500
test_data = simulate_data2(beta, nr, intercept=FALSE) #intercept = TRUE leads to "weird" results!

#DA algorithm
beta_0 = matrix(c(1,1,1), ncol=1)
niter=100000
res_DA = DA_algorithm(beta_0, V=test_data$V, Z = test_data$Z, niter)
beta_DA = res_DA$beta
apply(beta_DA, 1, summary)
#Diagnostic DA algorithm
diagnostics(c(-2,7,3), beta_DA)

#PX-DA Algorithm
beta_0 = matrix(c(1,1,1), ncol=1)
niter=10000
alpha = 5
delta = 1
beta_PXDA = PXDA_algorithm(beta_0, V=test_data$V, Z = test_data$Z, niter, alpha, delta)
apply(beta_PXDA, 1, summary)
#Diagnostic PX-DA Algorithm
diagnostics(c(-2,7,3), beta_PXDA)

#Frequentist approach
mydata <- data.frame(test_data$V, test_data$Z)
probit_freq <- glm(formula = test_data.Z ~ X1+X2+X3+0, data = mydata, family=binomial(link="probit"))
summary(probit_freq)

#DA-STEP
recap <- apply(beta_DA[,2:ncol(beta_DA)], 1, summary)
recap
C <- apply(beta_DA[,2:ncol(beta_DA)], 1, quantile, probs=c(.49,.51))[1,]
D <- apply(beta_DA[,2:ncol(beta_DA)], 1, quantile, probs=c(.49,.51))[2,]
beta_0 = matrix(c(1,1,1), ncol=1)
niter=100000
y_star <- apply(res_DA$Y, 2, mean)
set.seed(17)
step_chain <- DA_step_algorithm(beta_0, test_data$V, test_data$Z,
                                niter, C, D, y_star, VERBOSE = FALSE)
apply(step_chain$beta, 1, summary)

#M-H Algorithm
#DOES NOT WORK BECAUSE OF IMPROPER PRIOR!
beta_0 = matrix(c(-2,8,3), ncol=1)
niter=10000
beta_MH <- MH_algorithm(beta_0,test_data[[2]],test_data[[1]],niter)
diagnostics(c(-2,7,3), beta_MH)
