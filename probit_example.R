rm(list=ls())
source("probit_algorithms.R")

#Generate test data
set.seed(13)
beta = matrix(c(-2, 7, 3), ncol=1) # a column vector of regression coefficients
nr = 1000
test_data = simulate_data(beta, nr)

#DA algorithm
beta_0 = matrix(c(1,1,1), ncol=1)
niter=1000000
beta_DA = DA_algorithm(beta_0, V=test_data[[2]], Z = test_data[[1]], niter)
#Diagnostic DA algorithm
diagnostics(c(-2,7,3), beta_DA)

#PX-DA Algorithm
beta_0 = matrix(c(1,1,1), ncol=1)
niter=1000
alpha = 5
delta = 1
beta_PXDA = PXDA_algorithm(beta_0, V=test_data[[2]], Z = test_data[[1]], niter, alpha, delta)
#Diagnostic PX-DA Algorithm
diagnostics(c(-2,7,3), beta_PXDA)

#Frequentist approach
mydata <- data.frame(test_data$V, test_data$Z)
probit_freq <- glm(formula = test_data.Z ~ X1+X2+X3+0, data = mydata, family=binomial(link="probit"))

#M-H Algorithm
beta_0 = matrix(c(-2,8,3), ncol=1)
niter=10000
beta_MH <- MH_algorithm(beta_0,test_data[[2]],test_data[[1]],niter)
diagnostics(c(-2,7,3), beta_MH)
