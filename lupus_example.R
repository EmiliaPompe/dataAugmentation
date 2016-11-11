rm(list=ls())
source("probit_algorithms.R")
library(TruncatedNormal)
data("lupus")

#Probit
probit_freq <- glm(formula = response ~ const+x1+x2+0, data = as.data.frame(lupus), family=binomial(link="probit"))
mle_est <- probit_freq$coefficients

#Test likelihood
probit_loglikelihood(mle_est, V=lupus[,2:4], Z = lupus[,1])

#M-H algorithm
beta_0 = matrix(mle_est, ncol=1)
niter=100000
beta_MH <- MH_algorithm(beta_0, V=lupus[,2:4], Z = lupus[,1], niter)
diagnostics(mle_est, beta_MH)

#DA algorithm
beta_0 = matrix(mle_est, ncol=1)
niter=1000
beta_DA = DA_algorithm(beta_0, V=lupus[,2:4], Z = lupus[,1], niter)
diagnostics(mle_est, beta_DA)
par(mfrow=c(1,1))
plot(beta_DA[2,1:(niter-1)], beta_DA[2,2:niter])
acf(beta_DA[2,], lag.max = 100)


#PX-DA algorithm
beta_0 = matrix(mle_est, ncol=1)
niter=1000
alpha = 5
delta = 1
beta_PXDA = PXDA_algorithm(beta_0, V=lupus[,2:4], Z = lupus[,1], niter, alpha, delta)
diagnostics(mle_est, beta_PXDA)
par(mfrow=c(1,1))
plot(beta_PXDA[2,1:(niter-1)], beta_PXDA[2,2:niter])
acf(beta_PXDA[2,], lag.max = 100)
