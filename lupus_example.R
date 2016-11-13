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
niter=20000
set.seed(13)
res_DA =  DA_algorithm(beta_0, V=lupus[,2:4], Z = lupus[,1], niter)
beta_DA = res_DA$beta
diagnostics(mle_est, beta_DA)
par(mfrow=c(1,1))
plot(beta_DA[2,1:(niter-1)], beta_DA[2,2:niter])
acf(beta_DA[2,], lag.max = 100)

#STEP-DA Algorithm
C <- apply(beta_DA[,2:ncol(beta_DA)], 1, quantile, probs=c(.49,.51))[1,]
D <- apply(beta_DA[,2:ncol(beta_DA)], 1, quantile, probs=c(.49,.51))[2,]

#the interval that they used in the paper
#it is wider than the one Giulio proposes above so it gives better results
C_prim = apply(beta_DA[,2:ncol(beta_DA)], 1, mean) - 0.09*apply(beta_DA[,2:ncol(beta_DA)], 1, sd)
D_prim = apply(beta_DA[,2:ncol(beta_DA)], 1, mean) + 0.09*apply(beta_DA[,2:ncol(beta_DA)], 1, sd)

step_chain_prim <- DA_step_algorithm(beta_0, V=lupus[,2:4], Z = lupus[,1],
                                     niter, C_prim, D_prim, y_star, VERBOSE = FALSE)
apply(step_chain_prim$beta, 1, summary)



beta_0 = matrix(mle_est, ncol=1)
niter=1000000
y_star <- apply(res_DA$Y, 2, mean)
set.seed(17)
step_chain <- DA_step_algorithm(beta_0, V=lupus[,2:4], Z = lupus[,1],
                                niter, C, D, y_star, VERBOSE = FALSE)
apply(step_chain$beta, 1, summary)


#standard error sigma/sqrt(R)
DA_sderror(step_chain$beta, step_chain$delta, 100)

DA_sderror(step_chain_prim$beta, step_chain_prim$delta, 100)

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
