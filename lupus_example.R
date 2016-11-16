rm(list=ls())
source("probit_algorithms.R")
library(TruncatedNormal)
data("lupus")

#Probit
probit_freq <- glm(formula = response ~ const+x1+x2+0, data = as.data.frame(lupus), family=binomial(link="probit"))
mle_est <- probit_freq$coefficients
prob = pnorm(lupus[,2:4]%*%mle_est)
mean((prob-lupus[,1])^2) #Brier score

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

#SPLIT-DA Algorithm
# C <- apply(beta_DA[,2:ncol(beta_DA)], 1, quantile, probs=c(.49,.51))[1,]
# D <- apply(beta_DA[,2:ncol(beta_DA)], 1, quantile, probs=c(.49,.51))[2,]

#the interval that they used in the paper
#it is wider than the one Giulio proposes above so it gives better results
C_prim = apply(beta_DA[,2:ncol(beta_DA)], 1, mean) - 0.09*apply(beta_DA[,2:ncol(beta_DA)], 1, sd)
D_prim = apply(beta_DA[,2:ncol(beta_DA)], 1, mean) + 0.09*apply(beta_DA[,2:ncol(beta_DA)], 1, sd)
beta_0 = matrix(mle_est, ncol=1)
niter=1000000
y_star <- apply(res_DA$Y, 2, mean)

split_chain_prim <- DA_split_algorithm_original(beta_0, V=lupus[,2:4], Z = lupus[,1],
                                     niter, C_prim, D_prim, y_star, VERBOSE = FALSE)
apply(step_chain_prim$beta, 1, summary)

DA_sderror(split_chain_prim$beta, split_chain_prim$delta, 99)

Brier_score(beta_DA, V=lupus[,2:4], Z = lupus[,1])
# 
# #alter R
# set.seed(91)
# system.time({split_chain_prim <- DA_split_algorithm(beta_0, V=lupus[,2:4], Z = lupus[,1],
#                                        R=10, C=C_prim, D=D_prim, y_star=y_star, VERBOSE = FALSE)})
# 
# 
# set.seed(91)
# niter <- 151292
# system.time({step_chain <- DA_split_algorithm_original(beta_0, V=lupus[,2:4], Z = lupus[,1],
#                                 niter, C_prim, D_prim, y_star, VERBOSE = FALSE)})
# apply(step_chain$beta, 1, summary)


#standard error sigma/sqrt(R)
# DA_sderror(step_chain$beta, step_chain$delta, 100)





#PX-DA algorithm
beta_0 = matrix(mle_est, ncol=1)
niter=1000000
alpha = 5
delta = 1
beta_PXDA = PXDA_algorithm(beta_0, V=lupus[,2:4], Z = lupus[,1], niter, alpha, delta)
diagnostics(mle_est, beta_PXDA)
par(mfrow=c(1,1))
plot(beta_PXDA[2,1:(niter-1)], beta_PXDA[2,2:niter])
acf(beta_PXDA[2,], lag.max = 100)

Brier_score(beta_PXDA, V=lupus[,2:4], Z = lupus[,1])

Brier_DA = numeric(100)
Brier_PXDA = numeric(100)

for(i in 1:100){
  Brier_DA[i] = Brier_score(split_chain_prim$beta[,1:(i*1000)], V=lupus[,2:4], Z = lupus[,1])
  Brier_PXDA[i] = Brier_score(beta_PXDA[,1:(i*1000)], V=lupus[,2:4], Z = lupus[,1])
}
plot(seq(100,10000,by=100), Brier_DA, col="red", type="l")
lines(seq(100,10000,by=100), Brier_PXDA)

#Density plots
niter = dim(beta_DA)[2]
plot(density(beta_DA[1,floor(0.1*niter):niter]))
plot(density(beta_DA[2,floor(0.1*niter):niter]))
plot(density(beta_DA[3,floor(0.1*niter):niter]))

plot(density(beta_PXDA_12[1,floor(0.1*niter):niter]))
plot(density(beta_PXDA_12[2,floor(0.1*niter):niter]))
plot(density(beta_PXDA_12[3,floor(0.1*niter):niter]))


require(reshape)
require(ggplot2)
require(gridExtra)
niter <- 100001
plot_list <- vector("list", length = 3)
for(i in 1:3) {
  df = data.frame(method = c("DA","PX-DA 12"),
                  beta_DA = beta_DA[i,floor(0.1*niter):niter],
                  beta_PXDA_12 = beta_PXDA_12[i,floor(0.1*niter):niter])
  dfmelted = melt(df, id = 'method')
  plot_list[[i]] <- ggplot(dfmelted, aes(x=value)) + geom_density(aes(group=method, col=method, fill=method), alpha=0.3)
}
grid.arrange(plot_list[[1]], plot_list[[2]],plot_list[[3]], ncol=3)
