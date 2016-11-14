rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("probit_algorithms.R")

#Generate test data
set.seed(13)
beta = matrix(c(-2, 7, 3), ncol=1) # a column vector of regression coefficients
nr = 1000
test_data = simulate_data2(beta, nr, intercept=FALSE) #intercept = TRUE leads to "weird" results!

#DA algorithm
beta_0 = matrix(c(1,1,1), ncol=1)
niter=100000
res_DA = DA_algorithm(beta_0, V=test_data$V, Z = test_data$Z, niter)
beta_DA = res_DA$beta
apply(beta_DA, 1, summary)
#Diagnostic DA algorithm
diagnostics(c(-2,7,3), beta_DA)

Brier_score(beta_DA, V=test_data$V, Z = test_data$Z)

#PX-DA Algorithm
beta_0 = matrix(c(1,1,1), ncol=1)
niter=100000
alpha = 5
delta = 1
beta_PXDA = PXDA_algorithm(beta_0, V=test_data$V, Z = test_data$Z, niter, alpha, delta)
apply(beta_PXDA, 1, summary)
#Diagnostic PX-DA Algorithm
diagnostics(c(-2,7,3), beta_PXDA)
Brier_score(beta_PXDA, V=test_data$V, Z = test_data$Z)

#Brier score for the DA and PXDA chains 
Brier_DA = numeric(100)
Brier_PXDA = numeric(100)
for(i in 1:100){
  Brier_DA[i] = Brier_score(beta_DA[,1:(i*1000)], V=test_data$V, Z = test_data$Z)
  Brier_PXDA[i] = Brier_score(beta_PXDA[,1:(i*1000)], V=test_data$V, Z = test_data$Z)
}
plot(seq(1000,100000,by=1000), Brier_DA, col="red", type="l")
lines(seq(1000,100000,by=1000), Brier_PXDA)


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
split_chain_sim <- DA_split_algorithm_original(beta_0, test_data$V, test_data$Z,
                                niter, C, D, y_star, VERBOSE = FALSE)
apply(split_chain_sim$beta, 1, summary)



###run times####

all_iters = seq(10000,100000, by=10000)
run_times_matrix = matrix(NA,nrow=length(all_iters), ncol=2)
counter = 1

set.seed(13)
beta = matrix(c(-2, 7, 3), ncol=1) # a column vector of regression coefficients
nr = 1000
test_data = simulate_data2(beta, nr, intercept=FALSE) #intercept = TRUE leads to "weird" results!
V = test_data$V
Z = test_data$Z
alpha = 5
delta = 1

for(i in all_iters){
  aux = run_times(beta_0, V, Z, i, alpha, delta)
  run_times_matrix[counter,] = c(as.numeric(aux$run_time_DA), as.numeric(aux$run_time_PXDA))
  counter=counter+1
}

plot(all_iters, log(run_times_matrix[,1], base=10), col="red", ylim=c(0, max(log(run_times_matrix, base=10)))) #DA
points(all_iters, log(run_times_matrix[,2], base=10)) #PXDA
#about 15 times slower for PXDA


save(run_times_matrix, file="/data/tinamou/kindalov/Project3/dataAugmentation/run_times.RData")
