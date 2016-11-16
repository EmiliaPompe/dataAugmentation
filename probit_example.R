rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("probit_algorithms.R")

#Generate test data
set.seed(13)
beta = matrix(c(-2, 7, 3), ncol=1) # a column vector of regression coefficients
nr = 1000
test_data = simulate_data2(beta, nr, intercept=FALSE)

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


###Run times
set.seed(17)
size_observation <- c(200,500,1000,2000)
size_beta <- c(3,10,50,100)

list_run_times <- vector("list", length = length(size_observation)*length(size_beta))
counter_list <- 1

for (nr in size_observation) {
  for(nr_beta in size_beta) {
    beta = matrix(runif(nr_beta,min=-15, max=15), ncol=1)
    test_data = simulate_data2(beta, nr, intercept=FALSE) #intercept = TRUE leads to "weird" results!
    V = test_data$V
    Z = test_data$Z
    alpha = 5
    delta = 1
    niter <- 100000
    beta_0 = matrix(rep(1, nr_beta), ncol=1)
    time_snapshot <- seq(niter/1000,niter,by=niter/1000)
    chain_DA <- DA_algorithm(beta_0, V, Z, niter, time_snapshot)
    chain_PXDA <- PXDA_algorithm(beta_0, V, Z, niter, alpha, delta, time_snapshot)
    
    run_times_matrix = matrix(NA,nrow=length(time_snapshot), ncol=2)
    run_times_matrix[,1] <- chain_DA$time_elapsed
    run_times_matrix[,2] <- chain_PXDA$time_elapsed
    colnames(run_times_matrix) <- c("DA", "PX-DA")
    rownames(run_times_matrix) <- time_snapshot
    
    list_run_times[[counter_list]] <- list(nr=nr, size_beta=nr_beta, beta=beta,
                                           run_time_DA=chain_DA$time_elapsed, run_time_PXDA=chain_PXDA$time_elapsed)
    cat(counter_list,"/",length(size_observation)*length(size_beta),"\n",sep = "")
    save(list_run_times, file=paste("./all_run_times_",counter_list,".RData", sep = ""))
    counter_list <- counter_list+1
    
  }
}

save(list_run_times, file="./all_run_times.RData")


index <- 15
par(mfrow=c(4,4))
for(index in 1:16) {
plot(time_snapshot, log(list_run_times[[index]]$run_time_PXDA, base=10), type="l", 
     ylim=c(0,3.5),
     main = paste("Number of observations: ",list_run_times[[index]]$nr," - Length of beta: ", list_run_times[[index]]$size_beta,sep=""),
     xlab = "Number of iterations", ylab = "log(run time)")
lines(time_snapshot, log(list_run_times[[index]]$run_time_DA, base = 10), col = "red")
}


par(mfrow=c(4,4))
for(index in 1:16) {
  plot(time_snapshot, list_run_times[[index]]$run_time_PXDA, type="l", col="red",
       main = paste("Number of observations: ",list_run_times[[index]]$nr," - Length of beta: ", list_run_times[[index]]$size_beta,sep=""),
       xlab = "Number of iterations", ylab = "Run time (sec)")
  lines(time_snapshot, list_run_times[[index]]$run_time_DA)
}

index <- 16
list_run_times[[index]]$run_time_DA[length(list_run_times[[index]]$run_time_DA)]
list_run_times[[index]]$run_time_PXDA[length(list_run_times[[index]]$run_time_PXDA)]
list_run_times[[index]]$run_time_PXDA[length(list_run_times[[index]]$run_time_PXDA)]/list_run_times[[index]]$run_time_DA[length(list_run_times[[index]]$run_time_DA)]


#ggplot
plot_list <- vector("list", length = 4)
counter <- 1
for(index in c(1,4,13,16)) {
  df = data.frame(iteration = time_snapshot,
                  run_time_DA = list_run_times[[index]]$run_time_DA,
                  run_time_PXDA = list_run_times[[index]]$run_time_PXDA)
  colnames(df) <- c("Iterations","DA","PX-DA")
  library(reshape)
  library(ggplot2)
  library(gridExtra)
  dfmelted = melt(df, id = 'Iterations')
  colnames(dfmelted) <- c("Iterations","Method","value")
  plot_list[[counter]] <- ggplot(dfmelted, aes(x = Iterations, y = value, colour = Method)) + 
    geom_line() + 
    ylab(label="Run time (sec)") + 
    xlab("Number of iterations") + 
    scale_colour_manual(values=c("blue", "red")) +
    ggtitle(paste("m = ",list_run_times[[index]]$nr,", p = ",list_run_times[[index]]$size_beta, sep = ""))
  counter <- counter+1
}
grid.arrange(plot_list[[1]], plot_list[[2]],plot_list[[3]], plot_list[[4]], ncol=2)
