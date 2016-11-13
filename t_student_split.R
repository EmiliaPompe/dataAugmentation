rm(list=ls())

DA_split_algorithm_tstudent = function(x_0, niter, d1, d2, x_star, VERBOSE = FALSE){
  #C represent lower bound for the p-dimensional rectangle ; D upper bounds
  if(length(C) != length(D) || length(C) != length(x_0)) {
    stop("C,D and beta must have same length!")
  }
  if(length(x_star) != length(x_0)) {
    stop("x_star and x_0 must have same length!")
  }
  Y = numeric(niter)
  X = numeric(niter+1)
  X[1] = x_0
  delta = numeric(niter)
  num_in = 0
  all_eta = numeric(0)
  for(j in 1:niter){
    #Step 1
    Y[j] <- rgamma(1, shape = 5/2, scale = 1/(X[j]^2/2+2))
    #Step 2
    X[j+1] <- rnorm(1, 0, sqrt(Y[j]^{-1}))
    #Step 3: compute delta
    if(all(Y[j] >= d1) && (Y[j] <= d2)) {
      num_in = num_in+1
      #Compute probability eta
      if(X[j]^2 > x_star^2) {
        eta <- exp((X[j]^2-x_star^2)*(Y[j]/2-d2/2))
      } else {
        eta <- exp((X[j]^2-x_star^2)*(Y[j]/2-d1/2))
      }
      eta <- min(c(1,eta))
      if(VERBOSE)
        cat("INSIDE! j =",j,"eta =",eta,"\n")
      delta[j] = rbinom(1,1,eta)
      all_eta <- c(all_eta,eta)
    } else {
      delta[j] = 0
    }
    #Print progress
    if(j%%floor(niter/10)==0)
      cat(j/niter*100,"%\n",sep="")
  }
  cat("Number of times it was in: ",num_in,"(",num_in/niter*100,"%)\n",sep="")
  cat("Probability eta. Mean = ",mean(all_eta),", Min = ",min(all_eta),
      ", Max = ",max(all_eta),"\n",sep="")
  cat("Number of regenerations = ",sum(delta))
  return(list(X=X,Y=Y,delta=delta))
}

niter <- 100000
set.seed(17)
res <- DA_split_algorithm_tstudent(1,niter,0.5,2,0,VERBOSE = FALSE)
hist(res$X[1:niter], 100, col='yellow', las=1, border=T, freq=FALSE, xlim=c(c(-10, 10)), breaks = 200)
# comparing with the t distribution with 4 degrees of freedom
hist(rt(niter, df=4), 100, col='red', las=1, border=F, freq=FALSE,  xlim=c(-10, 10), add=TRUE, breaks=200)
hist(res$Y[1:niter], 100, col='yellow', las=1, border=T, freq=FALSE, xlim=c(c(-10, 10)))
