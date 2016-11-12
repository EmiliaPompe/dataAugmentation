library(MASS)
library(truncnorm)

simulate_data = function(beta, nr, intercept=FALSE){
  if(intercept) {
    stop("Sorry! Intercept = TRUE is still not implemented ^^'")
  }
  Z = numeric(nr)
  V = matrix(rnorm(nr*length(beta), 0, 4), nrow=nr, ncol=length(beta))
  prob = pnorm(V%*%beta)
  for (i in 1:length(prob))
    Z[i] = rbinom(1,1,prob[i])
  list(Z=Z, V=V, prob=prob)
}

simulate_data2= function(beta, nr, intercept=TRUE) {
  #If intercept equals TRUE, the first beta is treated as intercept
  V = matrix(NA, nrow=nr, ncol=length(beta))
  if(intercept) {
    V[,1] = 1
  }
  prob <- runif(nr,0,1)
  #Create observations that generates such result (in a probably stupid way)
  for (i in 1:nr) {
    if(intercept == FALSE) {
      #Choose a coordinate at random
      index <- sample(1:length(beta),1)
      #Choose a random value for all other coordinates
      V[i,-index] <- rnorm(length(beta)-1,0,4)
      #Find the value for the remaining value
      V[i,index] <- (prob[i]-beta[-index]%*%V[i,-index])/beta[index]
    } else {
      #Choose a coordinate at random
      index <- sample(2:length(beta),1)
      #Choose a random value for all other coordinates
      V[i,c(-1,-index)] <- rnorm(length(beta)-2,0,4)
      #Find the value for the remaining value
      V[i,index] <- (prob[i]-beta[-index]%*%V[i,-index])/beta[index]
    }
  }
  #Sample Z
  final_prob <- V%*%beta
  Z = numeric(nr)
  for(i in 1:nr) {
    Z[i] = rbinom(1,1,final_prob[i])
  }
  list(Z=Z, V=V, prob=final_prob)
}

DA_algorithm = function(beta_0, V, Z, niter){
  Y = matrix(nrow=niter, ncol=length(Z))
  beta = matrix(nrow=length(beta_0), ncol=niter+1)
  beta[,1] = beta_0
  cov_beta = solve(t(V)%*%V)
  
  for(j in 1:niter){
    
    Y[j,Z==0] = rtruncnorm(sum(Z==0), b=0, mean=V[Z==0,]%*%beta[,j], sd=1)
    Y[j,Z==1] = rtruncnorm(sum(Z==1), a=0, mean=V[Z==1,]%*%beta[,j], sd=1)
    # for(i in 1:length(Z)){
    #   if(Z[i]==0)
    #     Y[i] = rtruncnorm(1, b=0, mean=V[i,]%*%beta[,j], sd=1)
    #   if(Z[i]==1)
    #     Y[i] = rtruncnorm(1, a=0, mean=V[i,]%*%beta[,j], sd=1)
    # }
    #print(V[i,]%*%beta[,j])
    
    beta_hat = cov_beta%*%t(V)%*%Y[j,]
    
    beta[,j+1] = mvrnorm(1, beta_hat, cov_beta)
    #Print progress
    if(j%%floor(niter/10)==0)
      cat(j/niter*100,"%\n",sep="")
  }
  
  return(list(beta=beta,Y=Y))
}

DA_step_algorithm = function(beta_0, V, Z, niter, C, D, y_star, VERBOSE = FALSE){
  #C represent lower bound for the p-dimensional rectangle ; D upper bounds
  if(length(C) != length(D) || length(C) != length(beta_0)) {
    stop("C,D and beta must have the sime length!")
  }
  if(dim(V)[2] != length(beta_0)) {
    stop("Dimension of V and beta are not compatible")
  }
  if(dim(V)[1] != length(Z)) {
    stop("Dimension of V and Z are not compatible")
  }
  if(length(y_star) != length(Z)) {
    stop("Dimension of y_star and Z are not compatible")
  }
  Y = matrix(nrow=niter, ncol=length(Z))
  beta = matrix(nrow=length(beta_0), ncol=niter+1)
  beta[,1] = beta_0
  cov_beta = solve(t(V)%*%V)
  delta = numeric(niter)
  num_in = 0
  all_eta = numeric(0)
  for(j in 1:niter){
    #Step 1
    Y[j,Z==0] = rtruncnorm(sum(Z==0), b=0, mean=V[Z==0,]%*%beta[,j], sd=1)
    Y[j,Z==1] = rtruncnorm(sum(Z==1), a=0, mean=V[Z==1,]%*%beta[,j], sd=1)
    #Step 2
    beta_hat = cov_beta%*%t(V)%*%Y[j,]
    beta[,j+1] = mvrnorm(1, beta_hat, cov_beta)
    #Step 3: compute delta
    if(all(beta[,j+1] >= C) && (beta[,j+1] <= D)) {
      num_in = num_in+1
      #Compute probability eta
      t_vec = t(Y[j,]-y_star)%*%V
      aux_sum = 0
      for (i in 1:length(beta_0)) {
        if(t_vec[i]  >= 0) {
          aux_sum <- aux_sum + C[i]*t_vec[i]
        } else {
          aux_sum <- aux_sum + D[i]*t_vec[i]
        }
        aux_sum <- aux_sum - t_vec[i]*beta[i,j+1]
      }
      eta <- min(c(1,exp(aux_sum)))
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
  return(list(beta=beta,Y=Y,delta=delta))
}

diagnostics = function(true_beta, beta){
  n = ncol(beta)
  new_beta = beta[,floor(n*0.1):n]
  beta_ordered = matrix(ncol=dim(new_beta)[2] , nrow=dim(new_beta)[1])
  result = list()
  
  #create trace plots
  
  plot(1:n, beta[1,], main="beta 1")
  lines(1:n, rep(true_beta[1], n), type="l")
  hist(beta[1,])
  
  plot(1:n, beta[2,], main="beta 2")
  lines(1:n, rep(true_beta[2],n), type="l")
  hist(beta[2,])
  
  plot(1:n, beta[3,], main="beta 3")
  lines(1:n, rep(true_beta[3],n), type="l")
  hist(beta[3,])
  
  mean = numeric(3)
  lower = numeric(3)
  upper = numeric(3)
  
  #introduce a burnin period of 0.1*n
  
  for(i in 1:3){
    mean[i] = mean(new_beta[i,])
    beta_ordered[i,] = new_beta[i,order(new_beta[i,])]
    
    #quantile based credible interval
    lower[i] = beta_ordered[i, floor(0.05*0.9*n)]
    upper[i] = beta_ordered[i, floor(0.95*0.9*n)]
    
    result[[i]] = c(lower[i], mean[i], upper[i])
  }
  print(result)
}

PXDA_algorithm = function(beta_0, V, Z, niter, alpha, delta){
  m = length(Z)
  Y = numeric(m)
  beta = matrix(nrow=length(beta_0), ncol=niter+1)
  beta[,1] = beta_0
  
  cov_beta = solve(t(V)%*%V)
  H = V%*%cov_beta%*%t(V)
  I = diag(1, nrow=nrow(H))
  
  u = numeric(1)
  v = numeric(1)
  Y_tilda = numeric(m)
  
  all_u = rgamma(niter, alpha, delta)
  for(j in 1:niter){
    #step 1
    Y[Z==0] = rtruncnorm(sum(Z==0), b=0, mean=V[Z==0,]%*%beta[,j], sd=1)
    Y[Z==1] = rtruncnorm(sum(Z==1), a=0, mean=V[Z==1,]%*%beta[,j], sd=1)
    # for(i in 1:m){
    #   if(Z[i]==0)
    #     Y[i] = rtruncnorm(1, b=0, mean=V[i,]%*%beta[,j], sd=1)
    #   if(Z[i]==1)
    #     Y[i] = rtruncnorm(1, a=0, mean=V[i,]%*%beta[,j], sd=1)
    # }
    
    #step 2
    u = all_u[j]
    Y_tilda = Y / sqrt(u)
    v = rgamma(1, (m/2 + alpha), (t(Y_tilda)%*%(I-H)%*%Y_tilda)/2 + delta)
    
    Y_prim = sqrt(v)*Y_tilda
    
    #step 3
    beta_hat = cov_beta%*%t(V)%*%Y_prim
    
    beta[,j+1] = mvrnorm(1, beta_hat, cov_beta)
    if(j%%floor(niter/10)==0)
      cat(j/niter*100,"%\n",sep="")
  }
  
  return(beta)
}

probit_loglikelihood <- function(beta, V, Z) {
  #beta must be a column vector
  # 
  # 
  # 
  # result <- 0
  # for(i in 1:length(Z)) {
  #   if(Z[i] == 1) {
  #     aux <- log(pnorm(V[i,]%*%beta_prop))-log(pnorm(V[i,]%*%beta))
  #   } else {
  #     aux <- log((1-pnorm(V[i,]%*%beta_prop)))-log((1-pnorm(V[i,]%*%beta)))
  #   }
  #   if(is.nan(aux))
  #     aux <- 0
  #   result <- result + aux
  # }
  
  aux <- pnorm(V%*%beta)
  result <- sum(log(aux[Z==1]))+sum(log(1-aux[Z==0]))
    
    
  return(result)
}

MH_algorithm <- function(beta_0, V, Z, niter) {
  beta = matrix(nrow=length(beta_0), ncol=niter+1)
  beta[,1] = beta_0
  for(j in 1:niter) {
    beta_prop <- mvrnorm(1, beta[,j], 0.1*diag(length(beta_0)))
    alpha <- exp(probit_loglikelihood(beta_prop,V,Z)/probit_loglikelihood(beta[,j],V,Z))
    u <- runif(1,0,1)
    if(is.nan(alpha) == FALSE && u < alpha) {
      beta[,j+1] <- beta_prop
    } else {
      beta[,j+1] <- beta[,j]
    }
  }
  return(beta)
}