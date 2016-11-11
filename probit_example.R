library(MASS)
library(truncnorm)

simulate_data = function(beta, nr){
  Z = numeric(nr)
  V = matrix(rnorm(nr*length(beta), 0, 4), nrow=nr, ncol=length(beta))
  prob = pnorm(V%*%beta)
  for (i in 1:length(prob))
    Z[i] = rbinom(1,1,prob[i])
  list(Z=Z, V=V)
}

beta = matrix(c(-2, 7, 3), ncol=1) # a column vector of regression coefficients
nr = 1000
test_data = simulate_data(beta, nr)


DA_algorithm = function(beta_0, V, Z, niter){
  Y = numeric(length(Z))
  beta = matrix(nrow=length(beta), ncol=niter+1)
  beta[,1] = beta_0
  cov_beta = solve(t(V)%*%V)
  
  
  for(j in 1:niter){
    
    for(i in 1:length(Z)){
      if(Z[i]==0)
        Y[i] = rtruncnorm(1, b=0, mean=V[i,]%*%beta[,j], sd=1)
      if(Z[i]==1)
        Y[i] = rtruncnorm(1, a=0, mean=V[i,]%*%beta[,j], sd=1)
    }
    #print(V[i,]%*%beta[,j])
    
    beta_hat = cov_beta%*%t(V)%*%Y
    
    beta[,j+1] = mvrnorm(1, beta_hat, cov_beta)
    if(j%%10000==0)
      print(beta[,j])
  }

  return(beta)
}

beta_0 = matrix(c(1,1,1), ncol=1)
niter=100000
beta_DA = DA_algorithm(beta_0, V=test_data[[2]], Z = test_data[[1]], niter)


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

diagnostics(c(-2,7,3), beta_DA)


PXDA_algorithm = function(beta_0, V, Z, niter, alpha, delta){
  m = length(Z)
  Y = numeric(m)
  beta = matrix(nrow=length(beta), ncol=niter+1)
  beta[,1] = beta_0
  
  cov_beta = solve(t(V)%*%V)
  H = V%*%cov_beta%*%t(V)
  I = diag(1, nrow=nrow(H))
  
  u = numeric(1)
  v = numeric(1)
  Y_tilda = numeric(m)
  
  for(j in 1:niter){
    #step 1
    for(i in 1:m){
      if(Z[i]==0)
        Y[i] = rtruncnorm(1, b=0, mean=V[i,]%*%beta[,j], sd=1)
      if(Z[i]==1)
        Y[i] = rtruncnorm(1, a=0, mean=V[i,]%*%beta[,j], sd=1)
    }
    
    #step 2
    u = rgamma(1, alpha, delta)
    Y_tilda = Y / sqrt(u)
    v = rgamma(1, (m/2 + alpha), (t(Y_tilda)%*%(I-H)%*%Y_tilda)/2 + delta)
    
    Y_prim = sqrt(v)*Y_tilda
    
    #step 3
    beta_hat = cov_beta%*%t(V)%*%Y_prim
    
    beta[,j+1] = mvrnorm(1, beta_hat, cov_beta)
    if(j%%10000==0)
      print(beta[,j])
  }
  
  return(beta)
}

beta_0 = matrix(c(1,1,1), ncol=1)
niter=100000
alpha = 5
delta = 1
beta_PXDA = PXDA_algorithm(beta_0, V=test_data[[2]], Z = test_data[[1]], niter, alpha, delta)

diagnostics(c(-2,7,3), beta_PXDA)
