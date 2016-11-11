#THIS IS COMPLETELY RUBBISH! 
#DO NOT USE
#YOU ARE WARNED!

rm(list=ls())

qtrunc <- function(p, spec, a = -Inf, b = Inf, ...) {
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
  return(tt)
}
rtrunc <- function(n, spec, a = -Inf, b = Inf, ...) {
  x <- u <- runif(n, min = 0, max = 1)
  x <- qtrunc(u, spec, a = a, b = b,...)
  return(x)
}

DA_split_algorithm <- function(x_0,niter,d1,d2,x_star) {
  y <- numeric(length = niter)
  x <- numeric(length = niter+1)
  delta <- numeric(length = niter)
  x[1] <- x_0
  for (i in 1:niter) {
    # conditional distribution of Y given X
    y[i] <- rgamma(1, shape = 5/2, scale = 1/(x_star^2/2+2))
    #y[i] <- rtrunc(1, "gamma", d1, d2, shape=5/2, scale=1/(x_star^2/2+2))
    # conditional distribution of X given Y
    x[i+1] <- rnorm(1, 0, sqrt(y[i]^{-1}))
    #draw delta
    if(y[i] > d2 || y[i] < d1) {
      delta[i] = 0
    } else {
      if(x[i]^2 > x_star^2) {
        prob <- exp((x[i]^2-x_star^2)*(y[i]/2-(d2/2)))
      } else {
        prob <- exp((x[i]^2-x_star^2)*(y[i]/2-(d1/2)))
      }
      delta[i] <- rbinom(1,1,prob = prob)
    }
  }
  return(list(X=x,Y=y, delta=delta))
}

niter <- 1000
res <- DA_split_algorithm(1,niter,0.5,2,0)
hist(res$X[1:niter], 100, col='yellow', las=1, border=T, freq=FALSE, xlim=c(c(-10, 10)))
# comparing with the t distribution with 4 degrees of freedom
hist(rt(niter, df=4), 100, col='red', las=1, border=F, freq=FALSE,  xlim=c(-10, 10), add=TRUE)
hist(res$Y[1:niter], 100, col='yellow', las=1, border=T, freq=FALSE, xlim=c(c(-10, 10)))
