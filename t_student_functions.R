tStudentPXDA <- function(N,  x_0 =0,  alpha=1, beta=1){
  
# args: N - number of iterations
#       x_0 - starting point for x
#       alpha and beta - parameters; they specify which exactly PXDA we are implementing
#
# output: list of vectors x and y
  
  x <- rep(0, N+1)
  y <- rep(0, N)
  x[1] <- x_0
  
  
  for (i in seq_len(N)){
    # drawing Y given X
    u <- rf(1, 5, 2*alpha) #  f distribution
    y[i] <- u * 10*beta/(alpha*(x[i]^2+4))
    
    # drawing X given Y
    v <- rt(1, 2*alpha+4) #  t distribution
    x[i+1] <- v * sqrt(2*(y[i]+beta)/(y[i]*(alpha+2)))
  }
  
  return(list(x=x, y=y))
}

tStudentDA <- function(N,  x_0 = 0){
  
  # args: N - number of iterations
  #       x_0 - starting point for x
  #       
  # output: list of vectors x and y
  
  x <- rep(0, N+1)
  y <- rep(0, N)
  x[1] <- x_0
  
  
  for (i in seq_len(N)){
    # conditional distribution of Y given X
    y[i] <- rgamma(1, shape = 5/2, scale = 1/(x[i]^2/2+2))
    
    # conditional distribution of X given Y
    x[i+1] <- rnorm(1, 0, sqrt(y[i]^{-1}))
  }
  
  
  return(list(x=x, y=y))
}

# tests
result1 <- tStudentDA(1000, 0)
result2 <- tStudentPXDA(1000, 0, 1,2)
