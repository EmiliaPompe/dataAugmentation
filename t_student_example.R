# standard DA for a toy example - simulating from the t-student distribution
# N is the number of iterations; x_0 - starting value for x
N <- 100000
x_0 <- 1
burn_in <- 1000

x <- rep(0, N+1)
y <- rep(0, N)
x[1] <- x_0

set.seed(13)

for (i in seq_len(N)){
  # conditional distribution of Y given X
  y[i] <- rgamma(1, shape = 5/2, scale = 1/(x[i]^2/2+2))
  
  # conditional distribution of X given Y
  x[i+1] <- rnorm(1, 0, sqrt(y[i]^{-1}))
}


# plotting the histogram
par(mfrow=c(1,1))
set.seed(13)

hist(x[burn_in:N], 100, col='yellow', las=1, border=T, freq=FALSE, xlim=c(c(-30, 30)))
# comparing with the t distribution with 4 degrees of freedom
hist(rt(N-burn_in, df=4), 100, col='red', las=1, border=F, freq=FALSE,  xlim=c(-30, 30), add=TRUE)

# they are very similar

# PXDA for our toy example (page 29 in the Chapter)
# apart from N and x_0 we need to set alpha and beta
N <- 100000
x_0 <- 1
burn_in <- 1000
alpha <- 1
beta <- 1

x <- rep(0, N+1)
y <- rep(0, N)
x[1] <- x_0

set.seed(13)

for (i in seq_len(N)){
  # drawing Y given X
  u <- rf(1, 5, 2*alpha) #  f distribution
  y[i] <- u * 10*beta/(alpha*(x[i]^2+4))
  
  # drawing X given Y
  v <- rt(1, 2*alpha+4) #  t distribution
  x[i+1] <- v * sqrt(2*(y[i]+beta)/(y[i]*(alpha+2)))
}

hist(x[burn_in:N], 100, col='yellow', las=1, border=T, freq=FALSE, xlim=c(c(-30, 30)))
# comparing with the t distribution with 4 degrees of freedom
hist(rt(N-burn_in, df=4), 100, col='red', las=1, border=F, freq=FALSE,  xlim=c(-30, 30), add=TRUE)

