source('t_student_functions.R')

# example data frame
N <- 50000
df <- data.frame(iteration = rep(1:(N+1), times = 5),
                 x = c(tStudentDA(N)$x, tStudentPXDA(N, alpha=1, beta=1)$x,
                       tStudentPXDA(N, alpha=1, beta=0.1)$x, tStudentPXDA(N, alpha=1, beta=0.5)$x, 
                       rt(N+1, 4)),
                 method = as.factor(rep(c('DA', 'PXDA (1,1)', 'PXDA (1,10)', 'PXDA (1,0.5)', 'IID' ), each=(N+1))))
                 


PlotXFunction <- function(df, theoretical_result = NULL, size = 1,  ...){
  require(ggplot2)
   l <- list(...)
   p <- ggplot(df, aes(x=iteration, y=x, colour= method), ...) + geom_line(size=size) + l
   
   if(!is.null(theoretical_result)){
     p <- p + geom_hline(yintercept = theoretical_result, color = 'black')
   }
   p
}

# not very clear for comparing more than one method at the same time
PlotXFunction(df[1:4000,], NULL, 0.5, ggtitle('Value of X vs. the number of iteration'), 
              xlab('number of iteration'))


PlotXAverage <- function(df,burn_in = 0.1, theoretical_result = NULL, size = 1,  ...){
  require(ggplot2)
  require(plyr)
  
  av_df <- ddply(df, .(method), mutate, average = cumsum(x)/iteration)
  N <- nrow(df)/length(unique(df$method))
  burn_in_df <- av_df[which(av_df$iteration >= burn_in*N),]
  
  l <- list(...)
  p <- ggplot(burn_in_df, aes(x=iteration, y=average, colour = method), ...) + geom_line(size=size) + l
  
  if(!is.null(theoretical_result)){
    p <- p + geom_hline(yintercept = theoretical_result, color = 'black')
  }
  p
}

#average
PlotXAverage(df,0.3, 0, 1, ggtitle('Average vs. the number of iterations'), 
             xlab('number of iterations'))


PlotX2Average <- function(df,burn_in = 0.1, theoretical_result = NULL, size = 1,  ...){
  require(ggplot2)
  require(plyr)
  
  av_df <- ddply(df, .(method), mutate, average = cumsum(x^2)/iteration)
  N <- nrow(df)/length(unique(df$method))
  burn_in_df <- av_df[which(av_df$iteration >= burn_in*N),]
  
  l <- list(...)
  p <- ggplot(burn_in_df, aes(x=iteration, y=average, colour = method), ...) + geom_line(size=size) + l
  
  if(!is.null(theoretical_result)){
    p <- p + geom_hline(yintercept = theoretical_result, color = 'black')
  }
  p
}

#second moment
PlotX2Average(df,0.3, 2, 1, ggtitle('Second moment vs. the number of iterations'), 
             xlab('number of iterations'))

PlotHist <- function(df, burn_in = 0.1, ...){
  require(ggplot2)
  l <- list(...)
  
  N <- nrow(df)/length(unique(df$method))
  burn_in_df <- df[which(df$iteration >= burn_in*N),]
  
  p <- ggplot(burn_in_df, aes(x=x)) + geom_density() + facet_wrap(~method) + l
  p
}

PlotHist(df, 0.1)


ACFFunction  <- function(v, method, lag.max = NULL){
  # the next line creates autocorrelations with different lags
  d <- with(acf(v, lag.max, plot = FALSE), data.frame(lag, acf))
  
  colnames(d) <- c('lag', 'acf')
  d$method <- method
  d
}



PlotACF <-function(df, lag.max = NULL,  ...){
  require(ggplot2)
  require(plyr)
  l <- list(...)
  
  lag_list <- dlply(df, .(method), summarize, 
                    ACFFunction(x, method=method[1], lag.max=5))
  lag_df <- do.call("rbind", lapply(lag_list, function(x) x[,1]))
  lag_df$lag <- as.factor(lag_df$lag)
  
  
  q <- qplot(x = lag, y = acf, data = lag_df, geom = "bar", stat = "identity",
             position = "identity") + facet_wrap(~method) + l
  q
}
PlotACF(df, 10, ggtitle('Áutocorrelation function plot')) 

