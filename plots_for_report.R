setwd('C:\\Emilka\\Studies\\OxWaSP\\Project3\\DA_project\\dataAugmentation')
source('ggplot_functions.R')

library(gridExtra)

rm(df)

set.seed(13)
N <- 100000
chain <- tStudentDA(N, x_0=1)$x
iid_sample <- rt(N+1,4)
df <- data.frame(iteration = rep(1:(N+1), times = 4),
                 x = c(chain, iid_sample, chain^2, iid_sample^2),
                 method = as.factor(rep(c('DA', 'i.i.d.', 'DA', 'i.i.d.' ), each=(N+1))),
                 moment = as.factor(rep(c('First moment', 'First moment', 'Second moment', 'Second moment' ), each=(N+1))))



our_theme <-  theme(
  axis.text = element_text(size = 12),
  #legend.key = element_rect(fill = "navy"),
  #legend.background = element_rect(fill = "white"),
  #legend.position = c(0.14, 0.80),
  strip.text =element_text(size=14),
  legend.text =element_text(size=14),
  title = element_text(size=14))

# plots comparing means
# To DO: decide upon ylim ? (perhaps different ylims for the plots)    
p1 <- PlotXAverage(df[1:(2*N+2),], burn_in =0.1, theoretical_result=0, size =1.2, ggtitle('First moment'))
p1 <- p1 + scale_colour_brewer(palette= 'Set1') + our_theme

p2 <- PlotXAverage(df[(2*N+2+1):(4*N+4),], burn_in =0.1, theoretical_result=2, size =1.2, ggtitle('Second moment'))
p2 <- p2 + scale_colour_brewer(palette= 'Set1') + our_theme
p2

grid.arrange(p1,p2, ncol=2)


#plots for ACF
df_acf <- df[which(df$method=='DA'),]
df_acf$method <- df_acf$moment
p3 <- PlotACF(df_acf, lag.max=20, ggtitle('Autocorrelation functions'))
p3 + our_theme

# density after a certain number of iterations

df_iterations <- df[1:(N+1),]

theoretical_vector <- rt(100000,4)
df2 <- data.frame(x=theoretical_vector)
p4 <- PlotHistGif(df_iterations[1:3,], df2, ggtitle('Iteration 3'), xlim(-10,10), ylim(0,0.6)) + scale_colour_brewer(palette= 'Set1')+ our_theme
p5 <- PlotHistGif(df_iterations[1:50,], df2, ggtitle('Iteration 50'), xlim(-10,10), ylim(0,0.6)) + scale_colour_brewer(palette= 'Set1')+ our_theme
p6 <- PlotHistGif(df_iterations[1:1000,], df2, ggtitle('Iteration 1000'), xlim(-10,10), ylim(0,0.6)) + scale_colour_brewer(palette= 'Set1')+ our_theme

grid.arrange(p4,p5,p6, ncol=3)

# ACF for first moment only

d <- with(acf(df$x[10000:(N+1)], lag.max=10, plot = FALSE), data.frame(lag, acf))
q <- qplot(x = lag, y = acf, data = d, geom = "bar", stat = "identity",
           position = "identity") + ggtitle("Autocorrelation function for the t-student chain") + ylab("ACF")
q +our_theme
