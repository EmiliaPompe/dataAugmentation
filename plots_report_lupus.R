setwd('C:\\Emilka\\Studies\\OxWaSP\\Project3\\DA_project\\dataAugmentation')
source('ggplot_functions.R')

library(gridExtra)
load('lupus_chains_data2.RData')

n <- 99000
N <- 100000
df_beta1 <- data.frame( iteration = rep((n+1):N, times = 3),
                        x = c(beta_DA[1,(n+1):N], beta_PXDA_12[1, (n+1):N], beta_PXDA_51[1,(n+1):N]),
                        method = rep(c('DA', 'PXDA(1,2)',  'PXDA(5,1)'), each=N-n))

df_beta2 <- data.frame( iteration = rep((n+1):N, times = 3),
                        x = c(beta_DA[2,(n+1):N], beta_PXDA_12[2, (n+1):N], beta_PXDA_51[2,(n+1):N]),
                        method = rep(c('DA', 'PXDA(1,2)',  'PXDA(5,1)'), each=N-n))

df_beta3 <- data.frame( iteration = rep((n+1):N, times = 3),
                        x = c(beta_DA[3,(n+1):N], beta_PXDA_12[3, (n+1):N], beta_PXDA_51[3,(n+1):N]),
                        method = rep(c('DA', 'PXDA(1,2)',  'PXDA(5,1)'), each=N-n))


our_theme <-  theme(
  axis.text = element_text(size = 12),
  #legend.key = element_rect(fill = "navy"),
  #legend.background = element_rect(fill = "white"),
  #legend.position = c(0.14, 0.80),
  strip.text =element_text(size=14),
  legend.text =element_text(size=14),
  title = element_text(size=14))

p1<- PlotXFunction(df_beta1, NULL, 1, ylab(expression(paste( beta,"0"))), ggtitle(expression(paste("Trace plot, ", beta,"0"))))
p1 <- p1  + scale_colour_brewer(  palette="Set1") + our_theme
p1

p2<- PlotXFunction(df_beta2, NULL, 1, ylab(expression(paste( beta,"1"))), ggtitle(expression(paste("Trace plot, ", beta,"1"))))
p2 <- p2  + scale_colour_brewer( palette="Set1") + our_theme
p2

x11()
p3<- PlotXFunction(df_beta3, NULL, 1, ylab(expression(paste( beta_{2},"2"))), ggtitle(expression(paste("Trace plot, ", beta,"2"))))
#p3 <- p3 + scale_colour_brewer(palette= 'Set1') + our_theme 
p3 + scale_colour_brewer("method",
                      breaks = c("DA", "PXDA_12", "PXDA_51"),
                      labels = c("DA", "PXDA(1,2)", "PXDA(5,1)"), 
                      palette="Set1") + our_theme


# plots ACF

our_theme2 <-  theme(
  axis.text.y = element_text(size = 12),
  axis.text.x=element_blank(),
  axis.ticks.x = element_blank(),
  #legend.key = element_rect(fill = "navy"),
  #legend.background = element_rect(fill = "white"),
  #legend.position = c(0.14, 0.80),
  strip.text =element_text(size=14),
  legend.text =element_text(size=14),
  title = element_text(size=14))


n <- 120000
N <- 1200000
df_beta1 <- data.frame( iteration = rep((n+1):N, times = 3),
                        x = c(beta_DA[1,(n+1):N], beta_PXDA_12[1, (n+1):N], beta_PXDA_51[1,(n+1):N]),
                        method = rep(c('DA', 'PXDA(1,2)',  'PXDA(5,1)'), each=N-n))

df_beta2 <- data.frame( iteration = rep((n+1):N, times = 3),
                        x = c(beta_DA[2,(n+1):N], beta_PXDA_12[2, (n+1):N], beta_PXDA_51[2,(n+1):N]),
                        method = rep(c('DA', 'PXDA(1,2)',  'PXDA(5,1)'), each=N-n))

df_beta3 <- data.frame( iteration = rep((n+1):N, times = 3),
                        x = c(beta_DA[3,(n+1):N], beta_PXDA_12[3, (n+1):N], beta_PXDA_51[3,(n+1):N]),
                        method = rep(c('DA', 'PXDA(1,2)',  'PXDA(5,1)'), each=N-n))

p4 <- PlotACF(df_beta1, lag.max=100, ylab(expression(paste("ACF for ", beta,"0")))) + our_theme2
p5 <- PlotACF(df_beta2, lag.max=100, ylab(expression(paste("ACF for ", beta,"1")))) + our_theme2
p6 <- PlotACF(df_beta3, lag.max=100, ylab(expression(paste("ACF for ", beta,"2")))) + our_theme2
p4
p5
p6
