source("rdcm_fa.R")
source("rob_rdcm_yprime.R")
source("lspm_fa.R")

# Simulation Study

# 1. SETTING ################################################

# 1.1 Linear Gaussian: #######################################
set.seed(123)
n <- 20
xvect <- rnorm(n)
yvect <- 2 * xvect + rnorm(n)
xnew <- rnorm(1)

ysets=seq(min(yvect),max(yvect),0.1)

# True CDF
true_cdf <- function(y) pnorm(y, mean = 2 * xnew, sd = 1)
t_cdf=true_cdf(ysets)
#############################################################

# 1.2 Less isotonic: #########################################
set.seed(123)
n<-20
xvect <- runif(n,0,10)
ymean <- 2*xvect+5*sin(xvect)
yvar <- (xvect/5)^2
yvect <- rnorm(n,ymean,yvar)
xnew <- rnorm(1,mean(xvect))

ysets=seq(min(yvect),max(yvect),0.1)

# True CDF
true_cdf <- function(y, xnew) {
  mu <- 2 * xnew + 5 * sin(xnew)
  sigma <- xnew / 5       
  pnorm(y, mean = mu, sd = sigma)
}
t_cdf <- true_cdf(ysets, xnew)

#############################################################


# 2. PREDICTIVE BANDS #############################################

# RDCM Predictive Bands
rdcm_result=rdcm_fa(xvect,yvect,xnew,ysets)
PI_l=rdcm_result$PI_l
PI_u=rdcm_result$PI_u

# Robust RDCM
rob_rdcm_yprime_result=rob_rdcm_fa_yprime(xvect,yvect,xnew,ysets)
rob_PI_l=rob_rdcm_yprime_result$rob_PI_l
rob_PI_u=rob_rdcm_yprime_result$rob_PI_u

# LSPM Predictive Bands
lspm_result=lspm_fa(xvect,yvect,xnew,ysets)
PI_lpm=lspm_result$PI_l
PI_upm=lspm_result$PI_u


# 3. PLOT ####################################################

# plot of the real cdf
plot(ysets,t_cdf,type='l')

# rdcm
lines(ysets, PI_l, type = "l", col = "red", lwd = 1, lty = 1)
lines(ysets, PI_u, type = "l", col = "red", lwd = 1, lty = 1)

# robust rdcm
lines(ysets, rob_PI_l, type = "l", col = "green", lwd = 1, lty = 1)
lines(ysets, rob_PI_u, type = "l", col = "green", lwd = 1, lty = 1)

# lspm
lines(ysets, PI_lpm, type = "l", col = "blue", lwd = 1, lty = 1)
lines(ysets, PI_upm, type = "l", col = "blue", lwd = 1, lty = 1)

legend('topleft',legend = c("cdf", "rdcm", "rob_rdcm", "lspm"),
       col = c("black", "red", "green", "blue"),lty = 1, 
       lwd = 1, cex = 0.7, bty = "n")
