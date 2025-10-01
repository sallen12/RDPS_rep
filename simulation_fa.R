################################################################################
## Simulation Study

source("rdcm_fa.R")
source("rob_rdcm_yprime.R")
source("lspm_fa.R")

set.seed(75210322)
n <- 100
setting <- "nonlinear" # one of "linear" and "nonlinear"

# 1. DATA ################################################

if (setting == "linear") {
  # 1.1 Linear: ############################################
  
  x <- rnorm(n)
  y <- rnorm(n, x)
  x_new <- rnorm(n)
  
  y_sets <- seq(min(y), max(y), 0.1)
  
  true_cdf <- function(y, x) pnorm(y, mean = x, sd = 1)
  
} else if (setting == "nonlinear") {
  # 1.2 Nonlinear: #########################################
  
  x <- runif(n, 0, 10)
  y <- rnorm(n, 2*x + 5*sin(x), abs(x/5))
  x_new <- runif(1, 0, 10)
  
  y_sets <- seq(min(y), max(y), 0.1)

  true_cdf <- function(y, x) pnorm(y, mean = 2*x + 5*sin(x), sd = abs(x/5))
}


#############################################################


# 2. PREDICTIVE BANDS #############################################

# RDCM Predictive Bands
rdcm_result <- rdcm_fa(x, y, x_new, y_sets)
PI_l <- rdcm_result$PI_l
PI_u <- rdcm_result$PI_u

# Robust RDCM
rob_rdcm_yprime_result <- rob_rdcm_fa_yprime(x, y, x_new, y_sets)
rob_PI_l <- rob_rdcm_yprime_result$rob_PI_l
rob_PI_u <- rob_rdcm_yprime_result$rob_PI_u

# LSPM Predictive Bands
lspm_result <- lspm_fa(x, y, x_new, y_sets)
PI_lpm <- lspm_result$PI_l
PI_upm <- lspm_result$PI_u


# 3. PLOT ####################################################

library(ggplot2)

t_cdf <- true_cdf(y_sets, x_new)

df <- data.frame(x = y_sets,
                 l = c(t_cdf, PI_l, rob_PI_l, PI_lpm),
                 u = c(t_cdf, PI_u, rob_PI_u, PI_upm),
                 mth = rep(c(" Truth", "RDPS", "RDPS (rob.)", "LSPM"), each = length(y_sets)))
ggplot(df) + 
  geom_line(aes(x = x, y = l, col = mth)) +
  geom_line(aes(x = x, y = u, col = mth)) +
  scale_y_continuous(name = "G(x)") +
  scale_color_manual(values = c("black", scales::hue_pal()(3))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(), 
        legend.position = "bottom")
ggsave(paste0("Figures/ss_", setting, "_", n, ".png"), width = 5, height = 3.5)

