################################################################################
## Simulation Study

source("rdcm_fa.R")
source("rob_rdcm_yprime.R")
source("lspm_fa.R")

set.seed(75210322)
n_tr <- 100
n_ts <- 100
setting <- "linear" # one of "linear" and "nonlinear"

# 1. DATA ################################################

if (setting == "linear") {
  # 1.1 Linear: ############################################
  
  x_tr <- rnorm(n_tr)
  y_tr <- rnorm(n_tr, x_tr)
  x_ts <- rnorm(n_ts)
  y_ts <- rnorm(n_ts, x_ts)
  
  y_sets <- seq(min(y_tr), max(y_tr), 0.1)
  
  true_cdf <- function(y, x) pnorm(y, mean = x, sd = 1)
  
} else if (setting == "nonlinear") {
  # 1.2 Nonlinear: #########################################
  
  x_tr <- runif(n_tr, 0, 10)
  y_tr <- rnorm(n_tr, 2*x_tr + 5*sin(x_tr), abs(x_tr/5))
  x_ts <- runif(n_ts, 0, 10)
  y_ts <- rnorm(n_ts, 2*x_ts + 5*sin(x_ts), abs(x_ts/5))
  
  y_sets <- seq(min(y_tr), max(y_tr), 0.1)

  true_cdf <- function(y, x) pnorm(y, mean = 2*x + 5*sin(x), sd = abs(x/5))
}


#############################################################


# 2. PREDICTIVE BANDS #############################################

# RDCM Predictive Bands
rdcm_out <- rdcm_fa(x_tr, y_tr, x_ts, y_sets)

# Robust RDCM
rordcm_out <- rob_rdcm_fa_yprime(x_tr, y_tr, x_ts, y_sets)

# LSPM Predictive Bands
lspm_out <- lspm_fa(x_tr, y_tr, x_ts, y_sets)


# 3. PLOT ####################################################

library(ggplot2)

j <- sample(1:n_ts, 1)
t_cdf <- true_cdf(y_sets, x_ts[j])

df <- data.frame(x = y_sets,
                 l = c(t_cdf, rdcm_out[[j]]$PI_l, rordcm_out[[j]]$PI_l, lspm_out[[j]]$PI_l),
                 u = c(t_cdf, rdcm_out[[j]]$PI_u, rordcm_out[[j]]$PI_u, lspm_out[[j]]$PI_u),
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



# 4. EVALUATE #################################################

library(scoringRules)

get_predint <- function(int, alpha, ysets) {
  low <- ysets[int$PI_l <= alpha/2]
  upp <- ysets[int$PI_u >= 1 - alpha/2]
  if (length(low) == 0) {
    low <- min(ysets)
  } else {
    low <- max(low)
  }
  if (length(upp) == 0) {
    upp <- max(ysets)
  } else {
    upp <- min(upp)
  }  
  c(low, upp)
}

alpha_vec <- seq(0.5, 0.95, 0.01)

int_eval <- data.frame(alpha = alpha_vec, 
                       mth = rep(c("RDPS", "RDPS (rob.)", "LSPM"), each = length(alpha_vec)),
                       cov = NA, is = NA, wid = NA)
for (i in seq_along(alpha_vec)) {
  alpha <- 1 - alpha_vec[i]
  
  rdcm_int <- sapply(rdcm_out, get_predint, alpha, y_sets)
  rordcm_int <- sapply(rordcm_out, get_predint, alpha, y_sets)
  lspm_int <- sapply(lspm_out, get_predint, alpha, y_sets)
  
  int_eval$wid[int_eval$mth == "RDPS"][i] <- (rdcm_int[2, ] - rdcm_int[1, ]) |> mean()
  int_eval$wid[int_eval$mth == "RDPS (rob.)"][i] <- (rordcm_int[2, ] - rordcm_int[1, ]) |> mean()
  int_eval$wid[int_eval$mth == "LSPM"][i] <- (lspm_int[2, ] - lspm_int[1, ]) |> mean()
  
  int_eval$cov[int_eval$mth == "RDPS"][i] <- ((y_ts >= rdcm_int[1, ]) & (y_ts <= rdcm_int[2, ])) |> mean()
  int_eval$cov[int_eval$mth == "RDPS (rob.)"][i] <- ((y_ts >= rordcm_int[1, ]) & (y_ts <= rordcm_int[2, ])) |> mean()
  int_eval$cov[int_eval$mth == "LSPM"][i] <- ((y_ts >= lspm_int[1, ]) & (y_ts <= lspm_int[2, ])) |> mean()
  
  int_eval$is[int_eval$mth == "RDPS"][i] <- ints_sample(y_ts, t(rdcm_int), target_coverage = 1 - alpha) |> mean()
  int_eval$is[int_eval$mth == "RDPS (rob.)"][i] <- ints_sample(y_ts, t(rordcm_int), target_coverage = 1 - alpha) |> mean()
  int_eval$is[int_eval$mth == "LSPM"][i] <- ints_sample(y_ts, t(lspm_int), target_coverage = 1 - alpha) |> mean()
  
}

int_eval

## Interval score
ggplot(int_eval) +
  geom_line(aes(x = alpha, y = is, col = mth)) +
  scale_x_continuous(name = expression(alpha)) +
  scale_y_continuous(name = "Interval score", limits = c(2, 6)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

## Width
ggplot(int_eval) +
  geom_line(aes(x = alpha, y = wid, col = mth)) +
  scale_x_continuous(name = expression(alpha)) +
  scale_y_continuous(name = "Average width") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

## Coverage
ggplot(int_eval) +
  geom_line(aes(x = alpha, y = cov, col = mth)) +
  geom_abline(aes(slope = 1, intercept = 0), lty = "dotted") +
  scale_x_continuous(name = expression(alpha), limits = c(0.5, 1)) +
  scale_y_continuous(name = "Coverage", limits = c(0.2, 1)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")




