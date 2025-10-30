################################################################################
##################### Utility functions for analysis ###########################
################################################################################

library(scoringRules)

## function to get regularisation parameter for kernel ridge regression
get_krr_lambda <- function(n = 100, setting = "linear") {
  
  # range of lambda values to choose from
  l_vec <- 10^c(-6:2) 
  
  # generate independent data sets
  if (setting == "linear") {
    x_tr <- rnorm(n)
    y_tr <- rnorm(n, x_tr)
    x_ts <- rnorm(n)
    y_ts <- rnorm(n, x_ts)
  } else if (setting == "nonlinear") {
    x_tr <- runif(n, 0, 10)
    y_tr <- rnorm(n, 2*x_tr + 5*sin(x_tr), abs(x_tr/5))
    x_ts <- runif(n, 0, 10)
    y_ts <- rnorm(n, 2*x_ts + 5*sin(x_ts), abs(x_ts/5))
  }
  
  # calculate y_hat and calculate mean squared error
  K_tr <- exp(-abs(outer(x_tr, x_tr, "-"))) # Laplacian kernel
  K_ts <- exp(-abs(outer(x_ts, x_tr, "-"))) # Laplacian kernel
  mse <- sapply(l_vec, function(lambda) {
    H <- K_ts %*% solve(K_tr + lambda*diag(n))
    y_hat <- H %*% y_tr |> as.vector()
    mean((y_ts - y_hat)^2)
  })
  
  # choose the lambda values that minimises the MSE
  l_vec[which.min(mse)]
}

## function to get prediction intervals
get_predint <- function(int, alpha, y_eval) {
  low <- y_eval[int$PI_u <= alpha/2]
  upp <- y_eval[int$PI_l >= 1 - alpha/2]
  if (length(low) == 0) {
    low <- min(y_eval)
  } else {
    low <- max(low)
  }
  if (length(upp) == 0) {
    upp <- max(y_eval)
  } else {
    upp <- min(upp)
  }  
  c(low, upp)
}

## wrapper to evaluate prediction intervals
eval_predint <- function(ps_list, y_ts, alpha_vec, y_eval) {

  int_eval <- data.frame(alpha = alpha_vec, 
                         mth = rep(c("RDPS (LS)", "RDPS (KRR)", "LSPM", "KRRPM"), each = length(alpha_vec)),
                         cov = NA, is = NA, wid = NA)
  for (i in seq_along(alpha_vec)) {
    alpha <- 1 - alpha_vec[i]
    
    rdps_ols_int <- sapply(ps_list[['rdps_ols']], get_predint, alpha, y_eval)
    rdps_krr_int <- sapply(ps_list[['rdps_krr']], get_predint, alpha, y_eval)
    lspm_int <- sapply(ps_list[['lspm']], get_predint, alpha, y_eval)
    krrpm_int <- sapply(ps_list[['krrpm']], get_predint, alpha, y_eval)
    
    int_eval$wid[int_eval$mth == "RDPS (LS)"][i] <- (rdps_ols_int[2, ] - rdps_ols_int[1, ]) |> mean()
    int_eval$wid[int_eval$mth == "RDPS (KRR)"][i] <- (rdps_krr_int[2, ] - rdps_krr_int[1, ]) |> mean()
    int_eval$wid[int_eval$mth == "LSPM"][i] <- (lspm_int[2, ] - lspm_int[1, ]) |> mean()
    int_eval$wid[int_eval$mth == "KRRPM"][i] <- (krrpm_int[2, ] - krrpm_int[1, ]) |> mean()
    
    int_eval$cov[int_eval$mth == "RDPS (LS)"][i] <- ((y_ts >= rdps_ols_int[1, ]) & (y_ts <= rdps_ols_int[2, ])) |> mean()
    int_eval$cov[int_eval$mth == "RDPS (KRR)"][i] <- ((y_ts >= rdps_krr_int[1, ]) & (y_ts <= rdps_krr_int[2, ])) |> mean()
    int_eval$cov[int_eval$mth == "LSPM"][i] <- ((y_ts >= lspm_int[1, ]) & (y_ts <= lspm_int[2, ])) |> mean()
    int_eval$cov[int_eval$mth == "KRRPM"][i] <- ((y_ts >= krrpm_int[1, ]) & (y_ts <= krrpm_int[2, ])) |> mean()
    
    int_eval$is[int_eval$mth == "RDPS (LS)"][i] <- ints_sample(y_ts, t(rdps_ols_int), target_coverage = 1 - alpha) |> mean()
    int_eval$is[int_eval$mth == "RDPS (KRR)"][i] <- ints_sample(y_ts, t(rdps_krr_int), target_coverage = 1 - alpha) |> mean()
    int_eval$is[int_eval$mth == "LSPM"][i] <- ints_sample(y_ts, t(lspm_int), target_coverage = 1 - alpha) |> mean()
    int_eval$is[int_eval$mth == "KRRPM"][i] <- ints_sample(y_ts, t(krrpm_int), target_coverage = 1 - alpha) |> mean()
    
  }
  
  return(int_eval)
}
