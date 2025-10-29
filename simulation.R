################################################################################
## Simulation Study

source("rdcm.R")
source("rordcm.R")
source("lspm.R")

set.seed(75210322)
n_tr <- 100
n_ts <- 1000
setting <- "linear" # one of "linear" and "nonlinear"


# 1. DATA ################################################

get_krr_lambda <- function(n = 100, setting = "linear") {
  l_vec <- 10^c(-6:2)
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
  
  K_tr <- exp(-abs(outer(x_tr, x_tr, "-"))) # uses laplacian kernel
  K_ts <- exp(-abs(outer(x_ts, x_tr, "-"))) # uses laplacian kernel
  rmse <- sapply(l_vec, function(lambda) {
    H <- K_ts %*% solve(K_tr + lambda*diag(n))
    y_hat <- H %*% y_tr |> as.vector()
    mean((y_ts - y_hat)^2)
  })
  l_vec[which.min(rmse)]
}

if (setting == "linear") {
  # 1.1 Linear: ############################################
  
  x_tr <- rnorm(n_tr)
  y_tr <- rnorm(n_tr, x_tr)
  x_ts <- rnorm(n_ts)
  y_ts <- rnorm(n_ts, x_ts)
  
  y_sets <- seq(-10, 10, 0.1)
  
  true_cdf <- function(y, x) pnorm(y, mean = x, sd = 1)
  
  x_ind <- seq(-2, 2, 0.5)
  ind <- sapply(x_ind, function(x) which.min(abs(x_ts - x)))
  x_ind <- x_ts[ind]
  
  lambda <- get_krr_lambda(n = n_tr, setting = "linear")
  
} else if (setting == "nonlinear") {
  # 1.2 Nonlinear: #########################################
  
  x_tr <- runif(n_tr, 0, 10)
  y_tr <- rnorm(n_tr, 2*x_tr + 5*sin(x_tr), abs(x_tr/5))
  x_ts <- runif(n_ts, 0, 10)
  y_ts <- rnorm(n_ts, 2*x_ts + 5*sin(x_ts), abs(x_ts/5))
  
  y_sets <- seq(0, 50, 0.1)

  true_cdf <- function(y, x) pnorm(y, mean = 2*x + 5*sin(x), sd = abs(x/5))
  
  x_ind <- 1:9
  ind <- sapply(x_ind, function(x) which.min(abs(x_ts - x)))
  x_ind <- x_ts[ind]
  
  lambda <- get_krr_lambda(n = n_tr, setting = "nonlinear")
}

example_plot <- function(x_tr, y_tr) {
  
  ## fit models

  y_hat_ols <- lm(y_tr ~ x_tr) |> predict()
  y_hat_krr <- {
    K <- exp(-abs(outer(x_tr, x_tr, "-"))) # uses laplacian kernel
    H <- K %*% solve(K + lambda*diag(length(x_tr)))
    H %*% y_tr |> as.vector()
  }

  library(ggplot2)
  
  ## plot predictions
  df <- data.frame(x = x_tr, y = c(y_tr, y_hat_ols, y_hat_krr), mth = rep(c("Truth", "LS", "KRR"), each = length(x_tr)))
  ggplot(df) +
    geom_point(aes(x = x, y = y, col = mth, shape = mth)) +
    scale_y_continuous(name = "y", limits = c(0, 30)) +
    scale_color_manual(values = c("#00BFC4", "#C77CFF", "black")) +
    scale_shape_manual(values = c(15, 17, 16, 18)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.justification = c(0, 1),
          legend.position = c(0.01, 0.99))
  ggsave("Figures/ex_data.png", width = 3, height = 3)
  
  ## plot residuals
  df <- data.frame(r = y_tr - c(y_hat_ols, y_hat_krr), mth = rep(c("LS", "KRR"), each = length(y_tr)))
  ggplot(df) +
    geom_histogram(aes(x = r, y = after_stat(density), fill = mth), alpha = 0.3, position = "identity", bins = 20) +
    scale_x_continuous(name = "Residuals", limits = c(-11, 11)) +
    scale_y_continuous(name = "", limits = c(0, 0.6), expand = c(0, 0)) +
    scale_fill_manual(values = c("#00BFC4", "#C77CFF")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.justification = c(0, 1),
          legend.position = c(0.01, 0.99))
  ggsave("Figures/ex_resdist.png", width = 3, height = 3)

  ## plot predictive distributions (at x_{n+1} = 2.5)
  ind <- which.min(abs(x_tr - 2.5))
  G_ls <- ecdf(c(y_hat_ols[ind] + (y_tr - y_hat_ols)))
  G_krr <- ecdf(c(y_hat_krr[ind] + (y_tr - y_hat_krr)))
  df <- data.frame(x = y_sets, G = c(G_ls(y_sets), G_krr(y_sets)), mth = rep(c("LS", "KRR"), each = length(y_sets)))
  ggplot(df) + geom_step(aes(x = x, y = G, col = mth)) +
    scale_x_continuous(name = "y", limits = c(0, 15.5), expand = c(0, 0)) +
    scale_y_continuous(name = "G(y)", limits = c(0, 1)) +
    scale_color_manual(values = c("#00BFC4", "#C77CFF")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.justification = c(0, 1),
          legend.position = c(0.01, 0.99))
  ggsave("Figures/ex_preddist.png", width = 3, height = 3)
  
  ## add predictive systems
  PI_ols <- rordcm_fa(x_tr, y_tr, x_tr[ind], del = F, y_sets, mth = "ols")
  PI_krr <- rordcm_fa(x_tr, y_tr, x_tr[ind], del = T, y_sets, mth = "krr", lambda = lambda)
  df$PI_l <- c(PI_ols$PI_l, PI_krr$PI_l)
  df$PI_u <- c(PI_ols$PI_u, PI_krr$PI_u)
  ggplot(df) + 
    geom_step(aes(x = x, y = G, col = mth)) +
    geom_step(aes(x = x, y = PI_l, col = mth), linetype = "dotted") +
    geom_step(aes(x = x, y = PI_u, col = mth), linetype = "dotted") +
    scale_x_continuous(name = "y", limits = c(0, 15.5), expand = c(0, 0)) +
    scale_y_continuous(name = "G(y)", limits = c(0, 1)) +
    scale_color_manual(values = c("#00BFC4", "#C77CFF")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.justification = c(0, 1),
          legend.position = c(0.01, 0.99))
  ggsave("Figures/ex_predsyst.png", width = 3, height = 3)
}
if (setting == "noninear") example_plot(x_tr, y_tr)


#############################################################


# 2. PREDICTIVE BANDS #############################################

# RDCM Predictive Bands
#rdcm_out <- rdcm_fa(x_tr, y_tr, x_ts, y_sets)
#rdcm_out_knn <- rdcm_fa(x_tr, y_tr, x_ts, y_sets, method = "knn", k = 10)

# Robust RDCM
#rordcm_out <- rordcm_fa(x_tr, y_tr, x_ts, y_sets)
#rordcm_out_kw <- rordcm_fa(x_tr, y_tr, x_ts, y_sets, mth = "kw")
#rordcm_out_kwt <- rordcm_fa(x_tr, y_tr, x_ts, y_sets, mth = "kw_trim")
#rordcm_out_knn <- rordcm_fa(x_tr, y_tr, x_ts, y_sets, mth = "knn")
#rordcm_out_qrf <- rordcm_fa(x_tr, y_tr, x_ts, y_sets, mth = "qrf")
rordcm_out_ols <- rordcm_fa(x_tr, y_tr, x_ts, y_sets, mth = "ols")
rordcm_out_krr <- rordcm_fa(x_tr, y_tr, x_ts, y_sets, mth = "krr", lambda = lambda)
#rordcm_out_qrn <- rordcm_fa(x_tr, y_tr, x_ts, y_sets, mth = "qr_new")

# LSPM Predictive Bands
lspm_out <- lspm_fa(x_tr, y_tr, x_ts, y_sets)
krrpm_out <- krrpm_fa(x_tr, y_tr, x_ts, y_sets, a = lambda)

 
# saveRDS(list(rdps_ols = rordcm_out_ols,
#             rdps_krr = rordcm_out_krr,
#             krrpm = krrpm_out,
#             lspm = lspm_out),
#        file = "results_linear_100.RDS")
# results <- readRDS("results_linear_100.RDS")


# 3. PLOT ####################################################

library(ggplot2)

t_cdf <- sapply(x_ind, true_cdf, y = y_sets)
t_cdf <- rbind(t_cdf, t_cdf, t_cdf, t_cdf) |> c()
Pi_l <- sapply(ind, function(j) c(rordcm_out_ols[[j]]$PI_l, rordcm_out_krr[[j]]$PI_l, lspm_out[[j]]$PI_l, krrpm_out[[j]]$PI_l)) |> c()
Pi_u <- sapply(ind, function(j) c(rordcm_out_ols[[j]]$PI_u, rordcm_out_krr[[j]]$PI_u, lspm_out[[j]]$PI_u, krrpm_out[[j]]$PI_u)) |> c()
df <- data.frame(y = y_sets,
                 l = Pi_l,
                 u = Pi_u,
                 x = rep(round(x_ind, 1), each = length(y_sets)*4),
                 mth = rep(c("RDPS (MR)", "RDPS (KRR)", "LSPM", "KRRPM"), each = length(y_sets)))
ggplot(df) + 
  geom_step(aes(x = y, y = l, col = mth)) +
  geom_step(aes(x = y, y = u, col = mth)) +
  geom_line(aes(x = y, y = t_cdf), col = "black") +
  scale_y_continuous(name = "G(y)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(), 
        legend.position = "bottom") +
  facet_wrap(~x)

ggsave(paste0("Figures/ss_", setting, "_", n_tr, ".png"), width = 6, height = 6)

if (setting == "linear") {val <- 0} else {val <- 5}
t_cdf <- true_cdf(y = y_sets, x = val)
ggplot(df |> subset(x == val)) + 
  geom_step(aes(x = y, y = l, col = mth)) +
  geom_step(aes(x = y, y = u, col = mth)) +
  geom_line(aes(x = y, y = rep(t_cdf, 4)), col = "black") +
  scale_y_continuous(name = "G(y)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(), 
        legend.position = "bottom")
ggsave(paste0("Figures/ss_", setting, "_", n_tr, "_1.png"), width = 5, height = 3.5)


# 4. EVALUATE #################################################

library(scoringRules)

get_predint <- function(int, alpha, ysets) {
  low <- ysets[int$PI_u <= alpha/2]
  upp <- ysets[int$PI_l >= 1 - alpha/2]
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

alpha_vec <- seq(0.5, 0.98, 0.01)

int_eval <- data.frame(alpha = alpha_vec, 
                       mth = rep(c("RDPS (KRR)", "RDPS (LS)", "LSPM", "KRRPM"), each = length(alpha_vec)),
                       cov = NA, is = NA, wid = NA)
for (i in seq_along(alpha_vec)) {
  alpha <- 1 - alpha_vec[i]
  
  rordcm_krr_int <- sapply(rordcm_out_krr, get_predint, alpha, y_sets)
  rordcm_int <- sapply(rordcm_out_ols, get_predint, alpha, y_sets)
  lspm_int <- sapply(lspm_out, get_predint, alpha, y_sets)
  krrpm_int <- sapply(krrpm_out, get_predint, alpha, y_sets)
  
  int_eval$wid[int_eval$mth == "RDPS (KRR)"][i] <- (rordcm_krr_int[2, ] - rordcm_krr_int[1, ]) |> mean()
  int_eval$wid[int_eval$mth == "RDPS (LS)"][i] <- (rordcm_int[2, ] - rordcm_int[1, ]) |> mean()
  int_eval$wid[int_eval$mth == "LSPM"][i] <- (lspm_int[2, ] - lspm_int[1, ]) |> mean()
  int_eval$wid[int_eval$mth == "KRRPM"][i] <- (krrpm_int[2, ] - krrpm_int[1, ]) |> mean()
  
  int_eval$cov[int_eval$mth == "RDPS (KRR)"][i] <- ((y_ts >= rordcm_krr_int[1, ]) & (y_ts <= rordcm_krr_int[2, ])) |> mean()
  int_eval$cov[int_eval$mth == "RDPS (LS)"][i] <- ((y_ts >= rordcm_int[1, ]) & (y_ts <= rordcm_int[2, ])) |> mean()
  int_eval$cov[int_eval$mth == "LSPM"][i] <- ((y_ts >= lspm_int[1, ]) & (y_ts <= lspm_int[2, ])) |> mean()
  int_eval$cov[int_eval$mth == "KRRPM"][i] <- ((y_ts >= krrpm_int[1, ]) & (y_ts <= krrpm_int[2, ])) |> mean()
  
  int_eval$is[int_eval$mth == "RDPS (KRR)"][i] <- ints_sample(y_ts, t(rordcm_krr_int), target_coverage = 1 - alpha) |> mean()
  int_eval$is[int_eval$mth == "RDPS (LS)"][i] <- ints_sample(y_ts, t(rordcm_int), target_coverage = 1 - alpha) |> mean()
  int_eval$is[int_eval$mth == "LSPM"][i] <- ints_sample(y_ts, t(lspm_int), target_coverage = 1 - alpha) |> mean()
  int_eval$is[int_eval$mth == "KRRPM"][i] <- ints_sample(y_ts, t(krrpm_int), target_coverage = 1 - alpha) |> mean()
  
}


## Interval score
ggplot(int_eval) +
  geom_line(aes(x = alpha, y = is, col = mth)) +
  scale_x_continuous(name = "Level", limits = c(0.5, 1)) +
  #scale_y_continuous(name = "Interval score", limits = c(2, 22)) +
  scale_y_continuous(name = "Interval score", limits = c(2, 7)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.justification = c(0, 1),
        legend.position = c(0.01, 0.99))  +
  guides(col = guide_legend(nrow = 2))
ggsave(paste0("Figures/ss_", setting, "_", n_tr, "_is.png"), width = 3, height = 3)


## Width
ggplot(int_eval) +
  geom_line(aes(x = alpha, y = wid, col = mth)) +
  scale_x_continuous(name = "Level", limits = c(0.5, 1)) +
  #scale_y_continuous(name = "Average width", limits = c(0, 16)) +
  scale_y_continuous(name = "Average width", limits = c(0, 7)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.justification = c(0, 1),
        legend.position = c(0.01, 0.99)) +
  guides(col = guide_legend(nrow = 2))
ggsave(paste0("Figures/ss_", setting, "_", n_tr, "_wid.png"), width = 3, height = 3)


## Coverage
ggplot(int_eval) +
  geom_line(aes(x = alpha, y = cov, col = mth)) +
  geom_abline(aes(slope = 1, intercept = 0), lty = "dotted") +
  scale_x_continuous(name = "Level", limits = c(0.5, 1)) +
  scale_y_continuous(name = "Coverage", limits = c(0.4, 1)) +
  theme_bw() +
  theme(legend.title = element_blank(),,
        legend.justification = c(1, 0),
        legend.position = c(0.99, 0.01))
ggsave(paste0("Figures/ss_", setting, "_", n_tr, "_cov.png"), width = 3, height = 3)


## thickness
rordcm_th <- sapply(rordcm_out_ols, function(x) max(x$PI_u - x$PI_l))
rordcm_krr_th <- sapply(rordcm_out_krr, function(x) max(x$PI_u - x$PI_l))
lspm_th <- sapply(lspm_out, function(x) max(x$PI_u - x$PI_l))
krrpm_th <- sapply(krrpm_out, function(x) max(x$PI_u - x$PI_l))
c(mean(rordcm_th), mean(rordcm_krr_th), mean(lspm_th), mean(krrpm_th))

df <- data.frame(th = c(rordcm_krr_th, rordcm_th), mth = rep(c("RDPS (KRR)", "RDPS (LS)"), each = n_ts))
ggplot(df) +
  geom_histogram(aes(x = th, y = after_stat(density), fill = mth), alpha = 0.3, position = "identity", bins = 20) +
  scale_x_continuous(name = "Thickness") +
  scale_y_continuous(name = "", breaks = NULL, expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("#00BFC4", "#C77CFF")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1, 1),
        legend.position = c(0.99, 0.99))
ggsave(paste0("Figures/ss_", setting, "_", n_tr, "_th.png"), width = 3, height = 3)

