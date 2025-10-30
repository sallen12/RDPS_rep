################################################################################
##################### Utility functions for plotting ###########################
################################################################################

library(ggplot2)

## plot demonstration of the residual distribution forecasting procedure and predictive system
plot_demo <- function(x_tr, y_tr) {
  
  ## fit models
  
  y_hat_ols <- lm(y_tr ~ x_tr) |> predict()
  y_hat_krr <- {
    K <- exp(-abs(outer(x_tr, x_tr, "-"))) # uses laplacian kernel
    H <- K %*% solve(K + lambda*diag(length(x_tr)))
    H %*% y_tr |> as.vector()
  }
  
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

## plot predictive systems for different covariate values
plot_cov_vs_ps <- function(ps_list, y_eval, x_ts, setting) {
  if (setting == "linear") {
    # true conditional distribution function
    true_cdf <- function(y, x) pnorm(y, mean = x, sd = 1) 
    
    # covariate values at which to plot example predictive systems
    x_ind <- seq(-2, 2, 0.5)
    ind <- sapply(x_ind, function(x) which.min(abs(x_ts - x)))
    x_ind <- x_ts[ind]
    
  } else if (setting == "nonlinear") {
    # true conditional distribution function
    true_cdf <- function(y, x) pnorm(y, mean = 2*x + 5*sin(x), sd = abs(x/5))
    
    # covariate values at which to plot example predictive systems
    x_ind <- 1:9
    ind <- sapply(x_ind, function(x) which.min(abs(x_ts - x)))
    x_ind <- x_ts[ind]
  }
  
  t_cdf <- sapply(x_ind, true_cdf, y = y_sets)
  t_cdf <- rbind(t_cdf, t_cdf, t_cdf, t_cdf) |> c()
  Pi_l <- sapply(ind, function(j) c(ps_list[['rdps_ols']][[j]]$PI_l, 
                                    ps_list[['rdps_krr']][[j]]$PI_l, 
                                    ps_list[['lspm']][[j]]$PI_l, 
                                    ps_list[['krrpm']][[j]]$PI_l)) |> c()
  Pi_u <- sapply(ind, function(j) c(ps_list[['rdps_ols']][[j]]$PI_u, 
                                    ps_list[['rdps_krr']][[j]]$PI_u, 
                                    ps_list[['lspm']][[j]]$PI_u, 
                                    ps_list[['krrpm']][[j]]$PI_u)) |> c()
  df <- data.frame(y = y_eval,
                   l = Pi_l,
                   u = Pi_u,
                   x = rep(round(x_ind, 1), each = length(y_eval)*4),
                   mth = rep(c("RDPS (LS)", "RDPS (KRR)", "LSPM", "KRRPM"), each = length(y_eval)))
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
}

## plot example predictive systems at a particular covariate value
plot_example <- function(ps_list, y_eval, setting) {
  if (setting == "linear") {
    val <- 0
    true_cdf <- function(y, x) pnorm(y, mean = x, sd = 1) 
  } else if (setting == "nonlinear") {
    val <- 5
    true_cdf <- function(y, x) pnorm(y, mean = 2*x + 5*sin(x), sd = abs(x/5))
  }

  Pi_l <- sapply(ind, function(j) c(ps_list[['rdps_ols']][[j]]$PI_l, 
                                    ps_list[['rdps_krr']][[j]]$PI_l, 
                                    ps_list[['lspm']][[j]]$PI_l, 
                                    ps_list[['krrpm']][[j]]$PI_l)) |> c()
  Pi_u <- sapply(ind, function(j) c(ps_list[['rdps_ols']][[j]]$PI_u, 
                                    ps_list[['rdps_krr']][[j]]$PI_u, 
                                    ps_list[['lspm']][[j]]$PI_u, 
                                    ps_list[['krrpm']][[j]]$PI_u)) |> c()
  df <- data.frame(y = y_eval,
                   l = Pi_l,
                   u = Pi_u,
                   x = rep(round(x_ind, 1), each = length(y_eval)*4),
                   mth = rep(c("RDPS (LS)", "RDPS (KRR)", "LSPM", "KRRPM"), each = length(y_eval)))
  
  t_cdf <- true_cdf(y = y_eval, x = val)
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
}

## plot comparison of RDPS and deleted RDPS
plot_del_rdps <- function(ps_list, x_tr, y_tr, x_ts, y_eval) {
  rdcm_out <- rdps_lin(x_tr, y_tr, x_ts, y_eval)
  
  val <- 0
  t_cdf <- true_cdf(y = y_eval, x = val)
  Pi_l <- sapply(ind, function(j) c(ps_list[['rdps_ols']][[j]]$PI_l, 
                                    ps_list[['lspm']]$PI_l, 
                                    rdcm_out[[j]]$PI_l)) |> c()
  Pi_u <- sapply(ind, function(j) c(ps_list[['rdps_ols']]$PI_u, 
                                    ps_list[['lspm']][[j]]$PI_u, 
                                    rdcm_out[[j]]$PI_u)) |> c()
  df <- data.frame(y = y_eval,
                   l = Pi_l,
                   u = Pi_u,
                   x = rep(round(x_ind, 1), each = length(y_eval)*3),
                   mth = rep(c("RDPS (del)", "LSPM", "RDPS"), each = length(y_eval)))
  
  ggplot(df |> subset(x == val)) + 
    geom_line(aes(x = y, y = rep(t_cdf, 3)), col = "black") +
    geom_step(aes(x = y, y = l, col = mth)) +
    geom_step(aes(x = y, y = u, col = mth)) +
    scale_y_continuous(name = "G(y)") +
    scale_x_continuous(name = "y", limits = c(-3, 3)) +
    scale_color_manual(values = c("#7CAE00", "#E68613", "#C77CFF")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(), 
          legend.justification = c(1, 0),
          legend.position = c(0.99, 0.01))
  ggsave(paste0("Figures/ss_", setting, "_", n_tr, "_del.png"), width = 4, height = 3)
}

## plot interval score for prediction intervals
plot_is <- function(int_eval, setting) {
  if (setting == "linear") {lims <- c(2, 7)} else {lims <- c(2, 22)}
  ggplot(int_eval) +
    geom_line(aes(x = alpha, y = is, col = mth)) +
    scale_x_continuous(name = "Level", limits = c(0.5, 1)) +
    scale_y_continuous(name = "Interval score", limits = lims) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.justification = c(0, 1),
          legend.position = c(0.01, 0.99))  +
    guides(col = guide_legend(nrow = 2))
  ggsave(paste0("Figures/ss_", setting, "_", n_tr, "_is.png"), width = 3, height = 3)
}

## plot average width of prediction intervals
plot_width <- function(int_eval, setting) {
  if (setting == "linear") {lims <- c(0, 7)} else {lims <- c(0, 16)}
  ggplot(int_eval) +
    geom_line(aes(x = alpha, y = wid, col = mth)) +
    scale_x_continuous(name = "Level", limits = c(0.5, 1)) +
    scale_y_continuous(name = "Average width", limits = lims) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.justification = c(0, 1),
          legend.position = c(0.01, 0.99)) +
    guides(col = guide_legend(nrow = 2))
  ggsave(paste0("Figures/ss_", setting, "_", n_tr, "_wid.png"), width = 3, height = 3)
}

## plot unconditional coverage of prediction intervals
plot_cov <- function(int_eval, setting) {
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
}

## plot thickness of predictive systems
plot_thick <- function(ps_list, setting) {
  rdps_ols_th <- sapply(ps_list[['rdps_ols']], function(x) max(x$PI_u - x$PI_l))
  rdps_krr_th <- sapply(ps_list[['rdps_krr']], function(x) max(x$PI_u - x$PI_l))
  
  df <- data.frame(th = c(rdps_krr_th, rdps_ols_th), mth = rep(c("RDPS (KRR)", "RDPS (LS)"), each = n_ts))
  ggplot(df) +
    geom_histogram(aes(x = th, y = after_stat(density), fill = mth), alpha = 0.3, position = "identity", bins = 20) +
    scale_x_continuous(name = "Thickness") +
    scale_y_continuous(name = "", breaks = NULL, expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = c("#00BFC4", "#C77CFF")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.justification = c(0.5, 1),
          legend.position = c(0.5, 0.99)) +
    ggtitle("Nonlinear")
  ggsave(paste0("Figures/ss_", setting, "_", n_tr, "_th.png"), width = 3, height = 3.2)
}


