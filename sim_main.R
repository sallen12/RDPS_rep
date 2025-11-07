################################################################################
## Simulation Study

source("cps.R")         # conformal predictive systems
source("rdps.R")        # RDPS
source("rdps_lin.R")    # efficient non-robust RDPS
source("sim_utils.R")   # utility functions for the analysis
source("sim_plots.R")   # utility functions for plotting

set.seed(75210322)

n_tr <- 100             # training data size
n_ts <- 1000            # test data size

setting <- "nonlinear"  # one of "linear" and "nonlinear"

load <- F               # load or generate and save new results


# 1. DATA ######################################################################


# generate data
if (setting == "linear") {
  # 1.1 Linear: ############################################
  
  x_tr <- rnorm(n_tr)           # training covariates
  y_tr <- rnorm(n_tr, x_tr)     # training observations
  x_ts <- rnorm(n_ts)           # test covariates
  y_ts <- rnorm(n_ts, x_ts)     # test observations
  
  y_prime <- seq(-10, 10, 0.1)  # range of y_prime values to consider 
  y_eval <- seq(-5, 5, 0.1)     # range of values at which to evaluate the predictive distributions
  
  # get regularisation parameter for kernel ridge regression
  lambda <- get_krr_lambda(n = n_tr, setting = "linear")
  
} else if (setting == "nonlinear") {
  # 1.2 Nonlinear: #########################################
  
  x_tr <- runif(n_tr, 0, 10)                                # training covariates
  y_tr <- rnorm(n_tr, 2*x_tr + 5*sin(x_tr), abs(x_tr/5))    # training observations
  x_ts <- runif(n_ts, 0, 10)                                # test covariates
  y_ts <- rnorm(n_ts, 2*x_ts + 5*sin(x_ts), abs(x_ts/5))    # test observations
  
  y_prime <- seq(0, 50, 0.1)    # range of y_prime values to consider 
  y_eval <- seq(0, 20, 0.1)     # range of values at which to evaluate the predictive distributions

  # get regularisation parameter for kernel ridge regression
  lambda <- get_krr_lambda(n = n_tr, setting = "nonlinear")
  
  # plot RDPS demonstation
  plot_demo(x_tr, y_tr, y_eval, y_prime)
}


# 2. PREDICTIVE SYSTEMS ########################################################

ps_list <- get_predsys(x_tr, y_tr, x_ts, y_prime, y_eval, lambda, setting, load)


# 3. PLOT ####################################################

# plot examples of predictive systems
plot_example(ps_list, y_eval, setting)

# plot predictive systems for different covariate values
plot_cov_vs_ps(ps_list, y_eval, x_ts, setting)

# plot robust (deleted) vs non-robust RDPS bounds
if (setting == "linear") plot_del_rdps(ps_list, x_tr, y_tr, x_ts, y_eval)


# 4. EVALUATE #################################################

# range of prediction levels
alpha_vec <- seq(0.5, 0.98, 0.01)

# evaluate prediction intervals
int_eval <- eval_predint(ps_list, y_ts, alpha_vec, y_eval)

# plot performance

plot_thick(ps_list, setting)     # Thickness

plot_cov(int_eval, setting)      # Coverage

plot_width(int_eval, setting)    # Width

plot_is(int_eval, setting)       # Interval score

