# PROCEDURE:
# 1) I want to build lwr and upr bounds of G for a fixed y. And then extend that to all the different values of y in ysets
# 1.2) To do that I fix y=ysets[i] and I compute PI_l[i] and PI_u[i]. Then doing this for all the values in ysets, we get the bounds PL_l and PI_u
# 
# 2) To find the sup and inf of G(y) for a fixed y, we create a vector y_prime, and let it varies in R  
# 
# 3) Each value of y_prime, corresponds to a value of a_hat, that minimizes the loss (given that specific y_prime)
# 3.2) Hence we compute the vector of a_hat_vect (minimizing the loss)
#
# 4) For each a_hat in a_hat_vect (i.e, for each y_prime in y_prime_vect), we can compute the correspondent value of predictive distribution G(y) in the fixed point y
# 4.2) Then, letting a_hat varying in a_hat_vect (i.e, letting y_prime varying in R), we can find the lwr and upr values of the G(y) for that fixed y
# 4.3) We store these values that stands for PI_l(y) and PI_u(y)
#
# 5) We let the y varying in R and obtain teh whole predictive band

library(quantreg)

rordcm_fa <- function(x, y, x_new, y_sets, mth = "qr", del = T, lambda = NULL) {
  if (length(x_new) > 1) {
    lapply(seq_along(x_new), function(i) {print(i); rordcm_fa(x, y, x_new[i], y_sets, mth = mth, del = del, lambda = lambda)})
  } else {
    n <- length(y_sets)
    PI_l <- numeric(n)
    PI_u <- numeric(n)
    
    out <- lapply(seq_along(y_sets), function(i) find_y_hat(xx = c(x, x_new), yy = c(y, y_sets[i]), mth = mth, del = del, lambda = lambda))
    # sapply(1:100, function(i) all(diff(y_hat[101, ] - y_hat[i, ]) >= 0)) |> all() # check differences are increasing function of y'
    G <- sapply(seq_along(y_sets), function(i) get_G(out[[i]]$y_hat[length(out[[i]]$y_hat)], out[[i]]$res, y_sets))
    
    PI_l <- apply(G, 1, min)
    PI_u <- apply(G, 1, max)

    return(list(PI_l = PI_l, PI_u = PI_u))
  }
}

get_G <- function(y_hat, res, y_sets, lev = NULL) {
  n <- length(res)
  if (is.null(lev)) {
    G_y <- sapply(y_sets, function(yy) sum(y_hat + res <= yy)/n)
  } else {
    G_y <- sapply(y_sets, function(yy) sum(y_hat + (sqrt(1 - lev[length(lev)]) / sqrt(1 - lev[1:n])) * res <= yy)/n)
  }
  return(G_y)
}

find_y_hat <- function(xx, yy, mth = "qr", x_new = NULL, del = T, k = 10, rep = 100, lambda = NULL) {
  if (mth == "qr") {
    fit <- rq(yy ~ xx, tau = 0.5)
    y_hat <- predict(fit)
  } else if (mth == "qr_new") {
    y_hat <- qr(xx, yy, rep = rep)
  } else if (mth == "kw") {
    dif <- outer(xx, xx, "-") |> abs()
    dif <- exp(-5*dif)
    y_hat <- colSums(dif*yy) / colSums(dif)
  } else if (mth == "kw_trim") {
    yy <- pmax(pmin(yy, 5), -5)
    #yy <- pmax(pmin(yy, 30), 0)
    dif <- outer(xx, xx, "-") |> abs()
    dif <- exp(-5*dif)
    y_hat <- colSums(dif*yy) / colSums(dif)
  } else if (mth == "knn") {
    dif <- outer(xx, xx, "-") |> abs()
    neighs <- sapply(seq_along(xx), function(i) yy[which(dif[, i] %in% head(sort(dif[, i]), k))])
    y_hat <- colMeans(neighs)
  } else if (mth == "knn_med") {
    dif <- outer(xx, xx, "-") |> abs()
    neighs <- sapply(seq_along(xx), function(i) yy[which(dif[, i] %in% head(sort(dif[, i]), k))])
    y_hat <- apply(neighs, 1, median)
  } else if (mth == "qrf") {
    qrf <- quantregForest::quantregForest(x = as.matrix(xx), y = yy, ntree = 500, nodesize = 10)
    y_hat <- predict(qrf, as.matrix(xx), what = 0.5)
  } else if (mth == "ols") {
    fit <- lm(yy ~ xx)
    y_hat <- predict(fit)
  } else if (mth == "krr") {
    K <- exp(-abs(outer(xx, xx, "-"))) # uses laplacian kernel
    H <- K %*% solve(K + lambda*diag(length(xx)))
    y_hat <- H %*% yy |> as.vector()
    lev <- diag(H)
  }
  if (del) {
    ind <- order(abs(y_hat - yy), decreasing = TRUE)[1:5]
    x_new <- xx[length(xx)]
    xx <- xx[-ind]
    yy <- yy[-ind]
    out <- find_y_hat(xx, yy, mth = mth, del = F, k = k, rep = rep, lambda = lambda)
    res <- out$res
    lev <- out$lev
    out <- find_y_hat_predict(xx, yy, x_new, mth = mth, k = k, rep = rep, lambda = lambda)
    y_hat_new <- out$y_hat
    lev <- c(lev, out$lev)
    return(list(y_hat = y_hat_new, res = res, lev = lev))
  } else {
    if (!exists("lev")) lev <- NULL 
    return(list(y_hat = y_hat, res = yy - y_hat, lev = lev))
  }
  
}

find_y_hat_predict <- function(xx, yy, x_new, mth = "qr", k = 10, rep = 100, lambda = NULL) {
  if (mth == "qr") {
    fit <- rq(yy ~ xx, tau = 0.5)
    y_hat <- predict(fit, newdata = data.frame(xx = x_new))
  } else if (mth == "qr_new") {
    beta <- qr(xx, yy, rep = rep, all = F, predict = F)
    y_hat <- beta * x_new
  } else if (mth == "kw") {
    dif <- outer(x_new, xx, "-") |> abs() # check
    dif <- exp(-5*dif)
    y_hat <- colSums(dif*yy) / colSums(dif)
  } else if (mth == "kw_trim") {
    yy <- pmax(pmin(yy, 5), -5)
    #yy <- pmax(pmin(yy, 30), 0)
    dif <- outer(x_new, xx, "-") |> abs() # check
    dif <- exp(-5*dif)
    y_hat <- colSums(dif*yy) / colSums(dif)
  } else if (mth == "knn") {
    dif <- outer(x_new, xx, "-") |> abs() # check
    neighs <- sapply(seq_along(xx), function(i) yy[which(dif[, i] %in% head(sort(dif[, i]), k))])
    y_hat <- colMeans(neighs)
  } else if (mth == "knn_med") {
    dif <- outer(x_new, xx, "-") |> abs() # check
    neighs <- sapply(seq_along(xx), function(i) yy[which(dif[, i] %in% head(sort(dif[, i]), k))])
    y_hat <- apply(neighs, 1, median)
  } else if (mth == "qrf") {
    qrf <- quantregForest::quantregForest(x = as.matrix(xx), y = yy, ntree = 500, nodesize = 10)
    y_hat <- predict(qrf, as.matrix(x_new), what = 0.5)
  } else if (mth == "ols") {
    fit <- lm(yy ~ xx)
    y_hat <- predict(fit, newdata = data.frame(xx = x_new))
  } else if (mth == "krr") {
    K_tr <- exp(-abs(outer(xx, xx, "-")))
    K_ts <- exp(-abs(outer(x_new, xx, "-")))
    H <- K_ts %*% solve(K_tr + lambda*diag(length(xx)))
    y_hat <- H %*% yy |> as.vector()
    lev <-  H[length(H)]
  }
  if (!exists("lev")) lev <- NULL 
  return(list(y_hat = y_hat, lev = lev))
}

qr <- function(x_tr, y_tr, x_ts = x_tr, rep = 1, beta0 = 2, tau = 0.5, all = F, predict = T) {
  if (rep == 1) {
    n <- length(x_tr)
    ind <- sample(1:n)
    x_tr <- x_tr[ind]
    y_tr <- y_tr[ind]
    beta <- numeric(n)
    eta <- 1/n
    beta[1] <- beta0 - eta * (tau - (y_tr[1] > beta0*x_tr[1])) * x_tr[1]
    for (i in 2:n) beta[i] <- beta[i-1] - eta * (tau - (y_tr[i] > beta[i-1]*x_tr[i])) * x_tr[i]
    if (predict) {
      y_hat <- beta[n]*x_ts
      if (all) {
        return(list(pred = y_hat, par = beta))
      } else {
        return(y_hat)
      }
    } else {
      if (all) {
        return(beta)
      } else {
        return(beta[n])
      }
    }
  } else {
    out <- sapply(1:rep, qr, x_tr = x_tr, y_tr = y_tr, rep = 1, beta0 = beta0, tau = tau, all = T, predict = F)
    n <- length(x_tr)
    beta <- out[n, ] |> mean()
    if (predict) {
      y_hat <- beta*x_ts
      if (all) {
        return(list(pred = y_hat, par = beta))
      } else {
        return(y_hat)
      }
    } else {
      if (all) {
        return(out)
      } else {
        return(beta)
      }
    }
  }
}
