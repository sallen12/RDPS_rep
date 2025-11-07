################################################################################
## Brute force algorithm to calculate RDPS by fitting the model to a grid of y #
################################################################################

## wrapper to fit the RDPS
rdps <- function(x, y, x_new, y_prime, y_eval, mth = "ols", del = 1, ...) {
  if (length(x_new) > 1) {
    lapply(seq_along(x_new), function(i) {print(i); rdps(x, y, x_new[i], y_prime, y_eval, mth = mth, del = del, ...)})
  } else {
    out <- lapply(y_prime, function(yy) find_y_hat(x = c(x, x_new), y = c(y, yy), mth = mth, del = del, ...))
    G <- sapply(out, function(x) rd_G(x$y_hat, x$res, y_eval))
    PI_l <- apply(G, 1, min)
    PI_u <- apply(G, 1, max)
    return(list(PI_l = PI_l, PI_u = PI_u))
  }
}

## function to get residual distribution forecast distribution
rd_G <- function(y_hat, res, y_sets) {
  n <- length(res)
  G_y <- sapply(y_sets, function(y) sum(y_hat + res <= y)/n)
  return(G_y)
}

## function to get point-valued regression predictions
find_y_hat <- function(x, y, x_new = NULL, mth = "ols", del = 1, ...) {
  if (mth == "ols") {
    fit <- lm(y ~ x) 
    y_hat <- predict(fit)
    if (!is.null(x_new)) y_hat_new <- predict(fit, newdata = data.frame(x = x_new))
  } else if (mth == "krr") {
    K <- exp(-abs(outer(x, x, "-"))) # Laplacian kernel
    H <- K %*% solve(K + lambda*diag(length(x)))
    y_hat <- H %*% y |> as.vector()
    
    if (!is.null(x_new)) {
      K_new <- exp(-abs(outer(x_new, x, "-"))) # Laplacian kernel
      H <- K_new %*% solve(K + lambda*diag(length(x)))
      y_hat_new <- H %*% y |> as.vector()
    }
  }
  # deleted RDPS removes the del training points with the highest residuals and refits
  if (del) {
    ind <- order(abs(y_hat - y), decreasing = TRUE)[1:del]
    x_new <- x[length(x)]
    x <- x[-ind]
    y <- y[-ind]
    out <- find_y_hat(x, y, x_new = x_new, mth = mth, del = F, ...)
    return(out)
  } else {
    if (is.null(x_new)) y_hat_new <- y_hat[length(y_hat)]
    return(list(y_hat = y_hat_new, res = y - y_hat))
  }
}
