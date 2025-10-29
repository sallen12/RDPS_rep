rdcm_fa <- function(X, yvect, xnew, y, method = "lr", k = NULL) {
  if (length(xnew) > 1) {
    lapply(seq_along(xnew), function(i) {print(i); rdcm_fa(X, yvect, xnew[i], y, method = method, k = k)})
  } else {
    
    xvect <- c(X, xnew)
    n <- length(yvect)
    
    out <- quantities(xvect, yvect, n, method = method, k = k)
    Sp <- out$Sp
    Sm <- out$Sm
    S0 <- out$S0
    a <- out$a
    b <- out$b
    
    b <- outer(b, y, "+")
    if (length(S0) > 1) {
      Ky <- colSums(b[S0, ] >= 0)
    } else if (length(S0) == 1) {
      Ky <- as.numeric(b[S0, ] >= 0)
    } else {
      Ky <- numeric(length(y))
    }
    
    # Compute c_j = b_j / a_j
    c_y <- b / a
    c_y <- c_y[a != 0, ]
    c_sort <- rbind(rep(-1e16, length(y)), c_y, rep(1e16, length(y)), y)
    c_sort <- apply(c_sort, 2, sort)
      
    # Compute representative points
    y_prime <- (c_sort[-1, ] + c_sort[-nrow(c_sort), ]) / 2
    
    # Compute G
    G_values <- sapply(1:nrow(y_prime), function(i) sapply(1:ncol(y_prime), function(j) G_x(Sp, Sm, Ky[j], y_prime[i, j], c_y[, j], y[j], n)))
    
    PI_l <- apply(G_values, 1, min, na.rm = T)
    PI_u <- apply(G_values, 1, max, na.rm = T)
  
    return(list(PI_l = PI_l, PI_u = PI_u))
  }
}

G_x <- function(Sp, Sm, Ky, y_prime, c, y, n) {
  
  term_Sp <- if (length(Sp) > 0) sum(y_prime <= c[Sp], na.rm = TRUE) else 0
  term_Sm <- if (length(Sm) > 0) sum(y_prime >= c[Sm], na.rm = TRUE) else 0
  term_y <- as.numeric(y_prime <= y)
  
  # Compute G(x)
  G_value <- (1 / (n + 1)) * (term_Sp + term_Sm + Ky + term_y)
  
  return(G_value)
}

get_ab_lr <- function(xvect, yvect) {
  n <- length(xvect) - 1
  H <- xvect %*% solve(t(xvect) %*% xvect) %*% t(xvect)
  # Compute a_j = h_{n+1,n+1} - h_{j,n+1}
  a <- H[n + 1, n + 1] - H[1:n, n + 1]
  
  # Compute b_j = - y_j - sum((h_{n+1,k} - h_{j,k}) y_k)
  b <- sapply(1:n, function(j) - yvect[j] - (((H[n + 1, -(n + 1)] - H[j, -(n + 1)])*yvect) |> sum()))
  return(list(a = a, b = b))
}

get_ab_knn <- function(xvect, yvect, k) {
  n <- length(xvect) - 1
  dif_x <- outer(xvect, xvect, "-") |> abs()
  ind <- sapply(1:(n + 1), function(i) order(dif_x[i, ])[1:k])
  ind[ind == n + 1] <- NA
  
  # Compute a_j = (1{n + 1 in N_{n+1}} - 1{n + 1 in N_i}) / k
  a <- (1 - apply(ind, 2, function(vec) any(is.na(vec)))) / k
  a <- a[-(n + 1)]
  
  # Compute b_j = y_j - sum_{j in N_{n+1}}(y_j)/K + sum_{j in N_i}(y_j)/K
  b <- sapply(1:n, function(i) yvect[i] + sum(yvect[ind[, i]], na.rm = T) / k)
  b <- b - sum(yvect[ind[, n + 1]], na.rm = T) / k

  return(list(a = a, b = b))
}

quantities <- function(xvect, yvect, n, method = "lr", k = NULL) {

  if (method == "lr") {
    ab <- get_ab_lr(xvect, yvect)
  } else if (method == "knn") {
    if (is.null(k)) stop("'k' must be specified for method = 'knn'")
    ab <- get_ab_knn(xvect, yvect, k = k)
  } else {
    stop("'method' not recognised. must be one of 'lr' and 'knn'")
  }
  a <- ab$a
  b <- ab$b
  
  # Sp, Sm, S0
  Sp <- which(a > 0)
  Sm <- which(a < 0)
  S0 <- which(abs(a) < 1e-10)
  
  return(list(Sp = Sp, Sm = Sm, S0 = S0, a = a, b = b))
}
