rdcm_fa <- function(X, yvect, xnew, ysets) {
  if (length(xnew) > 1) {
    lapply(seq_along(xnew), function(i) {print(i); rdcm_fa(X, yvect, xnew[i], ysets)})
  } else {
    num <- length(ysets)
    PI_l <- numeric(num)
    PI_u <- numeric(num)
    
    for (i in 1:num) {
      y <- ysets[i]
      rp <- rdcm_point(X, yvect, xnew, y)
      PI_l[i] <- rp$PI_l
      PI_u[i] <- rp$PI_u
    }
    return(list(PI_l = PI_l, PI_u = PI_u))
  }
  

}

rdcm_point <- function(X, yvect, xnew, y) {
  Xf <- c(X,xnew)
  H <- Xf %*% solve(t(Xf) %*% Xf) %*% t(Xf)
  n <- length(yvect)
  
  quantities_result <- quantities(H, yvect, y, n)
  Sp <- quantities_result$Sp
  Sm <- quantities_result$Sm
  Ky <- quantities_result$Ky
  c <- quantities_result$c
  c_new=c(c,-1e16,1e16)
  
  # Compute c_j = b_j / a_j
  c_valid <- c_new[!is.na(c_new)]
  c_all<-c(c_valid,y)
  c_sorted <- sort(c_all)
  
  # Compute representative points
  rep <- numeric(length(c_sorted) - 1)
  if (length(c_sorted) > 1) {
    rep <- (c_sorted[-length(c_sorted)] + c_sorted[-1]) / 2
  }
  
  # Compute G
  G_values <- numeric(length(rep))
  if (length(rep) > 0) {
    for (i in 1:length(rep)) {
      G_values[i] <- G(Sp, Sm, Ky, rep[i], c, y, n)
    }
  }
  
  PI_l <- if (length(G_values) > 0) min(G_values) else NA
  PI_u <- if (length(G_values) > 0) max(G_values) else NA
  
  return(list(PI_l = PI_l, PI_u = PI_u))
}

G <- function(Sp, Sm, Ky, x, c, y, n) {
 
  term_Sp <- if (length(Sp) > 0) sum(x <= c[Sp], na.rm = TRUE) else 0
  term_Sm <- if (length(Sm) > 0) sum(x >= c[Sm], na.rm = TRUE) else 0
  term_y <- as.numeric(x <= y)
  
  # Compute G(x)
  G_value <- (1 / (n + 1)) * (term_Sp + term_Sm + Ky + term_y)
  
  return(G_value)
}

quantities <- function(H, yvect, y, n) {
  # Compute a_j = h_{n+1,n+1} - h_{j,n+1}
  a <- H[nrow(H), nrow(H)] - H[1:n, nrow(H)]
  
  # Compute b_j = y - y_j - sum((h_{n+1,k} - h_{j,k}) y_k)
  b=numeric(n)
  for (j in 1:n) {
    q=0
    for (k in 1:n) {
      q=q+(H[nrow(H),k]-H[j,k])*yvect[k]
    }
    b[j]=y-yvect[j]-q
  }

  # Compute c_j = b_j / a_j
  c <- ifelse(a == 0, NA, b / a)
  
  # Sp, Sm, S0
  Sp <- which(a > 0)
  Sm <- which(a < 0)
  S0 <- which(abs(a) < 1e-10)
  
  # K(y)
  Ky <- sum(b[S0] >= 0)
  
  return(list(Sp = Sp, Sm = Sm, Ky = Ky, c = c))
}