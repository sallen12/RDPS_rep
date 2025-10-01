lspm_fa=function(X,yvect,xnew,ysets) {
  num=length(ysets)
  PI_l=numeric(num)
  PI_u=numeric(num)
  
  Xf <- c(xvect,xnew)
  H <- Xf %*% solve(t(Xf) %*% Xf) %*% t(Xf)
  n <- length(yvect)
  a=numeric(n)
  b=numeric(n)
  c=numeric(n+2)
  for (i in 1:n) {
    b[i]=sqrt(1-H[n+1,n+1])+H[i,n+1]/sqrt(1-H[i,i])
    a[i]=(H[1:n,n+1]%*%yvect)/sqrt(1-H[n+1,n+1])+(yvect[i]-H[i,1:n]%*%yvect)/sqrt(1-H[i,i])
    c[i]=a[i]/b[i]
  }
  c[n+1]=-1e16
  c[n+2]=1e16
  c_valid <- c[!is.na(c)]
  c_sorted <- sort(c_valid)
  
  for (j in 1:num) {
    y=ysets[j]
    for (i in 1:(n+1)){
      if (c_sorted[i]<y && y<c_sorted[i+1]) {
        PI_l[j]=(i-1)/(n+1) # because the paper starts with c[0] but R with c[1]
        PI_u[j]=(i)/(n+1) # because the paper starts with c[0] but R with c[1]
      }
      else if (y==c_sorted[i]) {
        idx <- which(c_sorted == c_sorted[i])
        i_prime <- idx[1]
        i_secondo <- idx[length(idx)]
        PI_l[j]=(i_prime-2)/(n+1) # because the paper starts with c[0] but R with c[1]
        PI_u[j]=(i_secondo)/(n+1) # because the paper starts with c[0] but R with c[1]
      }
    }
    
  }
  
  return(list(PI_l=PI_l,PI_u=PI_u))
}