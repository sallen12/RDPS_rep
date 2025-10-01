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


rob_rdcm_fa_yprime <- function(xvect,yvect,xnew,ysets) {
  num=length(ysets)
  PI_l=numeric(num)
  PI_u=numeric(num)
  
  for (i in 1:num) {
    y=ysets[i]
    rp=rob_rdcm_point(xvect,yvect,xnew,y)
    PI_l[i]=rp$PI_l
    PI_u[i]=rp$PI_u
  }
  return(list(rob_PI_l=PI_l,rob_PI_u=PI_u))
}

rob_rdcm_point <- function(xvect,yvect,xnew,y) {
  y_prime_vect=seq(min(yvect)-100,max(yvect)+100,0.1)
  n=length(y_prime_vect)
  a_hat_vect=numeric(n) #-> I'll put the values of a_hat s.t (given the correspondent y_prime), minimize the loss
  
  # Now I start the minimization of the loss to fill the a_hat_vect
  for (i in 1:n) {
    y_prime=y_prime_vect[i]
    a_hat_vect[i]=find_a_hat(xvect,yvect,xnew,y_prime)
  }

  # For fixed y, in order to compute sup and inf of G(y), we compute G(y) for each different y_prime in R (i.e, for each values in a_hat_vect) and then take the max and min values reached as PI_l(y) and PI_u(y)
  G_values=numeric(n)
  for (i in 1:n) {
    a_hat=a_hat_vect[i] #-> the variation of y_prime is the variation of the a_hat
    G_values[i]=rob_G(a_hat,xvect,yvect,xnew,y)
  }
  
  PI_l=min(G_values)
  PI_u=max(G_values)
  
  return(list(PI_l=PI_l,PI_u=PI_u))
}

rob_G <- function(a_hat,xvect,yvect,xnew,y) {
  y_hat=a_hat*xnew
  res=yvect-a_hat*xvect
  
  # compute G(y)
  G_value=1/(n+1)*sum(y_hat+res<=y)
  return(G_value)
}

find_a_hat <- function(xvect,yvect,xnew,y_prime) {
  
  # full approach
  xx=c(xvect,xnew)
  yy=c(yvect,y_prime)
  nn=length(yy)
  
  a_range=compute_a_range(xx,yy,nn)
  nnn=length(a_range)
  res_loss=numeric(nnn)
  
  for (i in 1:nnn) {
    a=a_range[i]
    res_loss[i]=sum(abs(yy-a*xx))
  }
  
  # choose the a_hat which minimizes the loss
  j=which.min(res_loss)
  a_hat=a_range[j]
  return(a_hat)
}

compute_a_range <- function(xx,yy,nn) {
  vect_a=numeric(nn)
  
  for (i in 1:nn) {
    vect_a[i]=yy[i]/xx[i]
  }
  
  min_a=min(vect_a)
  max_a=max(vect_a)
  a_range=seq(min_a,max_a,0.1)
  return(a_range)
}