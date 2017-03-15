VAR = function(Z,p) {
  
  Z = as.matrix(Z)
  q = dim(Z)[2]
  n = dim(Z)[1]
  
  aux=embed(Z,p+1)
  Z=aux[,1:q]
  X=cbind(1,aux[,-c(1:q)])
  n=n-p
  
  beta = solve(t(X)%*%X,t(X)%*%Z)
  coef = list(C = t(beta)[,1])
  
  if (p>0) {
    AR = array(NA,dim = c(q,q,p))
    for (i in 1:p) AR[,,i] = t(beta)[,((i-1)*q+2):(i*q+1)]
    coef$AR = AR
  }
  
  e = Z - X%*%beta
  SSR = t(e)%*%e
  k = q + q^2*p
  V = SSR/(n-k)
  a = n*log(det(SSR/n)) + 2*k  
  b = n*log(det(SSR/n)) + k*log(n)
  LR = V
  
  if (p>0) {
    D = solve(diag(1,q)- apply(AR, c(1,2), sum))
    LR = t(D) %*% V %*% D
  }
  
  return(list(coef = coef, e = e, Cov = V, AIC = a, BIC = b, LR = LR, D=D))
}


neweywest = function(X,B=NULL, kernel.type = "QuadraticSpectral", prewhite = FALSE,VCOV.lag=1) {
  n=nrow(X)
  if (prewhite) {
    fit = VAR(X,VCOV.lag)
    X = fit$e
    J = fit$D
  }
  if (is.null(B)) B = autobw(X,kernel.type)
  T = dim(as.matrix(X))[1]
  V = cov(X,X)
  if(VCOV.lag!=0){
    for (k in 1:(VCOV.lag)) {
      D = kernel(k/B,kernel.type)*cov(X[(1 + k):nrow(X), ], X[1:(nrow(X)- k), ])
      V = V + D + t(D)
    }
  }
  
  if (prewhite) V = t(J)*V*J
  return(V)
}



kernel = function(x,type, normalize = TRUE) {
  if (normalize) ca = switch(type, Truncated = 2, Bartlett = 2/3, Parzen = 0.539285, 
                             TukeyHanning = 3/4, QuadraticSpectral = 1)
  else ca = 1
  
  switch(type, 
         Truncated = ifelse(abs(ca*x) > 1, 0, 1),
         Bartlett = ifelse(abs(ca*x) > 1, 0, 1 - abs(ca*x)),
         Parzen = ifelse(abs(ca*x) > 1, 0, ifelse(abs(ca*x) < 0.5, 1 - 6 * (ca* 
                                                                              x)^2 + 6 * abs(ca*x)^3, 2*(1 - abs(ca*x))^3)),
         TukeyHanning = ifelse(abs(ca*x) > 1, 0, (1 + cos(pi*ca*x))/2),
         QuadraticSpectral = {y =6*pi*x/5
         ifelse(abs(x) < 1e-04, 1, 3*(1/y)^2 * (sin(y)/y - cos(y)))})
}

autobw = function(X,type) {
  T = dim(as.matrix(X))[1]
  phi = 0.8
  alpha = 4*phi^2/(1+2*phi+phi^2)^2
  S = switch(type,
             Truncated = 0.6611*(alpha*T)^(1/5),
             Bartlett = 1.1447*(alpha*T)^(1/3),
             Parzen = 2.6614*(alpha*T)^(1/5), 
             TukeyHanning = 1.7462*(alpha*T)^(1/5),
             QuadraticSpectral = 1.3221*(alpha*T)^(1/5)
  )
}


VARHAC = function(X,lag_max) {
  M = 1e10
  for (i in 0:lag_max) {
    fit = VAR(X,i)
    if (fit$BIC<M) {
      M = fit$BIC
      b = i
      V = fit$LR
    }
  }
  return(V)
}
